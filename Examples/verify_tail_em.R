# Diagnose tail-mass changes under EM for right-censored data.
# Run from the repository root:
#   source("Examples/verify_tail_em.R")
#   Rscript Examples/verify_tail_em.R
# These commands recompile/load the current checkout with devtools.
# Each fit prints its diagnostics immediately; results are not written to files.
# No package implementation or default fitting settings are changed.
devtools::load_all(".", recompile = TRUE)

# Editable simulation and diagnostic settings.
n <- 500L
K <- 1
Lambda <- 0.12
q <- 0.3
max_depth <- 8L
iterations <- c(0:100, 200L, 500L, 1000L)
em_tolerance <- 1e-12
epsilon <- 1e-6

# At these defaults, seed 123 ends with an exact event and seed 1 with censoring.
# Administrative censoring supplies a counterexample with substantial tail mass.
set.seed(123)
event_times <- rweibull(n, shape = K, scale = 1 / Lambda)
cutoff <- 4
administrative <- data.frame(
  L = pmin(event_times, cutoff),
  R = ifelse(event_times <= cutoff, event_times, Inf),
  status = as.integer(event_times <= cutoff)
)
datasets <- list(
  random_seed123 = data_rc(123, n, K, Lambda, q),
  random_seed1 = data_rc(1, n, K, Lambda, q),
  administrative = administrative
)

clip <- function(x) pmin(1 - epsilon, pmax(epsilon, x))
rows <- summaries <- list()
for (case in names(datasets)) {
  data <- datasets[[case]]
  stopifnot(nrow(data) > 1, all(data$L == data$R | is.infinite(data$R)))
  censored <- is.infinite(data$R)
  tau <- max(data$L)
  censored_at_tau <- sum(censored & data$L == tau)
  tail_lower_bound <- clip(censored_at_tau / nrow(data))
  km <- survival::survfit(survival::Surv(L, !censored) ~ 1, data = data)
  km_tail <- tail(as.numeric(km$surv), 1)
  rate <- sum(!censored) / sum(data$L) # Exact exponential MLE for these data.
  bases <- list(uniform = list(prior = 1, params = NULL),
                exponential = list(prior = 2, params = rate))
  grid <- c(data$L[censored], tau)

  for (base in names(bases)) {
    spec <- bases[[base]]
    previous <- NULL
    trajectory <- vector("list", length(iterations))
    for (j in seq_along(iterations)) {
      k <- iterations[j]
      fit <- survPT(data[, c("L", "R")], test_seq = grid,
                    prior = spec$prior, prior_params = spec$params,
                    support = c(0, Inf), max_M = max_depth, allocated = 1,
                    max_EM_iter = k, EM_tol = em_tolerance, verbose = FALSE)
      p_tail <- as.numeric(fit$Posterior_Probs[3])
      n_tail <- as.numeric(fit$Allocations[3])
      total_mass <- sum(fit$Posterior_Probs[fit$Terminal_Nodes])
      survival_tau <- tail(as.numeric(fit$Survival_Probs), 1)
      cat(sprintf(paste0("[%s / %s] max_EM_iter=%d, actual=%d, ",
                         "body=%.10g, tail=%.10g, tail_n=%.10g, total=%.10g, ",
                         "logLik=%.10g, converged=%s\n"),
                  case, base, k, fit$EM_Iterations,
                  as.numeric(fit$Posterior_Probs[2]), p_tail, n_tail,
                  total_mass, fit$LogLik, fit$Converged))
      flush.console()
      # Beta(1,1) mode-style update at the added root:
      # p_tail = clip(n_tail / n), with epsilon bounds, also at k = 0.
      stopifnot(abs(p_tail - clip(n_tail / nrow(data))) < 1e-10,
                abs(total_mass - 1) < 1e-10,
                abs(sum(fit$Allocations[2:3]) - nrow(data)) < 1e-8,
                abs(survival_tau - p_tail) < 1e-10,
                p_tail >= tail_lower_bound - 1e-10,
                p_tail <= 1 - epsilon + 1e-10)

      # Independently verify initialization and consecutive E steps.
      # A censored row at C contributes p_tail / S(C) to the next tail count.
      e_step_error <- NA_real_
      if (k == 0L) {
        conditional_s <- if (spec$prior == 1) (tau - data$L[censored]) / tau else
          (exp(-rate * data$L[censored]) - exp(-rate * tau)) / (-expm1(-rate * tau))
        expected_n_tail <- sum(1 / (1 + conditional_s)) # Prior mean p = 1/2.
        e_step_error <- n_tail - expected_n_tail
      } else if (!is.null(previous) && k == iterations[j - 1L] + 1L) {
        old_survival <- head(as.numeric(previous$Survival_Probs), -1)
        stopifnot(all(old_survival > 0))
        expected_n_tail <- sum(as.numeric(previous$Posterior_Probs[3]) / old_survival)
        e_step_error <- n_tail - expected_n_tail
      }
      stopifnot(is.na(e_step_error) || abs(e_step_error) < 1e-8)
      trajectory[[j]] <- data.frame(
        case, base, iteration_limit = k, actual_iterations = fit$EM_Iterations,
        converged = fit$Converged, n = nrow(data), censoring_fraction = mean(censored),
        tau, censored_at_tau, tail_lower_bound, km_tail,
        body_mass = as.numeric(fit$Posterior_Probs[2]), tail_mass = p_tail,
        tail_allocations = n_tail, terminal_mass_sum = total_mass,
        survival_tau, loglik = fit$LogLik, e_step_error
      )
      previous <- fit
    }
    trajectory <- do.call(rbind, trajectory)
    first <- trajectory[1, ]
    last <- tail(trajectory, 1)
    floor_steps <- trajectory$actual_iterations[trajectory$tail_mass <= epsilon + 1e-10]
    summary_row <- data.frame(
      case, base, n = nrow(data), censoring_fraction = mean(censored), tau,
      censored_at_tau, initial_tail_mass = first$tail_mass,
      final_tail_mass = last$tail_mass,
      fraction_initial_tail_retained = last$tail_mass / first$tail_mass,
      final_body_mass = last$body_mass, tail_lower_bound, km_tail,
      true_survival_tau = pweibull(tau, K, 1 / Lambda, lower.tail = FALSE),
      final_at_floor = last$tail_mass <= epsilon + 1e-10,
      first_recorded_floor_iteration = if (length(floor_steps)) min(floor_steps) else NA_integer_,
      final_actual_iterations = last$actual_iterations, converged = last$converged
    )
    rows[[length(rows) + 1L]] <- trajectory
    summaries[[length(summaries) + 1L]] <- summary_row
    cat(sprintf("%s / %s: tail %.8f -> %.8f; censored at tau = %d; converged = %s\n",
                case, base, first$tail_mass, last$tail_mass, censored_at_tau, last$converged))
  }
}
trajectory <- do.call(rbind, rows)
summary_table <- do.call(rbind, summaries)
print(summary_table, row.names = FALSE, digits = 7)
cat("\nAll mass-conservation, root-update, and checked E-step identities passed.\n")
cat("Iteration 0 is AFTER initial allocation and its mode-style mass update.\n")
cat("Initialization uses prior mean 1/2; the returned initial tail mass need not be 1/2.\n")
cat("All split probabilities are bounded to [1e-6, 1 - 1e-6].\n")
cat("A small tail is not universal: censoring at tau enforces the reported lower bound.\n")
