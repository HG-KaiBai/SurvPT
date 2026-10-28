# Single-dataset left-truncated, interval-censored survival comparison.
# Run from the repository root: source("Examples/test_ltrunc.R")
# Requires devtools and the authors' ltic package (Rtools on Windows).
# Reference: Waller et al. (2026), Statistics and Computing, Section 5.1.
# https://doi.org/10.1007/s11222-026-10984-9
# https://github.com/ZacharyWaller/ltic
# Tested ltic revision: 33ccd31a5920f29d65ced2b6c840a405ce07b9f7
# One-time install, if needed (the script itself never installs packages):
# dir.create(".pt-validation/ltrunc-library", recursive = TRUE, showWarnings = FALSE)
# remotes::install_github(
#   "ZacharyWaller/ltic@33ccd31a5920f29d65ced2b6c840a405ce07b9f7",
#   lib = ".pt-validation/ltrunc-library", upgrade = "never")

local_library <- file.path(".pt-validation", "ltrunc-library")
if (dir.exists(local_library)) .libPaths(c(normalizePath(local_library), .libPaths()))
if (!requireNamespace("devtools", quietly = TRUE) ||
    !requireNamespace("ltic", quietly = TRUE)) {
  stop("Install devtools and ltic first; see the installation commands at the top of this file.")
}
devtools::load_all(".", recompile = TRUE)
if (!"entry" %in% names(formals(survPT))) stop("This example requires survPT(entry = ...).")

# Change these settings to explore another single dataset; no seed searching.
seed <- 123L
n_target <- 500L                  # Expected retained size, not a fixed-size sample
trunc_max <- 0.6                  # Paper settings: 0, 0.6, 1.6
inspection_max <- 4               # Published Section 5.1: U(0, 4)
inspection_gap <- 0.5
true_rate <- 1
max_M <- 8L
pt_max_iter <- 2000L
pt_tol <- 1e-8
npmle_max_iter <- 100000L
npmle_tol <- 1e-5                 # Paper's log-likelihood stopping tolerance
output_dir <- file.path("results", "ltrunc")

# The authors' analysis/simulation.R at the pinned revision uses U(0, 2)
# for the first assessment, whereas the published text specifies U(0, 4).
# This example follows the published text; set inspection_max = 2 to use
# that repository setting. All estimators below receive the same retained rows.
simulate_ltrunc <- function(n_target, trunc_max, inspection_max, gap, rate, seed) {
  stopifnot(n_target > 0, trunc_max >= 0, inspection_max > 0, gap > 0, rate > 0)
  keep_probability <- if (trunc_max == 0) 1 else
    -expm1(-rate * trunc_max) / (rate * trunc_max)
  n_source <- ceiling(n_target / keep_probability)
  set.seed(seed)
  entry <- runif(n_source, 0, trunc_max)
  event <- rexp(n_source, rate)
  first <- runif(n_source, 0, inspection_max)
  second <- first + gap
  # 0 is the finite lower support bound used by SurvPT; ltic accepts it too.
  left <- ifelse(event <= first, 0, ifelse(event <= second, first, second))
  right <- ifelse(event <= first, first, ifelse(event <= second, second, Inf))
  keep <- event > entry
  observations <- data.frame(L = left[keep], R = right[keep], entry = entry[keep])
  stopifnot(nrow(observations) > 1, all(observations$R > observations$entry))
  # Latent event times are used only to generate/filter observations, never to fit.
  list(data = observations, n_source = n_source, keep_probability = keep_probability)
}
generated <- simulate_ltrunc(n_target, trunc_max, inspection_max, inspection_gap,
                             true_rate, seed)
test_data <- generated$data
n <- nrow(test_data)
truncated_fraction <- 1 - n / generated$n_source
cat(sprintf("Seed %d: retained %d of %d subjects (%.1f%% truncated).\n",
            seed, n, generated$n_source, 100 * truncated_fraction))
cat(sprintf("Retained observations: %d left-, %d interval-, %d right-censored.\n",
            sum(test_data$L == 0), sum(test_data$L > 0 & is.finite(test_data$R)),
            sum(is.infinite(test_data$R))))

# Use the truncated interval likelihood to estimate the exponential base rate.
# estimate_exp(L, R) currently ignores entry and must not be used here.
fit_exp_ltrunc <- function(data) {
  left <- pmax(data$L, data$entry)
  width <- data$R - left
  exact <- data$L == data$R
  nll <- function(log_rate) {
    rate <- exp(log_rate)
    log_probability <- -rate * (left - data$entry)
    log_probability[!exact] <- log_probability[!exact] +
      log(-expm1(-rate * width[!exact]))
    log_probability[exact] <- log_probability[exact] + log_rate
    -sum(log_probability)
  }
  bounds <- log(c(1e-4, 1e4))
  fit <- optimize(nll, interval = bounds, tol = 1e-10)
  if (!is.finite(fit$objective) || min(abs(fit$minimum - bounds)) < 1e-5)
    stop("The left-truncated exponential fit failed or reached its search boundary.")
  exp(fit$minimum)
}
exp_par <- fit_exp_ltrunc(test_data)
cat(sprintf("Exponential base rate estimated with entry correction: %.6f\n", exp_par))

# Stop at tau because SurvPT leaves the shape beyond the finite body unspecified.
endpoints <- c(test_data$L, test_data$R, test_data$entry)
tau <- max(endpoints[is.finite(endpoints)])
# Fixed by the design, not chosen after seeing the fitted curves. The paper's
# Section 7 discusses early NPMLE drops and conditional survival estimates.
condition_time <- trunc_max
stopifnot(condition_time < tau)
grid <- sort(unique(c(seq(0, tau, length.out = 400L), condition_time)))
timed <- function(expr) {
  start <- proc.time()[["elapsed"]]
  fit <- force(expr)
  list(fit = fit, seconds = proc.time()[["elapsed"]] - start)
}
cat("Fitting SurvPT with uniform base...\n")
run_unif <- timed(survPT(test_data[, c("L", "R")], entry = test_data$entry,
  test_seq = grid, prior = 1, max_M = max_M, allocated = 1,
  max_EM_iter = pt_max_iter, EM_tol = pt_tol, verbose = FALSE))
result_unif <- run_unif$fit
cat("Fitting SurvPT with exponential base...\n")
run_exp <- timed(survPT(test_data[, c("L", "R")], entry = test_data$entry,
  test_seq = grid, prior = 2, prior_params = exp_par, max_M = max_M,
  allocated = 1, max_EM_iter = pt_max_iter, EM_tol = pt_tol, verbose = FALSE))
result_exp <- run_exp$fit

# Use the authors' native implementations. open_T = TRUE matches SurvPT's
# strict T > entry convention; the paper/package defaults to closed entry.
# With continuous generating distributions, equality has probability zero.
# Do not use interval::icfit here: it does not adjust for left truncation.
cat("Fitting the authors' PL-ICM NPMLE...\n")
run_pl_icm <- timed(ltic::ltic_np(test_data$L, test_data$R, trunc = test_data$entry,
  method = "PL-ICM", open_T = TRUE, tol = npmle_tol, max_it = npmle_max_iter))
result_pl_icm <- run_pl_icm$fit
cat("Fitting the authors' Product-Limit NPMLE...\n")
run_pl <- timed(ltic::ltic_np(test_data$L, test_data$R, trunc = test_data$entry,
  method = "Product-Limit", open_T = TRUE, tol = npmle_tol, max_it = npmle_max_iter))
result_pl <- run_pl$fit

# Evaluate the right-endpoint step convention used by the authors' plot.ltic.
# The last inner interval can end at Inf; its mass is retained at finite times.
ltic_survival <- function(fit, times, condition = NULL) {
  right <- fit$intervals$II$right
  origin <- if (is.null(condition)) 0L else findInterval(condition, right)
  steps <- switch(fit$type,
    expo = exp(-as.numeric(fit$res$lambda) + fit$res$lambda[origin + 1L]),
    prodlim = {
      hazards <- as.numeric(fit$res$lambda)
      if (origin > 0L) hazards[seq_len(origin)] <- 0
      c(1, cumprod(1 - hazards))
    },
    stop("Unexpected ltic output parameterization."))
  stopifnot(length(steps) == length(right) + 1L, !is.unsorted(right))
  as.numeric(steps[findInterval(times, right) + 1L])
}
curves <- data.frame(time = grid, truth = exp(-true_rate * grid),
  PT_uniform = as.numeric(result_unif$Survival_Probs),
  PT_exponential = as.numeric(result_exp$Survival_Probs),
  PL_ICM = ltic_survival(result_pl_icm, grid),
  PL = ltic_survival(result_pl, grid))
for (name in names(curves)[-1L]) {
  s <- curves[[name]]
  stopifnot(all(is.finite(s)), all(s >= -1e-10 & s <= 1 + 1e-10),
            all(diff(s) <= 1e-10))
}
stopifnot(is.finite(result_unif$LogLik), is.finite(result_exp$LogLik),
          is.finite(result_pl_icm$res$llike), is.finite(result_pl$res$llike))

# Independently check the PL-ICM endpoint extraction against its likelihood.
# Cumulative-hazard differences remain usable after a very large early drop.
log_s_icm <- function(t) -result_pl_icm$res$lambda[
  findInterval(t, result_pl_icm$intervals$II$right) + 1L]
log_left <- log_s_icm(pmax(test_data$L, test_data$entry))
checked_loglik <- sum(log_left - log_s_icm(test_data$entry) +
  log(-expm1(log_s_icm(test_data$R) - log_left)))
stopifnot(abs(checked_loglik - result_pl_icm$res$llike) < 1e-6)

# Condition EVERY estimator and the truth on the same, prespecified time.
# For PL, recompute the product after that time; do not divide two zeros.
after <- grid >= condition_time
at_origin <- match(condition_time, grid)
stopifnot(curves$PT_uniform[at_origin] > 0, curves$PT_exponential[at_origin] > 0)
conditional_curves <- data.frame(time = grid[after],
  truth = exp(-true_rate * (grid[after] - condition_time)),
  PT_uniform = curves$PT_uniform[after] / curves$PT_uniform[at_origin],
  PT_exponential = curves$PT_exponential[after] / curves$PT_exponential[at_origin],
  PL_ICM = ltic_survival(result_pl_icm, grid[after], condition_time),
  PL = ltic_survival(result_pl, grid[after], condition_time))
for (name in names(conditional_curves)[-1L]) {
  s <- conditional_curves[[name]]
  stopifnot(all(is.finite(s)), abs(s[1] - 1) < 1e-10,
            all(s >= -1e-10 & s <= 1 + 1e-10), all(diff(s) <= 1e-10))
}
if (any(diff(result_pl_icm$res$lambda)[-nrow(result_pl_icm$intervals$II)] > 10))
  cat("The NPMLE has a large early drop; raw and conditional survival plots are both saved.\n")

# Report distinct stopping rules. In ltic, res$conv is a score diagnostic,
# NOT a logical convergence flag. With finite output, stopping before the
# native loop cap indicates its likelihood-tolerance stopping condition.
# PL-ICM reports 11 steps for each native outer loop (10 PL + 1 ICM).
comparison <- data.frame(
  method = c("SurvPT uniform", "SurvPT exponential", "NPMLE PL-ICM", "NPMLE PL"),
  log_likelihood = c(result_unif$LogLik, result_exp$LogLik,
                     result_pl_icm$res$llike, result_pl$res$llike),
  iterations = c(result_unif$EM_Iterations, result_exp$EM_Iterations,
                 result_pl_icm$res$it, result_pl$res$it),
  converged = c(result_unif$Converged, result_exp$Converged,
                result_pl_icm$res$it / 11 < npmle_max_iter,
                result_pl$res$it < npmle_max_iter),
  stopping_rule = c(rep(sprintf("terminal-mass L2 <= %g", pt_tol), 2),
                    rep(sprintf("absolute logLik change < %g", npmle_tol), 2)),
  seconds = c(run_unif$seconds, run_exp$seconds, run_pl_icm$seconds, run_pl$seconds),
  grid_RMSE = vapply(curves[c("PT_uniform", "PT_exponential", "PL_ICM", "PL")],
    function(s) sqrt(mean((s - curves$truth)^2)), numeric(1)),
  conditional_grid_RMSE = vapply(
    conditional_curves[c("PT_uniform", "PT_exponential", "PL_ICM", "PL")],
    function(s) sqrt(mean((s - conditional_curves$truth)^2)), numeric(1)), row.names = NULL)
print(comparison, row.names = FALSE, digits = 6)
if (!all(comparison$converged)) warning("Some fits reached their iteration limit; see comparison.")
cat("RMSE is descriptive for this one dataset; timings include preprocessing, with different stopping rules.\n")

plot_survival_comparison <- function(conditional = FALSE) {
  shown <- if (conditional) conditional_curves else curves
  title <- if (conditional) sprintf("Conditional survival given T > %.1f", condition_time) else
    "Left-truncated, interval-censored survival"
  plot(shown$time, shown$truth, type = "l", col = "#C43C39", lwd = 2.5,
       xlim = range(shown$time), ylim = c(0, 1), xlab = "Event time",
       ylab = if (conditional) "Conditional survival probability" else "Survival probability",
       main = title,
       sub = sprintf("seed = %d | n = %d | truncation = %.1f%% | entry ~ U(0, %.1f)",
                     seed, n, 100 * truncated_fraction, trunc_max))
  abline(h = seq(0, 1, 0.2), col = "grey92", lty = 1)
  # Plot actual NPMLE jumps, not straight lines between grid predictions.
  origin_icm <- if (conditional) findInterval(condition_time, result_pl_icm$intervals$II$right) else NULL
  origin_pl <- if (conditional) findInterval(condition_time, result_pl$intervals$II$right) else NULL
  plot(result_pl_icm, cond = origin_icm, plot_type = "over", ends = "r", col = "#222222", lty = 1)
  plot(result_pl, cond = origin_pl, plot_type = "over", ends = "r", col = "#8B5FBF", lty = 3)
  lines(shown$time, shown$PT_uniform, col = "#2474B5", lwd = 2)
  lines(shown$time, shown$PT_exponential, col = "#E08B22", lwd = 2, lty = 2)
  lines(shown$time, shown$truth, col = "#C43C39", lwd = 2.5)
  legend("topright", bty = "n", lwd = c(2.5, 2, 2, 2, 2), lty = c(1, 1, 2, 1, 3),
         legend = c("True survival: Exp(1)", "SurvPT: uniform base",
                    "SurvPT: exponential base", "NPMLE: PL-ICM", "NPMLE: Product-Limit"),
         col = c("#C43C39", "#2474B5", "#E08B22", "#222222", "#8B5FBF"))
}
if (interactive()) {
  plot_survival_comparison()
  plot_survival_comparison(conditional = TRUE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
png_file <- file.path(output_dir, "survival_comparison.png")
save_plot <- function(path, conditional = FALSE) {
  png(path, width = 1500, height = 1050, res = 160)
  on.exit(dev.off())
  plot_survival_comparison(conditional)
}
save_plot(png_file)
save_plot(file.path(output_dir, "survival_conditional.png"), conditional = TRUE)
write.csv(test_data, file.path(output_dir, "observed_data.csv"), row.names = FALSE)
write.csv(curves, file.path(output_dir, "survival_curves.csv"), row.names = FALSE)
write.csv(conditional_curves, file.path(output_dir, "survival_conditional_curves.csv"), row.names = FALSE)
write.csv(comparison, file.path(output_dir, "fit_summary.csv"), row.names = FALSE)
saveRDS(list(settings = list(seed = seed, n_target = n_target, trunc_max = trunc_max,
  inspection_max = inspection_max, inspection_gap = inspection_gap, true_rate = true_rate,
  max_M = max_M, pt_max_iter = pt_max_iter, pt_tol = pt_tol,
  npmle_max_iter = npmle_max_iter, npmle_tol = npmle_tol, open_T = TRUE,
  condition_time = condition_time),
  n_source = generated$n_source, test_data = test_data, exponential_rate = exp_par,
  fits = list(PT_uniform = result_unif, PT_exponential = result_exp,
              PL_ICM = result_pl_icm, PL = result_pl),
  curves = curves, conditional_curves = conditional_curves, comparison = comparison),
  file.path(output_dir, "single_dataset.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("Saved comparison plot to", normalizePath(png_file, winslash = "/"), "\n")
