# Posterior survival intervals: SurvPT versus spBayesSurv MPT and TBP.
# Run in R after installing SurvPT >= 0.2.3:
#   source("Examples/posterior.R")
# Dependencies: SurvPT, survival, spBayesSurv 1.1.9, ggplot2.
# Results remain in memory. Tables are printed; no result or plot files are saved.
# In an interactive R session, the ggplot is printed to the graphics device.

required_packages <- c("SurvPT", "survival", "spBayesSurv", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) {
  stop("Install these packages before running the example: ",
       paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("SurvPT") < "0.2.3") {
  stop("Install SurvPT 0.2.3 or later, restart R, then source this example.")
}
if (as.character(utils::packageVersion("spBayesSurv")) != "1.1.9") {
  stop("The baseline initialization adapter below is validated for spBayesSurv 1.1.9.")
}
library(SurvPT)

# Same simulated data and PT settings as test_function.R.
K <- 1
Lambda <- 0.12
n <- 200L
data_seed <- 124
level <- 0.95
# Sampling controls below apply only to PT; spBayesSurv uses its own defaults.
n_samples <- 2000L
burn_in <- 3000L
thin <- 1L
max_M <- 8L
max_EM_iter <- 2000L
flat <- FALSE

test_data <- data_rc(seed = data_seed, n = n, K = K, Lambda = Lambda, q = 0.3)
exp_par <- estimate_exp(test_data$L, test_data$R)
endpoints <- as.matrix(test_data[, c("L", "R")])
x <- seq(0, max(endpoints[is.finite(endpoints)]), length.out = 201L)
true_survival <- stats::pweibull(x, shape = K, scale = 1/Lambda, lower.tail = FALSE)

# These are the same initialization adapters used in test_function.R.
# MPT 1.1.9 rejects ~1, so a zero design column represents a baseline model.
# Only the parametric initializer is adapted; the package's MCMC is unchanged.
bayes_survreg <- function(formula, dist, ...) {
  fit <- survival::survreg(formula, dist = dist, ...)
  # spBayesSurv inverts this covariance for its parametric MCMC initializer.
  # survreg can exhaust its iterations and return a singular covariance.
  invertible <- function(v) is.matrix(v) && all(is.finite(v)) &&
    !inherits(tryCatch(solve(v), error = identity), 'error')
  if (invertible(fit$var)) return(fit)
  # Only repair the intercept-only loglogistic initialization used here.
  # No random draws, jitter, changed observations, or replacement seeds.
  if (dist != 'loglogistic' || length(fit$coefficients) != 1L ||
      !all(is.finite(fit$coefficients))) stop('Unsupported singular survreg initializer')
  retry <- survival::survreg(formula, dist = dist, ...,
                             init = c(unname(fit$coefficients), 0),
                             control = survival::survreg.control(maxiter = 100L))
  if (!invertible(retry$var) || !is.finite(retry$scale) || retry$scale <= 0 ||
      !all(is.finite(retry$coefficients)) || retry$iter >= 100L ||
      !is.finite(tail(retry$loglik, 1L)) ||
      tail(retry$loglik, 1L) < tail(fit$loglik, 1L))
    stop('Deterministic survreg initializer retry failed')
  warning('Recovered singular loglogistic survreg initializer with deterministic scale=1 restart',
          call. = FALSE)
  retry
}

adapt_bayes_initializers <- function(f, initializer) {
  stopifnot(as.character(utils::packageVersion('spBayesSurv')) == '1.1.9')
  env <- new.env(parent = environment(f))
  env$baseline_survreg <- initializer
  replacements <- 0L
  adapt <- function(x) {
    if (!is.call(x)) return(x)
    if (identical(x[[1L]], quote(survival::survreg))) {
      x[[1L]] <- as.name('baseline_survreg')
      replacements <<- replacements + 1L
    }
    for (i in seq_along(x)[-1L]) if (is.call(x[[i]])) x[[i]] <- adapt(x[[i]])
    x
  }
  body(f) <- adapt(body(f))
  stopifnot(replacements == 7L)
  environment(f) <- env
  f
}

mpt_baseline_function <- function() {
  stopifnot(as.character(utils::packageVersion('spBayesSurv')) == '1.1.9')
  baseline_survreg <- function(formula, dist, ...) {
    frame <- environment(formula)
    X <- get('Xbase.scaled', envir = frame)
    stopifnot(ncol(X) == 1L, all(X == 0))
    # Fit the identifiable baseline, then pad the unused nuisance coefficient.
    baseline <- formula
    baseline[[3L]] <- 1
    fit <- bayes_survreg(baseline, dist = dist, ...)
    stopifnot(length(fit$coefficients) == 1L, all(dim(fit$var) == 2L))
    covariance <- matrix(0, 3L, 3L)
    covariance[c(1L, 3L), c(1L, 3L)] <- fit$var
    covariance[2L, 2L] <- 1
    fit$var <- covariance
    fit$coefficients <- c(fit$coefficients, empt_zero = 0)
    fit
  }
  adapt_bayes_initializers(spBayesSurv::survregbayes2, baseline_survreg)
}

# Omit mcmc so both spBayesSurv models use their package defaults.
# The package-default parametric initialization chain is also retained.
fit_mpt_posterior <- function(data) {
  data$empt_zero <- 0
  baseline_fit <- mpt_baseline_function()
  fit <- baseline_fit(survival::Surv(L, R, type = "interval2") ~ empt_zero,
                       data = data, scale.designX = FALSE)
  stopifnot(all(fit$X.scaled == 0), all(is.finite(fit$beta.scaled)),
            all(fit$X.scaled %*% fit$beta.scaled == 0),
            identical(fit$mcmc, eval(formals(spBayesSurv::survregbayes2)$mcmc)),
            fit$maxL == 4)
  fit
}
fit_tbp_posterior <- function(data) {
  baseline_fit <- adapt_bayes_initializers(spBayesSurv::survregbayes, bayes_survreg)
  fit <- baseline_fit(survival::Surv(L, R, type = "interval2") ~ 1,
                       data = data)
  stopifnot(identical(fit$mcmc, eval(formals(spBayesSurv::survregbayes)$mcmc)),
            fit$maxL == 15)
  fit
}

pt_specs <- list(
  "PT uniform" = list(prior = 1, prior_params = NULL),
  "PT exponential" = list(prior = 2, prior_params = exp_par),
  "PT uniform-exponential" = list(prior = 4, prior_params = exp_par)
)
pt_models <- pt_predictions <- setNames(vector("list", length(pt_specs)), names(pt_specs))
run_times <- numeric()
for (i in seq_along(pt_specs)) {
  method <- names(pt_specs)[i]
  cat("\nFitting and sampling", method, "...\n")
  started <- proc.time()[["elapsed"]]
  pt_models[[method]] <- do.call(survPT, c(list(
    X = test_data, max_M = max_M, max_EM_iter = max_EM_iter, EM_tol = 1e-6,
    flat = flat, verbose = FALSE), pt_specs[[method]]))
  set.seed(300000000L + i)
  pt_predictions[[method]] <- predict(
    pt_models[[method]], x, sample = TRUE, level = level,
    n_samples = n_samples, burn_in = burn_in, thin = thin)
  run_times[method] <- proc.time()[["elapsed"]] - started
}

# Encode censoring for Surv(interval2), leaving test_data unchanged.
bayes_data <- data.frame(L = test_data$L, R = test_data$R)
bayes_data$L[bayes_data$L == 0 & bayes_data$R > 0] <- NA_real_
bayes_data$R[is.infinite(bayes_data$R)] <- NA_real_
positive <- x > 0

cat("\nFitting spBayesSurv MPT ...\n")
started <- proc.time()[["elapsed"]]
set.seed(100000001L)
mpt_model <- fit_mpt_posterior(bayes_data)
mpt_curves <- spBayesSurv::plot.survregbayes2(
  mpt_model, xpred = matrix(0, 1L, 1L), tgrid = x[positive], CI = level, PLOT = FALSE)
run_times["spBayesSurv MPT"] <- proc.time()[["elapsed"]] - started

cat("\nFitting spBayesSurv TBP ...\n")
started <- proc.time()[["elapsed"]]
set.seed(200000001L)
tbp_model <- fit_tbp_posterior(bayes_data)
tbp_curves <- spBayesSurv::plot.survregbayes(
  tbp_model, tgrid = x[positive], CI = level, PLOT = FALSE)
run_times["spBayesSurv TBP"] <- proc.time()[["elapsed"]] - started

# Keep the PT MAP estimate distinct from its posterior mean and median.
# Shat is the survival summary returned by spBayesSurv; Shatlow/Shatup are
# its credible limits. No PT MAP or posterior median is imputed for that package.
pt_tables <- lapply(names(pt_predictions), function(method) {
  result <- pt_predictions[[method]]
  data.frame(method = method, x = x, truth = true_survival,
             center = result$mean, lower = result$lower, upper = result$upper,
             PT_MAP = result$estimate, PT_median = result$median)
})
bayes_table <- function(curves, method) {
  stopifnot(all(c("Shat", "Shatlow", "Shatup") %in% names(curves)))
  result <- data.frame(method = method, x = x, truth = true_survival,
                       center = 1, lower = 1, upper = 1,
                       PT_MAP = NA_real_, PT_median = NA_real_)
  result$center[positive] <- as.numeric(curves$Shat[, 1L])
  result$lower[positive] <- as.numeric(curves$Shatlow[, 1L])
  result$upper[positive] <- as.numeric(curves$Shatup[, 1L])
  result
}
posterior_intervals <- do.call(rbind, c(pt_tables, list(
  bayes_table(mpt_curves, "spBayesSurv MPT"),
  bayes_table(tbp_curves, "spBayesSurv TBP"))))
method_order <- c(names(pt_specs), "spBayesSurv MPT", "spBayesSurv TBP")
posterior_intervals$method <- factor(posterior_intervals$method, levels = method_order)
posterior_intervals$width <- posterior_intervals$upper - posterior_intervals$lower
stopifnot(all(is.finite(as.matrix(posterior_intervals[, c("center", "lower", "upper")]))),
          all(posterior_intervals$lower >= -1e-10),
          all(posterior_intervals$upper <= 1 + 1e-10),
          all(posterior_intervals$width >= -1e-10))

cat(sprintf("\nPointwise %.1f%% posterior intervals.\n", 100 * level))
sampling_settings <- data.frame(
  method = method_order,
  retained = c(rep(n_samples, length(pt_specs)), mpt_model$mcmc$nsave, tbp_model$mcmc$nsave),
  burn_in = c(rep(burn_in, length(pt_specs)), mpt_model$mcmc$nburn, tbp_model$mcmc$nburn),
  thin = c(rep(thin, length(pt_specs)), mpt_model$mcmc$nskip + 1, tbp_model$mcmc$nskip + 1))
cat("Main posterior chain settings (excluding spBayesSurv's parametric initialization):\n")
print(sampling_settings, row.names = FALSE)
cat("Center: PT posterior mean; spBayesSurv package Shat. PT_MAP is the fitted mode-based estimate.\n")
cat("Intervals condition on each model's settings; PT and spBayesSurv use different priors.\n")
cat("Versions: SurvPT", as.character(utils::packageVersion("SurvPT")),
    "; spBayesSurv", as.character(utils::packageVersion("spBayesSurv")), "\n")
# Print 11 common times per model; all 201 times remain in posterior_intervals.
print_times <- x[unique(round(seq(1, length(x), length.out = 11L)))]
print(posterior_intervals[posterior_intervals$x %in% print_times, ],
      row.names = FALSE, digits = 4)
interval_summary <- do.call(rbind, lapply(method_order, function(method) {
  result <- posterior_intervals[posterior_intervals$method == method, ]
  data.frame(method = method, mean_interval_width = mean(result$width),
             max_interval_width = max(result$width), elapsed_seconds = unname(run_times[method]))
}))
cat("\nInterval width summary on this grid (not repeated-sampling coverage):\n")
print(interval_summary, row.names = FALSE, digits = 4)

posterior_plot <- ggplot2::ggplot(posterior_intervals, ggplot2::aes(x = x)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = method), alpha = 0.25) +
  ggplot2::geom_line(ggplot2::aes(y = center, colour = method), linewidth = 0.8) +
  ggplot2::geom_line(ggplot2::aes(y = truth), colour = "black", linetype = "dashed", linewidth = 0.6) +
  ggplot2::facet_wrap(~ method, ncol = 2) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = sprintf("Pointwise %.0f%% posterior survival intervals", 100 * level),
    subtitle = "Colored line: posterior summary; dashed black line: simulated truth",
    x = "Time", y = "Survival probability",
    caption = "PT line: posterior mean; spBayesSurv line: Shat. PT MAP estimates are printed separately.") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "none")

if (interactive()) {
  print(posterior_plot)
} else {
  message("Plot available as posterior_plot. Source this script in interactive R to display it; no graphics file was written.")
}
