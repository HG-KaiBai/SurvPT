# Posterior survival intervals: SurvPT versus spBayesSurv MPT and TBP.
# Run in R after installing SurvPT >= 0.4.1:
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
if (utils::packageVersion("SurvPT") < "0.4.1") {
  stop("Install SurvPT 0.4.1 or later, restart R, then source this example.")
}
if (as.character(utils::packageVersion("spBayesSurv")) != "1.1.9") {
  stop("The baseline initialization adapter below is validated for spBayesSurv 1.1.9.")
}
library(SurvPT)

# Reproducible right-censored example; settings can be edited in R.
K <- 2
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
# Optional mixture: 0.5 * Weibull(shape = 5, scale = 1) +
#                   0.5 * Weibull(shape = 5, scale = 4.5); K/Lambda do not apply.
# Use a build containing data_mixture(), or run devtools::load_all(".") before
# sourcing this script. Uncomment ONE call below and the mixture truth below.
# All three types use non-informative censoring, with q = P(C < T).
# test_data <- data_mixture(seed = data_seed, n = n, type = "right", q = 0.3)
# test_data <- data_mixture(seed = data_seed, n = n, type = "interval", q = 0.3)
# test_data <- data_mixture(seed = data_seed, n = n, type = "current", q = 0.3)
# Fit using L/R: status indicates exact events, not event-by-inspection status.
exp_par <- estimate_exp(test_data$L, test_data$R)
endpoints <- as.matrix(test_data[, c("L", "R")])
x <- seq(0, max(endpoints[is.finite(endpoints)]), length.out = 201L)
true_survival <- stats::pweibull(x, shape = K, scale = 1/Lambda, lower.tail = FALSE)
# Matching truth for data_mixture(): uncomment both continuation lines.
# true_survival <- 0.5 * stats::pweibull(x, 5, 1, lower.tail = FALSE) +
#                  0.5 * stats::pweibull(x, 5, 4.5, lower.tail = FALSE)

source("Examples/bayes_helpers.R", local = TRUE)

pt_specs <- list(
  "PT uniform" = list(prior = 1, prior_params = NULL),
  "PT exponential" = list(prior = 2, prior_params = exp_par),
  "PT Weibull" = list(prior = 3, prior_params = estimate(test_data$L, test_data$R))
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
mpt_model <- fit_mpt(bayes_data)
mpt_curves <- spBayesSurv::plot.survregbayes2(
  mpt_model, xpred = matrix(0, 1L, 1L), tgrid = x[positive], CI = level, PLOT = FALSE)
run_times["spBayesSurv MPT"] <- proc.time()[["elapsed"]] - started

cat("\nFitting spBayesSurv TBP ...\n")
started <- proc.time()[["elapsed"]]
set.seed(200000001L)
tbp_model <- fit_tbp(bayes_data)
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

# Tail probability is available even if the target grid omits tau.
tail_probabilities <- do.call(rbind, lapply(names(pt_predictions), function(method) {
  data.frame(method = method, attr(pt_predictions[[method]], "tail"))
}))
cat("\nPT tail probabilities S(tau), with pointwise posterior intervals:\n")
print(tail_probabilities, row.names = FALSE)

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
