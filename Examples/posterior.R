# Compare PT, fiducial, MPT and TBP on the current test_function.R example.
# Run from the repository root: Rscript Examples/posterior.R
# Writes comparison tables and a survival plot under output/posterior-comparison.
devtools::load_all(".", quiet = TRUE)
source("./Examples/bayes_helpers.R", local = TRUE)
source("./Examples/LinInterpolation.R", local = TRUE)
source("./Examples/Fiducial.R", local = TRUE)

# Same data generator, parameters and seed as test_function.R.
K <- 1
Lambda <- 0.12
n <- 200
test_data <- data_rc(seed = 12, n = n, K = K, Lambda = Lambda, q = 0.3)
exp_par <- estimate_exp(test_data$L, test_data$R)
obs.mat <- c(test_data$L[is.finite(test_data$L)], test_data$R[is.finite(test_data$R)])
x <- seq(0, max(obs.mat), length.out = 201)
true.surv <- pweibull(x, shape = K, scale = 1 / Lambda, lower.tail = FALSE)
level <- 0.95  # The local fiducial helper returns fixed 95% pointwise intervals.
output_dir <- file.path("output", "posterior-comparison")

# Retain the deterministic MAP-EM curve alongside the posterior mean and interval.
pt_fits <- pt_predictions <- list()
elapsed <- numeric()
for (family in c("Uniform", "Exponential")) {
  message("Fitting and sampling PT: ", family)
  started <- proc.time()[["elapsed"]]
  prior <- if (family == "Uniform") 1 else 2
  fit <- survPT(test_data, prior = prior,
                prior_params = if (prior == 2) exp_par else NULL,
                max_M = 8, max_EM_iter = 2000, EM_tol = 1e-6, verbose = FALSE)
  set.seed(if (prior == 1) 2026L else 2027L)
  pred <- predict(fit, x, level, n_sweeps = 2000, burn_in = 500, keep_draws = TRUE)
  pt_fits[[family]] <- fit
  pt_predictions[[family]] <- pred
  elapsed[paste("PT", family)] <- proc.time()[["elapsed"]] - started
}

# The fiducial helper returns a CDF median and CDF limits. Transform and reverse
# the limits to obtain a survival interval: [1 - upper_F, 1 - lower_F].
message("Sampling fiducial curves...")
started <- proc.time()[["elapsed"]]
set.seed(300000001L)
fid_grid <- seq(0, max(x), length.out = 100)
fid <- fiducial(test_data$L, test_data$R, nsample = nrow(test_data),
                grid = fid_grid, testgrid = x, mfid = 1000, mburn = 100)
elapsed["Fiducial"] <- proc.time()[["elapsed"]] - started

# Match the censoring encoding and MPT/TBP defaults in test_function.R.
bayes_data <- data.frame(L = test_data$L, R = test_data$R)
bayes_data$L[!is.na(bayes_data$L) & bayes_data$L == 0 &
               !is.na(bayes_data$R) & bayes_data$R > 0] <- NA_real_
bayes_data$R[is.infinite(bayes_data$R)] <- NA_real_
positive <- x > 0

message("Fitting MPT with spBayesSurv defaults...")
started <- proc.time()[["elapsed"]]
set.seed(100000001L)
fit_mpt_example <- fit_mpt(bayes_data)
mpt_curves <- plot(fit_mpt_example, xpred = matrix(0, 1L, 1L),
                   tgrid = x[positive], CI = level, PLOT = FALSE)
elapsed["MPT"] <- proc.time()[["elapsed"]] - started

message("Fitting TBP with spBayesSurv defaults...")
started <- proc.time()[["elapsed"]]
set.seed(200000001L)
fit_tbp_example <- fit_tbp(bayes_data)
tbp_curves <- plot(fit_tbp_example, tgrid = x[positive], CI = level, PLOT = FALSE)
elapsed["TBP"] <- proc.time()[["elapsed"]] - started

# One row per target time. Fiducial intervals and Bayesian credible intervals
# retain their respective interpretations; none of these are simultaneous bands.
curves <- list()
for (family in names(pt_fits)) {
  pred <- pt_predictions[[family]]
  post <- pred$Posterior$Survival
  curves[[paste("PT", family)]] <- data.frame(
    time = x, estimate = post$mean, lower = post$lower, upper = post$upper,
    map = as.numeric(pred$Survival_Probs), estimator = "Posterior mean",
    interval = "Pointwise credible interval")
}
curves$Fiducial <- data.frame(
  time = x, estimate = 1 - fid$point_li, lower = 1 - fid$upper,
  upper = 1 - fid$lower, map = NA_real_, estimator = "Fiducial median",
  interval = "Pointwise fiducial interval")
for (method in c("MPT", "TBP")) {
  native <- if (method == "MPT") mpt_curves else tbp_curves
  curve <- data.frame(time = x, estimate = 1, lower = 1, upper = 1,
                      map = NA_real_, estimator = "Package survival estimate",
                      interval = "Pointwise credible interval")
  curve$estimate[positive] <- as.numeric(native$Shat[, 1L])
  curve$lower[positive] <- as.numeric(native$Shatlow[, 1L])
  curve$upper[positive] <- as.numeric(native$Shatup[, 1L])
  curves[[method]] <- curve
}

# All methods are evaluated on [0,tau]; PT does not specify survival beyond tau.
for (curve in curves) {
  stopifnot(identical(curve$time, x), all(is.finite(as.matrix(curve[, 2:4]))),
            all(curve$lower <= curve$upper + 1e-10),
            all(as.matrix(curve[, 2:4]) >= -1e-8),
            all(as.matrix(curve[, 2:4]) <= 1 + 1e-8),
            all(diff(curve$estimate) <= 1e-8))
}
for (family in names(pt_fits)) {
  fit <- pt_fits[[family]]
  pred <- pt_predictions[[family]]
  stopifnot(abs(sum(fit$Posterior_Probs[fit$Terminal_Nodes]) - 1) < 1e-10,
            max(abs(rowSums(pred$Posterior$Draws$Terminal_Probs) - 1)) < 1e-10,
            abs(tail(pred$Survival_Probs, 1) - fit$Tail_Prob) < 1e-10)
}

trapz_grid <- function(y) sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
metric_row <- function(method, estimator, estimate, lower = NULL, upper = NULL) {
  error <- estimate - true.surv
  data.frame(Method = method, Estimator = estimator,
             Survival_SSE = sum(error^2), Survival_RMSE = sqrt(mean(error^2)),
             Survival_ISE = trapz_grid(error^2),
             Average_Interval_Width = if (is.null(lower)) NA_real_ else
               trapz_grid(upper - lower) / diff(range(x)))
}
comparison <- do.call(rbind, lapply(names(curves), function(method) {
  curve <- curves[[method]]
  metric_row(method, curve$estimator[1], curve$estimate, curve$lower, curve$upper)
}))
for (family in names(pt_fits)) {
  method <- paste("PT", family)
  comparison <- rbind(comparison, metric_row(method, "MAP-EM", curves[[method]]$map))
}
rownames(comparison) <- NULL
timings <- data.frame(Method = names(elapsed), Seconds = as.numeric(elapsed))
convergence <- do.call(rbind, lapply(names(pt_fits), function(family) {
  do.call(rbind, lapply(names(pt_fits[[family]]$Convergence), function(stage) {
    s <- pt_fits[[family]]$Convergence[[stage]]
    data.frame(Method = paste("PT", family), Stage = stage, Iterations = s$iterations,
               Converged = s$converged, Difference = s$difference)
  }))
}))
pt_diagnostics <- do.call(rbind, lapply(names(pt_predictions), function(family) {
  data.frame(Method = paste("PT", family), pt_predictions[[family]]$Posterior$Survival)
}))
curve_table <- do.call(rbind, lapply(names(curves), function(method) {
  data.frame(Method = method, curves[[method]], truth = true.surv)
}))

print(comparison, row.names = FALSE)
print(convergence, row.names = FALSE)
print(timings, row.names = FALSE)
message("Errors describe this one simulated dataset; interval width is not coverage.")

plot_posterior_comparison <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  colors <- c("PT Uniform" = "#0072B2", "PT Exponential" = "#D55E00",
              "Fiducial" = "#009E73", "MPT" = "#CC79A7", "TBP" = "#7B5F25")
  plot(x, true.surv, type = "l", lwd = 2, lty = 3, ylim = c(0, 1),
       xlab = "Time", ylab = "Survival probability", main = "Point estimates")
  for (method in names(curves)) lines(x, curves[[method]]$estimate, col = colors[method], lwd = 2)
  legend("topright", c(names(curves), "Truth"), col = c(colors, "black"),
         lty = c(rep(1, length(curves)), 3), lwd = 2, bty = "n", cex = 0.85)
  for (method in names(curves)) {
    curve <- curves[[method]]
    plot(x, curve$estimate, type = "n", ylim = c(0, 1), xlab = "Time",
         ylab = "Survival probability", main = method)
    polygon(c(x, rev(x)), c(curve$lower, rev(curve$upper)), border = NA,
            col = adjustcolor(colors[method], alpha.f = 0.22))
    lines(x, curve$estimate, col = colors[method], lwd = 2)
    lines(x, true.surv, col = "black", lty = 3, lwd = 2)
    labels <- c(if (method == "Fiducial") "Fiducial median" else curve$estimator[1], "Truth")
    line_colors <- c(colors[method], "black"); line_types <- c(1, 3)
    if (any(is.finite(curve$map))) {
      lines(x, curve$map, col = colors[method], lty = 2, lwd = 2)
      labels <- c(labels, "MAP-EM")
      line_colors <- c(line_colors, colors[method]); line_types <- c(line_types, 2)
    }
    legend("topright", labels, col = line_colors, lty = line_types,
           lwd = 2, bty = "n", cex = 0.85)
  }
  mtext("Right-censored example: n = 200, seed = 12; shaded pointwise 95% intervals",
        outer = TRUE, cex = 1.1)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, file.path(output_dir, "comparison.csv"), row.names = FALSE)
write.csv(curve_table, file.path(output_dir, "survival-curves.csv"), row.names = FALSE)
write.csv(convergence, file.path(output_dir, "pt-convergence.csv"), row.names = FALSE)
write.csv(pt_diagnostics, file.path(output_dir, "pt-posterior-diagnostics.csv"), row.names = FALSE)
write.csv(timings, file.path(output_dir, "timings.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
png(file.path(output_dir, "survival-comparison.png"), width = 1800, height = 1100, res = 150)
tryCatch(plot_posterior_comparison(), finally = dev.off())
if (interactive()) plot_posterior_comparison()
message("Comparison saved to ", output_dir)
