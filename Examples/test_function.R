# Historical interactive comparison; see Examples/README.md before running.
# Requires devtools, interval and spBayesSurv 1.1.9.
# Optional fiducial plots require constructing the currently commented-out fid.
devtools::load_all(".", compile = TRUE)
library(interval)
library(spBayesSurv)

K = 2
Lambda = 0.12
n = 500
#test_data = data_rc(seed = 123, n = n, K = K, Lambda = Lambda, q = 0.3)
test_data = data_mixture(seed = 123, n = n, type = "right", q = 0.3)
# test_data = data_mixture(seed = 123, n = n, type = "interval", q = 0.3)
# test_data = data_mixture(seed = 123, n = n, type = "current", q = 0.3)
#test_data = fiducial_3(321, n)
exp_par = estimate_exp(test_data$L, test_data$R)
# Weibull parameters are c(lambda, shape), with scale = 1 / lambda.
weibull_par = estimate(test_data$L, test_data$R)

# PT predictions and the later comparisons use the same finite observation grid.
obs.mat = c(test_data$L[is.finite(test_data$L)], test_data$R[is.finite(test_data$R)])
x = seq(0, max(obs.mat), length.out = 201)
#true.surv = pweibull(x, shape = K, scale = 1/Lambda, lower.tail = FALSE)
#true.surv =  pgamma(x, 2, 1, lower.tail = FALSE)
true.density = dgamma(x, 2, 1)
true.surv = 0.5 * pweibull(x, 5, 1, lower.tail = FALSE) + 0.5 * pweibull(x, 5, 4.5, lower.tail = FALSE)
# true.density = 0.5 * dweibull(x, 5, 1) + 0.5 * dweibull(x, 5, 4.5)
cat("Calling SurvPT with uniform prior...\n")
fit_unif <- survPT(test_data, prior = 1,
                      max_M = 8, max_EM_iter = 2000, EM_tol = 1e-6)
result_unif <- predict(fit_unif, x, sample = FALSE)

cat("Calling SurvPT with exponential prior...\n")
fit_exp <- survPT(test_data, prior = 2, prior_params = exp_par,
                     max_M = 8, max_EM_iter = 2000, EM_tol = 1e-6)
result_exp <- predict(fit_exp, x, sample = FALSE)

cat("Calling SurvPT with Weibull prior...\n")
fit_weibull <- survPT(test_data, prior = 3, prior_params = weibull_par,
                         max_M = 8, max_EM_iter = 2000, EM_tol = 1e-6)
result_weibull <- predict(fit_weibull, x, sample = FALSE)

density_unif <- predict(fit_unif, x, type = "density", sample = FALSE)
density_exp <- predict(fit_exp, x, type = "density", sample = FALSE)
density_weibull <- predict(fit_weibull, x, type = "density", sample = FALSE)
tail_probabilities <- rbind(Unif = attr(result_unif, "tail"),
                            Exp = attr(result_exp, "tail"),
                            Weibull = attr(result_weibull, "tail"))
print(tail_probabilities)

# Retain the original survival SSE and probability-increment RMSE comparisons.
pt_results <- list(Unif = result_unif, Exp = result_exp, Weibull = result_weibull)
pt_fits <- list(Unif = fit_unif, Exp = fit_exp, Weibull = fit_weibull)
pt_checks <- do.call(rbind, lapply(names(pt_results), function(method) {
  fit <- pt_fits[[method]]
  survival <- as.numeric(pt_results[[method]]$estimate)
  data.frame(
    Method = method,
    Terminal_Mass = sum(fit$Posterior_Probs[fit$Terminal_Nodes]),
    Survival_SSE = sum((survival - true.surv)^2),
    Increment_RMSE = sqrt(mean((abs(diff(c(survival, 0))) -
                                abs(diff(c(true.surv, 0))))^2))
  )
}))
print(pt_checks, row.names = FALSE)

source("./Examples/bayes_helpers.R", local = TRUE)

result_opt <- NULL
# OPT's legacy interface requires exact events or right censoring.
if (all(test_data$L == test_data$R | is.infinite(test_data$R))) tryCatch({
  cat("Calling SurvOPT with uniform prior...\n")
  result_opt <- survOPT(test_data$L, as.integer(is.finite(test_data$R)),
                        test_seq = x, max_M = 8)

  cat("Success! Result structure:\n")
  print(str(result_opt))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

if (!is.null(result_opt)) {
for (l in 0:9){
  print(sum(result_opt$Posterior_Probs[result_opt$Level == l]))
}
print(sum(result_opt$Posterior_Probs[result_opt$Terminal_Nodes]))
}

grid = result_unif$x
#fid = fiducial(test_data$L, test_data$R, n = n, grid, grid)

# Encode censoring for Surv(interval2), leaving test_data unchanged.
bayes_data <- data.frame(L = test_data$L, R = test_data$R)
bayes_data$L[!is.na(bayes_data$L) & bayes_data$L == 0 &
               !is.na(bayes_data$R) & bayes_data$R > 0] <- NA_real_
bayes_data$R[is.infinite(bayes_data$R)] <- NA_real_
stopifnot(length(grid) > 1L, all(is.finite(grid)), all(grid >= 0),
          all(diff(grid) > 0))
positive_grid <- grid[grid > 0]

# Package-default MCMC/priors; only the agreed MPT baseline adapter and
# singular-initializer recovery are applied. Fixed seeds make this repeatable.
cat("Fitting MPT with package defaults...\n")
set.seed(100000001L)
fit_mpt_example <- fit_mpt(bayes_data)
mpt_curves <- plot(fit_mpt_example, xpred = matrix(0, 1L, 1L),
                   tgrid = positive_grid, PLOT = FALSE)

cat("Fitting TBP with package defaults...\n")
set.seed(200000001L)
fit_tbp_example <- fit_tbp(bayes_data)
tbp_curves <- plot(fit_tbp_example, tgrid = positive_grid, PLOT = FALSE)

# S(0)=1; package densities are evaluated only at positive times.
result_mpt <- data.frame(Test_seq = grid, Survival_Probs = 1, Density = NA_real_)
result_tbp <- result_mpt
result_mpt$Survival_Probs[grid > 0] <- as.numeric(mpt_curves$Shat[, 1L])
result_mpt$Density[grid > 0] <- as.numeric(mpt_curves$fhat[, 1L])
result_tbp$Survival_Probs[grid > 0] <- as.numeric(tbp_curves$Shat[, 1L])
result_tbp$Density[grid > 0] <- as.numeric(tbp_curves$fhat[, 1L])

kmf = icfit(test_data$L, test_data$R)
plot_comparison <- function() {
plot(kmf)
lines(result_unif$x, result_unif$estimate, col = "blue", lwd = 2)
lines(result_exp$x, result_exp$estimate, col = "orange", lwd = 2)
lines(result_weibull$x, result_weibull$estimate, col = "purple", lwd = 2)
#lines(result_opt$Test_seq, result_opt$Survival_Probs, col = "green", lwd = 2)
lines(result_unif$x, true.surv, col = "red", lwd = 2)
#lines(grid, 1 - fid$point_li, col = "purple", lwd = 2)
lines(result_mpt$Test_seq, result_mpt$Survival_Probs, col = "brown", lwd = 2)
lines(result_tbp$Test_seq, result_tbp$Survival_Probs, col = "cyan4", lwd = 2)
legend("topright", legend = c("NPMLE", "SurvPT Unif", "SurvPT Exp", "SurvPT Weibull",
                             "True Survival", "MPT", "TBP"),
       col = c("black", "blue", "orange", "purple", "red", "brown", "cyan4"), lwd = 2)


plot(result_unif$x, true.density, type = "l", col = "red", lwd = 2)
lines(result_unif$x, density_unif$estimate, col = "blue", lwd = 2)
lines(result_exp$x, density_exp$estimate, col = "orange", lwd = 2)
lines(result_weibull$x, density_weibull$estimate, col = "purple", lwd = 2)
if (!is.null(result_opt)) lines(result_opt$Test_seq, result_opt$Density, col = "green", lwd = 2)
#lines(grid, fid$density, col = "purple", lwd = 2)
lines(result_mpt$Test_seq, result_mpt$Density, col = "brown", lwd = 2)
lines(result_tbp$Test_seq, result_tbp$Density, col = "cyan4", lwd = 2)
legend("topright", legend = c("True Density", "SurvPT Unif", "SurvPT Exp", "SurvPT Weibull",
                             "MPT", "TBP",
                             if (!is.null(result_opt)) "OPT"),
       col = c("red", "blue", "orange", "purple", "brown", "cyan4",
               if (!is.null(result_opt)) "green"), lwd = 2)
}
if (interactive()) plot_comparison()







