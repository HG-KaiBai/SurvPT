# Historical interactive comparison; see Examples/README.md before running.
# Requires devtools and interval. View() calls need an interactive session;
# later fiducial plots also require constructing the currently commented-out fid.
devtools::load_all()
library(interval)
library(spBayesSurv)

test_data = data_rc(123, 200, 1, 0.12, q = 0.3)
# All MPT/TBP helpers are defined here; only installed R packages are used.
# Baseline mixed Polya tree using spBayesSurv 1.1.9.
# Its public wrapper rejects ~1. A zero design column contributes exactly zero
# to every linear predictor. The parametric INITIALIZER is adapted for this
# baseline design and singular-fit recovery; the native MCMC is unchanged.
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

fit_mpt <- function(data) {
  data$empt_zero <- 0
  f <- mpt_baseline_function()
  fit <- f(survival::Surv(L, R, type = 'interval2') ~ empt_zero,
           data = data, scale.designX = FALSE)
  stopifnot(all(fit$X.scaled == 0), all(is.finite(fit$beta.scaled)),
            all(fit$X.scaled %*% fit$beta.scaled == 0),
            identical(fit$mcmc, eval(formals(spBayesSurv::survregbayes2)$mcmc)),
            identical(unname(fit$prior$S0), matrix(1e10)),
            fit$maxL == 4, fit$prior$a0 == 1, fit$prior$b0 == 1)
  fit
}

fit_tbp <- function(data) {
  f <- adapt_bayes_initializers(spBayesSurv::survregbayes, bayes_survreg)
  fit <- f(survival::Surv(L, R, type = 'interval2') ~ 1, data = data)
  stopifnot(identical(fit$mcmc, eval(formals(spBayesSurv::survregbayes)$mcmc)),
            fit$maxL == 15, fit$prior$a0 == 1, fit$prior$b0 == 1)
  fit
}

tryCatch({
  cat("Calling SurvOPT with uniform prior...\n")
  result_opt <- survOPT(test_data$L, test_data$status, test_seq = x, max_M = 8)

  cat("Success! Result structure:\n")
  print(str(result_opt))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

for (l in 0:9){
  print(sum(result_opt$Posterior_Probs[result_opt$Level == l]))
}
print(sum(result_opt$Posterior_Probs[result_opt$Terminal_Nodes]))

grid = result_unif$Test_seq
fid = fiducial(test_data$L, test_data$R, n = n, grid, grid)

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
plot(kmf)
lines(result_unif$Test_seq, result_unif$Survival_Probs, col = "blue", lwd = 2)
lines(result_exp$Test_seq, result_exp$Survival_Probs, col = "orange", lwd = 2)
#lines(result_opt$Test_seq, result_opt$Survival_Probs, col = "green", lwd = 2)
lines(result_unif$Test_seq, 1-pweibull(result_unif$Test_seq, K, 1/Lambda), col = "red", lwd = 2)
#lines(grid, 1 - fid$point_li, col = "purple", lwd = 2)
lines(result_mpt$Test_seq, result_mpt$Survival_Probs, col = "brown", lwd = 2)
lines(result_tbp$Test_seq, result_tbp$Survival_Probs, col = "cyan4", lwd = 2)
legend("topright", legend = c("Kaplan-Meier", "SurvPT Unif", "SurvPT Exp", "OPT", "True Survival", "Fiducial", "MPT", "TBP"),
       col = c("black", "blue","orange", "green" ,"red", "purple", "brown", "cyan4"), lwd = 2)


plot(result_unif$Test_seq, dweibull(result_unif$Test_seq, K, 1/Lambda), type = "l", col = "red", lwd = 2)
lines(result_unif$Test_seq, result_unif$Density, col = "blue", lwd = 2)
lines(result_exp$Test_seq, result_exp$Density, col = "orange", lwd = 2)
lines(result_opt$Test_seq, result_opt$Density, col = "green", lwd = 2)
#lines(grid, fid$density, col = "purple", lwd = 2)
lines(result_mpt$Test_seq, result_mpt$Density, col = "brown", lwd = 2)
lines(result_tbp$Test_seq, result_tbp$Density, col = "cyan4", lwd = 2)
legend("topright", legend = c("True Density", "SurvPT Unif", "SurvPT Exp", "OPT", "Fiducial", "MPT", "TBP"),
       col = c("red", "blue", "orange", "green", "purple", "brown", "cyan4"), lwd = 2)







