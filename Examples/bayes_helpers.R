# Shared MPT/TBP fitting helpers for test_function.R and posterior.R.
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

