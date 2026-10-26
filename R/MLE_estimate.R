#' Weibull negative log-likelihood for censored observations
#'
#' Internal objective used by \code{\link{estimate}()}. Exact observations contribute a
#' density, finite intervals contribute interval probability, and right-
#' censored observations contribute a survival probability.
#'
#' @param theta Positive parameter vector `c(lambda, k)`, with inverse scale
#'   `lambda` and shape `k`. Weibull survival is `exp(-(lambda * t)^k)`.
#' @param LB Numeric vector of finite nonnegative left endpoints.
#' @param RB Numeric vector of right endpoints, the same length as `LB`.
#'   Use `RB == LB` for exact events and `Inf` for right censoring. Exact event
#'   times should be positive. Missing values are not supported.
#' @return A scalar negative log-likelihood. Extreme parameters can produce
#'   nonfinite results because direct exponential differences are used.
#' @details This helper is not exported. Use \code{\link{estimate}()} for fitting.
#' @seealso \code{\link{gradient}()}, \code{\link{estimate}()}
#' @keywords internal
#'
loglike <- function(theta, LB, RB){
  temp = 0
  for (i in 1:length(LB)){
    if (RB[i] == Inf) {
      temp = temp - (theta[1]*LB[i])^theta[2]
    } else if (LB[i] == RB[i]){
      temp = temp + log(theta[2]*theta[1]) + (theta[2]-1)*log(LB[i]*theta[1])-(theta[1]*LB[i])^theta[2]
    } else {
      temp = temp + log(exp(-(theta[1]*LB[i])^theta[2]) - exp(-(theta[1]*RB[i])^theta[2]))
    }
  }
  return(-temp)
}


#' Gradient of the Weibull negative log-likelihood
#'
#' Internal analytic gradient used by \code{\link{estimate}()}.
#'
#' @inheritParams loglike
#' @return A length-two gradient in the order `c(lambda, k)` for the negative
#'   log-likelihood, not the positive log-likelihood.
#' @seealso \code{\link{loglike}()}, \code{\link{estimate}()}
#' @keywords internal
#'
gradient <- function(theta, LB, RB){
  g = c(0, 0)
  for (i in 1:length(LB)){
    if (RB[i] == Inf){
      g[1] = g[1] + (-theta[2] * theta[1]^(theta[2]-1) * LB[i]^theta[2])
      if (LB[i] != 0){
        g[2] = g[2] - log(theta[1]*LB[i]) * (theta[1]*LB[i])^theta[2]
      }
    } else if (LB[i] == RB[i]){
      g[1] = g[1] + theta[2]/theta[1] - theta[2] * theta[1]^(theta[2]-1) * LB[i]^theta[2]
      if (LB[i] != 0){
        g[2] = g[2] + 1/theta[2] + log(LB[i]*theta[1]) - log(theta[1]*LB[i])*(theta[1]*LB[i])^theta[2]
      }
    } else {
      g[1] = g[1] +
        (-theta[2] * exp(-(theta[1]*LB[i])^theta[2])*theta[1]^(theta[2]-1)*LB[i]^theta[2]+
           theta[2] * exp(-(theta[1]*RB[i])^theta[2])*theta[1]^(theta[2]-1)*RB[i]^theta[2])/
        (exp(-(theta[1]*LB[i])^theta[2]) - exp(-(theta[1]*RB[i])^theta[2]))
      if (LB[i] != 0){
        g[2] = g[2] +
          (-log(theta[1]*LB[i])* (theta[1]*LB[i])^theta[2] * exp(-(theta[1]*LB[i])^theta[2])+
             log(theta[1]*RB[i])* (theta[1]*RB[i])^theta[2] * exp(-(theta[1]*RB[i])^theta[2]))/
          (exp(-(theta[1]*LB[i])^theta[2]) - exp(-(theta[1]*RB[i])^theta[2]))
      } else {
        g[2] = g[2] +
          (log(theta[1]*RB[i]) * (theta[1]*RB[i])^theta[2] * exp(-(theta[1]*RB[i])^theta[2]))/
          (1 - exp(-(theta[1]*RB[i])^theta[2]))
      }
    }
  }
  return(-g)
}

#' Fit a Weibull distribution to censored event times
#'
#' Minimize the negative log-likelihood of exact and censored observations
#' using `stats::optim()` with the `L-BFGS-B` method and an analytic gradient.
#'
#' @inheritParams loglike
#' @details The initial parameter vector is `c(0.01, 1)` and lower bounds are
#'   `c(0.001, 0.01)`. These controls are fixed by the current implementation.
#'   Weibull survival is `exp(-(lambda * t)^k)`, equivalent to R's Weibull with
#'   `shape = k, scale = 1 / lambda`. The result can be passed as `prior_params`
#'   to \code{\link{survPT}()} with `prior = 3`.
#'
#'   Only the optimizer's parameter estimates are returned: convergence codes,
#'   objective values, standard errors, and Hessians are not exposed. A returned
#'   parameter vector alone does not establish convergence or identifiability.
#' @return Numeric vector `c(lambda, k)` of estimated inverse scale and shape.
#' @seealso \code{\link{estimate_exp}()}, \code{\link{survPT}()}
#' @export
#'
#' @examples
#' d <- data_rc(seed = 123, n = 30, K = 2, Lambda = 0.4, q = 0.3)
#' estimate(d$L, d$R)
estimate <- function(LB, RB){
  return(optim(par = c(0.01, 1), fn = loglike, LB = LB, RB = RB,
               gr = gradient, method = "L-BFGS-B", lower = c(0.001,0.01))$par)
}



#' Exponential negative log-likelihood for censored observations
#'
#' Internal objective used by \code{\link{estimate_exp}()}, with exact, finite-interval,
#' and right-censored contributions.
#'
#' @param theta Positive exponential rate, supplied as a length-one vector.
#' @inheritParams loglike
#' @return Scalar negative log-likelihood. Extreme parameters can produce
#'   nonfinite results because direct exponential differences are used.
#' @details This helper is not exported. Use \code{\link{estimate_exp}()} for fitting.
#' @seealso \code{\link{gradient_exp}()}, \code{\link{estimate_exp}()}
#' @keywords internal
#'
loglike_exp <- function(theta, LB, RB){
  temp = 0
  for (i in 1:length(LB)){
    if (RB[i] == Inf) {
      temp = temp - (theta[1]*LB[i])
    } else if (LB[i] == RB[i]){
      temp = temp + log(theta[1]) -(theta[1]*LB[i])
    } else {
      temp = temp + log(exp(-(theta[1]*LB[i])) - exp(-(theta[1]*RB[i])))
    }
  }
  return(-temp)
}


#' Gradient of the exponential negative log-likelihood
#'
#' Internal analytic gradient used by \code{\link{estimate_exp}()}.
#'
#' @inheritParams loglike_exp
#' @return A scalar derivative with respect to the exponential rate for the
#'   negative log-likelihood, not the positive log-likelihood.
#' @seealso \code{\link{loglike_exp}()}, \code{\link{estimate_exp}()}
#' @keywords internal
#'
gradient_exp <- function(theta, LB, RB){
  g =0
  for (i in 1:length(LB)){
    if (RB[i] == Inf){
      g = g - LB[i]
    } else if (LB[i] == RB[i]){
      g[1] = g[1] + 1/theta[1] - LB[i]
    } else {
      g[1] = g[1] +
        (-exp(-(theta[1]*LB[i])) *LB[i] + exp(-(theta[1]*RB[i]))* RB[i])/
        (exp(-(theta[1]*LB[i])) - exp(-(theta[1]*RB[i])))
    }
  }
  return(-g)
}


#' Fit an exponential distribution to censored event times
#'
#' Minimize the negative log-likelihood with an analytic gradient using
#' `stats::optim()` and the `L-BFGS-B` method.
#'
#' @inheritParams loglike
#' @details The starting rate is `0.01` and the lower bound is `0.001`.
#'   The result can be passed to \code{\link{survPT}()} with `prior = 2`. Only the
#'   optimizer's parameter estimate is returned; convergence status, standard
#'   errors, and objective values are not exposed. A returned value alone does
#'   not establish convergence or identifiability.
#' @return A length-one numeric vector containing the estimated rate. Its
#'   reciprocal is the exponential mean.
#' @seealso \code{\link{estimate}()}, \code{\link{survPT}()}
#' @export
#'
#' @examples
#' d <- data_rc(seed = 123, n = 30, K = 1, Lambda = 0.4, q = 0.3)
#' estimate_exp(d$L, d$R)
estimate_exp <- function(LB, RB){
  return(optim(par = c(0.01), fn = loglike_exp, LB = LB, RB = RB,
               gr = gradient_exp, method = "L-BFGS-B", lower = c(0.001))$par)
}

