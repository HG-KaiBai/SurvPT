#' loglikelihood function for Weibull distribution
#'
#' @param theta parameter vector (rate, shape)
#' @param LB left bound of data
#' @param RB right bound of data
#'
#' @returns negative log-likelihood value
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


#' gradient of loglikelihood function for Weibull distribution
#'
#' @param theta parameter vector (rate, shape)
#' @param LB left bound of data
#' @param RB right bound of data
#'
#' @returns gradient vector
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

#' Estimate the parameters of Weibull distribution using MLE
#'
#' @param LB left bound of data
#' @param RB right bound of data
#'
#' @returns estimated parameter vector (rate, shape)
#' @export
#'
#' @examples
#' data = data_cs(seed = 123, n = 100, K = 2, Lambda = 0.1, q = 0.3)
#' estimate(data$L, data$R)
estimate <- function(LB, RB){
  return(optim(par = c(0.01, 1), fn = loglike, LB = LB, RB = RB,
               gr = gradient, method = "L-BFGS-B", lower = c(0.001,0.01))$par)
}



#' loglikelihood function for Exponential distribution
#'
#' @param theta parameter vector (rate)
#' @param LB left bound of data
#' @param RB right bound of data
#'
#' @returns negative log-likelihood value
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


#' gradient of loglikelihood function for Exponential distribution
#'
#' @param theta parameter vector (rate)
#' @param LB left bound of data
#' @param RB right bound of data
#'
#' @returns gradient vector
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


#' Estimate the parameter of Exponential distribution using MLE
#'
#' @param LB left bound of data
#' @param RB right bound of data
#'
#' @returns estimated parameter (rate)
#' @export
#'
#' @examples
#' data = data_cs(seed = 123, n = 100, K = 1, Lambda = 0.1, q = 0.3)
#' estimate_exp(data$L, data$R)
estimate_exp <- function(LB, RB){
  return(optim(par = c(0.01), fn = loglike_exp, LB = LB, RB = RB,
               gr = gradient_exp, method = "L-BFGS-B", lower = c(0.001))$par)
}

