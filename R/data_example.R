#' Create the Right Censored Test Data from Weibull Distribution
#'
#' @param seed Random seed
#' @param n Sample size
#' @param K Weibull shape parameter
#' @param Lambda Weibull rate parameter
#' @param q Right Censoring rate
#'
#' @returns A data frame with three columns: L (observed time), R (right censoring time), and status (1 if event observed, 0 if censored)
#' @export
#'
#' @examples
#' data_rc(seed = 123, n = 100, K = 2, Lambda = 0.1, q = 0.3)
#' data_ic(seed = 123, n = 100, K = 2, Lambda = 0.1, q = 0.3)
data_rc <- function(seed, n, K, Lambda, q = 0.3){
  death = c()
  censor = c()
  Lambda_c = ((Lambda^K) / (1-q) - Lambda^K)^(1/K)
  set.seed(seed)
  T = rweibull(n, K, 1/Lambda)
  C = rweibull(n, K, 1/Lambda_c)
  delta = as.numeric(T<C)
  observed = pmin(T, C)
  R = rep(Inf, n)
  R[delta == 1] = observed[delta==1]
  data = data.frame(L = observed, R = R, status = delta)
  return(data)
}

#' Create the Interval Censored Test Data from Weibull Distribution
#'
#' @param seed Random seed
#' @param n Sample size
#' @param K Weibull shape parameter
#' @param Lambda Weibull rate parameter
#' @param q Right Censoring rate
#'
#' @returns A data frame with three columns: L (left interval), R (right interval), and status (1 if event observed, 0 if censored)
#' @export
#'
#' @examples
#' data_ic(seed = 123, n = 100, K = 1, Lambda = 0.1, q = 0.3)
#' data_cs(seed = 123, n = 100, K = 2, Lambda = 0.1, q = 0.3)
data_ic <- function(seed, n, K, Lambda, q = 0.3){
  lb = c()
  rb = c()
  Lambda_c = ((Lambda^K) / (1-q) - Lambda^K)^(1/K)
  set.seed(seed)
  Time =  rweibull(n, K, 1/Lambda)
  Censor =  rweibull(n, K, 1/Lambda_c)
  for (s in 1:n) {
    if(Time[s] < Censor[s]){
      tau = Censor[s]
      Ri =  rnbinom(1, 3, 1/(1+3))
      s_lb = 0
      s_rb = Censor[s]
      for (i in 1:Ri){
        t = runif(1, 0, tau)
        if (t < Time[s] && t > s_lb){
          s_lb = t
        } else if (t > Time[s] && t < s_rb) {
          s_rb = t
        }
      }
      lb = c(lb, s_lb)
      rb = c(rb, s_rb)
    } else {
      lb = c(lb, Censor[s])
      rb = c(rb, Inf)
    }
  }
  data = data.frame(cbind(lb, rb))
  colnames(data) = c("L", "R")
  data$status = ifelse(data$L == data$R, 1, 0)
  return(data)
}


#' Create the Current Status Test Data from Weibull Distribution
#'
#' @param seed Random seed
#' @param n Sample size
#' @param K Weibull shape parameter
#' @param Lambda Weibull rate parameter
#' @param q Right Censoring rate
#'
#' @returns A data frame with four columns: time (observation time), status (1 if event observed, 0 if censored), L (left interval), and R (right interval)
#' @export
#'
#' @examples
#' data_cs(seed = 123, n = 100, K = 1, Lambda = 0.1, q = 0.3)
#' data_cs(seed = 123, n = 100, K = 2, Lambda = 0.1, q = 0.3)
data_cs <- function(seed, n, K, Lambda, q = 0.3){
  Lambda_c = ((Lambda^K) / (1-q) - Lambda^K)^(1/K)
  set.seed(seed)
  T = rweibull(n, K, 1/Lambda)
  C = rweibull(n, K, 1/Lambda_c)
  delta = as.numeric(T<C)
  data = data.frame(time = C, status = delta)
  data$lb = ifelse(data$status == 1, 0, data$time)
  data$rb = ifelse(data$status == 1, data$time, Inf)
  colnames(data) = c("time", "status", "L", "R")
  data$status = ifelse(data$L == data$R, 1, 0)
  return(data)
}
