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
  lb = ifelse(delta == 1, 0, C)
  rb = ifelse(delta == 1, C, Inf)
  data = data.frame(L = lb, R = rb, time = C, status = delta)
  colnames(data) = c("L", "R", "time", "status")
  data$status = ifelse(data$L == data$R, 1, 0)
  return(data)
}



#' Scenario 1 from Fiducial Approach paper (Cui et.al. 2024)
#'
#' @param t Random seed
#' @param n Sample size
#'
#' @returns A data frame with two columns: L (left interval) and R (right interval)
#' @export
#'
#' @examples
#' data = fiducial_1(t = 123, n = 100)
fiducial_1 <- function(t, n){
  set.seed(t)
  T = rexp(n, rate = 1)
  C = rexp(n, rate = 1)
  delta = as.numeric(T<C)
  l = ifelse(delta == 1, 0, C)
  r = ifelse(delta == 1, C, Inf)
  data = data.frame(L = l, R = r)
  colnames(data) = c("L", "R")
  return(data)
}

#' Scenario 2 from Fiducial Approach paper (Cui et.al. 2024)
#'
#' @param t Random seed
#' @param n Sample size
#'
#' @returns A data frame with two columns: L (left interval) and R (right interval)
#' @export
#'
#' @examples
#' data = fiducial_2(t = 123, n = 100)
fiducial_2 <- function(t, n){
  set.seed(t)
  T = rgamma(n, 3, 1)
  C = runif(n, 0, 10)
  delta = as.numeric(T<C)
  l = ifelse(delta == 1, 0, C)
  r = ifelse(delta == 1, C, Inf)
  data = data.frame(L = l, R = r)
  colnames(data) = c("L", "R")
  return(data)
}

#' Scenario 3 from Fiducial Approach paper (Cui et.al. 2024)
#'
#' @param t Random seed
#' @param n Sample size
#'
#' @returns A data frame with two columns: L (left interval) and R (right interval)
#' @export
#'
#' @examples
#' data = fiducial_3(t = 123, n = 100)
fiducial_3 <- function(t, n){
  set.seed(t)
  T = rgamma(n, 2, 1)
  C1 = runif(n, 0, 2)
  C2 = C1 + 0.5 + runif(n, 0, 2)
  lb = ifelse(T < C1, 0, ifelse(T >= C2, C2, C1))
  rb = ifelse(T < C1, C1, ifelse(T >= C2, Inf, C2))
  data = data.frame(L = lb, R = rb)
  colnames(data) = c("L", "R")
  return(data)
}

#' Scenario 4 from Fiducial Approach paper (Cui et.al. 2024)
#'
#' @param t Random seed
#' @param n Sample size
#'
#' @returns A data frame with two columns: L (left interval) and R (right interval)
#' @export
#'
#' @examples
#' data = fiducial_4(t = 123, n = 100)
fiducial_4 <- function(t, n){
  set.seed(t)
  T = rexp(n, 1)
  n_obs = sample(1:4, n, replace = TRUE, prob = rep(0.25, 4))
  bounds = matrix(0, nrow = n, ncol = 2)
  for (i in 1:n){
    times = sort(runif(n_obs[i], 0, 3))
    idx = which(times <= T[i])
    if (length(idx) == 0){
      bounds[i,1] <- 0
    } else {
      bounds[i,1] <- times[max(idx)]
    }
    idx = which(times > T[i])
    if (length(idx) == 0){
      bounds[i,2] <- Inf
    } else {
      bounds[i,2] <- times[min(idx)]
    }
  }
  data = data.frame(L = bounds[,1], R = bounds[,2])
  colnames(data) = c("L", "R")
  return(data)
}































