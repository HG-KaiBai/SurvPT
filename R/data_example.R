#' Simulate right-censored Weibull event times
#'
#' Generate independent Weibull event and censoring times with a common shape.
#' Event times have shape `K` and scale `1 / Lambda`.
#'
#' @param seed Integer random seed passed to `set.seed()`. The generator resets
#'   the global R random-number state.
#' @param n Positive integer sample size.
#' @param K Positive Weibull shape parameter.
#' @param Lambda Positive Weibull inverse-scale parameter. This is a rate only
#'   when `K = 1`; use `scale = 1 / Lambda` with [stats::rweibull()].
#' @param q Target probability that the censoring or inspection time precedes
#'   the event, with `0 < q < 1`. This determines the time distribution, not an
#'   exact realized sample proportion.
#' @details The censoring inverse scale is
#'   `(Lambda^K * q / (1 - q))^(1 / K)`. Arguments are not fully validated.
#' @return A data frame with `n` rows and columns `L` (observed follow-up time),
#'   `R` (equal to `L` for an exact event and `Inf` for right censoring), and
#'   `status` (`1` for an exact event, `0` for right censoring).
#' @seealso \code{\link{data_ic}()}, \code{\link{data_cs}()}, \code{\link{survPT}()}, \code{\link{survOPT}()}
#' @export
#'
#' @examples
#' d <- data_rc(seed = 123, n = 20, K = 2, Lambda = 0.4, q = 0.3)
#' head(d)
#' mean(d$status == 0)
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

#' Simulate interval-censored Weibull event times
#'
#' Generate Weibull event times, independent censoring times, and random
#' inspections before censoring. Observations are returned as interval endpoints.
#'
#' @inheritParams data_rc
#' @details Event and censoring distributions are as in \code{\link{data_rc}()}. For events
#'   before censoring, the inspection count is drawn with
#'   `rnbinom(1, size = 3, prob = 0.25)` and inspections are uniform before the
#'   censoring time. The legacy loop uses `1:Ri`; when `Ri` is zero this executes
#'   two iterations rather than none. This behavior is retained and should be
#'   considered when interpreting simulation results.
#' @return A data frame with `L`, `R`, and `status`. Right-censored rows have
#'   `R = Inf`; finite intervals bracket the event time. `status` tests `L == R`
#'   and does not indicate whether the event occurred before censoring. Pass
#'   `L` and `R` to \code{\link{survPT}()}.
#' @seealso \code{\link{data_rc}()}, \code{\link{data_cs}()}, \code{\link{survPT}()}
#' @export
#'
#' @examples
#' d <- data_ic(seed = 123, n = 20, K = 1.5, Lambda = 0.4, q = 0.3)
#' head(d[, c("L", "R")])
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


#' Simulate current-status observations of Weibull event times
#'
#' Observe only whether an event has occurred by an independent inspection time.
#'
#' @inheritParams data_rc
#' @details Event and inspection times use the Weibull distributions described
#'   in \code{\link{data_rc}()}. The returned `status` column is currently overwritten with
#'   `as.numeric(L == R)`: it is not the event-by-inspection indicator. For
#'   positive inspection times, recover that indicator with
#'   `as.integer(L == 0 & is.finite(R))`. Current-status data must be fitted via
#'   interval endpoints in \code{\link{survPT}()}, not as right-censored data in \code{\link{survOPT}()}.
#' @return A data frame with columns, in order, `L`, `R`, `time`, and `status`.
#'   If the event precedes inspection, the interval is `(0, time]`; otherwise
#'   it is `(time, Inf)`. See Details for the legacy `status` behavior.
#' @seealso \code{\link{data_rc}()}, \code{\link{data_ic}()}, \code{\link{survPT}()}
#' @export
#'
#' @examples
#' d <- data_cs(seed = 123, n = 20, K = 1.5, Lambda = 0.4, q = 0.3)
#' observed_by_inspection <- as.integer(d$L == 0 & is.finite(d$R))
#' head(data.frame(d[, c("L", "R", "time")], observed_by_inspection))
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



#' Generate fiducial-comparison scenario 1: exponential current status
#'
#' Generate current-status intervals from independent rate-one exponential
#' event and inspection times. This function generates data; it does not fit
#' a fiducial estimator.
#'
#' @param t Integer random seed passed to `set.seed()`. Resets the global
#'   R random-number state.
#' @param n Positive integer sample size.
#' @return A data frame with `n` rows and columns `L` and `R`, suitable for
#'   \code{\link{survPT}()}. Infinite right endpoints indicate right-open observation
#'   intervals, not an exact event at infinity.
#' @references
#' Cui, Y., Hannig, J., and Kosorok, M. R. (2024).
#' A Unified Nonparametric Fiducial Approach to Interval-Censored Data.
#' \emph{Journal of the American Statistical Association}, 119(547), 2230--2241.
#' \doi{10.1080/01621459.2023.2252143}.
#' @seealso \code{\link{fiducial_2}()}, \code{\link{fiducial_3}()}, \code{\link{fiducial_4}()}
#' @export
#'
#' @examples
#' head(fiducial_1(t = 123, n = 20))
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

#' Generate fiducial-comparison scenario 2: gamma current status
#'
#' Event times have a gamma distribution with shape 3 and rate 1. Independent
#' inspection times are uniform on `(0, 5)`. Only current status is observed.
#' This function generates data; it does not fit a fiducial estimator.
#'
#' @inheritParams fiducial_1
#' @inherit fiducial_1 return references
#' @seealso \code{\link{fiducial_1}()}, \code{\link{fiducial_3}()}, \code{\link{fiducial_4}()}
#' @export
#'
#' @examples
#' head(fiducial_2(t = 123, n = 20))
fiducial_2 <- function(t, n){
  set.seed(t)
  T = rgamma(n, 3, 1)
  C = runif(n, 0, 5)
  delta = as.numeric(T<C)
  l = ifelse(delta == 1, 0, C)
  r = ifelse(delta == 1, C, Inf)
  data = data.frame(L = l, R = r)
  colnames(data) = c("L", "R")
  return(data)
}

#' Generate fiducial-comparison scenario 3: two inspection times
#'
#' Event times are gamma with shape 2 and rate 1. The first inspection is
#' uniform on `(0, 2)`; the second is the first plus `0.5` plus an independent
#' uniform `(0, 2)` draw. Intervals represent events before, between, or after
#' the inspections. This function does not fit a fiducial estimator.
#'
#' @inheritParams fiducial_1
#' @inherit fiducial_1 return references
#' @seealso \code{\link{fiducial_1}()}, \code{\link{fiducial_2}()}, \code{\link{fiducial_4}()}
#' @export
#'
#' @examples
#' head(fiducial_3(t = 123, n = 20))
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

#' Generate fiducial-comparison scenario 4: random inspection counts
#'
#' Event times are exponential with rate 1. Each subject has one to four
#' inspections, selected with equal probabilities; inspection times are uniform
#' on `(0, 3)` and sorted. The returned interval brackets the event, with zero
#' or infinity used outside the inspection range. This function does not fit
#' a fiducial estimator.
#'
#' @inheritParams fiducial_1
#' @inherit fiducial_1 return references
#' @seealso \code{\link{fiducial_1}()}, \code{\link{fiducial_2}()}, \code{\link{fiducial_3}()}
#' @export
#'
#' @examples
#' head(fiducial_4(t = 123, n = 20))
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































