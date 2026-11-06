library(SurvPT)

# All three encodings use the same latent event and censoring times per seed.
for (q in c(0.1, 0.3, 0.5, 0.8)) {
  n <- 20000
  right <- data_mixture(2026, n, "right", q)
  interval <- data_mixture(2026, n, "interval", q)
  current <- data_mixture(2026, n, "current", q)
  censored <- is.infinite(right$R)
  exact <- !censored
  stopifnot(
    abs(mean(censored) - q) < 6 * sqrt(q * (1 - q) / n),
    identical(censored, is.infinite(interval$R)),
    identical(censored, is.infinite(current$R)),
    all(interval$L[exact] < right$L[exact]),
    all(interval$R[exact] >= right$L[exact]),
    all(interval$R[exact] <= current$time[exact]),
    identical(interval$L[censored], right$L[censored]),
    identical(current$L[censored], right$L[censored]),
    all(current$L[exact] == 0),
    identical(current$R[exact], current$time[exact]),
    all(right$status == as.numeric(exact)),
    all(interval$status == 0), all(current$status == 0)
  )
  for (type in c("right", "interval", "current")) {
    d <- data_mixture(123, 1, type, q)
    stopifnot(nrow(d) == 1L, !anyNA(d), d$L >= 0, d$L <= d$R,
              identical(d, data_mixture(123, 1, type, q)))
  }

  # With independent inspections, the event-by-inspection probability at C
  # must equal the specified mixture CDF, including within inspection bins.
  event_prob <- 0.5 * pweibull(current$time, 5, 1) +
                0.5 * pweibull(current$time, 5, 4.5)
  bins <- cut(current$time, c(0, 1, 2, 4, 5, Inf))
  for (bin in levels(bins)) {
    take <- bins == bin
    residual <- sum(exact[take] - event_prob[take])
    sd_sum <- sqrt(sum(event_prob[take] * (1 - event_prob[take])))
    stopifnot(abs(residual) < 6 * sd_sum + 1)
  }
}

# Almost uncensored samples should follow the mixture, not an average of
# Weibull draws or a single Weibull with the average component scale.
d <- data_mixture(739, 100000, "right", q = 1e-10)
for (x in c(0.8, 1.2, 3, 4.5, 5)) {
  prob <- 0.5 * pweibull(x, 5, 1) + 0.5 * pweibull(x, 5, 4.5)
  stopifnot(abs(mean(d$L <= x) - prob) < 6 * sqrt(prob * (1 - prob) / nrow(d)))
}

errors <- function(...) {
  inherits(tryCatch(data_mixture(...), error = identity), "error")
}
for (n in list(0, -1, 1.5, NA_real_, Inf, numeric(), c(1, 2), "10")) {
  stopifnot(errors(123, n))
}
for (q in list(0, 1, -0.1, 1.1, NA_real_, Inf, numeric(), c(0.2, 0.3), "0.3")) {
  stopifnot(errors(123, 10, q = q))
}
stopifnot(errors(123, 10, type = "unknown"))
cat("Weibull mixture censoring checks passed.\n")
