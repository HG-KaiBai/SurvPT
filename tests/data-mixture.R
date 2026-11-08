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
  event_prob <- 0.5 * pweibull(current$time, 3, 3) +
                0.5 * pweibull(current$time, 4, 10)
  bins <- cut(current$time, c(0, 2, 4, 6, 8, 10, Inf))
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
for (x in c(2, 4, 6, 8, 10, 12)) {
  prob <- 0.5 * pweibull(x, 3, 3) + 0.5 * pweibull(x, 4, 10)
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

# The underlying density retains two peaks separated by a valley.
mixture_density <- function(x) {
  0.5 * dweibull(x, 3, 3) + 0.5 * dweibull(x, 4, 10)
}
first_peak <- optimize(mixture_density, c(2, 5), maximum = TRUE)
second_peak <- optimize(mixture_density, c(7, 10), maximum = TRUE)
valley <- optimize(mixture_density, c(first_peak$maximum, second_peak$maximum))
stopifnot(first_peak$maximum < valley$minimum,
          valley$minimum < second_peak$maximum,
          first_peak$objective > valley$objective,
          second_peak$objective > valley$objective)

# This finite-seed regression checks the requested simulation settings, not a
# guarantee for arbitrary seeds. No failed datasets are discarded or redrawn.
mixture_survival <- function(x) {
  0.5 * pweibull(x, 3, 3, lower.tail = FALSE) +
    0.5 * pweibull(x, 4, 10, lower.tail = FALSE)
}
target_time <- uniroot(function(x) mixture_survival(x) - 0.25,
                       c(0, 30), tol = 1e-12)$root
coverage_checks <- expand.grid(n = c(200L, 500L), q = c(0.3, 0.5),
                               type = c("right", "interval", "current"),
                               stringsAsFactors = FALSE)
coverage_checks$failures <- integer(nrow(coverage_checks))
coverage_checks$min_tau <- numeric(nrow(coverage_checks))
for (j in seq_len(nrow(coverage_checks))) {
  setting <- coverage_checks[j, ]
  taus <- vapply(seq_len(1000L), function(seed) {
    d <- data_mixture(seed, setting$n, setting$type, setting$q)
    endpoints <- as.matrix(d[, c("L", "R")])
    max(endpoints[is.finite(endpoints)])
  }, numeric(1))
  coverage_checks$failures[j] <- sum(taus <= target_time)
  coverage_checks$min_tau[j] <- min(taus)
}
print(coverage_checks, row.names = FALSE)
stopifnot(all(coverage_checks$failures == 0L))
cat("Weibull mixture censoring checks passed.\n")
