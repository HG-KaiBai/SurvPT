library(SurvPT)

# Independent analytic survival formulas check the population truth, including
# inverse-scale Weibull parameters and the two distinct Gamma scenarios.
cases <- list(
  data_rc = list(draw = function(seed, n) data_rc(seed, n, K = 1.7, Lambda = 20),
                 survival = function(t) exp(-(20 * t)^1.7),
                 columns = c("L", "R", "status")),
  data_ic = list(draw = function(seed, n) data_ic(seed, n, K = 1.7, Lambda = 20),
                 survival = function(t) exp(-(20 * t)^1.7),
                 columns = c("L", "R", "status")),
  data_cs = list(draw = function(seed, n) data_cs(seed, n, K = 1.7, Lambda = 20),
                 survival = function(t) exp(-(20 * t)^1.7),
                 columns = c("L", "R", "time", "status")),
  mixture_right = list(draw = function(seed, n) data_mixture(seed, n, "right"),
                       survival = function(t) 0.5 * exp(-(t / 3)^3) +
                         0.5 * exp(-(t / 10)^4),
                       columns = c("L", "R", "status")),
  mixture_interval = list(draw = function(seed, n) data_mixture(seed, n, "interval"),
                          survival = function(t) 0.5 * exp(-(t / 3)^3) +
                            0.5 * exp(-(t / 10)^4),
                          columns = c("L", "R", "status")),
  mixture_current = list(draw = function(seed, n) data_mixture(seed, n, "current"),
                         survival = function(t) 0.5 * exp(-(t / 3)^3) +
                           0.5 * exp(-(t / 10)^4),
                         columns = c("L", "R", "time", "status")),
  fiducial_1 = list(draw = fiducial_1, survival = function(t) exp(-t),
                    columns = c("L", "R")),
  fiducial_2 = list(draw = fiducial_2,
                    survival = function(t) exp(-t) * (1 + t + t^2 / 2),
                    columns = c("L", "R")),
  fiducial_3 = list(draw = fiducial_3, survival = function(t) exp(-t) * (1 + t),
                    columns = c("L", "R")),
  fiducial_4 = list(draw = fiducial_4, survival = function(t) exp(-t),
                    columns = c("L", "R"))
)

for (case in cases) {
  for (n in c(1L, 40L)) {
    for (seed in c(1L, 90L, 123L)) {
      d <- case$draw(seed, n)
      endpoints <- c(d$L, d$R)
      expected_tau <- max(endpoints[is.finite(endpoints)])
      tau <- attr(d, "tau", exact = TRUE)
      truth <- attr(d, "true_S_tau", exact = TRUE)
      proportions <- attr(d, "censoring_proportions", exact = TRUE)
      stopifnot(is.data.frame(d), nrow(d) == n,
                identical(names(d), case$columns),
                is.numeric(tau), length(tau) == 1L, is.finite(tau),
                identical(tau, expected_tau),
                is.numeric(truth), length(truth) == 1L,
                is.finite(truth), truth >= 0, truth <= 1,
                isTRUE(all.equal(truth, case$survival(expected_tau),
                                 tolerance = 1e-12)),
                is.numeric(proportions),
                identical(names(proportions), c("Exact", "Left", "Interval", "Right")),
                all(is.finite(proportions)), all(proportions >= 0),
                all(proportions <= 1), abs(sum(proportions) - 1) < 1e-12,
                identical(d, case$draw(seed, n)))
    }
  }
}
# One exact event at zero, two left-censored, three interval-censored and four
# right-censored observations. An uninformative (0, Inf) row is right-censored;
# the deliberately incorrect status column must not affect classification.
fixture <- data.frame(L = c(0, 0, 0, 1, 2, 3, 0, 1, 4, 7),
                      R = c(0, 1, 3, 2, 4, 6, Inf, Inf, Inf, Inf),
                      status = rep(1, 10))
summarize <- get("with_simulation_truth", envir = environment(data_rc))
annotated <- summarize(fixture, function(t) exp(-t))
stopifnot(identical(attr(annotated, "censoring_proportions"),
                    c(Exact = 0.1, Left = 0.2, Interval = 0.3, Right = 0.4)))
for (rows in list(1L, 2:3, 4:6, 7:10)) {
  one_type <- summarize(fixture[rows, ], function(t) exp(-t))
  proportions <- attr(one_type, "censoring_proportions")
  stopifnot(sum(proportions == 1) == 1L, sum(proportions == 0) == 3L)
}
cat("All data-generator truth and censoring-proportion checks passed.\n")
