library(SurvPT)

assert_close <- function(actual, expected, tolerance = 1e-10) {
    stopifnot(length(actual) == length(expected),
              max(abs(actual - expected)) < tolerance)
}
expect_error <- function(expression, pattern) {
    error <- tryCatch(force(expression), error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}

# At one split, exact data give an analytically known Beta posterior.
times <- c(rep(0.25, 8), rep(0.75, 3))
model <- survPT(cbind(times, times), support = c(0, 1), max_M = 1,
                flat = TRUE, verbose = FALSE)
stopifnot(inherits(model, "survPT"), all(model$Alpha == 1),
          !"test_seq" %in% names(formals(survPT)),
          !any(c("Test_seq", "Survival_Probs", "Density") %in% names(model)))
expect_error(survPT(cbind(times, times), test_seq = 0.5), "unused argument")
set.seed(123)
posterior <- predict(model, c(0, 0.5, 1), n_samples = 12000, burn_in = 0,
                      keep_samples = TRUE)
assert_close(posterior$estimate, c(1, 3/11, 0))
assert_close(posterior$mean, c(1, 4/13, 0), 0.01)
assert_close(c(posterior$lower[2], posterior$median[2], posterior$upper[2]),
             qbeta(c(0.025, 0.5, 0.975), 4, 9), 0.015)
draws <- attr(posterior, "samples")
stopifnot(identical(dim(draws), c(12000L, 3L)), all(draws[, 3] == 0))
assert_close(draws[, 1], rep(1, nrow(draws)))
assert_close(posterior$mean, colMeans(draws))
assert_close(posterior$median, apply(draws, 2, median))

# At depth two the working prior is Beta(2,2). The occupied left sub-branch
# therefore uses mode 9/10, rather than posterior mean 10/12.
working <- survPT(cbind(times, times), support = c(0, 1), max_M = 2,
                  verbose = FALSE)
mode_prediction <- predict(working, 0.25, sample = FALSE)
assert_close(mode_prediction$estimate, 1 - (8/11) * (9/10))
stopifnot(abs(mode_prediction$estimate - (1 - (9/13) * (10/12))) > 0.05)

# Twelve completely uninformative intervals must NOT become twelve observed
# fractional counts. One partial interval has likelihood .2 + .6*S(.5).
# With one exact observation per side and flat priors the exact posterior is
# .4*Beta(2,2) + .6*Beta(3,2), independent of those twelve full intervals.
data <- rbind(c(0.1, 0.1), c(0.9, 0.9), c(0.4, 0.9),
              matrix(rep(c(0, 1), 12), ncol = 2, byrow = TRUE))
censored <- survPT(data, support = c(0, 1), max_M = 1, flat = TRUE,
                   max_EM_iter = 500, EM_tol = 1e-10, verbose = FALSE)
cdf <- function(q) 0.4 * pbeta(q, 2, 2) + 0.6 * pbeta(q, 3, 2)
expected_quantiles <- vapply(c(0.025, 0.5, 0.975), function(p)
    uniroot(function(q) cdf(q) - p, c(0, 1), tol = 1e-12)$root, numeric(1))
set.seed(456)
sampled <- predict(censored, 0.5, n_samples = 30000, burn_in = 1000,
                    keep_samples = TRUE)
assert_close(sampled$mean, 0.56, 0.02)
assert_close(c(sampled$lower, sampled$median, sampled$upper), expected_quantiles, 0.025)
assert_close(var(as.numeric(attr(sampled, "samples"))), 0.36 - 0.56^2, 0.007)

# Sampling must not be clipped to the MAP optimization box.
uninformative <- survPT(matrix(c(0, 1), 1, 2), support = c(0, 1),
                        max_M = 1, flat = TRUE, verbose = FALSE)
set.seed(321)
unrestricted <- predict(uninformative, 0.5, n_samples = 20000, burn_in = 100,
                         keep_samples = TRUE)
assert_close(unrestricted$mean, 0.5, 0.02)
assert_close(c(unrestricted$lower, unrestricted$upper), c(0.025, 0.975), 0.02)

# Reproducibility, thinning, input order, and serializing a fitted object.
points <- c(0.8, 0.2, 0.8, 0.5)
before <- serialize(censored, NULL)
set.seed(902)
first <- predict(censored, points, n_samples = 12, burn_in = 3, keep_samples = TRUE)
set.seed(902)
second <- predict(censored, points, n_samples = 12, burn_in = 3, keep_samples = TRUE)
stopifnot(identical(first, second), identical(serialize(censored, NULL), before),
          identical(first$x, points), first$estimate[1] == first$estimate[3],
          all(attr(first, "samples")[, 1] == attr(first, "samples")[, 3]))
set.seed(902)
thinned <- predict(censored, points, n_samples = 6, burn_in = 3, thin = 2,
                    keep_samples = TRUE)
stopifnot(identical(attr(thinned, "samples"), attr(first, "samples")[seq(2, 12, 2), ]))
saved <- tempfile(fileext = ".rds")
saveRDS(censored, saved)
restored <- readRDS(saved)
unlink(saved)
stopifnot(identical(predict(censored, points, sample = FALSE),
                    predict(restored, points, sample = FALSE)))
seed <- .Random.seed
point_only <- predict(censored, points, sample = FALSE)
stopifnot(identical(seed, .Random.seed), identical(names(point_only), c("x", "estimate")),
          identical(point_only$estimate, first$estimate), is.null(attr(point_only, "samples")))
stopifnot(identical(point_only, predict(censored, points, n_samples = 0)))
empty <- predict(censored, numeric(), keep_samples = TRUE)
stopifnot(nrow(empty) == 0L, identical(dim(attr(empty, "samples")), c(5000L, 0L)))

# All bases and both prior strengths: proper survival curves, ordered intervals,
# fixed root mass, mixed censoring, and agreement of repeated time points.
mixed <- rbind(c(0, 0), c(0.2, 0.2), c(0, 1), c(0.5, 2), c(1.5, Inf), c(3, Inf))
for (prior in 1:4) for (flat in c(FALSE, TRUE)) {
    fit <- survPT(mixed, prior = prior, prior_params = if (prior == 3) c(0.5, 1.5) else 0.5,
                  max_M = 4, flat = flat, verbose = FALSE)
    if (flat) stopifnot(all(fit$Alpha == 1))
    else assert_close(as.numeric(fit$Alpha), pmax(1, 0.5 * as.numeric(fit$Level)^2))
    set.seed(789)
    result <- predict(fit, c(-Inf, 0, 0.2, 1, 2, 4, Inf), n_samples = 150,
                       burn_in = 50, level = 0.8, keep_samples = TRUE)
    curves <- attr(result, "samples")
    stopifnot(all(is.finite(curves)), all(curves >= 0), all(curves <= 1),
              all(apply(curves, 1, function(s) all(diff(s) <= 1e-12))),
              all(result$lower <= result$median), all(result$median <= result$upper))
    assert_close(result$lower, apply(curves, 2, quantile, probs = 0.1))
    assert_close(result$upper, apply(curves, 2, quantile, probs = 0.9))
    root <- survPT(mixed, prior = prior, prior_params = fit$Prior_Params,
                   max_M = 0, verbose = FALSE)
    root_prediction <- predict(root, c(0, 1, Inf), n_samples = 4, burn_in = 0)
    assert_close(root_prediction$mean, root_prediction$estimate)
    assert_close(root_prediction$lower, root_prediction$upper)
}

for (bad in list(NA, NaN, NULL, c(0.8, 0.95), 0, 1, -1, Inf, "0.95")) {
    expect_error(predict(model, 0.5, level = bad), "level")
}
for (bad in list(NA, NaN, NULL, -1, 1.5, Inf, c(1, 2), "10")) {
    expect_error(predict(model, 0.5, n_samples = bad), "n_samples")
}
expect_error(predict(model, c(1, NA)), "x must")
expect_error(predict(model, matrix(1, 1, 1)), "x must")
expect_error(predict(model, 0.5, burn_in = -1), "burn_in")
expect_error(predict(model, 0.5, thin = 0), "thin")
expect_error(predict(model, 0.5, keep_samples = NA), "keep_samples")
for (bad in list(NA, NULL, c(TRUE, FALSE), 0, "FALSE")) {
    expect_error(predict(model, 0.5, sample = bad), "sample")
}
expect_error(predict(model, 0.5, typo = 1), "Unused arguments")
cat("Fit/predict separation, analytic posteriors, censoring Gibbs updates and prediction controls passed.\n")
