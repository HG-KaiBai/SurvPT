library(SurvPT)

assert_close <- function(actual, expected, label, tolerance = 1e-9) {
    error <- max(abs(actual - expected))
    if (!is.finite(error) || error > tolerance) stop(label, ': error = ', error)
}

# Independent reference base uses R's uniform and exponential distributions.
reference_base <- function(M, rate) {
    body <- pexp(M, rate)
    list(
        mass = function(left, right) {
            if (right <= left) return(0)
            body * (punif(min(right, M), 0, M) - punif(min(left, M), 0, M)) +
                pexp(max(left, M), rate, lower.tail = FALSE) -
                pexp(max(right, M), rate, lower.tail = FALSE)
        },
        density = function(t) {
            if (t <= M) body * dunif(t, 0, M) else dexp(t, rate)
        },
        cdf = function(t) ifelse(t <= M, body * punif(t, 0, M), pexp(t, rate)),
        quantile = function(p) if (p <= body) qunif(p/body, 0, M) else qexp(p, rate)
    )
}
mixture_mass <- function(fit, left, right, base) {
    sum(vapply(fit$Terminal_Nodes, function(j) {
        a <- fit$Left[j]; b <- fit$Right[j]
        fit$Posterior_Probs[j] * base$mass(max(left, a), min(right, b)) / base$mass(a, b)
    }, numeric(1)))
}
mixture_density <- function(fit, t, base) {
    if (!is.finite(t)) return(0)
    sum(vapply(fit$Terminal_Nodes, function(j) {
        a <- fit$Left[j]; b <- fit$Right[j]
        if (t < a || t > b || (t == a && a != 0)) return(0)
        fit$Posterior_Probs[j] * base$density(t) / base$mass(a, b)
    }, numeric(1)))
}

observations <- rbind(c(0, 0), c(0.2, 0.2), c(0, 1.2), c(1, 2.5),
                      c(2, 2), c(2, 4), c(3, Inf), c(4, Inf))
checks <- 0L
for (finite in c(FALSE, TRUE)) for (depth in c(0, 1, 5)) {
    data <- observations
    if (finite) data[is.infinite(data)] <- 6
    M <- max(data[is.finite(data)])
    for (rate in c(0.03, log(2)/M, 0.5, 5)) {
        base <- reference_base(M, rate)
        points <- sort(unique(c(-1, 0, seq(0, 9, length.out = 35), M - 1e-8,
                                M, M + 1e-8, M + 1/rate, Inf)))
        fit_at <- function(iterations, grid = points) {
            survPT(cbind(data, ignored = 9999), prior = 4, prior_params = rate,
                   support = c(0, if (finite) 6 else Inf), max_M = depth,
                   max_EM_iter = iterations, EM_tol = 0, test_seq = grid, verbose = FALSE)
        }
        fits <- lapply(c(0, 1, 2, 15), fit_at)
        for (step in 1:3) {
            current <- fits[[step]]
            previous <- if (step > 1) fits[[step - 1]] else NULL
            expected <- vapply(seq_along(current$Left), function(j) {
                a <- current$Left[j]; b <- current$Right[j]
                sum(vapply(seq_len(nrow(data)), function(i) {
                    left <- data[i, 1]; right <- data[i, 2]
                    if (left == right) return(as.numeric((left > a || (left == 0 && a == 0)) && right <= b))
                    if (is.null(previous)) {
                        return(base$mass(max(left, a), min(right, b)) / base$mass(left, right))
                    }
                    mixture_mass(previous, max(left, a), min(right, b), base) /
                        mixture_mass(previous, left, right, base)
                }, numeric(1)))
            }, numeric(1))
            assert_close(current$Allocations, expected, 'Initialization / EM allocations')
        }
        for (fit in fits) {
            assert_close(sum(fit$Posterior_Probs[fit$Terminal_Nodes]), 1, 'Total mass')
            assert_close(as.numeric(fit$Survival_Probs), vapply(points, function(t)
                mixture_mass(fit, t, Inf, base), numeric(1)), 'Survival')
            assert_close(as.numeric(fit$Density), vapply(points, function(t)
                mixture_density(fit, t, base), numeric(1)), 'Density')
            stopifnot(fit$Tail_Start == M, fit$Prior_Params == rate,
                      all(is.finite(fit$Posterior_Probs)), all(fit$Posterior_Probs > 0),
                      all(diff(as.numeric(fit$Survival_Probs)) <= 1e-12))
            for (parent in setdiff(seq_along(fit$Left), fit$Terminal_Nodes)) {
                children <- which(fit$Level == fit$Level[parent] + 1 &
                                  fit$Left >= fit$Left[parent] & fit$Right <= fit$Right[parent])
                stopifnot(length(children) == 2)
                split <- fit$Right[children[1]]
                assert_close(base$mass(fit$Left[parent], split) /
                             base$mass(fit$Left[parent], fit$Right[parent]), 0.5, 'Base median split')
                q <- fit$Posterior_Probs[children] / fit$Posterior_Probs[parent]
                stopifnot(all(q >= 1e-6 - 1e-15), all(q <= 1 - 1e-6 + 1e-15))
            }
        }
        if (depth == 1) {
            expected <- base$quantile((base$cdf(0) + base$cdf(if (finite) 6 else Inf))/2)
            assert_close(fits[[1]]$Right[2], expected, 'Inverse CDF')
        }
        if (!finite && depth == 0) {
            # At depth zero the output is exactly the specified spliced base.
            stopifnot(fits[[1]]$Survival_Probs[match(M + 1/rate, points)] > 0)
            assert_close(fits[[1]]$Survival_Probs[match(M, points)] / exp(-rate*M), 1, 'Base tail mass')
        }
        # Integrate the actual output, splitting at the density discontinuity.
        final <- fits[[4]]
        integral <- sum(vapply(final$Terminal_Nodes, function(j) {
            bounds <- sort(unique(c(final$Left[j], final$Right[j],
                                    if (M > final$Left[j] && M < final$Right[j]) M)))
            sum(vapply(seq_len(length(bounds) - 1L), function(k) {
                integrate(function(t) as.numeric(fit_at(15, t)$Density),
                          bounds[k], bounds[k + 1], rel.tol = 1e-7)$value
            }, numeric(1)))
        }, numeric(1)))
        assert_close(integral, 1, 'Density integral', 1e-6)
        checks <- checks + 1L
    }
}

# A root split at M assigns an exact event at the splice to the body child.
boundary <- rbind(c(0, 0), c(1, 1), c(2, 2), c(2, Inf))
rate <- log(2)/2
fit <- survPT(boundary, prior = 4, prior_params = rate, max_M = 1,
              max_EM_iter = 0, test_seq = c(2, 2 + 1e-8), verbose = FALSE)
assert_close(fit$Right[2], 2, 'Splice cut point')
assert_close(as.numeric(fit$Allocations[2:3]), c(3, 1), 'Splice allocation')

# Tail CDFs round to one here; survival-scale inversion must still split.
tail_data <- matrix(c(1, Inf), 1, 2)
tail <- survPT(tail_data, prior = 4, prior_params = 100, support = c(1, Inf),
               max_M = 1, max_EM_iter = 2, test_seq = c(1, 1.01, 1.1), verbose = FALSE)
assert_close(tail$Right[2], 1 + log(2)/100, 'Extreme-tail median')
stopifnot(length(tail$Terminal_Nodes) == 2L, all(is.finite(tail$Survival_Probs)))

# Legacy overlap mode, raw interface, and the unused uniform factor.
for (allocated in 0:1) {
    args <- list(X = observations, prior = 4, prior_params = 0.5,
                 allocated = allocated, max_M = 5, max_EM_iter = 5,
                 test_seq = c(0, 4, 10), verbose = FALSE)
    wrapped <- do.call(survPT, args)
    stopifnot(identical(wrapped, do.call(survPT, c(args, list(uniform_tail_factor = 2)))))
    raw <- do.call(SurvPT, c(args, list(support = c(0, Inf), EM_tol = 1e-6)))
    stopifnot(identical(wrapped$Allocations, raw$n_obs_allocated),
              identical(wrapped$Survival_Probs, raw$survival_probs))
    assert_close(as.numeric(raw$alpha), pmax(1, 0.5 * as.numeric(raw$node_levels)^2), 'Prior alpha')
}
for (bad in list(NULL, numeric(), c(0.5, 1), NA_real_, NaN, Inf, 0, -1, '0.5', 1i)) {
    error <- tryCatch(survPT(observations, prior = 4, prior_params = bad, verbose = FALSE), error = identity)
    stopifnot(inherits(error, 'error'))
}
for (bad_data in list(matrix(0, 1, 2), matrix(numeric(), 0, 2),
                      matrix(c(-1, 2), 1, 2), matrix(c(Inf, Inf), 1, 2))) {
    error <- tryCatch(survPT(bad_data, prior = 4, prior_params = 0.5, verbose = FALSE), error = identity)
    stopifnot(inherits(error, 'error'))
}
cat('All', checks, 'spliced-base configurations and boundary/interface checks passed.\n')
