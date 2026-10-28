library(SurvPT)

body_mass <- function(fit) if (is.na(fit$tau)) 1 else as.numeric(fit$Posterior_Probs[2])
tail_mass <- function(fit) if (is.na(fit$tau)) 0 else as.numeric(fit$Posterior_Probs[3])

near <- function(x, y, tolerance = 1e-8) {
    stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
              max(abs(x - y)) < tolerance)
}
fit_tail <- function(x, prior = 1, depth = 3, iterations = 0, points = NULL, ...) {
    survPT(x, prior = prior, prior_params = switch(prior, NULL, 0.7, c(0.7, 1.5)),
           max_M = depth, max_EM_iter = iterations, EM_tol = 1e-10,
           test_seq = points, verbose = FALSE, ...)
}

# Independently integrate finite leaves using R distribution functions.
body_probability <- function(fit, l, r, prior) {
    if (l >= r) return(0)
    cdf <- switch(prior, function(t) t, function(t) pexp(t, 0.7),
                  function(t) pweibull(t, shape = 1.5, scale = 1/0.7))
    sum(vapply(setdiff(fit$Terminal_Nodes, 3), function(j) {
        a <- fit$Left[j]; b <- fit$Right[j]
        lo <- max(a, l); hi <- min(b, r)
        if (lo >= hi) return(0)
        fit$Posterior_Probs[j]/body_mass(fit) * (cdf(hi)-cdf(lo))/(cdf(b)-cdf(a))
    }, numeric(1)))
}

# Includes exact tau, exact internal boundaries, finite intervals ending at tau,
# right censoring before tau and right censoring exactly at tau. Extra columns
# must not affect tau.
x <- cbind(rbind(c(0.25, 0.25), c(1, 1), c(2, 2), c(4, 4), c(0, 0.75),
                 c(1.5, 4), c(0.5, Inf), c(3, Inf), c(4, Inf)), ignored = 999)
for (prior in 1:3) for (depth in c(0, 3)) {
    points <- c(0, 0.25, 1, 2, 4, 4.01, 100, Inf)
    initial <- fit_tail(x, prior, depth, points = points)
    final <- fit_tail(x, prior, depth, iterations = 30, points = points)
    near(initial$tau, 4)

    near(body_mass(final) + tail_mass(final), 1)
    near(as.numeric(initial$Level[1:3]), c(0, 1, 1))
    stopifnot(3 %in% initial$Terminal_Nodes,
              max(initial$Level) <= depth + 1,
              initial$Right[2] == 4, is.infinite(initial$Right[3]))

    censored <- is.infinite(x[, 2])
    # Independent initial responsibilities: prior mean 1/2 and conditional base.
    base_cdf <- switch(prior, function(t) t, function(t) pexp(t, 0.7),
                       function(t) pweibull(t, shape = 1.5, scale = 1/0.7))
    s <- 1 - base_cdf(x[censored, 1])/base_cdf(4)
    expected_body <- sum(!censored) + sum(s/(1+s))
    near(body_mass(initial), (expected_body - 0.5)/(nrow(x)-1))

    # Follow two E/M steps using independent interval probabilities.
    previous <- initial
    for (step in 1:2) {
        current <- fit_tail(x, prior, depth, iterations = step, points = points)
        s <- vapply(x[censored, 1], function(c)
            body_probability(previous, c, 4, prior), numeric(1))
        expected_body <- sum(!censored) + sum(body_mass(previous)*s /
            (body_mass(previous)*s + tail_mass(previous)))
        near(current$Allocations[2], expected_body)
        near(body_mass(current), (expected_body - 0.5)/(nrow(x)-1))
        previous <- current
    }

    # Check every requested likelihood contribution against independently
    # computed finite-subtree probabilities and returned exact-event densities.
    for (fit in list(initial, final)) {
        likelihood <- vapply(seq_len(nrow(x)), function(i) {
            l <- x[i, 1]; r <- x[i, 2]
            if (l == r) return(as.numeric(fit$Density[match(l, points)]))
            if (is.infinite(r)) return(body_mass(fit) *
                body_probability(fit, l, fit$tau, prior) + tail_mass(fit))
            body_mass(fit) * body_probability(fit, l, r, prior)
        }, numeric(1))
        stopifnot(all(likelihood > 0))
        near(fit$LogLik, sum(log(likelihood)))
        near(as.numeric(fit$Survival_Probs[5]), tail_mass(fit))
        stopifnot(fit$Density[5] > 0, all(is.na(fit$Survival_Probs[6:7])),
                  all(is.na(fit$Density[6:7])), fit$Survival_Probs[8] == 0,
                  fit$Density[8] == 0)
        near(sum(fit$Posterior_Probs[fit$Terminal_Nodes]), 1)
        near(tail(as.numeric(fit$LogLik_Trace), 1), fit$LogLik)
    }
}

# With a single uniform body leaf, two finite observations and one censoring
# at tau give the Beta(1/2,1/2) mode (2-1/2)/(3-1) = 3/4.
simple <- rbind(c(1, 1), c(2, 2), c(3, Inf))
for (iterations in c(0, 1, 50)) {
    fit <- fit_tail(simple, depth = 0, iterations = iterations, points = c(1, 3))
    near(body_mass(fit), 3/4)
    near(as.numeric(fit$Density), rep(1/4, 2))
    near(as.numeric(fit$Survival_Probs), c(3/4, 1/4))
}

# Boundary mode updates and a default grid built from all finite endpoints.
for (prior in 1:3) {
    exact <- fit_tail(cbind(c(1, 2), c(1, 2)), prior, iterations = 10)
    near(body_mass(exact), 1)
    near(tail_mass(exact), 0)
    censored <- fit_tail(cbind(c(1, 4, 2), Inf), prior, iterations = 10)
    near(body_mass(censored), 0)
    near(tail_mass(censored), 1)
    near(range(censored$Test_seq), c(0, 4))
    near(as.numeric(censored$Survival_Probs), rep(1, 200))
    near(as.numeric(censored$Density), rep(0, 200))
    near(censored$LogLik, 0)
    stopifnot(censored$Converged)
}

# Bounded support has no added tail and preserves the original alpha update.
bounded <- fit_tail(rbind(c(1, 1), c(2, 2), c(2.5, 3)), depth = 1,
                    iterations = 1, support = c(0, 4))
stopifnot(is.na(bounded$tau), max(bounded$Level) == 1)
near(body_mass(bounded), 1)
near(tail_mass(bounded), 0)
near(as.numeric(bounded$Posterior_Probs[2:3]), c(0.75, 0.25))

# The extra root must not shift the finite tree's alpha schedule: first
# body split still has alpha = 1/2, hence (2 - 1/2)/(3 - 1) = 3/4.
unbounded <- fit_tail(cbind(c(1, 2, 4), c(1, 2, 4)), depth = 1, iterations = 2)
near(as.numeric(unbounded$Posterior_Probs[4:5]), c(0.75, 0.25))

# Reject a degenerate tau and invalid endpoints in both public entry points.
stopifnot(inherits(try(fit_tail(matrix(c(0, Inf), 1)), silent = TRUE), "try-error"),
          inherits(try(fit_tail(matrix(c(Inf, Inf), 1)), silent = TRUE), "try-error"),
          inherits(try(SurvPT(matrix(c(0, Inf), 1), c(0, 1), 1, c(0, Inf),
                              0L, c(0.12, 1), 1L, 1L, 1e-6, FALSE),
                       silent = TRUE), "try-error"))

# The first layer updates internally through the original public interface.
mixed <- rbind(c(1, 1), c(2, 2), c(1, Inf), c(3, Inf))
before <- fit_tail(mixed, depth = 0)
after <- fit_tail(mixed, depth = 0, iterations = 1)
initial_body <- 2 + (2/3)/(1 + 2/3)
near(body_mass(before), (initial_body - 0.5)/3)
expected_body <- 2 + (body_mass(before) * (2/3))/
    (body_mass(before) * (2/3) + tail_mass(before))
near(body_mass(after), (expected_body - 0.5)/3)
stopifnot(abs(body_mass(after) - body_mass(before)) > 1e-5)
native <- SurvPT(mixed, c(0, 3), 1, c(0, Inf), 0L, c(0.12, 1),
                 1L, 1L, 1e-10, FALSE)
near(as.numeric(native$alpha[2:3]), c(0.5, 0.5))
near(as.numeric(native$posterior_probs[2]), body_mass(after))
stopifnot(length(formals(survPT)) == 10L, length(formals(SurvPT)) == 10L,
          !any(c("q_0", "theta_0", "p_tail", "q0_prior") %in% names(after)),
          !any(c("q_0", "theta_0", "p_tail", "q0_prior") %in% names(native)))

# Known body/tail membership gives a directly checkable Beta posterior mode.
known <- fit_tail(simple, depth = 0, iterations = 3)
mode <- optimize(function(q) dbeta(q, 0.5+2, 0.5+1, log = TRUE),
                 c(0, 1), maximum = TRUE, tol = 1e-12)$maximum
near(body_mass(known), mode)
cat("Tail-root Beta initialization, EM updates, likelihood and boundary tests passed.\n")