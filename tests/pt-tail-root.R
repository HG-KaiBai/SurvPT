library(SurvPT)

epsilon <- 1e-6
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
    near(body_mass(initial), expected_body/nrow(x))

    # Follow two E/M steps using independent interval probabilities.
    previous <- initial
    for (step in 1:2) {
        current <- fit_tail(x, prior, depth, iterations = step, points = points)
        s <- vapply(x[censored, 1], function(c)
            body_probability(previous, c, 4, prior), numeric(1))
        expected_body <- sum(!censored) + sum(body_mass(previous)*s /
            (body_mass(previous)*s + tail_mass(previous)))
        near(current$Allocations[2], expected_body)
        near(body_mass(current), expected_body/nrow(x))
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
        if (prior == 1) {
            near(as.numeric(fit$Survival_Probs[6:7]), rep(tail_mass(fit), 2))
        } else if (prior == 2) {
            expected_tail <- tail_mass(fit) * pexp(points[6:7] - fit$tau,
                                                   rate = 0.7, lower.tail = FALSE)
            near(as.numeric(fit$Survival_Probs[6:7]), expected_tail)
            near(as.numeric(fit$Density[6:7]), 0.7 * expected_tail)
        } else {
            stopifnot(all(is.na(fit$Survival_Probs[6:7])))
        }
        if (prior != 2) stopifnot(all(is.na(fit$Density[6:7])))
        stopifnot(fit$Density[5] > 0,
                  fit$Survival_Probs[8] == 0,
                  fit$Density[8] == 0)
        near(sum(fit$Posterior_Probs[fit$Terminal_Nodes]), 1)
        near(tail(as.numeric(fit$LogLik_Trace), 1), fit$LogLik)
    }
}

# With a single uniform body leaf, two finite observations and one censoring
# at tau give the Beta(1,1) posterior mode 2/3.
simple <- rbind(c(1, 1), c(2, 2), c(3, Inf))
for (iterations in c(0, 1, 50)) {
    fit <- fit_tail(simple, depth = 0, iterations = iterations, points = c(1, 3))
    near(body_mass(fit), 2/3)
    near(as.numeric(fit$Density), rep(2/9, 2))
    near(as.numeric(fit$Survival_Probs), c(7/9, 1/3))
}

# Uniform tail predictions work without tau in the grid and preserve its order.
# Check both public interfaces and allocation modes, including missing values.
tail_points <- c(100, NA_real_, 3.01, NaN, Inf, -Inf, 1, 100)
for (allocated in 0:1) {
    fit <- fit_tail(simple, depth = 0, iterations = 3,
                    points = tail_points, allocated = allocated)
    native <- SurvPT(simple, tail_points, 1, c(0, Inf), 0L, c(0.12, 1),
                     allocated, 3L, 1e-10, FALSE)
    near(as.numeric(fit$Survival_Probs[c(1, 3, 8)]), rep(1/3, 3))
    near(as.numeric(native$survival_probs[c(1, 3, 8)]), rep(1/3, 3))
    for (predictions in list(fit$Survival_Probs, native$survival_probs)) {
        stopifnot(all(is.na(predictions[c(2, 4)])))
        near(as.numeric(predictions[c(5, 6, 7)]), c(0, 1, 7/9))
    }
    stopifnot(all(is.na(fit$Density[c(1, 2, 3, 4, 8)])),
              all(is.na(native$density_probs[c(1, 2, 3, 4, 8)])))
}

# Constant extension also handles tail masses close to the probability bounds.
for (data in list(cbind(c(1, 2), c(1, 2)), cbind(c(1, 2), Inf))) {
    fit <- fit_tail(data, iterations = 10, points = c(2, 2.01, 100))
    near(as.numeric(fit$Survival_Probs), rep(tail_mass(fit), 3))
}

# The conditional exponential tail retains its fitted mass and supplied rate.
# Use shifted R distribution functions as an independent reference, including
# large tau where the unconditional exponential survival underflows to zero.
for (tau in c(3, 2000)) for (rate in c(0.7, 2)) for (allocated in 0:1) {
    data <- rbind(c(1, 1), c(2, 2), c(tau, Inf))
    points <- c(tau + 1, NA_real_, tau, Inf, tau + 0.01, NaN, -Inf, tau + 1)
    fit <- survPT(data, test_seq = points, prior = 2, prior_params = rate,
                  max_M = 0, max_EM_iter = 3, allocated = allocated, verbose = FALSE)
    native <- SurvPT(data, points, 2, c(0, Inf), 0L, rate,
                     allocated, 3L, 1e-10, FALSE)
    near(tail_mass(fit), 1/3)
    valid <- c(1, 3, 4, 5, 8)
    expected_s <- 1/3 * pexp(points[valid] - tau, rate, lower.tail = FALSE)
    near(as.numeric(fit$Survival_Probs[valid]), expected_s)
    near(as.numeric(native$survival_probs[valid]), expected_s)
    beyond <- c(1, 4, 5, 8)
    expected_d <- 1/3 * dexp(points[beyond] - tau, rate)
    near(as.numeric(fit$Density[beyond]), expected_d)
    near(as.numeric(native$density_probs[beyond]), expected_d)
    # An exact event at tau still uses the finite body's density on its left.
    near(as.numeric(fit$Density[3]), 2/3 * dexp(tau, rate)/pexp(tau, rate))
    stopifnot(all(is.na(fit$Survival_Probs[c(2, 6)])),
              all(is.na(fit$Density[c(2, 6)])),
              all(is.na(native$survival_probs[c(2, 6)])),
              all(is.na(native$density_probs[c(2, 6)])),
              fit$Survival_Probs[7] == 1, fit$Density[7] == 0)
}

# Boundary mode updates and a default grid built from all finite endpoints.
for (prior in 1:3) {
    exact <- fit_tail(cbind(c(1, 2), c(1, 2)), prior, iterations = 10)
    near(body_mass(exact), 1 - epsilon)
    near(tail_mass(exact), epsilon)
    censored <- fit_tail(cbind(c(1, 4, 2), Inf), prior, iterations = 1000)
    near(body_mass(censored), epsilon)
    near(tail_mass(censored), 1 - epsilon)
    near(range(censored$Test_seq), c(0, 4))
    stopifnot(all(censored$Survival_Probs >= 1 - epsilon - 1e-12),
              all(censored$Survival_Probs <= 1 + 1e-12),
              all(censored$Density >= 0), any(censored$Density > 0),
              is.finite(censored$LogLik), censored$LogLik <= 1e-12,
              censored$LogLik >= 3 * log1p(-epsilon) - 1e-12,
              censored$Converged)
}

# Bounded support has no added tail and keeps first-level alpha = 1.
bounded <- fit_tail(rbind(c(1, 1), c(2, 2), c(2.5, 3)), depth = 1,
                    iterations = 1, support = c(0, 4))
stopifnot(is.na(bounded$tau), max(bounded$Level) == 1)
near(body_mass(bounded), 1)
near(tail_mass(bounded), 0)
near(as.numeric(bounded$Posterior_Probs[2:3]), c(2/3, 1/3))

# The extra root must not shift the finite tree's alpha schedule: first
# body split has alpha = 1, hence conditional masses 2/3 and 1/3.
unbounded <- fit_tail(cbind(c(1, 2, 4), c(1, 2, 4)), depth = 1, iterations = 2)
near(as.numeric(unbounded$Posterior_Probs[4:5]), (1 - epsilon) * c(2/3, 1/3))

# Clipping applies during initialization and EM, to both children of every
# split, for both allocation modes and bounded/unbounded support. Repeated
# exact observations force empty siblings, including first-level body splits.
for (prior in 1:3) for (finite in c(FALSE, TRUE)) for (allocated in 0:1) {
    for (iterations in c(0, 2)) for (time in c(0.25, 3)) {
        data <- matrix(rep(time, 12), ncol = 2)
        native <- SurvPT(data, c(0, 4), prior, c(0, if (finite) 4 else Inf),
                         3L, switch(prior, c(0.12, 1), 0.7, c(0.7, 1.5)),
                         allocated, iterations, 0, FALSE)
        levels <- as.numeric(native$node_levels)
        # Root/body-tail alpha = 1; finite-subtree depths 1, 2, 3 use 1, 2, 4.5.
        relative_depths <- pmax(0, levels - as.integer(!finite))
        expected_alpha <- c(1, 1, 2, 4.5)[relative_depths + 1]
        near(as.numeric(native$alpha), expected_alpha)
        stopifnot(max(expected_alpha) == 4.5)
        masses <- as.numeric(native$posterior_probs)
        near(sum(masses[native$terminal_nodes + 1]), 1)
        splits <- numeric()
        for (j in setdiff(seq_along(masses), native$terminal_nodes + 1)) {
            children <- which(levels == levels[j] + 1 &
                              native$left_bound >= native$left_bound[j] &
                              native$right_bound <= native$right_bound[j])
            stopifnot(length(children) == 2)
            near(sum(masses[children]), masses[j])
            splits <- c(splits, masses[children]/masses[j])
            # Six exact events in one child give conditional modes 7/8 at
            # depth 2 (alpha = 2) and 19/26 at depth 3 (alpha = 4.5).
            populated <- children[native$n_obs_allocated[children] == 6]
            depth <- relative_depths[children[1]]
            if (depth %in% c(2, 3)) {
                stopifnot(length(populated) == 1)
                near(masses[populated]/masses[j], c(7/8, 19/26)[depth - 1])
            }
        }
        stopifnot(all(splits >= epsilon - 1e-12),
                  all(splits <= 1 - epsilon + 1e-12),
                  any(abs(splits - epsilon) < 1e-12),
                  any(abs(splits - (1 - epsilon)) < 1e-12))
    }
}

# A single observation now has a positive root-update denominator. Exact
# membership and censoring at tau must reach opposite bounds, not fall back
# to 1/2. Test both initialization and EM for every base and allocation mode.
for (prior in 1:3) for (allocated in 0:1) for (iterations in c(0, 2)) {
    exact <- fit_tail(matrix(c(2, 2), 1), prior, depth = 0,
                      iterations = iterations, allocated = allocated)
    censored <- fit_tail(matrix(c(2, Inf), 1), prior, depth = 0,
                         iterations = iterations, allocated = allocated)
    near(c(body_mass(exact), tail_mass(exact)), c(1 - epsilon, epsilon), 1e-12)
    near(c(body_mass(censored), tail_mass(censored)), c(epsilon, 1 - epsilon), 1e-12)
}

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
near(body_mass(before), initial_body/4)
expected_body <- 2 + (body_mass(before) * (2/3))/
    (body_mass(before) * (2/3) + tail_mass(before))
near(body_mass(after), expected_body/4)
stopifnot(abs(body_mass(after) - body_mass(before)) > 1e-5)
native <- SurvPT(mixed, c(0, 3), 1, c(0, Inf), 0L, c(0.12, 1),
                 1L, 1L, 1e-10, FALSE)
near(as.numeric(native$alpha[2:3]), c(1, 1))
near(as.numeric(native$posterior_probs[2]), body_mass(after))
stopifnot(length(formals(survPT)) == 10L, length(formals(SurvPT)) == 10L,
          !any(c("q_0", "theta_0", "p_tail", "q0_prior") %in% names(after)),
          !any(c("q_0", "theta_0", "p_tail", "q0_prior") %in% names(native)))

# Known body/tail membership gives a directly checkable Beta posterior mode.
known <- fit_tail(simple, depth = 0, iterations = 3)
mode <- optimize(function(q) dbeta(q, 1+2, 1+1, log = TRUE),
                 c(epsilon, 1 - epsilon), maximum = TRUE, tol = 1e-12)$maximum
near(body_mass(known), mode)
cat("Tail-root Beta initialization, EM updates, likelihood and boundary tests passed.\n")
