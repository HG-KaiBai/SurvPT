library(SurvPT)

body_mass <- function(fit) if (is.na(fit$tau)) 1 else as.numeric(fit$Posterior_Probs[2])
tail_mass <- function(fit) if (is.na(fit$tau)) 0 else as.numeric(fit$Posterior_Probs[3])

assert_close <- function(actual, expected, label, tolerance = 1e-10) {
    stopifnot(identical(is.na(as.vector(actual)), is.na(as.vector(expected))))
    error <- max(c(0, abs(actual - expected)), na.rm = TRUE)
    if (!is.finite(error) || error > tolerance) {
        stop(sprintf("%s: maximum error %.12g exceeds %.12g", label, error, tolerance))
    }
    invisible(error)
}

times <- rep(c(0.1, 0.5, 1, 2), each = 4)
observations <- rbind(cbind(times, times), c(0, 0.3), c(0.1, 1.5),
                      c(1, 3), c(0.9, Inf), c(2.6, Inf))

fit_pt <- function(data, prior, support, depth, iterations, points) {
    survPT(data, prior = prior, support = support, max_M = depth,
           prior_params = switch(prior, NULL, 1, c(1, 2)),
           max_EM_iter = iterations, EM_tol = 0, test_seq = points,
           verbose = FALSE)
}

# Independent reference distributions use R's distribution functions.
base_survival <- function(x, prior) {
    if (prior == 2) pexp(x, rate = 1, lower.tail = FALSE)
    else pweibull(x, shape = 2, scale = 1, lower.tail = FALSE)
}
base_density <- function(x, prior) {
    if (prior == 2) dexp(x, rate = 1)
    else dweibull(x, shape = 2, scale = 1)
}
leaf_bounds <- function(fit, j, prior) {
    a <- fit$Left[j]
    b <- fit$Right[j]
    c(a, b)
}
reference_mass <- function(fit, left, right, prior) {
    sum(vapply(fit$Terminal_Nodes, function(j) {
        if (!is.na(fit$tau) && j == 3 && prior != 2)
            return(tail_mass(fit) * as.numeric(is.infinite(right) && left <= fit$tau))
        bounds <- leaf_bounds(fit, j, prior)
        a <- bounds[1]; b <- bounds[2]
        lo <- max(left, a); hi <- min(right, b)
        if (lo >= hi) return(0)
        fraction <- if (prior == 1) (hi-lo)/(b-a) else
            (base_survival(lo, prior)-base_survival(hi, prior)) /
            (base_survival(a, prior)-base_survival(b, prior))
        fit$Posterior_Probs[j] * fraction
    }, numeric(1)))
}
reference_density <- function(fit, t, prior) {
    if (!is.na(fit$tau) && prior != 2 && is.finite(t) && t > fit$tau)
        return(NA_real_)
    sum(vapply(fit$Terminal_Nodes, function(j) {
        if (!is.na(fit$tau) && j == 3 && prior != 2) return(0)
        bounds <- leaf_bounds(fit, j, prior)
        a <- bounds[1]; b <- bounds[2]
        if (t < a || t > b || (t == a && a > fit$Left[1])) return(0)
        density <- if (prior == 1) 1/(b-a) else
            base_density(t, prior)/(base_survival(a, prior)-base_survival(b, prior))
        fit$Posterior_Probs[j] * density
    }, numeric(1)))
}
reference_allocations <- function(fit, data, prior) {
    vapply(seq_along(fit$Left), function(j) {
        a <- fit$Left[j]; b <- fit$Right[j]
        sum(vapply(seq_len(nrow(data)), function(i) {
            l <- data[i, 1]; r <- data[i, 2]
            if (l == r) return(as.numeric(l > a && l <= b))
            numerator <- reference_mass(fit, max(l, a), min(r, b), prior)
            denominator <- reference_mass(fit, l, r, prior)
            if (denominator > 0) numerator/denominator else 0
        }, numeric(1)))
    }, numeric(1))
}

cases <- 0L
for (prior in 1:3) for (finite in c(FALSE, TRUE)) for (depth in c(2, 6)) {
    support <- c(0, if (finite) 4 else Inf)
    data <- observations
    if (finite) data[is.infinite(data[, 2]), 2] <- support[2]
    tau <- max(data[is.finite(data)])
    points <- sort(unique(c(-1, seq(0, 5, length.out = 101), tau, tau+1, 10, Inf)))
    fits <- lapply(c(0, 1, 2, 20), function(iterations)
        fit_pt(data, prior, support, depth, iterations, points))
    label <- sprintf("prior=%d finite=%s depth=%d", prior, finite, depth)

    for (step in 1:2) {
        previous <- fits[[step]]
        current <- fits[[step+1]]
        assert_close(current$Left, previous$Left, paste(label, "fixed partition"))
        expected <- reference_allocations(previous, data, prior)
        assert_close(current$Allocations, expected, paste(label, "EM step", step))
    }

    for (fit in fits) {
        reference_s <- vapply(points, function(t) {
            if (!is.na(fit$tau) && prior != 2 && is.finite(t) && t > fit$tau)
                return(if (prior == 1) tail_mass(fit) else NA_real_)
            reference_mass(fit, t, Inf, prior)
        }, numeric(1))
        reference_d <- vapply(points, function(t)
            reference_density(fit, t, prior), numeric(1))
        assert_close(as.vector(fit$Survival_Probs), reference_s, paste(label, "survival"))
        assert_close(as.vector(fit$Density), reference_d, paste(label, "density"))
        assert_close(sum(fit$Posterior_Probs[fit$Terminal_Nodes]), 1, paste(label, "leaf mass"))
        stopifnot(all(diff(na.omit(as.vector(fit$Survival_Probs))) <= 1e-12),
                  all(fit$Density >= 0, na.rm = TRUE), all(fit$Survival_Probs >= 0, na.rm = TRUE),
                  all(fit$Survival_Probs <= 1+1e-12, na.rm = TRUE))

        # Every internal node must conserve both counts and probability mass.
        for (j in setdiff(seq_along(fit$Left), fit$Terminal_Nodes)) {
            children <- which(fit$Level == fit$Level[j]+1 &
                              fit$Left >= fit$Left[j] & fit$Right <= fit$Right[j])
            stopifnot(length(children) == 2)
            assert_close(sum(fit$Allocations[children]), fit$Allocations[j],
                         paste(label, "child allocations"))
            assert_close(sum(fit$Posterior_Probs[children]), fit$Posterior_Probs[j],
                         paste(label, "child mass"))
        }
    }

    # Integrate actual C++ density output separately on each leaf; these calls
    # repeat the deterministic fit at the quadrature points.
    final <- fits[[4]]
    density_integral <- sum(vapply(final$Terminal_Nodes, function(j) {
        if (!finite && j == 3) return(0)
        bounds <- leaf_bounds(final, j, prior)
        integrate(function(t) as.vector(fit_pt(data, prior, support, depth, 20, t)$Density),
                  bounds[1], bounds[2], rel.tol = 1e-7)$value
    }, numeric(1)))
    assert_close(density_integral, body_mass(final), paste(label, "density integral"), 1e-6)
    if (!finite && prior == 2) {
        tail_integral <- integrate(function(t)
            as.vector(fit_pt(data, prior, support, depth, 20, t)$Density),
            final$tau, Inf, rel.tol = 1e-7)$value
        assert_close(tail_integral, tail_mass(final), paste(label, "tail integral"), 1e-6)
        assert_close(density_integral + tail_integral, 1,
                     paste(label, "total density integral"), 1e-6)
    }

    if (!finite) {
        s_far <- as.vector(final$Survival_Probs)[match(10, points)]
        if (prior == 1) assert_close(s_far, tail_mass(final), "uniform tail survival")
        else if (prior == 2)
            assert_close(s_far, tail_mass(final) * pexp(10 - final$tau, lower.tail = FALSE),
                         "exponential tail survival")
        else stopifnot(is.na(s_far))
    }
    if (depth == 6) stopifnot(any(final$Level[final$Terminal_Nodes] < depth))
    cases <- cases + 1L
    cat(label, "passed\n")
}

# A tree truncated at the root must retain its whole mass through EM.
for (prior in 1:3) {
    data <- observations
    data[is.infinite(data[, 2]), 2] <- 4
    points <- c(0, 0.5, 2, 10, Inf)
    before <- fit_pt(data, prior, c(0, 4), 0, 0, points)
    after <- fit_pt(data, prior, c(0, 4), 0, 5, points)
    assert_close(after$Survival_Probs, before$Survival_Probs, "root-only survival")
    assert_close(after$Density, before$Density, "root-only density")
    assert_close(after$Posterior_Probs, 1, "root-only mass")
}

# Exact observations require no redistribution, regardless of the base shape.
exact <- cbind(times, times)
for (prior in 1:3) {
    before <- fit_pt(exact, prior, c(0, Inf), 3, 0, c(0, 0.5, 2, 10))
    after <- fit_pt(exact, prior, c(0, Inf), 3, 10, c(0, 0.5, 2, 10))
    assert_close(after$Allocations, before$Allocations, "exact observations")
    assert_close(after$Posterior_Probs, before$Posterior_Probs, "exact observation masses")
}
cat(sprintf("All %d mixed-censoring cases, 3 root-only cases and 3 exact-data cases passed.\n", cases))
