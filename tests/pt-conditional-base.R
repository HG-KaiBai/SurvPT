library(SurvPT)

assert_close <- function(actual, expected, label, tolerance = 1e-10) {
    error <- max(abs(actual - expected))
    if (!is.finite(error) || error > tolerance) {
        stop(sprintf("%s: maximum error %.12g exceeds %.12g", label, error, tolerance))
    }
    invisible(error)
}

times <- rep(c(0, 0.1, 0.5, 1, 2), each = 4)
observations <- rbind(cbind(times, times), c(0, 0.3), c(0.1, 1.5),
                      c(1, 3), c(0.9, Inf), c(2.6, Inf))

fit_pt <- function(data, prior, support, depth, iterations, points, uniform_tail_factor = 1.1) {
    survPT(data, prior = prior, support = support, max_M = depth,
           prior_params = switch(prior, NULL, 1, c(1, 2)),
           max_EM_iter = iterations, EM_tol = 0, test_seq = points,
           verbose = FALSE, uniform_tail_factor = uniform_tail_factor)
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
leaf_bounds <- function(fit, j, prior, cap) {
    a <- fit$Left[j]
    b <- fit$Right[j]
    if (prior == 1 && is.infinite(b)) b <- cap
    c(a, b)
}
reference_mass <- function(fit, left, right, prior, cap) {
    sum(vapply(fit$Terminal_Nodes, function(j) {
        bounds <- leaf_bounds(fit, j, prior, cap)
        a <- bounds[1]; b <- bounds[2]
        lo <- max(left, a); hi <- min(right, b)
        if (lo >= hi) return(0)
        fraction <- if (prior == 1) (hi-lo)/(b-a) else
            (base_survival(lo, prior)-base_survival(hi, prior)) /
            (base_survival(a, prior)-base_survival(b, prior))
        fit$Posterior_Probs[j] * fraction
    }, numeric(1)))
}
reference_density <- function(fit, t, prior, cap) {
    sum(vapply(fit$Terminal_Nodes, function(j) {
        bounds <- leaf_bounds(fit, j, prior, cap)
        a <- bounds[1]; b <- bounds[2]
        if (!is.finite(t) || t < a || t > b || (t == a && a != 0)) return(0)
        density <- if (prior == 1) 1/(b-a) else
            base_density(t, prior)/(base_survival(a, prior)-base_survival(b, prior))
        fit$Posterior_Probs[j] * density
    }, numeric(1)))
}
reference_allocations <- function(fit, data, prior, cap) {
    vapply(seq_along(fit$Left), function(j) {
        a <- fit$Left[j]; b <- fit$Right[j]
        sum(vapply(seq_len(nrow(data)), function(i) {
            l <- data[i, 1]; r <- data[i, 2]
            if (l == r) return(as.numeric((l > a || (a == 0 && l == 0)) && l <= b))
            numerator <- reference_mass(fit, max(l, a), min(r, b), prior, cap)
            denominator <- reference_mass(fit, l, r, prior, cap)
            if (denominator > 0) numerator/denominator else 0
        }, numeric(1)))
    }, numeric(1))
}

cases <- 0L
for (prior in 1:3) for (finite in c(FALSE, TRUE)) for (depth in c(2, 6))
for (factor in if (prior == 1) c(1.1, 1.5, 2) else 1.1) {
    support <- c(0, if (finite) 4 else Inf)
    data <- observations
    if (finite) data[is.infinite(data[, 2]), 2] <- support[2]
    cap <- factor * max(data[is.finite(data)])
    points <- sort(unique(c(-1, seq(0, 5, length.out = 101), cap, cap+1, 10, Inf)))
    fits <- lapply(c(0, 1, 2, 20), function(iterations)
        fit_pt(data, prior, support, depth, iterations, points, factor))
    label <- sprintf("prior=%d finite=%s depth=%d factor=%g", prior, finite, depth, factor)

    for (step in 1:2) {
        previous <- fits[[step]]
        current <- fits[[step+1]]
        assert_close(current$Left, previous$Left, paste(label, "fixed partition"))
        expected <- reference_allocations(previous, data, prior, cap)
        assert_close(current$Allocations, expected, paste(label, "EM step", step))
    }

    for (fit in fits) {
        reference_s <- vapply(points, function(t)
            reference_mass(fit, t, Inf, prior, cap), numeric(1))
        reference_d <- vapply(points, function(t)
            reference_density(fit, t, prior, cap), numeric(1))
        assert_close(as.vector(fit$Survival_Probs), reference_s, paste(label, "survival"))
        assert_close(as.vector(fit$Density), reference_d, paste(label, "density"))
        assert_close(sum(fit$Posterior_Probs[fit$Terminal_Nodes]), 1, paste(label, "leaf mass"))
        stopifnot(all(is.finite(fit$Survival_Probs)), all(is.finite(fit$Density)),
                  all(diff(as.vector(fit$Survival_Probs)) <= 1e-12),
                  all(fit$Density >= 0), all(fit$Survival_Probs >= 0),
                  all(fit$Survival_Probs <= 1+1e-12))

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
        bounds <- leaf_bounds(final, j, prior, cap)
        integrate(function(t) as.vector(fit_pt(data, prior, support, depth, 20, t, factor)$Density),
                  bounds[1], bounds[2], rel.tol = 1e-7)$value
    }, numeric(1)))
    assert_close(density_integral, 1, paste(label, "density integral"), 1e-6)

    if (!finite) {
        s_far <- as.vector(final$Survival_Probs)[match(10, points)]
        if (prior == 1) stopifnot(s_far == 0)
        else stopifnot(s_far > 0)
    }
    if (depth == 6) stopifnot(any(final$Level[final$Terminal_Nodes] < depth))
    cases <- cases + 1L
    cat(label, "passed\n")
}

# A tree truncated at the root must retain its whole mass through EM.
for (prior in 1:3) {
    data <- observations
    points <- c(0, 0.5, 2, 10, Inf)
    before <- fit_pt(data, prior, c(0, Inf), 0, 0, points)
    after <- fit_pt(data, prior, c(0, Inf), 0, 5, points)
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

# Zero belongs to the first leaf. Positive split points and the finite upper
# endpoint belong to the leaf on their left, for allocation and density alike.
boundary_cases <- 0L
for (prior in 1:3) {
    params <- switch(prior, NULL, 1, c(1, 2))
    probe <- survPT(rbind(c(0, 0), c(2, 2)), prior = prior, prior_params = params,
                    support = c(0, 2), max_M = 1, max_EM_iter = 0,
                    test_seq = c(0, 2), verbose = FALSE)
    split <- as.numeric(probe$Right[2])
    exact_times <- c(0, 0, split, split, 2)
    data <- cbind(exact_times, exact_times)
    points <- c(0, split, 2, 2 + 1e-8, Inf)

    for (allocated in 0:1) for (depth in c(1, 3)) for (iterations in c(0, 5)) {
        fit <- survPT(data, prior = prior, prior_params = params,
                      support = c(0, 2), max_M = depth, allocated = allocated,
                      max_EM_iter = iterations, test_seq = points, verbose = FALSE)
        assert_close(as.vector(fit$Allocations[1:3]), c(5, 4, 1),
                     "zero and split-point allocations")
        assert_close(sum(fit$Allocations[fit$Terminal_Nodes]), 5,
                     "exact-event leaf count conservation")
        leaves <- fit$Terminal_Nodes
        boundary_leaves <- c(leaves[fit$Left[leaves] == 0],
                             leaves[fit$Right[leaves] == split],
                             leaves[fit$Right[leaves] == 2])
        stopifnot(length(boundary_leaves) == 3)
        expected <- vapply(seq_along(boundary_leaves), function(i) {
            j <- boundary_leaves[i]
            a <- fit$Left[j]; b <- fit$Right[j]
            conditional_density <- if (prior == 1) 1/(b-a) else
                base_density(points[i], prior)/(base_survival(a, prior)-base_survival(b, prior))
            fit$Posterior_Probs[j] * conditional_density
        }, numeric(1))
        assert_close(as.vector(fit$Density), c(expected, 0, 0),
                     "zero, split-point and upper-endpoint densities")
        stopifnot(fit$Density[2] > 0, fit$Density[3] > 0)
        boundary_cases <- boundary_cases + 1L
    }

    # The left-endpoint exception applies only at zero, not at a positive
    # lower support bound, including when the entire tree is one leaf.
    positive_lower <- survPT(rbind(c(1, 1), c(2, 2)), prior = prior,
                             prior_params = params, support = c(0.5, 2),
                             max_M = 0, test_seq = c(0.5, 2, Inf), verbose = FALSE)
    stopifnot(positive_lower$Density[1] == 0, positive_lower$Density[2] > 0,
              positive_lower$Density[3] == 0)
}
cat(sprintf("All %d endpoint cases and 3 positive-lower-bound cases passed.\n", boundary_cases))

# Default grids must include finite censoring times, ignore extra columns,
# and agree with predictions requested explicitly on the same grid.
grid_cases <- list(
    all_right = cbind(L = c(1, 2, 7), R = Inf, ignored = 999),
    mixed = cbind(rbind(c(1, 1), c(2, 2), c(7, Inf)), ignored = 999),
    finite_right = cbind(rbind(c(1, 5), c(2, Inf)), ignored = 999)
)
for (prior in 1:3) for (allocated in 0:1) for (data in grid_cases) {
    endpoints <- data[, 1:2]
    expected_grid <- seq(0, max(endpoints[is.finite(endpoints)]), length.out = 200)
    args <- list(X = data, prior = prior, prior_params = switch(prior, NULL, 1, c(1, 2)),
                 allocated = allocated, max_M = 3, max_EM_iter = 5, verbose = FALSE)
    automatic <- do.call(survPT, args)
    explicit <- do.call(survPT, c(args, list(test_seq = expected_grid)))
    assert_close(automatic$Test_seq, expected_grid, "default grid finite maximum")
    assert_close(automatic$Survival_Probs, explicit$Survival_Probs, "default grid survival")
    assert_close(automatic$Density, explicit$Density, "default grid density")
}
custom_grid <- c(0, 0.25, 10)
custom <- survPT(grid_cases$all_right, test_seq = custom_grid, max_M = 2, verbose = FALSE)
stopifnot(identical(custom$Test_seq, custom_grid))
invalid_grid <- tryCatch(survPT(matrix(c(Inf, Inf), 1, 2), verbose = FALSE),
                         error = identity)
stopifnot(inherits(invalid_grid, "error"),
          grepl("without finite observation endpoints", conditionMessage(invalid_grid), fixed = TRUE))
cat("All 18 default-grid cases, custom-grid and missing-finite-endpoint checks passed.\n")
