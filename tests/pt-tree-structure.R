library(SurvPT)

assert_close <- function(actual, expected, label, tolerance = 1e-10) {
    error <- max(abs(actual - expected))
    if (!is.finite(error) || error > tolerance)
        stop(sprintf("%s: maximum error %.12g", label, error))
}

native_fit <- function(data, prior, depth, allocated = 1, iterations = 0,
                       params = switch(prior, c(0, 1), 0.7, c(0.7, 1.5)),
                       points = c(0, 4, 5, Inf), support = c(0, Inf)) {
    SurvPT(data, points, prior, support, depth, params, allocated,
           iterations, 0, FALSE)
}

truncated_quantile <- function(p, prior, cap = 4) {
    switch(prior, cap*p,
           qexp(p*pexp(cap, rate = 0.7), rate = 0.7),
           qweibull(p*pweibull(cap, shape = 1.5, scale = 1/0.7),
                    shape = 1.5, scale = 1/0.7))
}

# Populate every ordinary leaf, so the expected partition is independent of
# count pruning. The extra column must not enter largest_obs.
cases <- 0L
for (prior in 1:3) for (depth in c(0, 1, 3)) for (allocated in 0:1) {
    cells <- 2^depth
    times <- rep(truncated_quantile((seq_len(cells)-0.5)/cells, prior), each = 4)
    data <- cbind(rbind(cbind(times, times), c(4, 4), c(4, Inf)), ignored = 999)
    for (iterations in c(0, 3)) {
        fit <- native_fit(data, prior, depth, allocated, iterations)
        leaves <- as.integer(fit$terminal_nodes)+1L
        tail <- leaves[is.infinite(fit$right_bound[leaves])]
        body <- setdiff(leaves, tail)
        stopifnot(length(fit$left_bound) == 2*cells+1,
                  fit$left_bound[1] == 0, fit$right_bound[1] == Inf,
                  length(tail) == 1, fit$left_bound[tail] == 4,
                  fit$node_levels[tail] == depth+1,
                  fit$node_levels[tail-1] == depth+1,
                  fit$right_bound[tail-1] == 4,
                  fit$n_obs_allocated[tail] == 1)
        expected <- c(truncated_quantile((0:(cells-1))/cells, prior), 4)
        assert_close(sort(fit$left_bound[leaves]), expected, "truncated quantile boundaries")
        assert_close(fit$alpha[-1], 0.5*fit$node_levels[-1]^2, "depth prior parameters")
        assert_close(fit$alpha[c(tail-1, tail)], rep(0.5*(depth+1)^2, 2), "tail Beta parameters")
        assert_close(sum(fit$posterior_probs[leaves]), 1, "total leaf mass")

        # Check the reported masses use the extra layer's prior in the same
        # mode update as every other branch, at initialization and after EM.
        for (node in setdiff(seq_along(fit$left_bound), leaves)) {
            children <- which(fit$node_levels == fit$node_levels[node]+1 &
                              fit$left_bound >= fit$left_bound[node] &
                              fit$right_bound <= fit$right_bound[node])
            stopifnot(length(children) == 2)
            a <- fit$alpha[children]
            counts <- fit$n_obs_allocated[children]
            denominator <- sum(a+counts)-2
            q <- if (denominator > 1e-300)
                min(max((a[1]+counts[1]-1)/denominator, 0), 1) else 0.5
            assert_close(fit$posterior_probs[children],
                         fit$posterior_probs[node]*c(q, 1-q), "branch mode update")
            assert_close(sum(counts), fit$n_obs_allocated[node], "branch count conservation")
        }
        if (prior == 1) {
            assert_close(as.vector(fit$survival_probs),
                         c(1, rep(fit$posterior_probs[tail], 2), 0), "infinite atom survival")
            assert_close(fit$density_probs[3:4], c(0, 0), "infinite atom density")
        }
        cases <- cases+1L
    }
}

# Uniform right-censored observations start wholly in the tail, regardless of
# censoring time, tree depth or allocation mode. Check every node against an
# independent observation-by-observation count, including finite intervals.
initial_data <- list(
    all_right = cbind(L = c(0, 1, 3, 8), R = Inf),
    mixed = rbind(cbind(c(0, 1, 2, 4, 8), c(0, 1, 2, 4, 8)),
                  c(0, 8), c(1, 5), c(6, 8),
                  cbind(L = c(0, 2, 7, 8), R = Inf))
)
initial_cases <- 0L
for (data in initial_data) for (depth in c(0, 2, 6)) for (allocated in 0:1) {
    fit <- native_fit(data, 1, depth, allocated = allocated)
    leaves <- as.integer(fit$terminal_nodes)+1L
    tail_node <- leaves[is.infinite(fit$right_bound[leaves])]
    expected <- vapply(seq_along(fit$left_bound), function(j) {
        a <- fit$left_bound[j]; b <- fit$right_bound[j]
        sum(vapply(seq_len(nrow(data)), function(i) {
            l <- data[i, 1]; r <- data[i, 2]
            if (r == Inf) return(as.numeric(b == Inf))
            if (l == r) return(as.numeric((l > a || (a == 0 && l == 0)) && l <= b))
            overlap <- max(min(r, b)-max(l, a), 0)
            if (allocated == 1) overlap/(r-l) else as.numeric(overlap > 0)
        }, numeric(1)))
    }, numeric(1))
    assert_close(fit$n_obs_allocated, expected, "uniform initial node allocations")
    assert_close(fit$n_obs_allocated[tail_node], sum(data[, 2] == Inf),
                 "all right-censored observations start in tail")
    if (all(data[, 2] == Inf)) {
        stopifnot(all(fit$n_obs_allocated[is.finite(fit$right_bound)] == 0),
                  length(fit$left_bound) == 2*depth+3,
                  fit$node_levels[tail_node] == depth+1)
    }
    initial_cases <- initial_cases+1L
}

# Zero initial weights must not remove finite nodes from subsequent EM
# allocation. Use exact finite events so any increase comes from censoring.
event_times <- rep(c(0, 2, 4, 8), each = 3)
em_data <- rbind(cbind(event_times, event_times), cbind(L = c(0, 2, 7, 8), R = Inf))
initial <- native_fit(em_data, 1, depth = 3)
updated <- native_fit(em_data, 1, depth = 3, iterations = 1)
leaves <- as.integer(initial$terminal_nodes)+1L
tail_node <- leaves[is.infinite(initial$right_bound[leaves])]
finite_leaves <- setdiff(leaves, tail_node)
stopifnot(identical(initial$left_bound, updated$left_bound),
          identical(initial$right_bound, updated$right_bound),
          updated$n_obs_allocated[tail_node] < initial$n_obs_allocated[tail_node],
          sum(updated$n_obs_allocated[finite_leaves]) > length(event_times))
assert_close(sum(updated$n_obs_allocated[leaves]), nrow(em_data), "EM count conservation")

# Exponential and Weibull initial allocations retain their conditional base
# probabilities, including right-censored observations below the cap.
for (prior in 2:3) for (depth in c(0, 3)) {
    data <- initial_data$mixed
    fit <- native_fit(data, prior, depth)
    leaves <- as.integer(fit$terminal_nodes)+1L
    tail_node <- leaves[is.infinite(fit$right_bound[leaves])]
    censor_times <- data[data[, 2] == Inf, 1]
    shape <- if (prior == 2) 1 else 1.5
    expected <- sum(exp((0.7*censor_times)^shape-(0.7*8)^shape))
    assert_close(fit$n_obs_allocated[tail_node], expected, "nonuniform initial tail counts")
}

# The final split is mandatory even if the rightmost node has fewer than one
# observation and stops before max_M. Its actual depth determines alpha.
for (prior in 1:3) {
    fit <- native_fit(matrix(c(0, 1), 1, 2), prior, depth = 8)
    leaves <- as.integer(fit$terminal_nodes)+1L
    tail <- leaves[is.infinite(fit$right_bound[leaves])]
    stopifnot(length(fit$left_bound) == 5, fit$left_bound[tail] == 1,
              fit$node_levels[tail] == 2, fit$alpha[tail] == 2,
              fit$alpha[tail-1] == 2)
}

# Small CDF values must still give interior truncated medians.
for (prior in 2:3) {
    params <- if (prior == 2) 1e-20 else c(1e-20, 1)
    fit <- native_fit(cbind(c(0.25, 0.75, 1.25, 2), c(0.25, 0.75, 1.25, 2)),
                      prior, depth = 2, params = params)
    assert_close(fit$right_bound[2], 1, "small-CDF root median")
    stopifnot(all(is.finite(fit$left_bound)), all(is.finite(fit$survival_probs)))
}

# Explicitly creating a tail at the largest endpoint must not lose its mass
# when the unconditioned base survival there underflows to zero.
for (prior in 2:3) for (iterations in c(0, 3)) {
    cap <- if (prior == 2) 1000 else 40
    params <- if (prior == 2) 1 else c(1, 2)
    times <- c(0.1, 0.2, 0.4, cap, cap)
    data <- rbind(cbind(times, times), matrix(rep(c(cap, Inf), 3), ncol = 2, byrow = TRUE))
    points <- c(cap, cap+0.01, Inf)
    fit <- native_fit(data, prior, depth = 3, iterations = iterations,
                      params = params, points = points)
    tail <- tail(as.integer(fit$terminal_nodes)+1L, 1)
    mass <- fit$posterior_probs[tail]
    delta <- if (prior == 2) 0.01 else (cap+0.01)^2-cap^2
    expected_density <- if (prior == 2) mass*exp(-delta) else
        mass*2*(cap+0.01)*exp(-delta)
    stopifnot(mass > 0, fit$n_obs_allocated[tail] == 3)
    assert_close(fit$survival_probs, c(mass, mass*exp(-delta), 0), "underflow-safe survival")
    assert_close(fit$density_probs[2], expected_density, "underflow-safe density")
    assert_close(sum(fit$n_obs_allocated[as.integer(fit$terminal_nodes)+1L]),
                 nrow(data), "underflow-safe allocations")
}

# Both interfaces require the common root support and a usable observed range.
for (support in list(c(0, 4), c(0.5, Inf))) {
    data <- cbind(1:2, 1:2)
    wrapper_error <- tryCatch(survPT(data, support = support, verbose = FALSE), error = identity)
    native_error <- tryCatch(native_fit(data, 1, 1, support = support), error = identity)
    stopifnot(inherits(wrapper_error, "error"), inherits(native_error, "error"),
              grepl("support must be c(0, Inf)", conditionMessage(wrapper_error), fixed = TRUE),
              grepl("support must be c(0, Inf)", conditionMessage(native_error), fixed = TRUE))
}
zero_error <- tryCatch(survPT(matrix(0, 1, 2), test_seq = 0, verbose = FALSE), error = identity)
stopifnot(inherits(zero_error, "error"),
          grepl("positive finite observation endpoint", conditionMessage(zero_error), fixed = TRUE))
cat(sprintf("All %d tree configurations, %d uniform initialization cases, EM redistribution, early-stop, numerical-tail and input checks passed.\n",
            cases, initial_cases))
