library(SurvPT)

assert_close <- function(actual, expected) {
    stopifnot(max(abs(actual - expected)) < 1e-10)
}
data <- cbind(rbind(c(0, 0.5), c(1, Inf), c(2, Inf)), ignored = 1000)
args <- list(X = data, prior = 1, support = c(0, Inf), max_M = 4,
             max_EM_iter = 5, test_seq = c(0, 1, 2, 3, 6), verbose = FALSE)

# Omitting the new trailing argument preserves the old public/native calls.
stopifnot(identical(do.call(survPT, args),
                    do.call(survPT, c(args, list(uniform_tail_factor = 1.1)))))
native <- list(X = data, test_seq = args$test_seq, prior = 1,
               support = args$support, max_M = 4, prior_params = 1,
               allocated = 1, max_EM_iter = 5, EM_tol = 1e-6, verbose = FALSE)
stopifnot(identical(do.call(SurvPT, native),
                    do.call(SurvPT, c(native, list(uniform_tail_factor = 1.1)))))

for (factor in c(1.01, 1.5, 2, 3)) {
    cutoff <- factor * 2 # The third column is not an observation endpoint.
    points <- c(0, cutoff/2, cutoff, cutoff + 1, Inf)
    root <- survPT(data, prior = 1, max_M = 0, max_EM_iter = 5,
                   test_seq = points, uniform_tail_factor = factor, verbose = FALSE)
    assert_close(as.numeric(root$Survival_Probs), c(1, 0.5, 0, 0, 0))
    assert_close(as.numeric(root$Density), c(rep(1/cutoff, 3), 0, 0))
    stopifnot(is.infinite(root$Right[1]))

    # Check the actual split and initial censoring allocations independently.
    initial <- survPT(data, prior = 1, max_M = 1, max_EM_iter = 0,
                      test_seq = points, uniform_tail_factor = factor, verbose = FALSE)
    split <- cutoff/2
    assert_close(initial$Right[2], split)
    left_counts <- sum(vapply(seq_len(nrow(data)), function(i) {
        left <- data[i, 1]
        right <- min(data[i, 2], cutoff)
        max(0, min(right, split) - left) / (right - left)
    }, numeric(1)))
    assert_close(as.numeric(initial$Allocations[2:3]), c(left_counts, 3 - left_counts))

    wrapped <- do.call(survPT, c(args, list(uniform_tail_factor = factor)))
    raw <- do.call(SurvPT, c(native, list(uniform_tail_factor = factor)))
    stopifnot(identical(wrapped$Survival_Probs, raw$survival_probs),
              identical(wrapped$Allocations, raw$n_obs_allocated))
}

# Finite compatible uniform support and the other bases do not use the factor.
finite_args <- args
finite_args$X <- data
finite_args$X[is.infinite(finite_args$X)] <- 6
finite_args$support <- c(0, 6)
stopifnot(identical(do.call(survPT, finite_args),
                    do.call(survPT, c(finite_args, list(uniform_tail_factor = 2)))))
for (prior in 2:3) {
    other_args <- args
    other_args$prior <- prior
    other_args$prior_params <- if (prior == 2) 0.5 else c(0.5, 1.5)
    stopifnot(identical(do.call(survPT, other_args),
                        do.call(survPT, c(other_args, list(uniform_tail_factor = 2)))))
}

for (bad in list(NULL, numeric(), c(1.1, 2), NA_real_, NaN, Inf, -Inf,
                 0, -1, 1, 0.9, '1.5', TRUE)) {
    error <- tryCatch(do.call(survPT, c(args, list(uniform_tail_factor = bad))),
                      error = identity)
    stopifnot(inherits(error, 'error'),
              grepl('uniform_tail_factor', conditionMessage(error), fixed = TRUE))
}
for (bad in c(NA_real_, NaN, Inf, -Inf, 0, -1, 1, 0.9)) {
    error <- tryCatch(do.call(SurvPT, c(native, list(uniform_tail_factor = bad))),
                      error = identity)
    stopifnot(inherits(error, 'error'),
              grepl('uniform_tail_factor', conditionMessage(error), fixed = TRUE))
}
cat('Uniform tail factor: cutoff, initialization, compatibility and validation passed.\n')
