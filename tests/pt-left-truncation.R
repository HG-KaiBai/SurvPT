library(SurvPT)

near <- function(x, y, label, tolerance = 1e-9) {
    error <- max(c(0, abs(x - y)))
    # Ghost sample sizes can be large when survival at entry is near zero.
    scale <- max(1, abs(y))
    if (!is.finite(error) || error > tolerance * scale)
        stop(sprintf("%s: error %.12g", label, error))
}
base_cdf <- function(t, prior) switch(prior, t, pexp(t, 0.7),
                                      pweibull(t, shape = 1.5, scale = 1/0.7))
base_density <- function(t, prior) switch(prior, rep(1, length(t)),
    dexp(t, 0.7), dweibull(t, shape = 1.5, scale = 1/0.7))

# Independent probabilities use R distribution functions, not the C++ cache.
mass <- function(fit, l, r, prior) {
    if (l >= r) return(0)
    sum(vapply(fit$Terminal_Nodes, function(j) {
        if (!is.na(fit$tau) && j == 3)
            return(fit$Posterior_Probs[j] * as.numeric(is.infinite(r) && l <= fit$tau))
        a <- fit$Left[j]; b <- fit$Right[j]
        lo <- max(l, a); hi <- min(r, b)
        if (lo >= hi) return(0)
        fit$Posterior_Probs[j] * (base_cdf(hi, prior) - base_cdf(lo, prior)) /
            (base_cdf(b, prior) - base_cdf(a, prior))
    }, numeric(1)))
}
density <- function(fit, t, prior) {
    sum(vapply(fit$Terminal_Nodes, function(j) {
        if (!is.na(fit$tau) && j == 3) return(0)
        a <- fit$Left[j]; b <- fit$Right[j]
        if (t <= a || t > b) return(0)
        fit$Posterior_Probs[j] * base_density(t, prior) /
            (base_cdf(b, prior) - base_cdf(a, prior))
    }, numeric(1)))
}
reference_counts <- function(bounds, probability, x, entry) {
    vapply(seq_along(bounds$Left), function(j) {
        a <- bounds$Left[j]; b <- bounds$Right[j]
        sum(vapply(seq_len(nrow(x)), function(i) {
            l <- max(x[i, 1], entry[i]); r <- x[i, 2]
            observed <- if (x[i, 1] == r) as.numeric(r > a && r <= b) else
                probability(max(a, l), min(b, r)) / probability(l, r)
            observed + probability(a, min(b, entry[i])) / probability(entry[i], Inf)
        }, numeric(1)))
    }, numeric(1))
}

times <- rep(c(0.4, 1, 2, 4), each = 3)
x <- rbind(cbind(times, times), c(0, 1.5), c(1, 3.5),
           c(1, Inf), c(4, Inf), c(2, 4))
entry <- c(times * rep(c(0, 0.2, 0.7), 4), 0.8, 2, 0.5, 4, 3.8)
cases <- 0L
for (prior in 1:3) for (bounded in c(FALSE, TRUE)) for (depth in c(0, 3, 6)) {
    data <- x
    if (bounded) data[is.infinite(data[, 2]), 2] <- 6
    support <- c(0, if (bounded) 6 else Inf)
    args <- list(X = data, prior = prior,
        prior_params = switch(prior, NULL, 0.7, c(0.7, 1.5)),
        support = support, max_M = depth, EM_tol = 0, verbose = FALSE,
        test_seq = c(0, 0.4, 1, 2, 4, 6, Inf))
    label <- paste(prior, bounded, depth, sep = "/")
    fit_at <- function(iter, a) do.call(survPT, c(args, list(max_EM_iter = iter, entry = a)))

    # Exact, whole-list equality, including iteration traces and convergence.
    for (iter in c(0, 1, 20)) for (allocation in 0:1) {
        old_args <- c(args, list(max_EM_iter = iter, allocated = allocation))
        omitted <- do.call(survPT, old_args)
        stopifnot(identical(omitted, do.call(survPT, c(old_args, list(entry = 0)))),
                  identical(omitted, do.call(survPT, c(old_args, list(entry = rep(0, nrow(data)))))))
    }

    fits <- lapply(c(0, 1, 2, 30), function(iter) fit_at(iter, entry))
    initial <- fits[[1]]
    body_right <- if (bounded) support[2] else initial$tau
    body_weight <- if (bounded) 1 else 0.5
    initial_mass <- function(l, r) {
        if (l >= r) return(0)
        lo <- max(l, support[1]); hi <- min(r, body_right)
        body <- if (lo < hi) body_weight *
            (base_cdf(hi, prior) - base_cdf(lo, prior)) /
            (base_cdf(body_right, prior) - base_cdf(support[1], prior)) else 0
        body + if (!bounded && is.infinite(r) && l <= body_right) 1-body_weight else 0
    }
    near(initial$Allocations, reference_counts(initial, initial_mass, data, entry),
         paste(label, "initial ghost counts"))

    for (step in 1:2) {
        previous <- fits[[step]]; current <- fits[[step + 1]]
        stopifnot(identical(previous$Left, current$Left),
                  identical(previous$Right, current$Right))
        probability <- function(l, r) mass(previous, l, r, prior)
        near(current$Allocations, reference_counts(current, probability, data, entry),
             paste(label, "E-step", step))
        near(current$Allocations[1], sum(vapply(entry, function(a)
            1/probability(a, Inf), numeric(1))), paste(label, "root count"))

        # Confirm the unchanged alpha + expected-count - 1 branch update.
        for (j in setdiff(seq_along(current$Left), current$Terminal_Nodes)) {
            children <- which(current$Level == current$Level[j] + 1 &
                current$Left >= current$Left[j] & current$Right <= current$Right[j])
            branch_depth <- if (!bounded && j != 1) current$Level[children[1]]-1 else
                current$Level[children[1]]
            alpha <- 0.5 * branch_depth^2
            total <- 2*alpha + sum(current$Allocations[children])-2
            p <- if (total > 1e-300) max(0, min(1,
                (alpha + current$Allocations[children[1]]-1)/total)) else 0.5
            near(current$Posterior_Probs[children], current$Posterior_Probs[j]*c(p, 1-p),
                 paste(label, "M-step"))
        }
    }
    for (fit in fits) {
        contributions <- vapply(seq_len(nrow(data)), function(i) {
            numerator <- if (data[i, 1] == data[i, 2]) density(fit, data[i, 1], prior) else
                mass(fit, max(data[i, 1], entry[i]), data[i, 2], prior)
            log(numerator) - log(mass(fit, entry[i], Inf, prior))
        }, numeric(1))
        near(fit$LogLik, sum(contributions), paste(label, "conditional likelihood"))
        near(tail(as.numeric(fit$LogLik_Trace), 1), fit$LogLik, "likelihood trace")
        near(sum(fit$Posterior_Probs[fit$Terminal_Nodes]), 1, "leaf mass")
        for (j in setdiff(seq_along(fit$Left), fit$Terminal_Nodes)) {
            children <- which(fit$Level == fit$Level[j]+1 &
                fit$Left >= fit$Left[j] & fit$Right <= fit$Right[j])
            near(sum(fit$Allocations[children]), fit$Allocations[j], "count conservation")
        }
    }
    clipped <- data
    clipped[, 1] <- pmax(clipped[, 1], entry)
    clipped_args <- args; clipped_args$X <- clipped
    stopifnot(identical(fits[[3]], do.call(survPT,
        c(clipped_args, list(max_EM_iter = 2, entry = entry)))))
    cases <- cases + 1L
}

# Scalar entry, lower support, exact boundary rejection, and direct native API.
small <- rbind(c(1, 1), c(2, 3), c(3, Inf))
stopifnot(identical(survPT(small, entry = 0.5, verbose = FALSE),
                    survPT(small, entry = rep(0.5, 3), verbose = FALSE)))
bounded_data <- rbind(c(2, 2), c(3, 4))
stopifnot(identical(survPT(bounded_data, support = c(1, 5), verbose = FALSE),
                    survPT(bounded_data, support = c(1, 5), entry = 1, verbose = FALSE)))
native_args <- list(X = small, test_seq = c(0, 1, 3), prior = 1L,
    support = c(0, Inf), max_M = 3L, prior_params = c(0.12, 1),
    allocated = 1L, max_EM_iter = 2L, EM_tol = 0, verbose = FALSE)
native <- do.call(SurvPT, native_args)
stopifnot(identical(native, do.call(SurvPT, c(native_args, list(entry = 0)))),
          identical(native, do.call(SurvPT, c(native_args, list(entry = rep(0, 3))))))
native_lt <- do.call(SurvPT, c(native_args, list(entry = 0.5)))
public_lt <- do.call(survPT, c(native_args, list(entry = 0.5)))
near(native_lt$n_obs_allocated, public_lt$Allocations, "native counts")
near(native_lt$log_likelihood, public_lt$LogLik, "native likelihood")
bad <- list(-1, NA_real_, Inf, numeric(0), c(0, 0), c(0, 0, 0, 0), c(1, 0, 0),
            c(0, 3, 0))
for (a in bad) {
    stopifnot(inherits(try(survPT(small, entry = a, verbose = FALSE), silent = TRUE), "try-error"),
              inherits(try(do.call(SurvPT, c(native_args, list(entry = a))), silent = TRUE), "try-error"))
}
stopifnot(inherits(try(survPT(small, entry = "0", verbose = FALSE), silent = TRUE), "try-error"),
          inherits(try(survPT(small, entry = matrix(0, 3), verbose = FALSE), silent = TRUE), "try-error"),
          inherits(try(survPT(small, entry = 0.5, allocated = 0, verbose = FALSE), silent = TRUE), "try-error"))
native_args$allocated <- 0L
stopifnot(inherits(try(do.call(SurvPT, c(native_args, list(entry = 0.5))), silent = TRUE), "try-error"))

# Entry can refine the finite endpoint when an interval begins before entry.
late <- survPT(rbind(c(4, 4), c(3, Inf)), entry = c(0, 5), max_M = 0,
               max_EM_iter = 2, verbose = FALSE)
near(late$tau, 5, "entry defines tau")
near(max(late$Test_seq), 5, "entry defines default grid")
stopifnot(is.finite(late$LogLik))
cat(sprintf("All %d truncation cases, zero-entry comparisons and validation checks passed.\n", cases))
