#' Fit a finite Polya tree with a Beta-prior tail split
#'
#' @param X Numeric matrix or data frame. The first two columns are left and
#'   right observation endpoints. Use equal endpoints for exact events and
#'   a right endpoint of `Inf` for right censoring. Additional columns are ignored.
#' @param prior Prior type, (1: Uniform, 2: Exponential, 3: Weibull)
#' @param support Numeric vector of lower and upper support bounds, default
#'   `c(0, Inf)`. Left endpoints must be finite and within nonnegative support;
#'   exact events must exceed the lower bound.
#' @param max_M Maximum depth within the finite subtree. The added body/tail
#'   split does not consume this budget; reported levels include it.
#' @param prior_params `NULL` for uniform, a positive rate for exponential,
#'   or positive `c(lambda, k)` for Weibull (shape `k`, scale `1/lambda`).
#' @param allocated `1` for fractional EM allocation; `0` for legacy overlap
#'   initialization without EM.
#' @param max_EM_iter Maximum number of EM iterations
#' @param EM_tol Tolerance for EM convergence
#' @param verbose Logical, whether to print progress messages
#' @param test_seq Evaluation grid; by default, 200 points from the lower
#'   support bound to the largest finite observed endpoint, including censoring times.
#'
#' @details
#' With infinite upper support, `tau` is the largest finite observed endpoint.
#' A new root separates the finite subtree from an unresolved `[tau, Inf)`
#' tail. Its left probability has the internal first-level prior
#' `q_0 ~ Beta(0.5, 0.5)`; `p_tail = 1 - q_0`.
#' Initial censored allocations use the Beta prior mean `0.5` and the
#' selected base conditioned on the finite body. The first layer then uses
#' the same clipped mode-style update as other PT branches, at initialization
#' and in every EM iteration: `q_0 = clip((n_body - 0.5)/(n - 1))`,
#' with a one-half fallback when the denominator is nonpositive or tiny.
#' This is not a posterior-mean update; the fallback need not be an exact MAP
#' for Beta shapes below one. Internal PT alpha values retain their original
#' finite-subtree depths. No separate scalar MLE is performed.
#'
#' If `f` and `P` are the conditional finite-tree density and probability,
#' exact events contribute `q_0 * f(t)`, finite intervals contribute
#' `q_0 * P(L < T <= R)`, and right censoring at `C < tau` contributes
#' `q_0 * P(C < T <= tau) + p_tail`. Right censoring at `tau` contributes
#' `p_tail`. Exact events at split points, including `tau`, use the interval
#' on the left. The first-layer probability is updated along with the tree.
#'
#' Survival at `tau` equals `p_tail`. Tail shape is unspecified, so survival
#' and density at finite points greater than `tau` return `NA`; at positive
#' infinity they return zero. Bounded support adds no tail node. The observed
#' log-likelihood need not increase under the prior-adjusted mode-style updates.
#'
#' @returns A list containing `Left`, `Right`, `Level`, `Allocations`,
#'   `Posterior_Probs`, and one-based `Terminal_Nodes`; terminal masses sum to
#'   one including the tail. `Test_seq`, `Survival_Probs`, and `Density` contain
#'   grid predictions. `tau` is `NA` for bounded support. Body and tail masses
#'   are included in the existing node arrays without separate return fields.
#'   `LogLik` and `LogLik_Trace` report observed
#'   log-likelihood (initialization and each EM step; empty trace without EM).
#'   `EM_Iterations` counts updates; `Converged` indicates whether terminal-mass
#'   tolerance was reached, and is `FALSE` when no updates were run.
#' @export
#'
#' @examples
#' data <- data_ic(seed = 123, n = 100, K = 1, Lambda = 0.1, q = 0.3)
#' out <- survPT(X = data)
survPT <- function(X, test_seq = NULL, prior = 1, support = c(0, Inf), max_M = 8, prior_params = NULL, allocated = 1, max_EM_iter = 100, EM_tol = 1e-6, verbose = TRUE) {
    X = as.matrix(X)
    # The first column and second column of X are the observed left bound and right bound, respectively.
    if (!is.numeric(X) || nrow(X) == 0 || ncol(X) < 2 ||
        anyNA(X[, 1:2]) || any(!is.finite(X[, 1])) || any(X[,2] < X[,1])) {
        stop("Invalid input data X.")
    }
    if (! is.numeric(prior) || length(prior) != 1 || is.na(prior) || prior %in% c(1, 2, 3) == FALSE){
        stop("Invalid input prior.")
    }
    if (prior == 2 && (length(prior_params) != 1 || any(prior_params <= 0))) {
        stop("Please provide a single prior parameter for Exponential prior.")
    }
    if (prior == 3 && (length(prior_params) != 2 || any(prior_params <= 0))) {
        stop("Please provide two positive prior parameters for Weibull prior.")
    }
    if (! is.numeric(support) || length(support) != 2 || anyNA(support) ||
        !is.finite(support[1]) || support[1] < 0 || support[1] >= support[2]) {
        stop("Invalid input support.")
    }

    for (control in list(max_M, max_EM_iter)) {
        if (!is.numeric(control) || length(control) != 1 || !is.finite(control) ||
            control < 0 || control != floor(control) || control > .Machine$integer.max - 1)
            stop("Tree depth and iteration limits must be nonnegative integers.")
    }
    if (is.null(test_seq)){
        endpoints = X[, 1:2]
        test_seq = seq(support[1], max(endpoints[is.finite(endpoints)]), length.out = 200)
    }

    # Handle NULL prior_params for uniform prior
    if (prior == 1 && is.null(prior_params)) {
        prior_params = c(0.12, 1)  # Default values, not used for uniform prior
    }

    out = SurvPT(X, test_seq, prior, support, max_M, prior_params, allocated, max_EM_iter, EM_tol, verbose)
    out_PT = list("Left" = out$left_bound,
                  "Right" = out$right_bound,
                  "Level" = out$node_levels,
                  "Allocations" = out$n_obs_allocated,
                  "Posterior_Probs" = out$posterior_probs,
                  "Terminal_Nodes" = out$terminal_nodes + 1,
                  "Test_seq" = test_seq,
                  "Survival_Probs" = out$survival_probs,
                  "Density" = out$density_probs,
                  "tau" = out$tau,
                  "LogLik" = out$log_likelihood,
                  "LogLik_Trace" = out$loglik_trace,
                  "EM_Iterations" = out$em_iterations,
                  "Converged" = out$converged
                  )
    return(out_PT)
}
