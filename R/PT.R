#' Fit the Polya tree for Interval censored Data
#'
#' @param X Data frame with three columns: L (left interval), R (right interval), and status (1 if event observed, 0 if censored)
#' @param prior Prior type, (1: Uniform, 2: Exponential, 3: Weibull)
#' @param support Support type, "Finite" or "Infinite"
#' @param max_M Maximum level of the Polya tree
#' @param prior_params Prior parameters, a vector of length 2
#' @param allocated Logical, whether to use allocated version
#' @param max_EM_iter Maximum number of EM iterations
#' @param EM_tol Tolerance for EM convergence
#' @param verbose Logical, whether to print progress messages
#' @param test_seq Sequence of test points for survival probability estimation (optional)
#'
#' @returns A list with the following components: Left (left interval), Right (right interval), Alpha (prior parameters), Allocations (number of observations allocated to each partition), Posterior_Probs (posterior probabilities for each partition)
#' @export
#'
#' @examples
#' data <- data_ic(seed = 123, n = 100, K = 1, Lambda = 0.1, q = 0.3)
#' out <- survPT(X = data)
survPT <- function(X, test_seq = NULL, prior = 1, support = c(0, Inf), max_M = 8, prior_params = NULL, allocated = 1, max_EM_iter = 100, EM_tol = 1e-6, verbose = TRUE) {
    X = as.matrix(X)
    # The first column and second column of X are the observed left bound and right bound, respectively.
    if (!is.matrix(X) || ncol(X) < 2 || (any(X[,2] < X[,1]))) {
        stop("Invalid input data X.")
    }
    if (! is.numeric(prior) || length(prior) != 1 || prior %in% c(1, 2, 3) == FALSE){
        stop("Invalid input prior.")
    }
    if (prior == 2 && (length(prior_params) != 1 || any(prior_params <= 0))) {
        stop("Please provide a single prior parameter for Exponential prior.")
    }
    if (prior == 3 && (length(prior_params) != 2 || any(prior_params <= 0))) {
        stop("Please provide two positive prior parameters for Weibull prior.")
    }
    if (! is.numeric(support) || length(support) != 2 || support[1] >= support[2]) {
        stop("Invalid input support.")
    }

    if (is.null(test_seq)){
        test_seq = seq(support[1], max(X[is.finite(X[,2]),1:2]), length.out = 200)
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
                  "Density" = out$density_probs
                  )
    return(out_PT)
}
