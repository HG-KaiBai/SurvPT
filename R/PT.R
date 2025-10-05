#' Fit the Polya tree for Interval censored Data
#'
#' @param X Data frame with three columns: L (left interval), R (right interval), and status (1 if event observed, 0 if censored)
#' @param prior Prior type, "Dirichlet" or "Beta"
#' @param support Support type, "Finite" or "Infinite"
#' @param max_M Maximum level of the Polya tree
#' @param prior_params Prior parameters, a vector of length 2
#' @param allocated Logical, whether to use allocated version
#'
#' @returns A list with the following components: Left (left interval), Right (right interval), Alpha (prior parameters), Allocations (number of observations allocated to each partition), Posterior_Probs (posterior probabilities for each partition)
#' @export
#'
#' @examples
#' data <- data_ic(seed = 123, n = 100, K = 1, Lambda = 0.1, q = 0.3)
#' out <- survPT(X = data)
survPT <- function(X, prior = 1, support = c(0, Inf), max_M = 8, prior_params = NULL, allocated = 1) {
    X = as.matrix(X)
    # The first column and second column of X are the observed left bound and right bound, respectively.
    if (!is.matrix(X) || ncol(X) < 2 || (any(X[,2] < X[,1]))) {
        stop("Invalid input data X.")
    }
    if (! is.numeric(prior) || length(prior) != 1 || prior %in% c(1, 2, 3) == FALSE){
        stop("Invalid input prior.")
    }
    if (prior %in% c(2, 3) && is.null(prior_params)) {
        stop("Please provide prior parameters.")
    }
    if (! is.numeric(support) || length(support) != 2 || support[1] >= support[2]) {
        stop("Invalid input support.")
    }

    # Handle NULL prior_params for uniform prior
    if (prior == 1 && is.null(prior_params)) {
        prior_params = c(0.12, 1)  # Default values, not used for uniform prior
    }

    out = SurvPT(X, prior, support, max_M, prior_params, allocated)
    out_PT = list("Left" = out$left_bound,
                  "Right" = out$right_bound,
                  "Level" = out$node_levels,
                  "Alpha" = out$alpha,
                  "Allocations" = out$n_obs_allocated,
                  "Posterior_Probs" = out$posterior_probs,
                  "Terminal_Nodes" = out$terminal_nodes
                  )
    return(out_PT)
}
