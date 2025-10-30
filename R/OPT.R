#' Fit the Polya tree for Interval censored Data
#'
#' @param X vector of observed data
#' @param support Support type, "Finite" or "Infinite"
#' @param max_M Maximum level of the Polya tree
#' @param verbose Logical, whether to print progress messages
#' @param test_seq Sequence of test points for survival probability estimation (optional)
#' @param delta censoring indicator (1: event time, 0: right censored)
#' @param Rho prior stop probability
#' @param max_iter max number of iteration
#' @param iter_tol tolerance of iteration
#'
#' @returns A list with the following components: Left (left interval), Right (right interval), Alpha (prior parameters), Allocations (number of observations allocated to each partition), Posterior_Probs (posterior probabilities for each partition)
#' @export
#'
#' @examples
#' data <- data_ic(seed = 123, n = 100, K = 1, Lambda = 0.1, q = 0.3)
#' out <- survPT(X = data)
survOPT <- function(X, delta, test_seq = NULL, Rho = 0.5, support = c(0, Inf), max_M = 8, max_iter = 100, iter_tol = 1e-6, verbose = TRUE) {
    # The first column and second column of X are the observed left bound and right bound, respectively.
    X = as.numeric(X)
    delta = as.numeric(delta)
    if (!is.numeric(X) || length(X) < 1) {
        stop("Invalid input data X.")
    }
    if (!is.numeric(delta) || length(delta) != length(X) || any(!(delta %in% c(0,1)))) {
        stop("Invalid input delta.")
    }
    if (! is.numeric(support) || length(support) != 2 || support[1] >= support[2]) {
        stop("Invalid input support.")
    }
    if (! is.numeric(Rho) || length(Rho) != 1 || Rho <= 0 || Rho >=1){
        stop("Invalid input Rho.")
    }
    if (is.null(test_seq)){
        test_seq = seq(support[1], max(X[is.finite(X)]), length.out = 200)
    }
    out = OPT_right(X, delta, test_seq, Rho, support, max_M, max_iter, iter_tol, verbose)
    out_PT = list("Left" = out$left_bound,
                  "Right" = out$right_bound,
                  "Level" = out$node_levels,
                  "Rho" = out$rho,
                  "Allocations" = out$n_obs_allocated,
                  "Posterior_Probs" = out$posterior_probs,
                  "Terminal_Nodes" = out$terminal_nodes + 1,
                  "Test_seq" = test_seq,
                  "Survival_Probs" = out$survival_probs
                  )
    return(out_PT)
}
