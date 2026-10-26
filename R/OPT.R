#' @export
#' @noRd
survOPT <- function(X, delta, test_seq = NULL, Rho = 0.5, support = c(0, Inf), max_M = 8, max_iter = 100, iter_tol = 1e-6, verbose = TRUE) {
    # One observed follow-up time and one event indicator per subject.
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
                  "Survival_Probs" = out$survival_probs,
                  "Density" = out$density_probs
                  )
    return(out_PT)
}
