#' Fit an optional Polya tree to univariate right-censored data
#'
#' Estimate a univariate survival function and density using an optional Polya
#' tree with uniform finite terminal intervals. The censored-data OPT
#' methodology is described by Seok, Tian, and Wong (2014).
#'
#' @param X Numeric vector of finite observed event or censoring times, one per
#'   subject. This is not a matrix of interval endpoints.
#' @param delta Numeric or logical vector of the same length as `X`: `1` for an
#'   exact event and `0` for right censoring. A current-status indicator is not
#'   a right-censoring indicator.
#' @param test_seq Numeric evaluation points. If `NULL`, uses 200 points from
#'   `support[1]` to the largest finite observed time. A sorted finite grid is
#'   convenient for plotting and comparison.
#' @param Rho Prior stopping probability, strictly between zero and one.
#'   This is a stopping probability, not a splitting probability.
#' @param support Numeric vector `c(lower, upper)`, default `c(0, Inf)`.
#'   Use nonnegative support compatible with the observations. A finite upper
#'   bound must be justified by the event-time model.
#' @param max_M Nonnegative integer maximum depth; the root is level zero.
#'   Initial allocated counts below one also stop splitting. Use modest depths.
#' @param max_iter Nonnegative integer maximum number of allocation updates.
#'   Zero returns the initialized fit.
#' @param iter_tol Nonnegative tolerance for the Euclidean norm of the change
#'   in terminal masses. No convergence flag is returned.
#' @param verbose Logical; print iteration progress if `TRUE`.
#'
#' @details
#' The cited paper studies multivariate censored data. This interface implements
#' a univariate right-censoring procedure with weighted Kaplan-Meier
#' initialization and a binary partition; it does not expose the complete
#' multivariate method or claim exact reproduction of the paper's experiments.
#'
#' Finite terminal intervals use uniform conditional distributions. For an
#' infinite terminal interval, output retains the entire terminal mass in
#' survival at every finite time beyond its left endpoint and adds zero finite
#' density. A positive survival plateau can therefore occur. Some internal
#' width calculations instead use `1.1 * max(X)`; these tail conventions have
#' not been harmonized. The PT conditional-base correction in version 0.2.0
#' does not modify OPT. Do not interpret the plateau as a validated cure model.
#'
#' Exact events must lie strictly above the root's lower endpoint. Input
#' validation is incomplete; check missing values, support compatibility, and
#' iteration controls before fitting. No covariate model, interval-censoring
#' interface, posterior draws, or uncertainty intervals are returned.
#'
#' @return A list containing `Left`, `Right`, `Level`, `Allocations`,
#'   `Posterior_Probs`, `Terminal_Nodes`, `Test_seq`, `Survival_Probs`, and
#'   `Density`, as described for \code{\link{survPT}()}, plus `Rho`, the node stopping
#'   probabilities. Only terminal entries of `Posterior_Probs` are populated;
#'   internal entries are zero. `Terminal_Nodes` uses one-based R indexing.
#'   `Rho` can include unused trailing storage: use its first `length(Left)`
#'   elements. Numeric arrays can be one-column matrices. No `Alpha` is
#'   returned by this wrapper.
#' @references
#' Seok, J., Tian, L., and Wong, W. H. (2014).
#' Density estimation on multivariate censored data with optional Polya tree.
#' \emph{Biostatistics}, 15(1), 182--195.
#' \doi{10.1093/biostatistics/kxt025}.
#' @seealso \code{\link{survPT}()}, \code{\link{data_rc}()}, \code{\link{citation}()} for package references.
#' @export
#'
#' @examples
#' d <- data_rc(seed = 123, n = 30, K = 1.5, Lambda = 0.4, q = 0.3)
#' fit <- survOPT(d$L, d$status, test_seq = seq(0, 12, length.out = 31),
#'                max_M = 4, max_iter = 30, verbose = FALSE)
#' sum(fit$Posterior_Probs[fit$Terminal_Nodes])
#' citation("SurvPT")
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
