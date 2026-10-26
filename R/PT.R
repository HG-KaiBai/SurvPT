#' Fit a finite Polya tree to censored event times
#'
#' Estimate a univariate survival function and density from exact, interval-
#' censored, right-censored, or current-status observations. The partition is
#' fixed after initialization; the default fit iteratively reallocates censored
#' observations using conditional terminal-node base distributions.
#'
#' @param X Numeric matrix or data frame whose first two columns contain the
#'   left and right observation endpoints, in that order. Additional columns
#'   are ignored. Use `L == R` for an exact event, `R = Inf` for right censoring,
#'   and `L = 0` for an event observed before a positive inspection time.
#'   A status column is not used. Missing values are not supported.
#' @param prior Base distribution: `1` for uniform, `2` for exponential,
#'   or `3` for Weibull.
#' @param support Numeric vector `c(lower, upper)` with `lower < upper`.
#'   The default is `c(0, Inf)`. Use nonnegative time support and observation
#'   intervals compatible with it. A finite bound is a modeling assumption.
#' @param max_M Nonnegative integer maximum tree depth; the root is level zero.
#'   A node also stops splitting if its initial allocated count is below one.
#'   Large depths can require excessive memory; start with a modest value.
#' @param prior_params `NULL` for uniform (ignored), one positive rate `lambda`
#'   for exponential, or two positive values `c(lambda, k)` for Weibull.
#'   The Weibull shape is `k` and scale is `1 / lambda`.
#' @param allocated `1` (default) for fractional allocation with EM updates;
#'   `0` for the legacy overlap-count initialization without EM.
#' @param max_EM_iter Nonnegative integer maximum number of EM updates.
#'   Zero returns the initialized fit.
#' @param EM_tol Nonnegative stopping tolerance for the Euclidean norm of the
#'   change in terminal-node masses. No convergence flag is returned.
#' @param verbose Logical; print iteration progress if `TRUE`.
#' @param test_seq Numeric evaluation points, preferably a sorted finite grid.
#'   If `NULL`, uses 200 points from `support[1]` to the largest endpoint among
#'   rows with finite right endpoints. Supply a grid explicitly when no such
#'   rows exist or when comparing fits.
#'
#' @details
#' Since version 0.2.0, EM allocation probabilities and the returned survival
#' function and density use the same selected base distribution conditioned on
#' each terminal interval. The exponential base survival is `exp(-lambda * t)`;
#' the Weibull base survival is `exp(-(lambda * t)^k)`.
#' A uniform terminal interval with an infinite right endpoint is capped at
#' 1.1 times the largest finite endpoint in the first two columns of `X`.
#' Exponential and Weibull terminal intervals retain their infinite tails.
#' This cutoff does not replace `Inf` in the returned `Right` boundaries.
#'
#' Intervals represent `(L, R]`, with exact events treated separately. Exact
#' events must lie strictly above the root's lower endpoint. Input validation
#' is incomplete: check numeric types, missing values, endpoint ordering,
#' support compatibility, and integer iteration controls before fitting.
#'
#' Branch masses use a clipped, mode-style update with depth-dependent prior
#' parameters, not posterior draws or the usual posterior-mean update.
#' The function does not estimate uncertainty intervals, covariate effects,
#' multivariate distributions, or a sampling-truncation likelihood.
#'
#' @return A list with the following components. Armadillo vectors can appear
#'   as one-column matrices; use `as.numeric()` if a plain vector is needed.
#'   \describe{
#'     \item{Left, Right}{Node interval boundaries.}
#'     \item{Level}{Node depths, with the root at zero.}
#'     \item{Allocations}{Allocated observation counts, possibly fractional.}
#'     \item{Posterior_Probs}{Estimated masses for internal and terminal nodes.}
#'     \item{Terminal_Nodes}{One-based R indices of the terminal nodes. Sum
#'       masses at these indices to check total probability.}
#'     \item{Test_seq}{The evaluation grid.}
#'     \item{Survival_Probs}{Estimated survival probabilities at the grid.}
#'     \item{Density}{Estimated density at the grid.}
#'   }
#'   The wrapper does not return `Alpha`, posterior samples, or convergence
#'   diagnostics.
#' @seealso \code{\link{estimate}()}, \code{\link{estimate_exp}()}, \code{\link{data_rc}()},
#'   \code{\link{SurvPT}()} for the low-level interface.
#' @export
#'
#' @examples
#' x <- data.frame(L = c(0.3, 0.7, 1.2, 1.8, 2.4, 3),
#'                 R = c(0.3, 1.1, 1.2, 2.5, Inf, Inf))
#' fit <- survPT(x, prior = 2, prior_params = 0.5,
#'               test_seq = seq(0, 6, length.out = 31),
#'               max_M = 4, max_EM_iter = 30, verbose = FALSE)
#' sum(fit$Posterior_Probs[fit$Terminal_Nodes])
#' head(as.numeric(fit$Survival_Probs))
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
