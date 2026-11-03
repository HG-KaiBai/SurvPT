#' Fit a finite Polya tree to censored event times
#'
#' Fit a univariate survival model from exact, interval-
#' censored, right-censored, or current-status observations. The partition is
#' fixed after initialization; the default fit iteratively reallocates censored
#' observations using conditional terminal-node base distributions.
#' Evaluation times and posterior sampling are handled by \code{\link{predict.survPT}()}.
#'
#' @param X Numeric matrix or data frame whose first two columns contain the
#'   left and right observation endpoints, in that order. Additional columns
#'   are ignored. Use `L == R` for an exact event, `R = Inf` for right censoring,
#'   and `L = 0` for an event observed before a positive inspection time.
#'   A status column is not used. Missing values are not supported.
#' @param prior Base distribution: `1` for uniform, `2` for exponential,
#'   `3` for Weibull, or `4` for a uniform body with an exponential tail.
#' @param support Numeric vector `c(lower, upper)` with `lower < upper`.
#'   The default is `c(0, Inf)`. Use nonnegative time support and observation
#'   intervals compatible with it. A finite bound is a modeling assumption.
#' @param max_M Nonnegative integer maximum tree depth; the root is level zero.
#'   A node also stops splitting if its initial allocated count is below one.
#'   Large depths can require excessive memory; start with a modest value.
#' @param prior_params `NULL` for uniform (ignored), one positive rate `lambda`
#'   for exponential or the uniform/exponential splice (`prior = 4`), or two
#'   positive values `c(lambda, k)` for Weibull. The Weibull shape is `k` and
#'   scale is `1 / lambda`. For `prior = 4`, a finite positive rate is required;
#'   pass an estimate from `estimate_exp()` to use an MLE. It is held fixed.
#' @param allocated `1` (default) for fractional allocation with EM updates;
#'   `0` for the legacy overlap-count initialization without EM.
#' @param max_EM_iter Nonnegative integer maximum number of EM updates.
#'   Zero returns the initialized fit.
#' @param EM_tol Nonnegative stopping tolerance for the Euclidean norm of the
#'   change in terminal-node masses. No convergence flag is returned.
#' @param verbose Logical; print iteration progress if `TRUE`.
#' @param uniform_tail_factor Finite numeric scalar greater than 1, default
#'   `1.1`. For a uniform base, multiplies the largest finite observation
#'   endpoint when replacing an infinite bound in partitioning, initialization,
#'   EM allocation, and prediction. Used only for `prior = 1`.
#'   With finite support, keep observation endpoints within that support.
#' @param flat Logical scalar, default `FALSE`. If `TRUE`, all node prior
#'   parameters are 1. If `FALSE`, children at depth `m` use
#'   `alpha_m = max(1, 0.5 * m^2)`. The root parameter is always 1.
#'
#' @details
#' EM allocation probabilities and predictions use the same selected base
#' distribution conditioned on
#' each terminal interval. The exponential base survival is `exp(-lambda * t)`;
#' the Weibull base survival is `exp(-(lambda * t)^k)`.
#' A uniform terminal interval with an infinite right endpoint is capped at
#' `uniform_tail_factor` times the largest finite endpoint in the first two
#' columns of `X` (default multiplier `1.1`).
#' Exponential and Weibull terminal intervals retain their infinite tails.
#' This cutoff does not replace `Inf` in the returned `Right` boundaries.
#'
#' For `prior = 4`, let `M` be the largest finite endpoint in the first two
#' columns of `X` and `lambda = prior_params[1]`. The base CDF is
#' `(1 - exp(-lambda * M)) * t / M` on `[0, M]`, and
#' `1 - exp(-lambda * t)` for `t > M`. Its survival at `M` is
#' `exp(-lambda * M)`. The CDF is continuous; the density may jump at `M`,
#' where the left-hand density is returned. Splits bisect this base probability,
#' and allocation and prediction condition the same base on each terminal
#' interval. The exponential tail is unbounded with infinite upper support.
#' `M` must be positive and finite, and left endpoints finite and nonnegative.
#' The tree updates interval masses, so its fitted tail mass is not fixed to
#' the base tail mass.
#'
#' Intervals represent `(L, R]`, including zero when `L = 0`, with exact events
#' treated separately. Node allocation and density use `(a, b]`, or `[0, b]`
#' when `a = 0`; exact events at positive split points belong to the left node.
#' Input validation is incomplete: check numeric types, missing values,
#' endpoint ordering, support compatibility, and integer iteration controls.
#'
#' Branch masses use a clipped, mode-style update with prior parameters
#' `alpha_m = max(1, 0.5 * m^2)` for children at depth `m` by default,
#' or `alpha_m = 1` at every depth when `flat = TRUE`.
#' During initialization and EM, conditional split probabilities are bounded
#' by `[epsilon, 1 - epsilon]`, with `epsilon = 1e-6`. These are not posterior
#' draws or the usual posterior-mean update.
#' The fit does not generate predictions or posterior draws. Use `predict()`
#' for point estimates and pointwise posterior summaries. Covariate effects,
#' multivariate distributions and sampling-truncation likelihoods are not supported.
#'
#' @return An S3 object of class `survPT` with the following components.
#'   Armadillo vectors can appear
#'   as one-column matrices; use `as.numeric()` if a plain vector is needed.
#'   \describe{
#'     \item{Left, Right}{Node interval boundaries.}
#'     \item{Level}{Node depths, with the root at zero.}
#'     \item{Alpha}{Node prior Beta parameters, including the root placeholder.}
#'     \item{Parent}{One-based parent indices; the root has parent zero.}
#'     \item{Allocations}{Allocated observation counts, possibly fractional.}
#'     \item{Posterior_Probs}{Estimated masses for internal and terminal nodes.}
#'     \item{Terminal_Nodes}{One-based R indices of the terminal nodes. Sum
#'       masses at these indices to check total probability.}
#'     \item{X}{The two observation endpoint columns, retained for posterior sampling.}
#'     \item{Prior, Support, Max_M, Prior_Params, Allocated, Uniform_Tail_Factor, Flat}{
#'       The fitting configuration needed for prediction.}
#'     \item{Largest_Obs}{Largest finite observation endpoint.}
#'     \item{Tail_Start}{For `prior = 4` only: the splice point `M`.}
#'   }
#'   No evaluation grid is constructed. The object contains ordinary R data
#'   and can be saved with `saveRDS()` and restored with `readRDS()`.
#' @seealso \code{\link{estimate}()}, \code{\link{estimate_exp}()}, \code{\link{data_rc}()},
#'   \code{\link{SurvPT}()} for the low-level interface.
#' @export
#'
#' @examples
#' x <- data.frame(L = c(0.3, 0.7, 1.2, 1.8, 2.4, 3),
#'                 R = c(0.3, 1.1, 1.2, 2.5, Inf, Inf))
#' fit <- survPT(x, prior = 2, prior_params = 0.5,
#'               max_M = 4, max_EM_iter = 30, verbose = FALSE)
#' sum(fit$Posterior_Probs[fit$Terminal_Nodes])
#' set.seed(123)
#' head(predict(fit, seq(0, 6, length.out = 31), level = 0.95))
survPT <- function(X, prior = 1, support = c(0, Inf), max_M = 8, prior_params = NULL, allocated = 1, max_EM_iter = 100, EM_tol = 1e-6, verbose = TRUE, uniform_tail_factor = 1.1, flat = FALSE) {
    if (!is.logical(flat) || length(flat) != 1L || is.na(flat)) {
        stop("flat must be TRUE or FALSE.")
    }
    if (!is.numeric(uniform_tail_factor) || length(uniform_tail_factor) != 1L ||
        !is.finite(uniform_tail_factor) || uniform_tail_factor <= 1) {
        stop("uniform_tail_factor must be a finite number greater than 1.")
    }
    X = as.matrix(X)
    # The first column and second column of X are the observed left bound and right bound, respectively.
    if (!is.numeric(X) || is.complex(X) || nrow(X) == 0L || ncol(X) < 2L ||
        anyNA(X[, 1:2]) || any(!is.finite(X[, 1])) || any(X[, 1] < 0) ||
        any(X[, 2] < X[, 1])) {
        stop("X must contain nonempty, nonmissing, ordered endpoints with finite nonnegative left endpoints.")
    }
    X <- X[, 1:2, drop = FALSE]
    if (!is.numeric(prior) || length(prior) != 1L || is.na(prior) || !prior %in% 1:4) {
        stop("Invalid input prior.")
    }
    if (prior != 1 && (!is.numeric(prior_params) || is.complex(prior_params) ||
        length(prior_params) != (if (prior == 3) 2L else 1L) ||
        any(!is.finite(prior_params)) || any(prior_params <= 0))) {
        stop("Provide finite positive prior_params: c(rate, shape) for Weibull, or one rate otherwise.")
    }
    if (!is.numeric(support) || is.complex(support) || length(support) != 2L ||
        anyNA(support) || !is.finite(support[1]) || support[1] < 0 || support[1] >= support[2]) {
        stop("Invalid input support.")
    }
    .survpt_integer(max_M, "max_M", 0, 30)
    .survpt_integer(max_EM_iter, "max_EM_iter", 0)
    .survpt_integer(allocated, "allocated", 0, 1)
    if (!is.numeric(EM_tol) || length(EM_tol) != 1L || !is.finite(EM_tol) || EM_tol < 0) {
        stop("EM_tol must be a finite nonnegative scalar.")
    }
    exact <- X[, 1] == X[, 2]
    if (any(exact & (X[, 1] < support[1] | X[, 1] > support[2] |
                    (support[1] > 0 & X[, 1] == support[1]))) ||
        any(!exact & (X[, 2] <= support[1] | X[, 1] >= support[2]))) {
        stop("Every observation must have positive overlap with the support.")
    }

    if (prior == 4) {
        if (!is.numeric(prior_params) || is.complex(prior_params) ||
            length(prior_params) != 1L || !is.finite(prior_params) || prior_params <= 0) {
            stop("prior = 4 requires one finite positive exponential rate in prior_params.")
        }
        endpoints <- X[, 1:2, drop = FALSE]
        if (!is.numeric(endpoints) || is.complex(endpoints) || nrow(endpoints) == 0L ||
            anyNA(endpoints) || any(!is.finite(endpoints[, 1])) || any(endpoints[, 1] < 0)) {
            stop("prior = 4 requires finite nonnegative left endpoints and nonmissing right endpoints.")
        }
        tail_start <- max(endpoints[is.finite(endpoints)])
        if (!is.finite(tail_start) || tail_start <= 0) {
            stop("prior = 4 requires a positive finite largest observation endpoint.")
        }
        if (!is.finite(support[1]) || support[1] < 0) {
            stop("prior = 4 requires nonnegative, ordered support bounds.")
        }
    }

    # Handle NULL prior_params for uniform prior
    if (prior == 1 && is.null(prior_params)) {
        prior_params = c(0.12, 1)  # Default values, not used for uniform prior
    }

    out = SurvPT(X, prior, support, max_M, prior_params, allocated, max_EM_iter, EM_tol, verbose, uniform_tail_factor, flat)
    out_PT = list("Left" = out$left_bound,
                  "Right" = out$right_bound,
                  "Level" = out$node_levels,
                  "Alpha" = out$alpha,
                  "Parent" = out$parent_nodes + 1,
                  "Allocations" = out$n_obs_allocated,
                  "Posterior_Probs" = out$posterior_probs,
                  "Terminal_Nodes" = out$terminal_nodes + 1,
                  "X" = X,
                  "Prior" = prior,
                  "Support" = support,
                  "Max_M" = max_M,
                  "Prior_Params" = prior_params,
                  "Allocated" = allocated,
                  "Uniform_Tail_Factor" = uniform_tail_factor,
                  "Flat" = flat,
                  "Largest_Obs" = out$largest_obs
                  )
    if (prior == 4) {
        out_PT$Prior_Params <- as.numeric(prior_params)
        out_PT$Tail_Start <- tail_start
    }
    class(out_PT) <- "survPT"
    return(out_PT)
}
