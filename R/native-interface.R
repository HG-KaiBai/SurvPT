#' Low-level native finite Polya tree interface
#'
#' Direct Rcpp entry point for the PT implementation. Most users should call
#' \code{\link{survPT}()}, which validates some inputs, provides defaults, and converts
#' terminal-node indices to R indexing. This interface returns native tree arrays.
#' All arguments except `uniform_tail_factor` and `flat` are required.
#'
#' @name SurvPT-native
#' @aliases SurvPT
#' @usage SurvPT(X, prior, support, max_M, prior_params, allocated,
#'   max_EM_iter, EM_tol, verbose, uniform_tail_factor = 1.1, flat = FALSE)
#' @param X Numeric matrix with left and right endpoints in its first two
#'   columns; see \code{\link{survPT}()} for supported observation encodings.
#' @param prior Integer base selector: `1` uniform, `2` exponential, `3` Weibull,
#'   or `4` uniform body with exponential tail.
#' @param support Numeric vector of lower and upper support bounds.
#' @param max_M Nonnegative integer maximum tree depth.
#' @param prior_params Numeric base parameters: one rate for exponential or
#'   the uniform/exponential splice (`prior = 4`),
#'   `c(lambda, k)` for Weibull, or an ignored numeric placeholder for uniform.
#' @param allocated `1` for EM allocation or `0` for legacy overlap counts.
#' @param max_EM_iter Nonnegative integer maximum number of EM updates.
#' @param EM_tol Nonnegative tolerance for terminal-mass changes.
#' @param verbose Logical progress-output flag.
#' @param uniform_tail_factor Finite numeric scalar greater than 1, default
#'   `1.1`. Multiplier for the uniform base's finite tail cutoff; see
#'   \code{\link{survPT}()}. Used only for `prior = 1`.
#' @param flat Logical scalar, default `FALSE`. If `TRUE`, all node prior
#'   parameters are 1. If `FALSE`, children at depth `m` use
#'   `alpha_m = max(1, 0.5 * m^2)`. The root parameter is always 1.
#' @details
#' Pass a numeric matrix as `X` and numeric vectors for vector arguments.
#' For a uniform base, pass a numeric placeholder such as `c(0.12, 1)` for
#' `prior_params`; it is ignored. The wrapper's input validation is bypassed.
#' Use the endpoint, support, and parameter conventions in \code{\link{survPT}()}.
#' For `prior = 4`, supply a finite positive exponential rate in `prior_params`.
#' No rate estimation is performed by this interface.
#'
#' @return A list with native names: `left_bound`, `right_bound`, `node_levels`,
#'   `alpha`, `n_obs_allocated`, `posterior_probs`, `terminal_nodes`,
#'   `parent_nodes`, and `largest_obs`. No prediction grid or probabilities
#'   are evaluated. The `alpha` component contains
#'   branch prior parameters. Native `terminal_nodes` indices are **zero-based**;
#'   add one before using them to index R arrays. `parent_nodes` also uses
#'   zero-based indices, with -1 for the root. Other components correspond
#'   to those documented for \code{\link{survPT}()}.
#' @seealso \code{\link{survPT}()}
#' @export SurvPT
NULL
