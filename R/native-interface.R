#' Low-level native finite Polya tree interface
#'
#' Direct Rcpp entry point for the PT implementation. Most users should call
#' \code{\link{survPT}()}, which validates some inputs, provides defaults, and converts
#' terminal-node indices to R indexing. This interface is retained for
#' compatibility and requires all arguments explicitly.
#'
#' @name SurvPT-native
#' @aliases SurvPT
#' @usage SurvPT(X, test_seq, prior, support, max_M, prior_params, allocated,
#'   max_EM_iter, EM_tol, verbose)
#' @param X Numeric matrix with left and right endpoints in its first two
#'   columns; see \code{\link{survPT}()} for supported observation encodings.
#' @param test_seq Numeric evaluation points.
#' @param prior Integer base selector: `1` uniform, `2` exponential, `3` Weibull.
#' @param support Numeric vector of lower and upper support bounds.
#' @param max_M Nonnegative integer maximum tree depth.
#' @param prior_params Numeric base parameters: one rate for exponential,
#'   `c(lambda, k)` for Weibull, or an ignored numeric placeholder for uniform.
#' @param allocated `1` for EM allocation or `0` for legacy overlap counts.
#' @param max_EM_iter Nonnegative integer maximum number of EM updates.
#' @param EM_tol Nonnegative tolerance for terminal-mass changes.
#' @param verbose Logical progress-output flag.
#' @details
#' Pass a numeric matrix as `X` and numeric vectors for vector arguments.
#' For a uniform base, pass a numeric placeholder such as `c(0.12, 1)` for
#' `prior_params`; it is ignored. The wrapper's input validation is bypassed.
#' Use the endpoint, support, and parameter conventions in \code{\link{survPT}()}.
#'
#' @return A list with native names: `left_bound`, `right_bound`, `node_levels`,
#'   `alpha`, `n_obs_allocated`, `posterior_probs`, `terminal_nodes`,
#'   `survival_probs`, and `density_probs`. The `alpha` component contains
#'   branch prior parameters. Native `terminal_nodes` indices are **zero-based**;
#'   add one before using them to index R arrays. Other components correspond
#'   to those documented for \code{\link{survPT}()}.
#' @seealso \code{\link{survPT}()}
#' @export SurvPT
NULL
