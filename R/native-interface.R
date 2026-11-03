#' Low-level native self-consistent Polya tree interface
#'
#' Native entry point for the retained tree and deterministic MAP-EM fit.
#' Prefer \code{\link{survPT}()} for validation, automatic rate estimation,
#' one-based indices and a serializable S3 model.
#' @name SurvPT-native
#' @aliases SurvPT
#' @usage SurvPT(X, prior, max_M, prior_params, allocated,
#'   max_EM_iter, EM_tol, verbose, gamma = 0.5, kappa = 0.1, s0 = -1,
#'   min_count = 1, split_epsilon = 1e-6, flat_init = TRUE)
#' @param X Numeric matrix with left and right endpoints in the first two columns.
#' @param prior Body family: uniform (1), exponential (2), or Weibull (3).
#' @param max_M Maximum ordinary depth, from 0 to 30.
#' @param prior_params Exponential rate or Weibull `c(rate, shape)`; use
#'   `numeric()` for uniform. Automatic MLE fitting is wrapper-only.
#' @param allocated MAP-EM (1) or historical overlap counts (0).
#' @param max_EM_iter Maximum iterations per EM stage.
#' @param EM_tol Leaf-mass change tolerance.
#' @param verbose Print progress.
#' @param gamma Nonnegative Beta-prior strength.
#' @param kappa Positive boundary-density calibration factor.
#' @param s0 Reference tail mass; `-1` requests automatic calibration.
#' @param min_count Reference-count threshold for splitting.
#' @param split_epsilon MAP-only branch probability bound.
#' @param flat_init Fit a flat-prior initialization before the working prior.
#' @return Native node arrays and diagnostics. All node indices are zero-based;
#'   absent links are `-1`. `alpha` stores each child's prior parameter.
#'   `model_state` stores the fitted state for prediction and sampling.
#'   Support is fixed at `[0, Inf]`. See \code{\link{survPT}()} for interpretation.
#' @export SurvPT
NULL
