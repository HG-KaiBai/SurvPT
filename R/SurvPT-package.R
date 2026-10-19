#' Survival estimation with finite and optional Polya trees
#'
#' SurvPT provides univariate finite Polya tree estimation for censored event
#' times, an optional Polya tree interface for right-censored observations,
#' parametric likelihood helpers, and simulation generators.
#'
#' Start with \code{\link{survPT}()} for endpoint data or \code{\link{survOPT}()} for observed times
#' and event indicators. See their help pages for support and tail conventions.
#' The optional Polya tree interface is univariate; the cited article has a
#' broader multivariate scope. Use `citation("SurvPT")` for references.
#'
#' This is a private research package. No software license has been selected;
#' obtain the maintainer's permission for use and redistribution. The package
#' provides point estimates, not posterior samples or uncertainty intervals.
#'
#' @references
#' Seok, J., Tian, L., and Wong, W. H. (2014).
#' Density estimation on multivariate censored data with optional Polya tree.
#' \emph{Biostatistics}, 15(1), 182--195.
#' \doi{10.1093/biostatistics/kxt025}.
#' @name SurvPT-package
#' @aliases NULL SurvPT-package
#' @useDynLib SurvPT, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim rexp rgamma rnbinom runif rweibull
"_PACKAGE"
