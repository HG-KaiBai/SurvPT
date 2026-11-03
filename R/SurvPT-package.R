#' Survival estimation with finite Polya trees
#'
#' SurvPT provides univariate finite Polya tree estimation for censored event
#' times, with uniform, exponential, or Weibull conditional base distributions.
#' Parametric likelihood helpers and simulation generators support fitting
#' and evaluating the estimator.
#'
#' Start with \code{\link{survPT}()} and encode observations as left and right
#' interval endpoints. Exact events, interval censoring, right censoring,
#' and current status share this interface. See its help page for base
#' parameters, support, tail conventions, and returned survival and density
#' estimates. Use `citation("SurvPT")` for the software reference.
#' Model fitting is independent of the prediction grid: call
#' `predict(fit, times, 0.95)` to run posterior sampling and prediction at the
#' requested times. `samplePT()` also exposes standalone sampling.
#'
#' SurvPT is distributed under the MIT license. It provides deterministic
#' MAP-EM fits and optional DA-Gibbs posterior means and pointwise credible
#' intervals. The tail is represented only by its probability mass for every
#' body centering family; predictions beyond the observed range are undefined.
#'
#' @seealso \code{\link{survPT}()}, \code{\link{estimate}()}, \code{\link{estimate_exp}()},
#'   \code{\link{data_rc}()}, \code{\link{data_ic}()}, \code{\link{data_cs}()}
#' @name SurvPT-package
#' @aliases NULL SurvPT-package
#' @useDynLib SurvPT, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim rexp rgamma rnbinom runif rweibull
#' @rawNamespace export(survOPT)
"_PACKAGE"
