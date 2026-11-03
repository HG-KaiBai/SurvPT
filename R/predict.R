.survpt_integer <- function(value, name, minimum, maximum = .Machine$integer.max) {
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value != floor(value) || value < minimum || value > maximum) {
        stop(name, " must be an integer between ", minimum, " and ", maximum, ".", call. = FALSE)
    }
    invisible(value)
}

#' Predict from a fitted Polya tree with posterior uncertainty
#'
#' Evaluate the fitted MAP-EM distribution and summarize posterior draws at
#' specified times. By default, returns survival probabilities `P(T > x)`.
#'
#' @param object A fitted object returned by \code{\link{survPT}()}.
#' @param x Numeric vector of evaluation times, in any order. Duplicates and
#'   infinite endpoints are allowed; missing values are not.
#' @param level Credible interval probability, strictly between 0 and 1.
#' @param sample Logical, default `TRUE`. If `FALSE`, skip posterior sampling
#'   and return only `x` and the fitted `estimate`.
#' @param n_samples Number of retained posterior draws, default 5000. Use zero
#'   as an alternative to `sample = FALSE`.
#' @param burn_in Number of initial Gibbs sweeps to discard, default 3000.
#' @param thin Retain every `thin`-th sweep after burn-in, default 1.
#' @param keep_samples Logical; attach the retained draws as a matrix in the
#'   `samples` attribute. Rows are draws and columns correspond to `x`.
#' @param type Either `"survival"` (default) or `"density"`.
#' @param ... Reserved for future arguments. Unrecognized arguments are rejected.
#'
#' @details
#' `estimate` evaluates the distribution at the fitted tree parameters. It is
#' the MAP-EM plug-in estimate, not the mode of the marginal posterior of
#' `S(x)`. At each split the fitted update uses the Beta posterior mode
#' `(alpha_left + N_left - 1) / (alpha_left + alpha_right + N_left + N_right - 2)`,
#' with the existing clipping and zero-denominator fallback. It does not use
#' the Beta posterior mean. With `allocated = 0` or `max_EM_iter = 0`, it evaluates the
#' initialized fit, as no EM optimization was performed.
#'
#' Posterior draws use data-augmentation Gibbs sampling on the fixed fitted
#' tree. Each sweep assigns every censored observation to a terminal node with
#' probability proportional to that node's current mass times its conditional
#' base probability over the observation interval. Exact-event labels are fixed.
#' Integer counts are aggregated up the tree and each left split is drawn from
#' `Beta(alpha_left + count_left, alpha_right + count_right)`. Counts and
#' allocations are sampled afresh each sweep, rather than treating fitted
#' fractional EM counts as observed counts. Sampling starts at the fitted
#' masses and uses unrestricted Beta draws, without the MAP clipping bound.
#'
#' Intervals are pointwise equal-tailed credible intervals, conditional on the
#' retained tree and the supplied base parameters. They are not simultaneous
#' bands. `mean` and `median` summarize the posterior of the predicted survival
#' probability (or density), not the distribution of event times. The total
#' number of sweeps is `burn_in + n_samples * thin`. Use `set.seed()` for
#' reproducibility, and `keep_samples = TRUE` to inspect Monte Carlo behavior.
#' Finite runs have Monte Carlo error; mixing is not checked automatically.
#' All original partition and within-node base conventions are retained.
#'
#' @return A data frame with columns `x`, `estimate`, `lower`, `upper`, `mean`
#'   and `median`, preserving the order of `x`. If `sample = FALSE` (or
#'   `n_samples = 0`), only `x` and `estimate` are returned. Attributes record
#'   `type` and `sample`, plus `level`, `n_samples`, `burn_in` and `thin` when
#'   sampling is enabled. If requested, `samples` contains the
#'   retained draws. An empty `x` returns a zero-row data frame.
#' @seealso \code{\link{survPT}()}
#' @export
#' @importFrom stats predict
#' @examples
#' observations <- rbind(c(0.2, 0.2), c(0.4, 0.8), c(1, Inf))
#' model <- survPT(observations, prior = 2, prior_params = 1, verbose = FALSE)
#' set.seed(123)
#' predict(model, c(0, 0.5, 1, 2), level = 0.95)
#' predict(model, c(0, 0.5, 1, 2), sample = FALSE)
predict.survPT <- function(object, x, level = 0.95, sample = TRUE, n_samples = 5000L,
                           burn_in = 3000L, thin = 1L, keep_samples = FALSE,
                           type = c("survival", "density"), ...) {
    if (length(list(...))) stop("Unused arguments in predict.survPT().")
    if (!inherits(object, "survPT")) stop("object must be a fitted survPT model.")
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || anyNA(x)) {
        stop("x must be a numeric vector of times without missing values.")
    }
    if (!is.logical(sample) || length(sample) != 1L || is.na(sample)) {
        stop("sample must be TRUE or FALSE.")
    }
    if (!is.numeric(level) || is.complex(level) || length(level) != 1L ||
        !is.finite(level) || level <= 0 || level >= 1) {
        stop("level must be a number strictly between 0 and 1.")
    }
    .survpt_integer(n_samples, "n_samples", 0)
    .survpt_integer(burn_in, "burn_in", 0)
    .survpt_integer(thin, "thin", 1)
    if (!is.logical(keep_samples) || length(keep_samples) != 1L || is.na(keep_samples)) {
        stop("keep_samples must be TRUE or FALSE.")
    }
    type <- match.arg(type)
    sample <- sample && n_samples > 0L
    prediction <- predict_survPT_cpp(object, as.numeric(x), if (sample) as.integer(n_samples) else 0L,
                                    as.integer(burn_in), as.integer(thin), type == "density")
    result <- data.frame(x = x, estimate = as.numeric(prediction$estimate))
    attr(result, "type") <- type
    attr(result, "sample") <- sample
    if (!sample) return(result)
    result$lower <- rep(NA_real_, length(x))
    result$upper <- rep(NA_real_, length(x))
    result$mean <- rep(NA_real_, length(x))
    result$median <- rep(NA_real_, length(x))
    if (n_samples > 0L && length(x)) {
        probabilities <- c((1 - level)/2, 0.5, (1 + level)/2)
        intervals <- vapply(seq_along(x), function(j)
            stats::quantile(prediction$samples[, j], probabilities, names = FALSE), numeric(3))
        result$lower <- intervals[1, ]
        result$upper <- intervals[3, ]
        result$mean <- colMeans(prediction$samples)
        result$median <- intervals[2, ]
    }
    attr(result, "level") <- level
    attr(result, "n_samples") <- n_samples
    attr(result, "burn_in") <- burn_in
    attr(result, "thin") <- thin
    if (keep_samples) attr(result, "samples") <- prediction$samples
    result
}
