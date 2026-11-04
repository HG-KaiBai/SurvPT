#' Predict from a fitted Polya tree with posterior uncertainty
#'
#' Evaluate the fitted MAP-EM distribution and optionally sample its posterior
#' at specified times. The interface follows version 0.2.3, with a separate
#' tail probability summary for the unspecified tail in version 0.4.0.
#'
#' @param object A fitted model returned by \code{\link{survPT}()}.
#' @param x Numeric vector of nonnegative evaluation times, in any order.
#'   Duplicates and `Inf` are allowed; missing values are not. An empty vector
#'   returns no curve rows but still provides the tail summary.
#' @param level Pointwise credible interval probability, strictly between 0 and 1.
#' @param sample Logical, default `TRUE`. Set `FALSE` for fitted estimates only.
#' @param n_samples Number of posterior draws used for summaries, default 2000. Zero also
#'   disables sampling.
#' @param burn_in Number of initial Gibbs sweeps to discard, default 3000.
#' @param thin Retain every `thin`-th sweep after burn-in, default 1.
#' @param type Either `"survival"` (default) or `"density"`.
#' @param ... Reserved for future arguments. Unrecognized arguments are rejected.
#'
#' @details
#' `estimate` is the fitted MAP-EM plug-in prediction. It is not the mode of the
#' marginal posterior of survival or density. No refitting or tree rebuilding
#' takes place during prediction. With `max_EM_iter=0`, the estimate evaluates
#' the reference initialization. With `allocated=0`, only `sample=FALSE` (or
#' `n_samples=0`) is supported.
#'
#' Each Gibbs sweep allocates censored observations to leaves with probabilities
#' proportional to current leaf mass times conditional interval coverage.
#' Exact-event labels are fixed. Integer counts are aggregated up the tree,
#' then each left split is drawn from `Beta(alpha_left + count_left,
#' alpha_right + count_right)`. Fractional fitted counts are not treated as
#' observed counts. Draws are unrestricted, without MAP clipping. The stored
#' prior is used, including Beta(1,1) when the fit used `flat=TRUE`.
#'
#' Each call starts from the fitted masses and runs `burn_in + n_samples*thin`
#' sweeps. Retention starts at `burn_in + thin`. `set.seed()` controls sampling;
#' `sample=FALSE` consumes no random draws. The input model is unchanged and
#' remains usable after `saveRDS()` and `readRDS()`.
#'
#' Intervals are equal-tailed and pointwise, conditional on the retained tree
#' and fixed centering parameters. Finite chains have Monte Carlo error;
#' compare summaries across independent runs and longer chains to assess stability.
#'
#' At times above `object$Tau`, including `Inf`, curve estimates and posterior
#' summaries are `NA`: the model specifies only total tail probability, with
#' no shape inside the tail. `attr(result, "tail")` always reports this
#' probability `S(Tau)`, even when `x` omits `Tau` or `type="density"`.
#'
#' @return A data frame with columns `x`, `estimate`, `lower`, `upper`, `mean`
#'   and `median`, preserving the order of `x`. Without sampling, only `x` and
#'   `estimate` are returned. Attributes record `type` and `sample`, plus
#'   `level`, `n_samples`, `burn_in` and `thin` when sampling is enabled.
#'   The `tail` attribute is a one-row data frame with the same columns:
#'   `x` is `Tau`, `estimate` is `object$Tail_Prob`, and the remaining columns
#'   summarize the posterior tail probability. Draws are used internally to
#'   compute summaries and are not retained in the returned object.
#' @seealso \code{\link{survPT}()}
#' @export
#' @importFrom stats predict
#' @examples
#' d <- data.frame(L=c(0,1,2,3), R=c(1,2,Inf,Inf))
#' fit <- survPT(d, max_M=3, verbose=FALSE)
#' set.seed(42)
#' p <- predict(fit, c(0,0.5,1,2,3), n_samples=200, burn_in=50)
#' p
#' attr(p, "tail")
#' predict(fit, c(1,2), sample=FALSE)
predict.survPT <- function(object, x, level = 0.95, sample = TRUE,
                           n_samples = 2000L, burn_in = 3000L, thin = 1L,
                           type = c("survival", "density"), ...) {
    pt_model(object)
    if (length(list(...))) stop('Unused arguments in predict.survPT().')
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || anyNA(x) || any(x < 0))
        stop('x must be a numeric vector of nonnegative times without missing values.')
    pt_scalar(level, 'level', 0, 1)
    if (level <= 0 || level >= 1) stop('level must lie in (0,1).')
    if (!is.logical(sample) || length(sample) != 1L || is.na(sample))
        stop('sample must be TRUE or FALSE.')
    pt_scalar(n_samples, 'n_samples', 0, .Machine$integer.max, TRUE)
    pt_scalar(burn_in, 'burn_in', 0, .Machine$integer.max, TRUE)
    pt_scalar(thin, 'thin', 1, .Machine$integer.max, TRUE)
    type <- match.arg(type)
    sample <- sample && n_samples > 0L
    if (sample && object$Model$allocated != 1L)
        stop('DA-Gibbs requires allocated=1; use sample=FALSE for overlap-count fits.')

    leaf_draws <- matrix(numeric(), 0L, 0L)
    if (sample) {
        draws <- SurvPTSample(object$Model, n_samples, burn_in, thin, FALSE)
        leaf_draws <- draws$terminal_probs
    }
    out <- SurvPTPredict(object$Model, as.numeric(x), leaf_draws)
    estimate <- if (type == 'survival') out$survival_probs else out$density_probs
    result <- data.frame(x=x, estimate=as.numeric(estimate))
    tail <- data.frame(x=object$Tau, estimate=object$Tail_Prob)
    if (sample) {
        curve_draws <- if (type == 'survival') out$survival_draws else out$density_draws
        tail_draws <- leaf_draws[, match(object$Tail_Node, object$Terminal_Nodes)]
        result <- cbind(result, pt_prediction_summary(curve_draws, level))
        tail <- cbind(tail, pt_prediction_summary(matrix(tail_draws, ncol=1L), level))
        attr(result, 'level') <- level
        attr(result, 'n_samples') <- n_samples
        attr(result, 'burn_in') <- burn_in
        attr(result, 'thin') <- thin
    }
    attr(result, 'type') <- type
    attr(result, 'sample') <- sample
    attr(result, 'tail') <- tail
    result
}

pt_prediction_summary <- function(draws, level) {
    values <- vapply(seq_len(ncol(draws)), function(j) {
        x <- draws[, j]
        if (anyNA(x)) return(rep(NA_real_, 4L))
        qs <- stats::quantile(x, c((1-level)/2, (1+level)/2, 0.5), names=FALSE)
        c(qs[1:2], mean(x), qs[3])
    }, numeric(4L))
    result <- as.data.frame(t(values))
    names(result) <- c('lower', 'upper', 'mean', 'median')
    result
}
