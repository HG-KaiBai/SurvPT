#' Predict survival at target times, with optional posterior sampling
#'
#' `predict(fit, x, 0.95)` accepts a numeric vector of target times `x` and
#' returns deterministic survival probabilities and posterior summaries at
#' exactly those times. By default it runs DA-Gibbs during this call; prior
#' sampling is not required. It never rebuilds the tree or refits MAP-EM.
#'
#' @param object A model returned by \code{\link{survPT}()}.
#' @param newdata Numeric vector of target event times, not a covariate matrix.
#'   Unsorted and repeated times are allowed and their order is preserved.
#'   If omitted, use 200 points on `[0,Tau]`.
#' @param level Equal-tailed pointwise credible level, supplied as the third
#'   positional argument in `predict(fit, x, 0.95)`.
#' @param ... Reserved for future use; unused arguments raise an error.
#' @param posterior Run posterior sampling and return posterior summaries
#'   (default `TRUE`). Set `FALSE` for deterministic MAP predictions only.
#' @param n_sweeps,burn_in,thin Sampler controls passed to \code{\link{samplePT}()}
#'   when `posterior=TRUE` and no explicit `samples` are supplied.
#' @param samples Optional result of \code{\link{samplePT}(object)}. When supplied,
#'   reuse these draws without further sampling. Sampler controls are ignored.
#' @param keep_draws Include evaluated curve draws and the underlying tree draws
#'   in `Posterior$Draws` (default `FALSE`).
#' @param verbose Print sampler progress.
#'
#' @details
#' Model fitting, posterior sampling and evaluation are independent operations.
#' The default prediction call orchestrates sampling and evaluation on the
#' requested time vector. The posterior on tree probabilities is independent of
#' the prediction grid; each sampled tree is evaluated at the supplied times.
#' There is no interpolation from a previously stored curve or implicit extra
#' prediction grid. `set.seed()` controls fresh sampling.
#'
#' Every call with `posterior=TRUE, samples=NULL` runs a fresh chain starting
#' from the stored MAP fit. To compare grids or credible levels using the same
#' Monte Carlo draws, obtain `samples <- samplePT(fit)` and pass it explicitly.
#' `posterior=FALSE` consumes no random draws. Neither mode mutates `object`.
#'
#' At times above `object$Tau`, including `Inf`, survival and density outputs
#' and summaries are `NA`. The model specifies only the total tail probability.
#' Credible intervals are pointwise and condition on the retained tree and
#' fixed base parameters. The separate overlap-count fit (`allocated=0`)
#' supports only `posterior=FALSE`.
#'
#' @return A list with `Test_seq`, deterministic `Survival_Probs`, `Density`,
#'   `Tail_Prob` and `Posterior`. With `posterior=FALSE`, `Posterior` is `NULL`.
#'   Otherwise it contains `Survival` and `Density` data frames with one row per
#'   target time and columns `time`, `mean`, `sd`, `lower`, `median`, `upper`,
#'   `ess`, `mcse`; `Terminal_Probs` and `Tail` summaries; sampler `Settings`;
#'   and optional `Draws`. These means are posterior means of probabilities
#'   and densities, not mean event time. ESS/MCSE are approximate single-chain
#'   diagnostics. The deterministic result is retained alongside the posterior.
#' @importFrom stats predict
#' @export
#' @examples
#' d <- data.frame(L=c(0,1,2,3), R=c(1,2,Inf,Inf))
#' fit <- survPT(d, prior=1, max_M=3, verbose=FALSE)
#' x <- c(0,0.5,1,2,3)
#' set.seed(42)
#' p <- predict(fit, x, 0.95, n_sweeps=200, burn_in=50)
#' p$Survival_Probs
#' p$Posterior$Survival
#' predict(fit, x, posterior=FALSE)$Survival_Probs
predict.survPT <- function(object, newdata = NULL, level = 0.95, ...,
                           posterior = TRUE, n_sweeps = 2000, burn_in = 500,
                           thin = 1, samples = NULL, keep_draws = FALSE,
                           verbose = FALSE) {
    pt_model(object)
    if (length(list(...))) stop('Unused prediction arguments; supply target times as newdata and the credible level as level.')
    if (is.null(newdata)) newdata <- seq(0,object$Tau,length.out=200)
    pt_times(newdata)
    pt_scalar(level,'level',0,1)
    if (level<=0 || level>=1) stop('level must lie in (0,1).')
    for (nm in c('posterior','keep_draws','verbose'))
        if (!is.logical(get(nm)) || length(get(nm))!=1L || is.na(get(nm))) stop(nm,' must be TRUE or FALSE.')
    if (!posterior && !is.null(samples)) stop('Use posterior=TRUE when supplying samples.')
    if (posterior) {
        if (is.null(samples)) samples <- samplePT(object,n_sweeps=n_sweeps,burn_in=burn_in,thin=thin,verbose=verbose)
        else if (!inherits(samples,'survPT_samples') || !identical(samples$Model,object$Model))
            stop('samples must come from samplePT() for this fitted model.')
        leaf_draws <- samples$Terminal_Probs
    } else leaf_draws <- matrix(numeric(),0L,0L)
    out <- SurvPTPredict(object$Model,as.numeric(newdata),leaf_draws)
    result <- list(Test_seq=newdata,Survival_Probs=out$survival_probs,
                   Density=out$density_probs,Tail_Prob=object$Tail_Prob,Posterior=NULL)
    if (posterior) {
        result$Posterior <- list(
            Survival=cbind(time=newdata,pt_draw_summary(out$survival_draws,level)),
            Density=cbind(time=newdata,pt_draw_summary(out$density_draws,level)),
            Terminal_Probs=cbind(node=object$Terminal_Nodes,pt_draw_summary(samples$Terminal_Probs,level)),
            Tail=pt_draw_summary(matrix(samples$Tail,ncol=1),level),
            Settings=c(samples$Settings,list(credible_level=level)))
        if (keep_draws) result$Posterior$Draws <- list(
            Survival=out$survival_draws,Density=out$density_draws,
            Terminal_Probs=samples$Terminal_Probs,Branch_Probs=samples$Branch_Probs,
            Tail=samples$Tail,Iterations=samples$Iterations)
    }
    result
}

pt_times <- function(x) {
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x) || anyNA(x) || any(x<0))
        stop('Prediction times must be a nonempty numeric vector of nonnegative, nonmissing target times.')
    invisible(x)
}

pt_draw_summary <- function(draws, level) {
    result <- lapply(seq_len(ncol(draws)),function(j) {
        x <- draws[,j]; n <- length(x)
        if (anyNA(x)) return(rep(NA_real_,7))
        qs <- stats::quantile(x,c((1-level)/2,0.5,(1+level)/2),names=FALSE)
        ess <- NA_real_
        constant <- all(x == x[1])
        s <- if (constant) 0 else stats::sd(x)
        if (all(is.finite(x)) && n >= 4L) {
            if (constant) ess <- n else {
                ac <- as.numeric(stats::acf(x,lag.max=min(n-1L,1000L),plot=FALSE)$acf)[-1]
                pairs <- ac[seq(1,length(ac)-1,by=2)] + ac[seq(2,length(ac),by=2)]
                stop_at <- which(pairs <= 0)[1]
                if (!is.na(stop_at)) pairs <- if (stop_at==1L) numeric() else pairs[seq_len(stop_at-1L)]
                ess <- min(n,n/(1+2*sum(pairs)))
            }
        }
        c(mean(x),s,qs,ess,s/sqrt(ess))
    })
    result <- as.data.frame(do.call(rbind,result))
    names(result) <- c('mean','sd','lower','median','upper','ess','mcse')
    result
}
