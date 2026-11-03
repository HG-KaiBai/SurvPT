#' Sample the posterior of an already fitted Polya tree
#'
#' Run DA-Gibbs on a stored model without refitting, rebuilding its partition,
#' estimating base parameters, or evaluating a prediction grid. This is the
#' sampling function used by \code{\link{predict.survPT}()}.
#'
#' @param object A model returned by \code{\link{survPT}()}, with `allocated=1`.
#' @param n_sweeps Total Gibbs sweeps including burn-in, default 2000.
#' @param burn_in Number of initial sweeps to discard, default 500.
#' @param thin Retain every `thin`th sweep starting at `burn_in+1`, default 1.
#' @param verbose Print sampler progress, default `FALSE`.
#' @details Each call starts from the stored deterministic fit and runs a fresh
#'   chain. Use `set.seed()` for reproducibility. No continuous event times are
#'   drawn: labels and unrestricted conjugate Beta branches are sampled. The
#'   posterior on tree masses does not depend on target prediction times.
#'   Pass the result to `predict(fit, times, level, samples=draws)` to evaluate
#'   the same draws on different grids or at different credible levels.
#' @return A serializable object of class `survPT_samples`, containing matrices
#'   `Terminal_Probs` and `Branch_Probs`, vector `Tail`, retained `Iterations`,
#'   and sampler `Settings`. Rows index retained sweeps; columns follow the
#'   model's `Terminal_Nodes` or `Split_Nodes`. `Model` records the fitted state
#'   for checking compatibility when reusing draws. The input model is unchanged.
#' @seealso \code{\link{survPT}()}, \code{\link{predict.survPT}()}
#' @export
#' @examples
#' d <- data.frame(L=c(0,1,2,3), R=c(1,2,Inf,Inf))
#' fit <- survPT(d, max_M=3, verbose=FALSE)
#' set.seed(42)
#' draws <- samplePT(fit, n_sweeps=200, burn_in=50)
#' p <- predict(fit, c(0,1,2,3), 0.95, samples=draws)
#' p$Posterior$Survival
samplePT <- function(object, n_sweeps = 2000, burn_in = 500,
                     thin = 1, verbose = FALSE) {
    pt_model(object)
    for (nm in c('n_sweeps','thin')) pt_scalar(get(nm),nm,1,.Machine$integer.max,TRUE)
    pt_scalar(burn_in,'burn_in',0,.Machine$integer.max,TRUE)
    if (burn_in>=n_sweeps) stop('n_sweeps must exceed burn_in.')
    if (!is.logical(verbose) || length(verbose)!=1L || is.na(verbose)) stop('verbose must be TRUE or FALSE.')
    if (object$Model$allocated!=1L) stop('DA-Gibbs requires allocated=1; use predict(..., posterior=FALSE) for overlap-count fits.')
    draws <- SurvPTSample(object$Model,n_sweeps,burn_in,thin,verbose)
    tail_column <- match(object$Tail_Node,object$Terminal_Nodes)
    result <- list(Terminal_Probs=draws$terminal_probs,Branch_Probs=draws$branch_probs,
        Tail=as.numeric(draws$terminal_probs[,tail_column]),Iterations=draws$iterations,
        Settings=list(n_sweeps=n_sweeps,burn_in=burn_in,thin=thin,
            retained=nrow(draws$terminal_probs),
            target='unrestricted posterior conditional on the retained tree and base parameters'),
        Model=object$Model)
    class(result) <- c('survPT_samples','list')
    result
}
