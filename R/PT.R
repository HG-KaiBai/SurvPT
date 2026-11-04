#' Fit a self-consistent Polya tree model
#'
#' Fit the deterministic MAP-EM model independently of prediction times and
#' posterior sampling. Use `predict(fit, times, 0.95)` to run sampling and obtain
#' survival predictions, credible intervals and tail probability summaries.
#'
#' @param X Numeric matrix or data frame with left and right endpoints in its
#'   first two columns. Use `L == R` for exact events, `R = Inf` for right
#'   censoring, and `L = 0` for left censoring. Extra columns are ignored.
#' @param prior Body centering family: `1` uniform, `2` truncated exponential,
#'   or `3` truncated Weibull. Weibull is an extension of the manuscript.
#' @param max_M Maximum ordinary tree depth, from 0 to 30; root depth is zero.
#'   The separate tail split can add one more level.
#' @param prior_params Exponential rate, or Weibull `c(rate, shape)` with scale
#'   `1/rate`. If omitted for exponential, estimate the censored-data MLE once.
#'   Uniform ignores this argument; Weibull requires explicit parameters.
#' @param allocated `1` for MAP-EM; `0` retains the historical overlap-count
#'   estimator, with the new partition and priors, without EM or sampling.
#' @param max_EM_iter Iteration cap for the single MAP-EM fit (default 2000).
#'   Zero returns the reference distribution and its allocations without EM
#'   when `allocated=1`.
#' @param EM_tol Stopping tolerance for the Euclidean change in leaf masses.
#' @param verbose Print fitting progress.
#' @param gamma Nonnegative prior strength, default `0.5`. At child depth `m`
#'   and reference left share `r`, the Beta parameters are
#'   `max(1, 2*gamma*m^2*r)` and `max(1, 2*gamma*m^2*(1-r))`.
#'   Zero gives the flat Beta(1,1) prior.
#' @param kappa Positive boundary-extension calibration factor, default `0.1`.
#'   Used only to choose `s0`; it imposes no finite tail cutoff.
#' @param s0 Optional reference tail mass in `(0,1)`. Defaults to
#'   `kappa*tau*p0(tau)/(1 + kappa*tau*p0(tau))`, where `p0` is the normalized
#'   body density. For uniform centering this is `kappa/(1+kappa)`.
#' @param min_count Minimum reference expected count for ordinary splitting,
#'   default 1. Zero disables count-based pruning. The tail is always separated.
#' @param split_epsilon Base scale for the MAP-only probability bound, default
#'   `1e-6`. Each child at level `m` uses bounds
#'   `[split_epsilon * 2^(-m), 1 - split_epsilon * 2^(-m)]`, including the
#'   final tail split. Root level is zero. The Gibbs sampler uses unrestricted
#'   Beta draws.
#' @param flat Use Beta(1,1) at every split for both fitting and posterior
#'   sampling (default `FALSE`). Overrides `gamma` when `TRUE`, without changing
#'   the centering distribution or partition. There is no preliminary fit.
#'
#' @details
#' The support is fixed at `[0, Inf]`. Ordinary cuts are medians of the centering
#' distribution on the node's intersection with `[0,tau]`. The rightmost leaf is split at `tau`, producing
#' a separate tail leaf. Initial expected counts under the reference measure
#' `(1-s0)*P0 + s0*tail` determine pruning. The resulting tree stays fixed.
#' Allocation starts from this reference distribution and a single MAP-EM fit
#' uses the selected prior throughout. `flat=TRUE` is equivalent to `gamma=0`.
#'
#' Each EM E-step uses leaf masses times conditional body probabilities, or a
#' tail factor of one for right censoring and zero otherwise. Each M-step uses
#' the clipped Beta mode with fractional expected counts. The initial reference
#' distribution is distinct from the actual Beta-prior mean when the lower
#' bound of one on a Beta parameter is active. Both tail quantities are returned.
#'
#' Fitting leaves the tail shape unspecified. `Tail_Prob` is the MAP estimate
#' of `S(tau)`. Prediction adds an exponential tail whose rate is determined
#' by this mass, or by each posterior draw's tail mass; see
#' \code{\link{predict.survPT}()}. Unrestricted mean event time is not identified
#' by the fit without such an extrapolation assumption.
#' Exact events at a split belong to its left interval; zero is included in
#' intervals starting at zero. `Posterior_Probs` remains the historical name
#' for deterministic MAP node masses, not posterior mean masses.
#'
#' @return A serializable S3 model of class `survPT`, containing `Left`, `Right`,
#'   `Level`, `Allocations`, `Posterior_Probs`, and `Terminal_Nodes`.
#'   Predictions and posterior summaries are returned separately by `predict()`.
#'   Additional fields are:
#'   \describe{
#'     \item{Alpha, Parent, Left_Child, Right_Child, Split_Nodes, Branch_Probs}{
#'       Prior parameters, one-based tree links, and fitted conditional branch
#'       probabilities. Root parents and absent children are `NA`.}
#'     \item{Reference_Probs, Reference_Branch_Probs, Initial_Allocations}{
#'       Reference masses, shares and expected counts used to build the tree.}
#'     \item{Tau, Tail_Node, Tail_Prob}{Largest finite endpoint, tail index and
#'       fitted tail probability.}
#'     \item{Prior}{Centering parameters, gamma, flat, kappa, reference `s0`, and
#'       `tail_mean`, the actual induced Beta-prior mean after parameter flooring.}
#'     \item{Convergence}{The single `map` fit's iteration count, convergence flag,
#'       final change, and log posterior kernel omitting constants independent
#'       of the branch probabilities. A cap does not imply convergence.}
#'     \item{Model}{Internal plain-R fitted state, including training intervals,
#'       retained tree and fitted probabilities. It supports prediction and
#'       sampling after `saveRDS()` and `readRDS()`; do not edit it manually.}
#'   }
#' @seealso \code{\link{predict.survPT}()}, \code{\link{estimate_exp}()}
#' @export
#' @examples
#' x <- data.frame(L=c(0.3,0.7,1.2,1.8,2.4,3),
#'                 R=c(0.3,1.1,1.2,2.5,Inf,Inf))
#' fit <- survPT(x, prior=2, prior_params=0.5, max_M=4, verbose=FALSE)
#' sum(fit$Posterior_Probs[fit$Terminal_Nodes])
#' fit$Tail_Prob
#' set.seed(42)
#' predicted <- predict(fit, c(0,0.5,1,2,3), 0.95,
#'                      n_samples=200, burn_in=50)
#' predicted
#' attr(predicted, "tail")
survPT <- function(X, prior = 1, max_M = 8, prior_params = NULL, allocated = 1,
                   max_EM_iter = 2000, EM_tol = 1e-6, verbose = TRUE,
                   gamma = 0.5, kappa = 0.1, s0 = NULL, min_count = 1,
                   split_epsilon = 1e-6, flat = FALSE) {
    if (length(dim(X)) != 2L || ncol(X) < 2L || nrow(X) == 0L)
        stop('X must have rows and two endpoint columns.')
    X <- as.matrix(X[, 1:2, drop = FALSE])
    if (!is.numeric(X) || anyNA(X) || any(!is.finite(X[,1])) ||
        any(X[,1] < 0) || any(X[,2] < X[,1])) stop('Require numeric 0 <= L <= R, finite L, and no missing endpoints.')
    scalar <- pt_scalar
    scalar(prior, 'prior', 1, 3, TRUE)
    scalar(allocated, 'allocated', 0, 1, TRUE)
    scalar(max_M, 'max_M', 0, 30, TRUE)
    scalar(max_EM_iter,'max_EM_iter',0,.Machine$integer.max,TRUE)
    for (nm in c('EM_tol','gamma','min_count')) scalar(get(nm),nm,0)
    scalar(kappa,'kappa',.Machine$double.xmin)
    scalar(split_epsilon,'split_epsilon',.Machine$double.xmin,0.5)
    if (split_epsilon >= 0.5) stop('split_epsilon must be less than 0.5.')
    if (!is.null(s0)) {
        scalar(s0,'s0',0,1)
        if (s0 <= 0 || s0 >= 1) stop('s0 must lie in (0,1).')
    }
    for (nm in c('verbose','flat'))
        if (!is.logical(get(nm)) || length(get(nm)) != 1L || is.na(get(nm))) stop(nm,' must be TRUE or FALSE.')
    tau <- max(X[is.finite(X)])
    if (tau <= 0) stop('At least one positive finite endpoint is required.')
    if (prior == 2 && is.null(prior_params)) prior_params <- estimate_exp(X[,1],X[,2])
    if (prior == 1) prior_params <- numeric()
    if (prior != 1 && (!is.numeric(prior_params) ||
        length(prior_params) != (if (prior == 2) 1L else 2L) ||
        any(!is.finite(prior_params)) || any(prior_params <= 0)))
        stop('Supply a positive exponential rate, or Weibull c(rate, shape).')
    out <- SurvPT(X,prior,max_M,prior_params,allocated,
                  max_EM_iter,EM_tol,verbose,gamma,kappa,if (is.null(s0)) -1 else s0,
                  min_count,split_epsilon,flat)
    link <- function(index) ifelse(index < 0, NA_integer_, as.integer(index)+1L)
    fit <- list(Left=out$left_bound, Right=out$right_bound, Level=out$node_levels,
        Allocations=out$n_obs_allocated, Posterior_Probs=out$posterior_probs,
        Terminal_Nodes=as.integer(out$terminal_nodes)+1L,
        Alpha=out$alpha, Parent=link(out$parent_nodes), Left_Child=link(out$left_children),
        Right_Child=link(out$right_children), Split_Nodes=as.integer(out$split_nodes)+1L,
        Branch_Probs=out$branch_probs, Reference_Probs=out$reference_probs,
        Reference_Branch_Probs=out$reference_branch_probs,
        Initial_Allocations=out$initial_allocations, Tau=out$largest_obs,
        Tail_Node=out$tail_node+1L, Tail_Prob=out$tail_prob,
        Prior=list(family=c('uniform','exponential','weibull')[prior],
                   parameters=prior_params,gamma=gamma,flat=flat,kappa=kappa,
                   s0=out$tail_prior,tail_mean=out$prior_tail_mean),
        Flat=flat, Convergence=list(map=out$map_fit),
        Model=out$model_state)
    class(fit) <- c('survPT','list')
    fit
}
