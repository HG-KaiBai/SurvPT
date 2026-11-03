# User guide

[Back to README](../README.md)

## Encode observations

Supply `L,R` as the first two columns of a matrix or data frame. Exact events
use `L=R`, interval censoring uses finite `L<R`, right censoring uses `R=Inf`,
and left censoring uses `L=0`. The same encoding covers current status. Require
finite nonnegative left endpoints, ordered endpoints, no missing values, and
at least one positive finite endpoint. Extra columns are ignored.

```r
library(SurvPT)
dat <- data.frame(L=c(0,0.4,1.1,2,3), R=c(0.8,0.4,2,Inf,Inf))
```

## Fit the deterministic estimator

```r
fit <- survPT(dat, prior=2, max_M=8, verbose=FALSE)
fit$Prior$parameters  # exponential MLE, computed once before fitting
fit$Convergence$map[c("iterations", "converged", "difference")]
fit$Tail_Prob         # S(tau), with tau = largest finite endpoint
x <- seq(0,fit$Tau,length.out=101)  # numeric vector of target times t
pred <- predict(fit,x,posterior=FALSE)
plot(pred$Test_seq,as.numeric(pred$Survival_Probs),type="l",ylim=c(0,1))
```

`survPT()` fits only: it does not evaluate a grid or run a
sampler. It returns an S3 model with the retained partition, prior, fitted node
masses, training intervals and convergence information. Prediction and sampling
reuse this state without rebuilding the tree or running EM again.

Use `prior=1` for uniform body centering; no parameters are required. For
`prior=2`, supply a positive rate via `prior_params` or leave it `NULL` for the
censored-data exponential MLE. All-right-censored samples do not have a finite
positive exponential MLE and require a supplied rate. `prior=3` retains Weibull
body centering with explicit `prior_params=c(rate,shape)` and scale `1/rate`.
Weibull is an extension of the supplied manuscript.

Every family uses root support `c(0,Inf)`, a body ending at `tau`, and one
unspecified tail. Predictions above `tau` are `NA`. Use a grid within `[0,tau]`
or the default 200-point grid when `predict()` is called without times. No cure fraction or parametric extrapolation is
implied by the tail mass.

## Control the reference prior

```r
custom <- survPT(dat, prior=1, s0=0.15, gamma=0.5,
                 min_count=1, max_M=8, verbose=FALSE)
custom$Prior[c("s0", "tail_mean")]
```

`s0` is the reference tail mass. If omitted it is calibrated from the body
boundary density with `kappa=0.1`; for uniform centering this gives `1/11`.
`gamma` sets depth-dependent Beta strength. `gamma=0` gives a flat prior.
Because each Beta parameter is floored at one, `s0` can differ from the actual
induced prior mean `tail_mean`. See [Methods](methods.md).

The tree grows using expected counts under the reference measure. `min_count`
controls pruning (zero disables count pruning); `max_M` limits ordinary depth.
The tail split is always retained and can add one level. The tree remains fixed
after construction.

By default a flat-prior fit initializes the working-prior MAP-EM. Each stage
has cap `max_EM_iter=2000` and tolerance `EM_tol=1e-6`. Inspect both stages in
`Convergence`; reaching the cap is not convergence. `max_EM_iter=0` returns the
reference initialization. `flat_init=FALSE` skips the preliminary fit.

## Add posterior means and uncertainty

```r
set.seed(2026)
pred <- predict(fit,x,0.95,n_sweeps=2000,burn_in=500,keep_draws=TRUE)

map_curve <- as.numeric(pred$Survival_Probs)  # deterministic MAP-EM evaluated at x
post <- pred$Posterior$Survival
post[, c("time", "mean", "lower", "upper", "ess", "mcse")]
pred$Posterior$Tail                         # posterior S(tau) summary
pred$Posterior$Draws$Survival                # retained sweeps x target times
```

`x` is a numeric vector of target times, not observation intervals or covariates.
`predict(fit,x,0.95)` runs sampling during the call and evaluates survival at
exactly those times. No preliminary sampler call is needed. Repeated, unsorted
and single times retain their supplied order. The input model is not changed.

DA-Gibbs draws leaf labels then Beta branches with integer counts. Its target
is the unrestricted posterior: `split_epsilon` applies only to MAP estimation.
`n_sweeps` includes burn-in. Thinning starts at `burn_in+1`; the default keeps
all 1500 post-burn-in sweeps. `keep_draws=FALSE` is the prediction default;
set it to `TRUE` to retain evaluated curve draws as well as summaries.

The `mean` column is the posterior mean of the curve, not mean event time.
`lower` and `upper` form pointwise equal-tailed credible intervals, not a
simultaneous confidence band. ESS and MCSE are approximate single-chain
diagnostics. Use larger runs, trace plots and independent seeds to investigate
mixing. The posterior conditions on the retained tree and fixed base parameters.

## Reuse samples or save the model

Every ordinary `predict(..., posterior=TRUE)` call runs a fresh chain from the
stored MAP fit. To compare different grids or levels with the same draws,
sample explicitly and pass the result:

```r
set.seed(2026)
draws <- samplePT(fit,n_sweeps=2000,burn_in=500)
p95 <- predict(fit,c(0.5,1,2),0.95,samples=draws)
p90 <- predict(fit,c(2,1,0.5,1),0.90,samples=draws)
saveRDS(fit,"fit.rds")
restored <- readRDS("fit.rds")
predict(restored,c(1,2),posterior=FALSE)
```

`samplePT()` samples tree probabilities independently of an evaluation grid;
`predict()` maps those sampled trees to the supplied target times. When `samples`
is supplied, no additional random draws are consumed. Models and samples are
plain R objects and can be restored in a new R session. `posterior=FALSE` always
gives deterministic prediction without sampling.

## Inspect nodes and deterministic outputs

The model fields `Left`, `Right`, `Level`, `Allocations`, `Posterior_Probs`, and
`Terminal_Nodes` retain their roles. `Posterior_Probs` contains MAP node masses.
`Allocations` are refreshed at the final masses.

The prediction result contains `Test_seq`, `Survival_Probs`, `Density`,
`Tail_Prob`, and `Posterior`. Only `Tail_Prob`, a fitted node mass, also belongs
to the model. The model contains no prediction or sampler output fields.

As of 0.3.2, `survPT()` no longer accepts `test_seq`, `n_sweeps`, `burn_in`,
`thin`, `credible_level`, `keep_draws`, or the fixed `support` argument.
Supply target times and interval level to `predict(fit, times, level)`;
sampler controls belong to `predict()` or `samplePT()`.

`Alpha` stores each child's Beta parameter; sibling entries form one Beta
distribution. `Branch_Probs` are fitted conditional branch probabilities.
`Parent`, `Left_Child`, `Right_Child` and `Split_Nodes` use one-based R indices;
absent links are `NA`. `Reference_Probs`, `Reference_Branch_Probs` and
`Initial_Allocations` expose construction quantities.

```r
stopifnot(abs(sum(fit$Posterior_Probs[fit$Terminal_Nodes])-1) < 1e-10)
nodes <- data.frame(left=fit$Left, right=fit$Right, depth=fit$Level,
                    alpha=fit$Alpha, mass=fit$Posterior_Probs)
```

Only terminal masses should sum to one. Raw depth-level sums need not equal one
on a pruned tree. Armadillo vectors may appear as one-column matrices; use
`as.numeric()` when needed. Native `SurvPT()` uses zero-based indices instead.

## Legacy options and generators

`allocated=0` remains a separate overlap-count estimator without EM or Gibbs
sampling. It now uses the revised partition, priors and unspecified-tail output.
It is not an alternative implementation of the same MAP estimator.

The simulation generators `data_rc`, `data_ic`, `data_cs` and `fiducial_1` through
`fiducial_4` are unchanged. See [Limitations](limitations.md) before using the
legacy current-status `status` field or inspection-count logic. Historical
bootstrap scripts are separate from the new posterior sampler.
