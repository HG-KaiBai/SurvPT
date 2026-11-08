# User guide

[Back to README](../README.md)

## Encode observations and fit

Supply `L,R` as the first two columns of a numeric matrix or data frame. Exact
events use `L=R`, interval censoring uses finite `L<R`, right censoring uses
`R=Inf`, and left censoring uses `L=0`. Require finite nonnegative left endpoints,
ordered endpoints, no missing values, and a positive finite endpoint.

```r
library(SurvPT)
dat <- data.frame(L=c(0,0.4,1.1,2,3), R=c(0.8,0.4,2,Inf,Inf))
fit <- survPT(dat, prior=2, max_M=8, verbose=FALSE)
fit$Prior$parameters
fit$Convergence$map
fit$Tail_Prob
x <- seq(0,fit$Tau,length.out=101)
pred <- predict(fit,x,sample=FALSE)
pred                        # x and estimate
attr(pred,"tail")           # x = Tau; estimate = fitted S(Tau)
if (interactive()) plot(pred$x,pred$estimate,type="l",ylim=c(0,1))
```

`survPT()` returns a serializable fitted model without evaluating a grid or
sampling. `predict.survPT()` handles both estimation and posterior sampling.
There is no standalone `samplePT()` interface.

Use `prior=1` for uniform body centering. For `prior=2`, supply a positive
exponential rate via `prior_params`, or leave it `NULL` to estimate the
censored-data MLE once. All-right-censored data require an explicit rate.
For Weibull centering, use `prior=3, prior_params=c(rate,shape)` with scale `1/rate`.

Every family uses `[0,Inf]`, a body ending at the largest finite endpoint `Tau`,
and one tail whose shape is unspecified. For times above `Tau`, including `Inf`,
survival and density estimates and all posterior summary columns are `NA`.
At `Tau`, survival equals the fitted or sampled tail mass; density retains
the body endpoint convention.

## Generated data and population truth

`data_rc()`, `data_ic()`, `data_cs()`, `data_mixture()` and `fiducial_1()` through
`fiducial_4()` return their usual data frames with two scalar attributes and a
named proportion vector:

```r
d <- data_mixture(seed=123, n=200, type="interval", q=0.5)
c(tau=attr(d,"tau"), true_S_tau=attr(d,"true_S_tau"))
attr(d,"censoring_proportions")
```

`tau` is the largest finite endpoint among `L` and `R`. `true_S_tau` is the
survival probability from the generating population distribution at that time;
it is not a fitted estimate. `censoring_proportions` contains the observed
proportions in the fixed order `Exact`, `Left`, `Interval`, `Right`:

- `Exact`: finite `L == R`.
- `Left`: `L == 0 < R < Inf`.
- `Interval`: `0 < L < R < Inf`.
- `Right`: `R == Inf`, including `(0, Inf)`.

Classification uses only the endpoints and the four proportions sum to one.
The attributes describe the original generated dataset and are not recomputed
when rows or endpoints are changed.

## Select the prior

```r
regular <- survPT(dat, prior=2, gamma=0.5, verbose=FALSE)
flat <- survPT(dat, prior=2, flat=TRUE, verbose=FALSE)
```

The default `flat=FALSE` uses the 0.4.0 depth- and reference-dependent alpha
parameters. `flat=TRUE` selects Beta(1,1) on every branch for both EM and Gibbs,
equivalent to `gamma=0`. It does not change the centering distribution or tree.
Both models begin allocation at the centering reference distribution and run
one MAP-EM fit. There is no preliminary flat-prior optimization.

`s0` is the reference tail mass, calibrated from the body boundary density with
`kappa=0.1` if omitted. Uniform centering gives `s0=1/11`. Because Beta parameters
are floored at one, `Prior$s0` can differ from the induced prior tail mean
`Prior$tail_mean`. `min_count` controls pruning using reference expected counts;
`max_M` limits ordinary depth. The final tail split is always retained.

The fit uses at most `max_EM_iter=2000` updates and stops when the leaf-mass
change is at most `EM_tol=1e-6`. Inspect `Convergence$map`; reaching the iteration
cap does not imply convergence. With `allocated=1`, `max_EM_iter=0` returns
the reference initialization and its allocations.

At child level `m`, fitted branch probabilities use bounds
`[split_epsilon * 2^(-m), 1 - split_epsilon * 2^(-m)]`, with base scale `1e-6`.
The root has level zero. Its children use `5e-7`; level-two children use `2.5e-7`.
This also applies to the final tail split and the overlap-count fit.

## Estimate curves and sample the posterior

```r
set.seed(2026)
pred <- predict(fit,x,level=0.95,n_samples=2000,burn_in=3000,thin=1)
head(pred)                     # x, estimate, lower, upper, mean, median
attr(pred,"tail")              # same columns, for S(Tau)
density <- predict(fit,x,type="density",sample=FALSE)
tail_curve <- predict(fit,c(fit$Tau,2*fit$Tau,Inf),sample=FALSE)
# estimates: fit$Tail_Prob, NA, NA
```

The prediction interface follows 0.2.3. `estimate` is the fitted MAP-EM plug-in
curve. `mean` and `median` summarize its posterior; they are not mean or median
event times. `lower` and `upper` are pointwise equal-tailed credible limits.
For density predictions, the separate `tail` attribute still reports a probability.
It is available even when the target grid omits `Tau`, or when `x=numeric()`.
No curve is evaluated beyond `Tau`. Rows for these times remain in the output,
with `NA` in `estimate`, `lower`, `upper`, `mean` and `median`. The separate
tail summary remains available even when all requested times exceed `Tau`.

The defaults use `n_samples=2000` draws after `burn_in=3000` sweeps, with
`thin=1`. The total is `burn_in+n_samples*thin`; the first retained sweep is
`burn_in+thin`. Each call starts at the stored fitted masses. Censored labels
are drawn afresh, integer counts are aggregated, and conjugate Beta branches
are sampled without the MAP clipping bound. The prior stored in the fit is used.
Only curve and tail summaries are returned; raw draws are not retained.

The input model is unchanged. Repeated and unsorted times keep their order;
posterior tree sampling does not depend on the prediction grid. Reset the same
seed to compare grids or credible levels with the same draws. `sample=FALSE`
or `n_samples=0` consumes no random draws and returns only `x` and `estimate`,
plus the deterministic `tail` attribute. Saved models can still be restored
with `readRDS()`; examples do not automatically save models or results.

Intervals condition on the retained tree and fixed centering parameters. They
are not simultaneous bands. Compare summaries across independent runs and
longer chains to assess Monte Carlo stability.

## Inspect fitted state

`Left`, `Right`, `Level`, `Allocations`, `Posterior_Probs` and `Terminal_Nodes`
retain their roles. `Posterior_Probs` contains fitted node masses, not posterior
mean masses. Final allocations are refreshed at the returned masses. `Alpha`
stores each child's Beta parameter; siblings form a pair. `Branch_Probs`
contains fitted conditional probabilities. Tree links use one-based indices.

`Reference_Probs`, `Reference_Branch_Probs` and `Initial_Allocations` describe
the centering reference used to construct the partition. Only terminal masses
sum to one on a pruned tree. Some native vectors appear as one-column matrices;
use `as.numeric()` if needed. Prediction columns are ordinary numeric vectors.

`allocated=0` retains the historical overlap-count fit without EM or Gibbs;
use `sample=FALSE` for its predictions. Simulation and parametric helpers are
unchanged. See [Methods](methods.md) and [Limitations](limitations.md).
