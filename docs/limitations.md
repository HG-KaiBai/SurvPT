# Known limitations and troubleshooting

[Back to README](../README.md)

## Interpretation

- Fitting estimates only tail mass `S(tau)`. Prediction adds an exponential
  shape with rate `-log(p)/tau` for each fitted or sampled mass `p`. Outputs
  beyond `tau` depend on this assumption; the data do not identify the tail
  shape. Rate uncertainty is induced by the mass draws, not separately sampled.
- DA-Gibbs targets the posterior on the retained tree with fixed base
  parameters. It does not integrate over tree selection or uncertainty in an
  estimated body-centering rate. The Weibull body center is a software extension
  of the manuscript's uniform/exponential choices.
- Reported credible intervals are pointwise, not simultaneous bands. Finite
  Gibbs runs have Monte Carlo error; compare independent runs and longer chains.
  Summary stability alone does not certify convergence.
- MAP-EM may have multiple stationary points. Allocation starts at the centering
  reference; monotone objective ascent alone does not prove global optimality.
  Inspect `Convergence$map` and its iteration cap.
- Beta parameters are floored at one. `Prior$s0` is the reference tail target;
  `Prior$tail_mean` is the actual mean induced by the floored prior, and can differ.
- The model is univariate with no covariates, delayed entry or sampling-truncation
  likelihood. All bases use support `c(0,Inf)`.
- Finite precision still limits extreme rates, shapes and observation scales.
  Unrepresentable automatic tail calibration raises an error requesting explicit
  `s0`. At zero, Weibull body density with shape below one is mathematically
  infinite; survival and allocation remain well-defined.

## Compatibility with 0.2.1

Version 0.3.0 intentionally changes the reference construction, tail-path priors,
initialization and tails. Old numerical results are not expected to be reproduced.
It preserves deterministic result names and provides posterior summaries
separately. Exponential and Weibull no longer extrapolate a parametric tail;
uniform no longer treats the tail as an atom at infinity.
Version 0.4.1 adds a prediction-only exponential tail matched to the fitted or
sampled mass, while retaining the unspecified tail in fitting.

Version 0.3.2 makes `survPT()` exclusively fit-only. Use
`predict(fit,times,level)` for curves and uncertainty: it runs sampling during
prediction unless `sample=FALSE` or `n_samples=0`. The
previous combined arguments and fixed `support` argument have been removed. Model
lists saved by versions before 0.3.1 must be refitted to obtain serializable
state for the new prediction and sampling methods.

`allocated=0` is retained as a separate historical overlap-count estimator with
the new priors and partition. It does not support posterior sampling.

## Numerical and package use

**No finite exponential MLE:** all-right-censored data, or data whose likelihood
is maximized at infinite rate, require an explicit positive `prior_params` when
using exponential centering. The automatic estimator no longer uses an arbitrary
rate floor of 0.001.

**Predictions beyond the data:** exponential extrapolation preserves the fitted
or sampled tail mass. For `x>tau`, it uses `S(x)=p^(x/tau)`. Survival and density
at `Inf` are zero; `attr(prediction,"tail")` still reports `S(tau)`. A mass of
exactly one cannot determine a positive matching rate and raises an error.

**Slow or large runs:** lower tree depth, raise `min_count`, reduce grid size, or
use fewer draws. Prediction returns summaries only, but temporarily stores
draws internally to calculate them. Coverage and prediction caches each have
a 64 MiB cap; native sample arrays have a 512 MiB limit before R overhead.

**Loaded old version:** restart R and use `devtools::load_all()` for the checkout,
or install it and verify `find.package("SurvPT")` and `packageVersion("SurvPT")`.
Build dependencies are described in [Installation](installation.md).

## Unchanged historical research scripts

- `data_cs()` sets `status` using `L==R`; fit its interval columns rather than
  treating that field as the event-by-inspection indicator.
- `data_ic()` retains its historical `1:Ri` inspection loop, which needs care
  when `Ri=0`.
- The Weibull `estimate()` helper retains its existing numerical optimization;
  it returns only parameters and requires separate diagnostic assessment.
- Several scripts launch large experiments. Start from `Examples/quickstart.R`
  and `Examples/posterior.R`. Results stay in R; plots display interactively.
  Grid points beyond each fit's `tau` use the exponential extrapolation.

Run `tests/pt-exponential-tail.R` for tail prediction checks.
