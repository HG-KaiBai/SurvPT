# User guide

[Back to README](../README.md)

## Encode observations

`survPT()` expects a numeric matrix or a data frame coercible to a numeric matrix. Its first two columns are `L` and `R`; later columns are ignored by the C++ estimator.

| Observation | `L` | `R` |
| --- | --- | --- |
| Exact event at `t` | `t` | `t` |
| Event in an interval | lower endpoint | upper endpoint |
| Right censored at `c` | `c` | `Inf` |
| Event occurred by inspection at `c` | `0` | `c` |
| Event not yet observed at inspection at `c` | `c` | `Inf` |

Intervals are interpreted as `(L, R]`, including zero when `L = 0`, with `L == R` handled separately as an exact event. Endpoint conventions have no effect on continuous interval probabilities. Both exact-event allocation and density at a positive partition boundary use the interval on its left; nodes starting at zero include zero.

Use finite nonnegative left endpoints and finite right endpoints or positive `Inf`. Require `L <= R`, no missing values, at least one row, and nonnegative exact event times. With a positive lower support bound, exact events must lie strictly above that bound. Input checking is incomplete, so validate these conditions before fitting. There is no delayed-entry or sampling-truncation interface.

```r
x <- data.frame(L = c(0, 0.4, 1.1, 2.0),
                R = c(0.8, 0.4, 2.0, Inf))
stopifnot(all(is.finite(x$L)), all(x$L >= 0),
          !anyNA(x), all(x$R >= x$L))
```

## Choose the PT base distribution

| `prior` | `prior_params` | Base survival function |
| --- | --- | --- |
| `1` | `NULL` (ignored) | Uniform within each effective terminal interval |
| `2` | positive rate `lambda` | `exp(-lambda * t)` |
| `3` | `c(lambda, k)`, both positive | `exp(-(lambda * t)^k)` |
| `4` | one finite positive rate `lambda` | Uniform body on `[0, M]`, exponential survival `exp(-lambda * t)` after `M = largest_obs` |

For Weibull calculations with R, use `shape = k, scale = 1 / lambda`. The base distribution controls partition locations and the shape inside a terminal interval; the tree estimates the interval masses.

```r
library(SurvPT)
d <- data_rc(seed = 123, n = 60, K = 1.5, Lambda = 0.4, q = 0.3)
grid <- seq(0, 12, length.out = 121)
rate <- estimate_exp(d$L, d$R)
weibull <- estimate(d$L, d$R)  # c(lambda, k)

fit_exp <- survPT(d[, c("L", "R")], prior = 2, prior_params = rate,
                  max_M = 6, verbose = FALSE)
fit_weibull <- survPT(d[, c("L", "R")], prior = 3,
                      prior_params = weibull,
                      max_M = 6, verbose = FALSE)
```

These plug-in parameters are held fixed by the tree fit. Their estimation uncertainty is not propagated into the returned curve. The optimization helpers return parameter estimates only; they do not return convergence diagnostics or standard errors.

For `prior = 4`, use `survPT(d[, c("L", "R")], prior = 4, prior_params = rate)`.
The rate must be supplied; the function does not estimate it automatically.
The body receives probability `1 - exp(-rate * M)` and is uniform conditional on `[0, M]`.
The tail receives probability `exp(-rate * M)` and, conditional on exceeding `M`,
has distribution `M + Exp(rate)`. The base CDF is continuous, although the density
can jump at `M`. A finite upper support conditions and truncates this base as usual.
The result includes `Prior_Params` for every prior and `Tail_Start` for this prior only.

## Choose support and evaluation points

The default `support = c(0, Inf)` permits an unbounded time axis, but the **uniform PT output has an effective finite tail cutoff** at `uniform_tail_factor * max(finite observed endpoints)`. The multiplier defaults to `1.1` and must be a finite number greater than one. It controls uniform partitioning, initialization, EM allocation, and prediction. Priors `2`, `3`, and `4` retain infinite tails and ignore this multiplier. See [Methods](methods.md) for the conditional distributions within terminal intervals.

```r
fit_uniform <- survPT(x, prior = 1, uniform_tail_factor = 1.5,
                      verbose = FALSE)
```

Use a finite upper support only when a finite event-time bound is part of the model. In that case, keep observation intervals within that support, including replacing an infinite observation endpoint with the known bound. Do not discard observations merely to force a finite support. Initial allocation can fail to conserve counts for intervals extending beyond the declared support.

`survPT()` does not accept `test_seq` or construct a grid. Supply the vector of evaluation times to `predict(fit, x)`. The same fitted model can be evaluated on multiple grids without refitting. Unsorted times and duplicates are preserved; use a sorted finite grid for plotting.

## Control tree size and iteration

- `max_M` is the maximum depth, with the root at level zero. Start with a modest nonnegative integer, such as 6 or the default 8. Large depths can exhaust memory or exceed integer-shift limits.
- A node also stops splitting when its **initial allocated count is below one**. The partition is fixed during subsequent updates.
- `max_EM_iter` and `EM_tol` control allocation iterations.
- `flat = TRUE` uses unit Beta prior parameters at all depths; `FALSE` (default) uses `max(1, 0.5 * m^2)` at depth `m`.
- The stopping criterion is the Euclidean norm of the change in terminal masses. It is not a likelihood or parameter-distance criterion.
- `verbose = TRUE` prints progress. The returned list contains no convergence flag, iteration count, or objective history.
- PT `allocated = 0` selects the legacy overlap-count initialization and skips EM. It is a different estimator, not a speed setting for the default fit.

## Interpret the result

| Field | Meaning |
| --- | --- |
| `Left`, `Right` | Node interval boundaries; the last right boundary may still be `Inf` even when uniform PT uses a finite effective cutoff |
| `Level` | Depth of each node; root is zero |
| `Alpha`, `Parent` | Prior Beta parameters and one-based parent indices (zero for the root) |
| `Allocations` | Fractional counts assigned to nodes at the returned iteration |
| `Posterior_Probs` | Estimated node masses for both internal and terminal nodes |
| `Terminal_Nodes` | **One-based R indices** of terminal nodes |
| `X` | Observation endpoints retained for posterior sampling |
| `Prior`, `Prior_Params`, `Support`, `Largest_Obs` | Base distribution settings needed for prediction |

The fit is an S3 object of class `survPT` and can be saved with `saveRDS()`. Numeric arrays from RcppArmadillo may have one column; `as.numeric()` removes those dimensions. Sum only terminal masses when checking total probability:

```r
idx <- fit_exp$Terminal_Nodes
stopifnot(abs(sum(fit_exp$Posterior_Probs[idx]) - 1) < 1e-8)
set.seed(123)
summary_table <- predict(fit_exp, grid, level = 0.95)
head(summary_table)
# Point estimates only, with no posterior simulation:
predict(fit_exp, grid, sample = FALSE)
# Density predictions are also available:
predict(fit_exp, grid, type = "density", sample = FALSE)
```

The prediction columns are `x`, `estimate`, `lower`, `upper`, `mean`, and `median`. `estimate` is the fitted MAP-EM distribution evaluated at `x`, using the same formula as the former `Survival_Probs` output. It is not a marginal posterior mode or a Monte Carlo mean. The other columns summarize the posterior of `S(x)` (or density for `type = "density"`), not event-time means or medians. Intervals are pointwise and equal-tailed, conditional on the fitted partition and base parameters.

Each call uses 3,000 burn-in sweeps and 5,000 retained draws by default (8,000 total sweeps with the default thinning). Change `burn_in`, `n_samples`, and `thin` as needed. Set a seed for reproducibility; `keep_samples = TRUE` attaches the draws as `attr(result, "samples")`, with draws in rows and times in columns. Inspect mixing and Monte Carlo stability for substantive analyses. For `sample = FALSE`, only `x` and `estimate` are returned. The model is not modified by prediction.

The field `Posterior_Probs` stores the fitted node masses; posterior samples are generated separately by `predict()`. [Methods](methods.md) describes both updates.

## Simulation helpers

`data_rc()`, `data_ic()`, and `data_cs()` simulate Weibull event times with inverse-scale parameter `Lambda` and shape `K`. Use `0 < q < 1`; `q` determines a censoring or inspection-time distribution, not a guaranteed realized censoring proportion. Each function resets R's random seed.

**Current-status caveat:** `data_cs()` currently overwrites `status` with `as.numeric(L == R)`. It is not the event-by-inspection indicator. Use its `L` and `R` columns for PT; for positive inspection times the inspection indicator can be recovered as `as.integer(L == 0 & is.finite(R))`. Do not use that indicator as a right-censoring `delta`.

`fiducial_1()` through `fiducial_4()` generate fixed simulation scenarios; they do not fit a fiducial estimator. Their distributions and the associated reference are documented in their help pages. The experimental estimator in `Examples/Fiducial.R` is separate from the package API.
