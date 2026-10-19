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

Intervals are interpreted as `(L, R]`, with `L == R` handled separately as an exact event. Endpoint conventions have no effect on continuous interval probabilities, but exact events at a partition boundary belong to the interval on its left.

Use finite nonnegative left endpoints and finite right endpoints or positive `Inf`. Require `L <= R`, no missing values, at least one row, and positive exact event times. Input checking is incomplete, so validate these conditions before fitting. There is no delayed-entry or sampling-truncation interface.

```r
x <- data.frame(L = c(0, 0.4, 1.1, 2.0),
                R = c(0.8, 0.4, 2.0, Inf))
stopifnot(all(is.finite(x$L)), all(x$L >= 0),
          !anyNA(x), all(x$R >= x$L))
```

`survOPT()` accepts a vector `X` of finite observed follow-up times and an equally long `delta` vector, with `1` for an exact event and `0` for right censoring. Do not pass an interval matrix or current-status indicator to this interface.

## Choose the PT base distribution

| `prior` | `prior_params` | Base survival function |
| --- | --- | --- |
| `1` | `NULL` (ignored) | Uniform within each effective terminal interval |
| `2` | positive rate `lambda` | `exp(-lambda * t)` |
| `3` | `c(lambda, k)`, both positive | `exp(-(lambda * t)^k)` |

For Weibull calculations with R, use `shape = k, scale = 1 / lambda`. The base distribution controls partition locations and the shape inside a terminal interval; the tree estimates the interval masses.

```r
library(SurvPT)
d <- data_rc(seed = 123, n = 60, K = 1.5, Lambda = 0.4, q = 0.3)
grid <- seq(0, 12, length.out = 121)
rate <- estimate_exp(d$L, d$R)
weibull <- estimate(d$L, d$R)  # c(lambda, k)

fit_exp <- survPT(d[, c("L", "R")], prior = 2, prior_params = rate,
                  test_seq = grid, max_M = 6, verbose = FALSE)
fit_weibull <- survPT(d[, c("L", "R")], prior = 3,
                      prior_params = weibull, test_seq = grid,
                      max_M = 6, verbose = FALSE)
```

These plug-in parameters are held fixed by the tree fit. Their estimation uncertainty is not propagated into the returned curve. The optimization helpers return parameter estimates only; they do not return convergence diagnostics or standard errors.

## Choose support and evaluation points

The default `support = c(0, Inf)` permits an unbounded time axis, but the **uniform PT output has an effective finite tail cutoff** at `1.1 * max(finite observed endpoints)`. Exponential and Weibull PT outputs retain infinite tails. See [Methods](methods.md) for OPT's different infinite-tail behavior.

Use a finite upper support only when a finite event-time bound is part of the model. In that case, keep observation intervals within that support, including replacing an infinite observation endpoint with the known bound. Do not discard observations merely to force a finite support. Initial allocation can fail to conserve counts for intervals extending beyond the declared support.

Supply `test_seq` explicitly for reproducible comparisons. The PT default grid uses only rows with a finite right endpoint, so it may omit late right-censoring times and cannot be constructed from an entirely right-censored sample. For plotting, choose a sorted finite grid. Fitting is unchanged by the choice of evaluation grid.

## Control tree size and iteration

- `max_M` is the maximum depth, with the root at level zero. Start with a modest nonnegative integer, such as 6 or the default 8. Large depths can exhaust memory or exceed integer-shift limits.
- A node also stops splitting when its **initial allocated count is below one**. The partition is fixed during subsequent updates.
- PT `max_EM_iter` and `EM_tol` control allocation iterations; OPT uses `max_iter` and `iter_tol`.
- The stopping criterion is the Euclidean norm of the change in terminal masses. It is not a likelihood or parameter-distance criterion.
- `verbose = TRUE` prints progress. The returned list contains no convergence flag, iteration count, or objective history.
- PT `allocated = 0` selects the legacy overlap-count initialization and skips EM. It is a different estimator, not a speed setting for the default fit.

## Interpret the result

| Field | Meaning |
| --- | --- |
| `Left`, `Right` | Node interval boundaries; the last right boundary may still be `Inf` even when uniform PT uses a finite effective cutoff |
| `Level` | Depth of each node; root is zero |
| `Allocations` | Fractional counts assigned to nodes at the returned iteration |
| `Posterior_Probs` | Estimated node masses; for PT, internal and terminal nodes are populated; for OPT, only terminal nodes are populated |
| `Terminal_Nodes` | **One-based R indices** of terminal nodes |
| `Test_seq` | Evaluation points supplied to the fit |
| `Survival_Probs` | Estimated `P(T > t)` at those points |
| `Density` | Estimated density at those points |
| `Rho` | OPT only: node stopping probabilities; storage can include unused trailing entries, so use the first `length(Left)` entries |

The user-facing PT wrapper does not return `Alpha`. Numeric arrays from RcppArmadillo may have one column; `as.numeric()` removes those dimensions. Sum only terminal masses when checking total probability:

```r
idx <- fit_exp$Terminal_Nodes
stopifnot(abs(sum(fit_exp$Posterior_Probs[idx]) - 1) < 1e-8)
summary_table <- data.frame(
  time = fit_exp$Test_seq,
  survival = as.numeric(fit_exp$Survival_Probs),
  density = as.numeric(fit_exp$Density)
)
head(summary_table)
```

The field name `Posterior_Probs` does not imply posterior simulation or uncertainty intervals. PT currently uses a clipped, mode-style branch-mass update; [Methods](methods.md) describes the formula.

## Right-censored OPT example

```r
opt <- survOPT(d$L, d$status, test_seq = grid,
               max_M = 6, max_iter = 100, verbose = FALSE)
sum(opt$Posterior_Probs[opt$Terminal_Nodes])
```

Here `d` comes from `data_rc()`, whose `status` is an event indicator. Cite Seok, Tian, and Wong (2014), [doi:10.1093/biostatistics/kxt025](https://doi.org/10.1093/biostatistics/kxt025), when using OPT. The interface is univariate and right-censored; the cited paper has a broader multivariate scope.

## Simulation helpers

`data_rc()`, `data_ic()`, and `data_cs()` simulate Weibull event times with inverse-scale parameter `Lambda` and shape `K`. Use `0 < q < 1`; `q` determines a censoring or inspection-time distribution, not a guaranteed realized censoring proportion. Each function resets R's random seed.

**Current-status caveat:** `data_cs()` currently overwrites `status` with `as.numeric(L == R)`. It is not the event-by-inspection indicator. Use its `L` and `R` columns for PT; for positive inspection times the inspection indicator can be recovered as `as.integer(L == 0 & is.finite(R))`. Do not use that indicator as a right-censoring `delta`.

`fiducial_1()` through `fiducial_4()` generate fixed simulation scenarios; they do not fit a fiducial estimator. Their distributions and the associated reference are documented in their help pages. The experimental estimator in `Examples/Fiducial.R` is separate from the package API.
