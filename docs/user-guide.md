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

Use finite nonnegative left endpoints and finite right endpoints or positive `Inf`. Require `L <= R`, no missing values, at least one row, and at least one positive finite endpoint. Input checking is incomplete, so validate these conditions before fitting. There is no delayed-entry or sampling-truncation interface.

```r
x <- data.frame(L = c(0, 0.4, 1.1, 2.0),
                R = c(0.8, 0.4, 2.0, Inf))
stopifnot(all(is.finite(x$L)), all(x$L >= 0),
          !anyNA(x), all(x$R >= x$L))
```

## Choose the PT base distribution

| `prior` | `prior_params` | Base survival function |
| --- | --- | --- |
| `1` | `NULL` (ignored) | Uniform within finite terminal intervals, plus a point mass at infinity |
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

All bases require `support = c(0, Inf)`. Ordinary splits use the base distribution truncated to the observed range, and a final split at the largest finite endpoint separates the tail. The **uniform tail is a point mass at infinity**: finite-time density is zero beyond that endpoint, and survival stays at the fitted tail mass. Exponential and Weibull tails retain their continuous conditional distributions. See [Methods](methods.md) for the partition and terminal distributions.

Finite upper support and positive lower support bounds are not accepted by this tree construction. Keep the original censoring endpoints in the data.

Supply `test_seq` explicitly for reproducible comparisons. The PT default grid uses 200 points from the lower support bound to the largest finite endpoint in the first two columns of the data, including right-censoring times. It can therefore be constructed from an entirely right-censored sample with finite censoring times. For plotting, choose a sorted finite grid. Fitting is unchanged by the choice of evaluation grid.

## Control tree size and iteration

- `max_M` is the maximum ordinary depth, with the root at level zero. The final tail split adds children one level deeper than its parent, even if `max_M = 0`. Start with a modest nonnegative integer, such as 6 or the default 8. Large depths can exhaust memory.
- A node also stops splitting when its **initial allocated count is below one**. The partition is fixed during subsequent updates.
- `max_EM_iter` and `EM_tol` control allocation iterations.
- The stopping criterion is the Euclidean norm of the change in terminal masses. It is not a likelihood or parameter-distance criterion.
- `verbose = TRUE` prints progress. The returned list contains no convergence flag, iteration count, or objective history.
- With a uniform base, every right-censored observation starts entirely in the tail in either allocation mode. In fractional mode, finite interval-censored observations start with overlap-length weights; EM can subsequently move right-censored weight into compatible finite nodes.
- PT `allocated = 0` uses overlap counts for finite intervals (and all censored intervals for exponential/Weibull bases), and skips EM. It is a different estimator, not a speed setting for the default fit.

## Interpret the result

| Field | Meaning |
| --- | --- |
| `Left`, `Right` | Node interval boundaries; the final tail starts at the largest finite endpoint and ends at `Inf` |
| `Level` | Depth of each node; root is zero |
| `Allocations` | Fractional counts assigned to nodes at the returned iteration |
| `Posterior_Probs` | Estimated node masses for both internal and terminal nodes |
| `Terminal_Nodes` | **One-based R indices** of terminal nodes |
| `Test_seq` | Evaluation points supplied to the fit |
| `Survival_Probs` | Estimated `P(T > t)` at those points |
| `Density` | Estimated density at those points |

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

## Simulation helpers

`data_rc()`, `data_ic()`, and `data_cs()` simulate Weibull event times with inverse-scale parameter `Lambda` and shape `K`. Use `0 < q < 1`; `q` determines a censoring or inspection-time distribution, not a guaranteed realized censoring proportion. Each function resets R's random seed.

**Current-status caveat:** `data_cs()` currently overwrites `status` with `as.numeric(L == R)`. It is not the event-by-inspection indicator. Use its `L` and `R` columns for PT; for positive inspection times the inspection indicator can be recovered as `as.integer(L == 0 & is.finite(R))`. Do not use that indicator as a right-censoring `delta`.

`fiducial_1()` through `fiducial_4()` generate fixed simulation scenarios; they do not fit a fiducial estimator. Their distributions and the associated reference are documented in their help pages. The experimental estimator in `Examples/Fiducial.R` is separate from the package API.
