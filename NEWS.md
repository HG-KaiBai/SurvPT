# SurvPT 0.3.2

- Bound all PT split probabilities to `[1e-6, 1 - 1e-6]` during
  initialization and every EM update.
- Use `Beta(1, 1)` for the body/tail split and the finite-subtree schedule
  `alpha_m = max(1, 0.5 * m^2)`, including bounded support. The first finite
  split uses alpha 1; deeper levels retain the original `0.5 * m^2` schedule.
- Update regression checks and the tail diagnostic for the new priors and
  probability bounds.

# SurvPT 0.3.1

- Uniform-base survival at finite times beyond tau now stays at S(tau)
  instead of returning NA. Tail density remains NA.
- Exponential-base tail survival and density now use the conditional
  exponential distribution beyond tau with the supplied rate. Tail mass
  continues to update during EM. Shifted-time calculations avoid dividing
  by an underflowed base survival probability at tau.
- Regression tests cover both public interfaces, allocation modes, missing
  evaluation times, tail boundaries, density integrals and large tau.
  Weibull tail predictions and bounded-support fits remain unchanged.

# SurvPT 0.3.0

- Unbounded-support PT fits now split first at the largest finite observed
  endpoint into a finite subtree and an unresolved tail. The first split has
  an internal Beta(0.5, 0.5) prior and is updated during EM.
- The finite subtree keeps its existing depth budget and alpha schedule.
  Exact events at partition boundaries, including the upper finite endpoint,
  use the interval on the left.
- All finite leaves use the selected conditional base distribution for both
  censored allocation and prediction. The tail retains probability mass
  without imposing a density; predictions at finite times beyond the largest
  endpoint are NA.
- EM caches observation/leaf probabilities and aggregates them up the fixed
  tree, reducing repeated probability calculations without changing the
  branch update rules.
- Fits report the finite endpoint, observed log-likelihood history and
  terminal-mass convergence diagnostics. Default grids include right-censoring
  times. Existing fitting arguments remain unchanged.
- Examples include uniform and exponential PT comparisons; regression tests
  cover conditional bases, the first-layer Beta update and tail boundaries.
