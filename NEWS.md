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
