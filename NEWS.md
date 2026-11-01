# SurvPT 0.2.1

- All PT bases now use root support `[0, Inf]`. Ordinary splits use the base
  median within the node's intersection with `[0, largest_obs]`.
- The rightmost terminal interval receives a final split at `largest_obs`.
  For a parent at depth `m`, the new branch prior is
  `Beta(0.5*(m+1)^2, 0.5*(m+1)^2)`, including after early stopping.
- Uniform tails now hold a point mass at infinity. All uniform right-censored
  observations start entirely in the tail; subsequent EM steps can redistribute
  them into compatible finite nodes. Exponential and Weibull tails retain
  their conditional base distributions.
- Median splits and conditional probabilities use stable cumulative-hazard
  calculations for small probabilities and distant tails. Regression tests
  cover tree structure, allocation, endpoint conventions and probability mass.
