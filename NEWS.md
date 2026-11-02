# SurvPT 0.2.1

- Both Beta prior parameters at child depth `m` are `max(1, 0.5*m^2)`
  for every base and allocation mode, including the final tail split.
- All conditional PT split probabilities are constrained to `[1e-6, 1 - 1e-6]`
  in initialization and EM updates, including the final tail split and both
  allocation modes.
- All PT bases now use root support `[0, Inf]`. Ordinary splits use the base
  median within the node's intersection with `[0, largest_obs]`.
- The rightmost terminal interval receives a final split at `largest_obs`.
  For a parent at depth `m`, the new branch prior is `Beta(alpha, alpha)`
  with `alpha = max(1, 0.5*(m+1)^2)`, including after early stopping.
- Fractional allocation initializes every conditional branch probability to
  one half, including the final tail split. Initial and subsequent E-steps use
  the current leaf masses and conditional terminal-base distributions. The
  first count-based mode update occurs in the first EM iteration.
- Uniform tails hold a point mass at infinity, with observation factor one for
  right censoring and zero for finite intervals. Exponential and Weibull tails
  retain their conditional base distributions. The separate overlap-count
  option (`allocated = 0`) retains its previous initialization.
- Median splits and conditional probabilities use stable cumulative-hazard
  calculations for small probabilities and distant tails. Regression tests
  cover tree structure, allocation, endpoint conventions and probability mass.
