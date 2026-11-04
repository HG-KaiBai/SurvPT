# SurvPT 0.4.1

- Keeps the infinite support, separate unspecified tail, reference-based
  partition/pruning and asymmetric alpha construction from 0.4.0.
- Starts allocation from the centering reference and runs one MAP-EM fit.
  Replaces automatic flat-prior initialization with `flat=FALSE/TRUE`; a flat
  model uses Beta(1,1) throughout both fitting and posterior sampling.
- Uses two PT interfaces: `survPT()` fits and `predict.survPT()` evaluates and
  samples. Removes the standalone `samplePT()` and public native fit entry point.
- Uses the 0.2.3 prediction arguments (`sample`, `n_samples`, `burn_in`,
  `thin`, `type`) and data-frame output, without `keep_samples` or raw draw
  attributes. Tail summaries appear in `attr(prediction, "tail")`.
- Defaults to `n_samples=2000` and `burn_in=3000` during prediction.
- Uses level-dependent MAP bounds: at child level `m`, the lower bound is
  `split_epsilon * 2^(-m)` (base scale `1e-6`) and the upper bound is one minus
  that value, including the final tail split. Gibbs draws remain unrestricted.
- Posterior retention uses `burn_in + n_samples*thin` sweeps and keeps exactly
  `n_samples` draws. Predictions beyond `Tau` remain `NA`.
- Updates examples to the 0.2.3 workflow; all results and plots stay in the R
  session. No CSV, RData or image files are automatically saved.
- Adds regression checks for centering initialization, one-stage EM, flat-prior
  sampling, analytic and numerical posterior targets, tail outputs and interfaces.

# SurvPT 0.4.0

- Releases the separated PT workflow: `survPT()` fits a serializable model,
  `predict(fit, times, level)` samples and predicts at target times, and
  `samplePT()` provides reusable posterior tree draws.
- Retains deterministic MAP-EM curves alongside posterior means and pointwise
  credible intervals, with prior, tail, EM-convergence and sampling diagnostics.
- Updates research examples to use the fitting and prediction interfaces.
  `Examples/posterior.R` compares uniform/exponential PT, fiducial, MPT and TBP
  on the same simulated dataset and writes curves, error tables and diagnostics.
- Shares MPT/TBP fitting helpers between the comparison examples and updates
  version metadata, citation information and local installation instructions.

# SurvPT 0.3.2

- Makes `survPT()` exclusively a fitting function. Removes `test_seq`,
  `n_sweeps`, `burn_in`, `thin`, `credible_level`, `keep_draws`, and `support`
  (the support remains fixed at `[0, Inf]`). Supply times, level and sampling
  controls to `predict()`; use `samplePT()` for standalone draws.
- Removes prediction placeholders and combined-interface delegation from fitted
  models. MAP node masses, priors, tail probability and convergence stay available.
- Consolidates the native fitting entry points into fit-only `SurvPT()` and
  removes grid evaluation from fitting and sampling. Prediction alone evaluates
  deterministic and sampled curves at requested times.
- Updates examples, documentation and regression checks to use separate calls.

# SurvPT 0.3.1

- Separates model fitting, sampling and prediction. `survPT()` now returns a
  serializable S3 model without evaluating a grid by default.
- `predict(fit, times, 0.95)` treats `times` as a numeric target-time vector,
  runs fresh posterior sampling during the call, and returns deterministic
  `Survival_Probs` alongside `Posterior` at exactly those times. It never refits.
- `posterior=FALSE` gives deterministic prediction without sampling. Unsorted,
  repeated and single target times retain their supplied order and dimensions.
- `samplePT(fit)` exposes standalone grid-independent DA-Gibbs draws.
  `predict(..., samples=draws)` explicitly reuses them across grids and levels.
- Saved models and draws work after `saveRDS()`/`readRDS()`. Prediction and sampling
  reconstruct the retained state without rebuilding the tree or changing MAP fits.
- Previous `survPT(test_seq=..., n_sweeps=...)` calls remain convenience wrappers
  around the separated operations. Default prediction draws are not retained;
  use `keep_draws=TRUE` to return them.

# SurvPT 0.3.0

- Implements Algorithms 1 and 2 of the supplied SCPT manuscript: reference-based
  pruning, flat-prior initialization, working-prior MAP-EM, and DA-Gibbs sampling.
- Beta parameters now depend on the extended body/tail reference measure:
  `max(1, 2*gamma*m^2*r)` and `max(1, 2*gamma*m^2*(1-r))`. Defaults are
  `gamma=0.5`, `kappa=0.1`, and `min_count=1`. The user can specify `s0` directly.
- All families now leave the tail shape unspecified. `Tail_Prob` estimates
  survival at the largest finite endpoint; survival and density beyond it are
  `NA`. The former uniform atom and parametric tail extrapolations are removed.
- Deterministic result fields are preserved. Set `n_sweeps=2000` (with default
  `burn_in=500`) for separate posterior means, pointwise credible intervals,
  raw draws, and approximate single-chain ESS/MCSE in `Posterior`.
- Sampling draws discrete labels and unrestricted conjugate Beta branches;
  it does not reuse fractional EM counts as fixed posterior counts.
- `Convergence` reports both EM stages and objective histories. Final allocations
  are refreshed at the returned fitted masses. Default cap is 2000 per stage.
- Exponential centering can estimate its rate automatically using a stable
  likelihood-score root, without the former artificial rate floor of 0.001.
- Weibull body centering remains available as an extension of the manuscript.
  `allocated=0` remains a separate overlap-count estimator without sampling.
- Adds analytic-posterior, numerical-integration, conservation, tail, pruning,
  reproducibility, and MAP-ascent regression checks in `tests/scpt.R`.

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
