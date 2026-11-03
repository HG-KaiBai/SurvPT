# Implementation and statistical interpretation

[Back to README](../README.md)

Version 0.4.0 implements Sections 2.2-2.5 and Algorithms 1-2 of the supplied
manuscript *Self-Consistent Polya Trees for Arbitrarily Censored Survival Data*.
The manuscript's uniform and exponential body centers are supported. Weibull
body centering is retained as an explicit extension, using the same construction.

## Body, tail and reference measure

Let `tau` be the largest finite endpoint in the first two input columns, including
right-censoring times. The root is `[0, Inf]`. The body is `(0,tau]`; the tail is
the single set `(tau,Inf)`. No distributional shape is specified inside the tail,
for any centering family. This is not an atom at infinity or a parametric tail.

`P0` is uniform, exponential or Weibull restricted and normalized to `(0,tau]`.
Ordinary splits use its conditional median on each node's body portion. Its
right child retains the true parent endpoint, including infinity. A node splits
only when its depth is below `max_M` and its expected count under the reference
measure is at least `min_count`. Even if pruned early, the rightmost leaf is
always split at `tau`. There is exactly one tail leaf. The retained tree stays
fixed during both optimization and sampling.

The reference measure puts mass `1-s0` in the body according to `P0` and `s0`
in the undivided tail. By default, equation (4) is used:

\[
s_0=\frac{\kappa\tau p_0(\tau)}{1+\kappa\tau p_0(\tau)}.
\]

`kappa=0.1`. Uniform centering gives `s0=1/11`. For exponential centering,
`tau*p0(tau) = lambda*tau/expm1(lambda*tau)`. For the Weibull extension, with
`h=(lambda*tau)^shape`, it is `shape*h/expm1(h)`. Computation uses log scale to
avoid cancellation and overflow. If automatic calibration is too small to
represent as a positive double, an explicit representable `s0` is required.
This calibration does not set a finite upper endpoint or a tail density.

## Beta prior and initialization

For children at depth `m`, let `r` be the reference left-child mass divided by
its parent's reference mass. Equation (2) gives

\[
\alpha_L=\max(1,2\gamma m^2r),\qquad
\alpha_R=\max(1,2\gamma m^2(1-r)).
\]

The default is `gamma=0.5`. Body-only median splits have `r=1/2` and recover
the old symmetric `max(1,m^2/2)` parameters. Splits on the rightmost path are
generally asymmetric, including the final body/tail split. `gamma=0` gives
the flat prior on every branch.

The initial measure is the reference measure itself: branches start at `r`,
not necessarily one half or the floored Beta means. Initial expected counts
use the reference probability of the intersection with an observation divided
by the probability of that observation. These counts determine early stopping.
They are not counts propagated by equal splitting of overlapping observations.

Flooring a Beta parameter at one can change its mean. Accordingly `Prior$s0`
is the reference tail target, while `Prior$tail_mean` is the actual induced
prior mean, computed as the product of right-branch Beta means. They need not
coincide, including on pruned trees. Both are exposed rather than silently
calling the reference mass the exact prior mean.

## Deterministic MAP-EM (Algorithm 1)

For leaf `j` with current mass `p_j`, the coverage weight `h_ij` is the
conditional body probability of observation `i` within that leaf. A tail leaf
has weight one for every right-censored observation and zero for finite
observations. For a censored observation,

\[
z_{ij}=\frac{p_jh_{ij}}{\sum_k p_kh_{ik}}.
\]

Exact observations have known labels. Their conditional density factors are
constant in the branch probabilities, so the implementation uses the known
label directly. Positive boundary events belong to the interval on the left;
zero belongs to intervals beginning at zero.

Leaf expected counts are summed over observations, then aggregated up the tree.
Each M-step uses

\[
q_L=\operatorname{clip}_{[\epsilon,1-\epsilon]}
\frac{\alpha_L+N_L-1}{\alpha_L+\alpha_R+N_L+N_R-2}.
\]

All prior parameters are at least one. If the denominator is zero, both
posterior shapes equal one and `q_L=1/2` is used. `epsilon=1e-6` is configurable.
Child masses are parent mass times conditional branch probabilities.

By default a flat-prior EM fit is run first and initializes the working-prior
fit. Each stage stops when the Euclidean leaf-mass change is at most `EM_tol`
(default `1e-6`), or after `max_EM_iter` updates (default 2000 per stage).
`flat_init=FALSE` skips the first stage. `max_EM_iter=0` returns the reference
measure and reference allocations without either optimization stage.

`Convergence$flat` and `$map` return iteration counts, convergence flags, final
changes and objective histories. The recorded log posterior kernels omit Beta
normalization constants and exact-observation conditional-density constants;
these are independent of the branch probabilities on the retained tree.
Compare monotonicity within a stage, not across the two different priors.
Final `Allocations` are recomputed from the returned masses.

## Exact-target DA-Gibbs (Algorithm 2)

Starting from the deterministic fit, each sweep:

1. Samples each censored observation's leaf label from `z_ij`; exact labels stay fixed.
2. Aggregates integer leaf counts through the retained tree.
3. Samples each left branch independently from
   `Beta(alpha_L + n_L, alpha_R + n_R)` and rebuilds node masses.
4. Records leaf masses, branch probabilities, survival and body density.

Beta draws are unrestricted: the MAP epsilon clipping is not applied. The
chain targets the posterior conditional on this tree and fixed base parameters.
This is not an independent Beta approximation based on fractional MAP counts.
No continuous event times are imputed.

Model fitting is deterministic and independent of target times. Prediction
with `predict(fit,times,0.95)` runs DA-Gibbs by default during that call; use
`posterior=FALSE` for deterministic predictions alone. The manuscript's setup is
`n_sweeps=2000, burn_in=500, thin=1`, retaining 1500 sweeps. `set.seed()` controls
R's random stream. Deterministic predictions remain in `Survival_Probs`;
summaries and optional raw draws appear separately in the prediction's `Posterior`.

`samplePT(fit)` exposes grid-independent tree draws. The ordinary prediction
call samples these trees internally and evaluates them only at the supplied
target times. Passing `samples=draws` explicitly reuses a previous chain without
sampling again. Neither prediction nor sampling changes the input model,
re-estimates base parameters, reruns EM or rebuilds the partition. Fitted state
is stored as plain R arrays, so saved models and samples work in a fresh session.

Reported intervals are equal-tailed pointwise credible intervals. Approximate
single-chain ESS uses the initial positive sequence of paired autocorrelations,
with at most 1000 lags; MCSE is the sample SD divided by the square root of ESS.
These summaries do not prove mixing. Inspect trace plots and compare independent
chains when using inference substantively. Tree and base-parameter uncertainty
are not integrated out.

## Evaluation and implementation

On `[0,tau]`, survival and density are mixtures of the same conditional leaf
distributions used in allocation. `S(tau)` equals `Tail_Prob`; body density
integrates to `1-Tail_Prob`. At times above `tau`, including `Inf`, both outputs
are `NA`, since the tail shape has not been supplied. Full mean event time is
therefore not identified by this model. Posterior means refer to means of
survival probabilities, densities, leaf masses or tail mass across draws.

Coverage probabilities and allocation normalization use log scale. Exact
labels are handled separately. Expected/integer counts are aggregated once
up the tree, so the allocation work is `O(n * number_of_leaves)` per sweep.
Coverage and prediction caches each have a 64 MiB cap, with uncached fallbacks.
Raw native sample matrices are limited to 512 MiB before R summary overhead.

`allocated=0` retains the separate historical overlap-count estimator, using
the revised reference tree and Beta prior. It has no EM or posterior sampling;
it is not Algorithm 1 and is not a faster setting for the same estimator.

| Source | Responsibility |
| --- | --- |
| `R/PT.R` | Fit-only S3 model, validation and rate selection |
| `R/predict.R` | Target-time prediction, sampling orchestration and posterior summaries |
| `R/sample.R` | Standalone posterior sampling from a stored model |
| `src/PT.cpp` | Reference tree, asymmetric priors, MAP-EM, DA-Gibbs, prediction |
| `R/MLE_estimate.R` | Parametric fits; exponential score solved without a rate floor |
| `tests/scpt.R` | Reference/count checks, MAP ascent, analytic and quadrature posterior checks |
| `tests/interfaces.R` | Fresh prediction sampling, time ordering, reusable draws, serialization and unchanged model state |
