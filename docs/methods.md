# Implementation and statistical interpretation

[Back to README](../README.md)

Version 0.4.1 retains the reference tree and tail construction from the supplied
manuscript *Self-Consistent Polya Trees for Arbitrarily Censored Survival Data*,
with a single fit initialized at the centering reference distribution.
The manuscript's uniform and exponential body centers are supported. Weibull
body centering is retained as an explicit extension, using the same construction.

## Body, tail and reference measure

Let `tau` be the largest finite endpoint in the first two input columns, including
right-censoring times. The root is `[0, Inf]`. The body is `(0,tau]`; the tail is
the single set `(tau,Inf)`. During fitting and Gibbs allocation, no distributional
shape is specified inside the tail, for any centering family. Prediction alone
adds an exponential extrapolation matched to each fitted or sampled tail mass.

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
the flat prior on every branch. The explicit option `flat=TRUE` has the same
effect, overriding `gamma` for both fitting and posterior sampling. It does
not change the reference measure, splits or pruning.

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
q_L=\operatorname{clip}_{[\epsilon_m,1-\epsilon_m]}
\frac{\alpha_L+N_L-1}{\alpha_L+\alpha_R+N_L+N_R-2}.
\]

All prior parameters are at least one. If the denominator is zero, both
posterior shapes equal one and `q_L=1/2` is used. At child level `m`,
`epsilon_m = split_epsilon * 2^(-m)`, with base scale `split_epsilon=1e-6`.
The root is level zero, so its children use a bound of `5e-7`. The final
body/tail split follows its children's actual level, including after pruning.
Child masses are parent mass times conditional branch probabilities.

The first E-step uses the centering reference masses. A single MAP-EM fit uses
the selected alpha parameters throughout and stops when the Euclidean leaf-mass
change is at most `EM_tol` (default `1e-6`), or after `max_EM_iter` updates
(default 2000). There is no flat-prior warm start. `flat=TRUE` selects a flat
model for the entire fit. With `allocated=1`, `max_EM_iter=0` returns the
reference measure and reference allocations without optimization.

`Convergence$map` returns the iteration count, convergence flag, final change
and objective history. The recorded log posterior kernels omit Beta
normalization constants and exact-observation conditional-density constants;
these are independent of the branch probabilities on the retained tree.
The history refers to the single selected prior.
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
`sample=FALSE` for deterministic predictions alone. The defaults are
`n_samples=2000, burn_in=3000, thin=1`. There are `burn_in+n_samples*thin`
sweeps; retention begins at `burn_in+thin`. `set.seed()` controls R's random
stream. Deterministic predictions appear in `estimate`; `lower`, `upper`,
`mean` and `median` summarize the posterior at each supplied `x`.

Prediction samples tree masses internally and evaluates them at the supplied
times. Draws are used temporarily to compute summaries and are not returned.
The one-row `tail` attribute reports
`S(tau)` and its posterior summary independently of the target grid or curve type.
Neither prediction nor sampling changes the input model,
re-estimates base parameters, reruns EM or rebuilds the partition. Fitted state
is stored as plain R arrays, so saved models work in a fresh session.

Reported intervals are equal-tailed pointwise credible intervals.
Compare summaries across independent runs and longer chains to assess stability.
Tree and base-parameter uncertainty
are not integrated out.

## Evaluation and implementation

On `[0,tau]`, survival and density are mixtures of the same conditional leaf
distributions used in allocation. `S(tau)` equals `Tail_Prob`; body density
integrates to `1-Tail_Prob`. For prediction above `tau`, let `p` denote the
fitted tail mass or one Gibbs draw's tail mass. Impose the additional matching
condition `exp(-lambda*tau)=p`, which gives `lambda=-log(p)/tau`. The output is

\[
S(x)=p\exp\{-\lambda(x-\tau)\}=p^{x/\tau},\qquad
f(x)=\lambda S(x),\quad x>\tau.
\]

The tail density integrates to `p`; survival is continuous at `tau`. Density
uses the body endpoint convention at `tau` and need not be continuous there.
At `Inf`, survival and density are zero. Zero mass gives a zero contribution;
mass exactly one cannot yield a positive matching rate and raises an error
during extrapolation. Computation uses log scale for tiny masses and distant times.

Each posterior draw uses its own `p` and `lambda` before quantiles and moments
are calculated. This is not extrapolation of posterior summary columns with a
single fitted rate. The shape is an output assumption, not part of the fitting
likelihood or Gibbs allocation. A normalized truncated exponential alone would
not identify a rate; the matching condition above supplies that extra assumption.
Full mean event time is not identified by the fit without a tail-shape assumption.
Posterior means reported by prediction are means of survival probabilities or
densities, not mean event times.

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
| `src/PT.cpp` | Reference tree, asymmetric priors, MAP-EM, DA-Gibbs, prediction |
| `R/MLE_estimate.R` | Parametric fits; exponential score solved without a rate floor |
| `tests/pt-exponential-tail.R` | Prediction-only exponential extension, per-draw rates, density normalization and boundary behavior |
