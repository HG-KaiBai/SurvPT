# Implementation and statistical interpretation

[Back to README](../README.md)

This document describes the finite Polya tree estimator implemented by `survPT()` in version 0.2.3, including its partition, allocation updates, and tail conventions.

## Finite PT partition

`survPT()` partitions the declared support into a binary tree. For a finite uniform interval the split is its midpoint. Exponential and Weibull splits bisect the base probability of the parent interval. For a uniform interval with infinite right endpoint, the split uses the data-dependent cutoff `C = uniform_tail_factor * largest_obs` and the current depth, where `largest_obs` is the largest finite value in the first two columns of `X`. The multiplier defaults to `1.1` and must be a finite number greater than one. The same multiplier is used for infinite observation bounds during uniform initialization and for infinite terminal bounds during EM and prediction.

Splitting ends at `max_M` or when the initial allocated count is below one. Consequently, terminal intervals can occur at different depths. The tree is initialized once and does not grow again during EM.

For children at depth `m`, the branch prior parameter is `alpha_m = max(1, 0.5 * m^2)` by default (`flat = FALSE`). Set `flat = TRUE` to use `alpha_m = 1` at every depth. The root parameter is always 1. Both initialization and EM use the branch-mass update

\[
q_L=\operatorname{clip}_{[\epsilon,1-\epsilon]}
\frac{\alpha_m+n_L-1}{2\alpha_m+n_L+n_R-2},\qquad q_R=1-q_L,
\]

with `epsilon = 1e-6` and a one-half fallback when the denominator is too small. Both split probabilities lie in `[epsilon, 1-epsilon]`. Child masses are parent mass times these branch proportions; the bounds apply to conditional split probabilities, not node masses. This is a mode-style update with safeguards, not the usual posterior-mean formula `(alpha + n_L) / (2 * alpha + n_L + n_R)`.

## Uniform body with exponential tail

For `prior = 4`, let `M = largest_obs > 0`, `lambda = prior_params[1] > 0`,
and `w = 1 - exp(-lambda * M)`. The rate is supplied directly and held fixed;
an exponential MLE from `estimate_exp()` can be passed by the caller. The base is

\[
F_0(t) = \begin{cases}
0, & t < 0,\\
w\,t/M, & 0 \le t \le M,\\
1-\exp(-\lambda t), & t > M,
\end{cases}
\qquad
f_0(t) = \begin{cases}
w/M, & 0 \le t \le M,\\
\lambda\exp(-\lambda t), & t > M,\\
0, & t < 0.
\end{cases}
\]

The CDF is continuous at `M`; the density need not be. The implementation uses
the body density at exactly `M`. The base tail mass is `exp(-lambda * M)`, and
the tail conditional on exceeding `M` is `M + Exp(lambda)`.
The inverse CDF, for `0 < p < 1`, is

\[
F_0^{-1}(p) = \begin{cases}
M p/w, & p \le w,\\
-\log(1-p)/\lambda, & p > w.
\end{cases}
\]

Splits bisect the parent interval's base probability. In the upper half of
the distribution, the implementation averages survival probabilities and uses
the inverse survival function, avoiding loss of the exponential tail when the
CDF rounds to one. Unrepresentable interior splits terminate the node.
Interval probabilities add the uniform and exponential contributions when an
interval crosses `M`; `expm1` avoids subtracting nearly equal exponentials.
Initialization, EM, density, and survival all use this same base, conditioned
on each node or terminal interval. A finite upper support truncates the base.
`uniform_tail_factor` is not used by this prior.

The specified tail mass belongs to the base distribution. Fitting updates
the tree masses, so fitted survival at `M` need not equal `exp(-lambda * M)`.

## Conditional base distribution inside terminal intervals

Let `A = (a,b]` be a terminal interval, using `[0,b]` when `a = 0`, `p_A` its fitted mass, and `F_0`, `f_0`, `S_0` the base CDF, density, and survival function. For exponential, Weibull, and spliced uniform-exponential bases,

\[
f_A(t)=\frac{f_0(t)\,1\{t\in A\}}{F_0(b)-F_0(a)},
\qquad
P(T\in(u,v)\mid T\in A)
=\frac{F_0(\min(v,b))-F_0(\max(u,a))}{F_0(b)-F_0(a)},
\]

where the probability is zero if the intersection is empty. For `a < t < b`,

\[
S_A(t)=\frac{S_0(t)-S_0(b)}{S_0(a)-S_0(b)}.
\]

The fitted density and survival function mix these conditional distributions using `p_A`. In version 0.2.0, EM integrates the same terminal distributions when allocating an interval-censored observation to child nodes:

\[
w_{i,B}=\frac{P(T_i\in B\cap I_i)}{P(T_i\in I_i)}.
\]

Exact observations are assigned directly by their event times. At a positive split point, both allocation and density use the interval on the left; zero is included in intervals starting at zero. Density is zero at infinite evaluation times. Probability computations use direct exponential differences and a small-denominator guard; extreme tail parameters can still suffer numerical cancellation or underflow.

During EM, the partition, observation intervals and base parameters stay fixed. The implementation caches observation/leaf conditional probabilities and child links, and computes each observation's denominator once per iteration. A node's numerator sums its descendant leaves in the original `terminal_nodes` order. Observation accumulation, branch updates and the convergence calculation retain their original order; no child-total regrouping or prefix-sum subtraction is used. Initialization and prediction are unchanged. The additional dense probability cache is limited to 64 MiB per fit; larger problems and degenerate numerical partitions use the original allocation calculation.

## Posterior prediction and DA-Gibbs sampling

`survPT()` only fits and retains the tree, its prior parameters, fitted masses, and observation endpoints. `predict(model, x, level = 0.95)` evaluates the same fitted distribution as the previous survival output in its `estimate` column. The prediction grid never enters fitting.

For posterior uncertainty, prediction applies Algorithm 2's two-block data-augmentation Gibbs scheme to this existing tree. For each censored observation, compute its conditional base coverage `w[i, leaf]` in every compatible leaf. At each sweep, draw one leaf label with probability proportional to `mass[leaf] * w[i, leaf]`; exact-event labels are fixed. Aggregate the resulting integer leaf counts up the tree and draw each split from `Beta(alpha_left + count_left, alpha_right + count_right)`. Recompute node masses along the paths and evaluate the resulting survival curve. Allocations change each sweep; fractional EM counts are not frozen and substituted into Beta distributions.

The chain starts at the fitted masses and samples unrestricted Beta probabilities, without the clipping used for MAP optimization. The partition, prior parameters and within-leaf base distributions remain fixed. No special tail-node construction is introduced. Nonzero coverage weights are cached per observation, and label probabilities are normalized on the log scale.

The defaults discard 3,000 sweeps and retain 5,000 draws (8,000 total sweeps with the default thinning). `n_samples`, `burn_in` and `thin` control the run; `set.seed()` controls reproducibility. The result contains posterior means, medians and pointwise equal-tailed credible intervals. These are conditional on the retained tree and base parameters, and are Monte Carlo estimates rather than a guarantee of convergence. `keep_samples = TRUE` retains the sampled curves for diagnostics. `sample = FALSE` evaluates only the fitted point estimate. `type = "density"` applies the same calculation to density.

## Three meanings of truncation

1. **Tree truncation:** stopping recursive partitioning at a finite depth or low initial count.
2. **Conditional base truncation:** restricting and renormalizing a base distribution inside one terminal interval.
3. **Uniform time cutoff:** replacing the infinite right endpoint of a uniform PT leaf by `C = uniform_tail_factor * largest_obs` for both EM and output (default multiplier `1.1`).

The third is a modeling choice in this implementation. Exponential and Weibull PT bases do not impose that finite cutoff: an exponential tail leaf has distribution `a + Exp(lambda)`, while a Weibull tail leaf follows the original Weibull law conditional on `T > a`. None of these operations implements delayed entry or a likelihood for sampling-truncated data.

## Implementation map

| Source | Responsibility |
| --- | --- |
| `R/PT.R` | Input wrapper, defaults, result names, conversion of terminal indices to R indexing |
| `R/predict.R` | S3 prediction interface, sampling controls and posterior summaries |
| `src/PT.cpp` | PT tree, initialization, EM, shared conditional-base probabilities, output |
| `src/PT_predict.cpp` | Restore fitted trees, evaluate predictions, and sample posterior leaf labels and Beta splits |
| `R/MLE_estimate.R` | Parametric likelihoods, gradients, and optimization |
| `R/data_example.R` | Simulation generators |
| `tests/pt-conditional-base.R` | Conditional allocation, output, mass, and tail regression checks |
| `tests/pt-uniform-exp.R` | Spliced-base CDF/quantile, allocation, density, tail, and interface checks |
| `tests/pt-predict.R` | Fit/predict separation, analytic posterior checks, Gibbs sampling and reproducibility |

Neither fitting nor prediction performs automatic tuning selection or certifies convergence.
