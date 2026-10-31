# Implementation and statistical interpretation

[Back to README](../README.md)

This document describes the finite Polya tree estimator implemented by `survPT()` in version 0.3.2, including its partition, allocation updates, and tail conventions.

## Finite PT partition

With infinite upper support, `survPT()` first splits at `tau`, the largest finite endpoint in the first two columns of `X`, including right-censoring times. The root's left child is the finite PT subtree with mass `theta_0`; its right child is a `[tau, Inf)` tail with mass `p_tail = 1 - theta_0`. Only the finite subtree is subdivided. For a uniform interval the split is its midpoint. Exponential and Weibull splits bisect the base probability of the finite parent interval. A finite upper support retains the bounded tree without an extra root split.

Splitting ends at `max_M` within the finite subtree or when the initial allocated count is below one. The added root split does not consume this depth budget: with unbounded support, the finite subtree root is reported at level one. Consequently, terminal intervals can occur at different depths. The tree is initialized once and does not grow again during EM.

For finite-subtree children at relative depth `l`, the branch prior parameter is `alpha_l = gamma * l^2`, with `gamma = 1`. All split probabilities are bounded to `[epsilon, 1 - epsilon]`, with `epsilon = 1e-6`, at initialization and in every EM iteration. The branch-mass update computes

\[
q_L=\operatorname{clip}_{[\epsilon,1-\epsilon]}
\frac{\alpha_l+n_L-1}{2\alpha_l+n_L+n_R-2},\qquad q_R=1-q_L,
\]

with a one-half fallback when the denominator is too small. Child masses are parent mass times these branch proportions. The bounds apply to conditional split probabilities; a descendant's absolute mass can be smaller than `epsilon`. This is a mode-style update with safeguards, not the usual posterior-mean formula `(alpha + n_L) / (2 * alpha + n_L + n_R)`.

## Beta prior for the first layer

Let `f` and `P_f` denote the conditional finite-subtree density and probabilities. The observed likelihood contributions are:

| Observation | Contribution |
| --- | --- |
| Exact `t <= tau` | `theta_0 * f(t)` |
| Finite interval `(L, R]`, `R <= tau` | `theta_0 * P_f(L < T <= R)` |
| Right censored at `C < tau` | `theta_0 * P_f(C < T <= tau) + p_tail` |
| Right censored at `tau` | `p_tail` |

The first layer has the internal prior `q_0 ~ Beta(a,b)` with `a = b = 1`. Here `q_0 = theta_0` and `p_tail = 1 - q_0` are mathematical notation, not separate return fields; masses remain in the existing node arrays. Initialization uses the prior mean `0.5` to allocate censored observations between the finite conditional base and the tail. Once allocated counts are available, the first layer uses the same clipped mode-style rule as the PT branches:

\[
q_0^{\mathrm{new}}=\operatorname{clip}_{[\epsilon,1-\epsilon]}
\frac{a+n_{\mathrm{body}}-1}{a+b+n_{\mathrm{body}}+n_{\mathrm{tail}}-2}.
\]

When a branch denominator is nonpositive or at most `1e-300`, the existing one-half fallback is used. With fractional allocations, `n_body + n_tail = n`, so the first-layer update simplifies to `q_0 = clip(n_body / n)` and `p_tail = clip(n_tail / n)`. Both body and tail masses therefore remain between `epsilon` and `1 - epsilon`, including when every observation is exact or every observation is right censored at `tau`. The corresponding complete-data Beta parameters are `1 + n_body` and `1 + n_tail`; this is not their posterior-mean update.

In every EM iteration the E-step includes current tail mass in each right-censoring denominator, and **updates `q_0` and `p_tail` along with the finite tree**. There is no separate scalar MLE or frozen first-layer probability. The finite subtree uses `alpha_l = l^2` relative to its own root; adding this layer does not shift those alpha values or consume `max_M`. `LogLik_Trace` reports observed log-likelihood, which need not increase under the prior-adjusted mode-style updates. Bounded support has no extra split and uses the same finite-tree priors and split bounds.

## Conditional base distribution inside terminal intervals

Let `A = (a,b)` be a terminal interval, `p_A` its fitted mass, and `F_0`, `f_0`, `S_0` the base CDF, density, and survival function. For exponential and Weibull bases,

\[
f_A(t)=\frac{f_0(t)\,1\{a<t<b\}}{F_0(b)-F_0(a)},
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

Exact observations are assigned directly by their event times. At internal boundaries and at `tau`, exact-event densities use the interval on the left, matching the allocation convention. Probability computations use direct exponential differences and a small-denominator guard; extreme tail parameters can still suffer numerical cancellation or underflow.

The fixed partition allows conditional observation probabilities (and exact-event densities) to be computed once and cached as nonzero observation/leaf pairs. Each E-step weights these entries by the previous leaf masses and sums them up the tree once per observation. This replaces repeated integration over every leaf for every child node. With `n` observations, `K` nodes and `L` leaves, the worst-case allocation cost per iteration drops from `O(n K L)` to `O(n K)`; the additional cache stores only nonzero pairs (at most `n L`). Child links are cached too, and the log-likelihood reuses the same observation probabilities. The partition, alpha updates, containment shortcuts and stopping criterion are unchanged; reordered sums can introduce floating-point rounding differences.

## Tail interpretation

1. **Tree truncation:** stopping recursive partitioning at a finite depth or low initial count.
2. **Conditional base truncation:** restricting and renormalizing a base distribution inside one terminal interval.
3. **Tail mass:** retaining `p_tail` beyond the largest finite observed endpoint. The exponential base supplies a conditional tail distribution; uniform and Weibull tail densities remain unspecified.

There is no longer a factor-1.1 uniform cutoff. For all three bases, `S(tau) = p_tail`, and the finite density integrates to `theta_0`. At finite times strictly greater than `tau`, uniform survival is carried forward as `S(tau)`; this is an output convention and does not specify a tail density, which remains `NA`. Weibull survival and density remain `NA` beyond `tau`.

For the exponential base with fixed rate `lambda = prior_params`, the terminal tail uses the base conditioned on `T > tau`:

\[
S(t)=p_{\mathrm{tail}}\exp\{-\lambda(t-\tau)\},\qquad
f(t)=p_{\mathrm{tail}}\lambda\exp\{-\lambda(t-\tau)\},\quad t>\tau.
\]

The exponential tail density integrates to `p_tail`, so the full density integrates to one. Tail computations use the shifted time `t - tau` to avoid dividing by a potentially underflowed `exp(-lambda * tau)`. The rate is not re-estimated during EM; the tail mass is updated with the other node masses. Because every finite observed endpoint is at most `tau`, conditional tail probabilities for right-censored observations are one, preserving the likelihood contributions above. Exact events at `tau` still use the finite interval on the left.

At positive infinity, survival and density return zero for all bases. The tail is included in `Terminal_Nodes`, so all terminal masses sum to one. None of these operations implements delayed entry or a likelihood for sampling-truncated data.

## Implementation map

| Source | Responsibility |
| --- | --- |
| `R/PT.R` | Input wrapper, defaults, result names, conversion of terminal indices to R indexing |
| `src/PT.cpp` | PT tree, initialization, EM, shared conditional-base probabilities, output |
| `R/MLE_estimate.R` | Parametric likelihoods, gradients, and optimization |
| `R/data_example.R` | Simulation generators |
| `tests/pt-conditional-base.R` | Conditional allocation, output, mass, and tail regression checks |
| `tests/pt-tail-root.R` | Beta-prior initialization, first-layer EM updates, likelihood and endpoint checks |

`survPT()` provides point estimates, an observed log-likelihood trace, and a terminal-mass convergence flag. It does not return posterior samples, credible intervals, automatic tuning selection, or proof of convergence to a global optimum.
