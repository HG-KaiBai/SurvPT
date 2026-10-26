# Implementation and statistical interpretation

[Back to README](../README.md)

This document describes the finite Polya tree estimator implemented by `survPT()` in version 0.2.0, including its partition, allocation updates, and tail conventions.

## Finite PT partition

`survPT()` partitions the declared support into a binary tree. For a finite uniform interval the split is its midpoint. Exponential and Weibull splits bisect the base probability of the parent interval. For a uniform interval with infinite right endpoint, the split uses the data-dependent cutoff `C = 1.1 * largest_obs` and the current depth, where `largest_obs` is the largest finite value in the first two columns of `X`.

Splitting ends at `max_M` or when the initial allocated count is below one. Consequently, terminal intervals can occur at different depths. The tree is initialized once and does not grow again during EM.

For children at depth `l`, the branch prior parameter is `alpha_l = l^2 / 2`. The current branch-mass update computes

\[
q_L=\operatorname{clip}_{[0,1]}
\frac{\alpha_l+n_L-1}{2\alpha_l+n_L+n_R-2},\qquad q_R=1-q_L,
\]

with a one-half fallback when the denominator is too small. Child masses are parent mass times these branch proportions. This is a mode-style update with safeguards, not the usual posterior-mean formula `(alpha + n_L) / (2 * alpha + n_L + n_R)`.

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

Exact observations are assigned directly by their event times. Probability computations use direct exponential differences and a small-denominator guard; extreme tail parameters can still suffer numerical cancellation or underflow.

## Three meanings of truncation

1. **Tree truncation:** stopping recursive partitioning at a finite depth or low initial count.
2. **Conditional base truncation:** restricting and renormalizing a base distribution inside one terminal interval.
3. **Uniform time cutoff:** replacing the infinite right endpoint of a uniform PT leaf by `C = 1.1 * largest_obs` for both EM and output.

The third is a modeling choice in this implementation. Exponential and Weibull PT bases do not impose that finite cutoff: an exponential tail leaf has distribution `a + Exp(lambda)`, while a Weibull tail leaf follows the original Weibull law conditional on `T > a`. None of these operations implements delayed entry or a likelihood for sampling-truncated data.

## Implementation map

| Source | Responsibility |
| --- | --- |
| `R/PT.R` | Input wrapper, defaults, result names, conversion of terminal indices to R indexing |
| `src/PT.cpp` | PT tree, initialization, EM, shared conditional-base probabilities, output |
| `R/MLE_estimate.R` | Parametric likelihoods, gradients, and optimization |
| `R/data_example.R` | Simulation generators |
| `tests/pt-conditional-base.R` | Conditional allocation, output, mass, and tail regression checks |

`survPT()` provides point estimates. It does not return posterior samples, credible intervals, automatic tuning selection, or proof of convergence.
