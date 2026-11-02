# Implementation and statistical interpretation

[Back to README](../README.md)

This document describes the finite Polya tree estimator implemented by `survPT()` in version 0.2.1, including its partition, allocation updates, and tail conventions.

## Finite PT partition

Every base starts with the root interval `[0, Inf]`; `support` must be `c(0, Inf)`. Let `T = largest_obs`, the largest finite endpoint in the first two columns of `X`, including right-censoring times and ignoring extra columns. `T` must be positive. For an ordinary node `(a,b]`, the split is the median of the base restricted to `[a, min(b,T)]`. Uniform splits are interval midpoints. For exponential and Weibull bases the split `s` satisfies

\[
F_0(s)=\frac{F_0(a)+F_0(\min(b,T))}{2}.
\]

The right child retains its parent's actual right endpoint, including `Inf`. Thus truncation determines the split locations without truncating the root support.

Ordinary splitting ends at `max_M` or when the construction weight is below one. During fractional tree construction, an observation's current weight is passed to its compatible children in proportion to their branch probabilities: equally when both children overlap its interval, and entirely to the sole compatible child otherwise. No base probability is used for these construction weights. Consequently, terminal intervals can occur at different depths. Once these splits are complete, the rightmost terminal interval `(a, Inf)` is always split once more at `T`, into `(a,T]` and `(T,Inf)`. If its depth is `m`, these two children have depth `m+1`, and their branch prior is `Beta(alpha, alpha)` with `alpha = max(1, 0.5*(m+1)^2)`. This final split also applies when `max_M = 0` or the rightmost node stopped early. The tree is initialized once and does not grow again during EM.

With fractional allocation, every branch initializes to `q_L = q_R = 1/2`, including the final tail split, so a depth-`d` leaf starts with mass `2^(-d)`. Initialization computes observation allocations under that complete tree but does not perform a count-based mass update. For children at depth `l`, the branch prior parameter is `alpha_l = max(1, l^2 / 2)` in both allocation modes and for all bases. Each EM M-step computes

\[
q_L=\operatorname{clip}_{[\epsilon,1-\epsilon]}
\frac{\alpha_l+n_L-1}{2\alpha_l+n_L+n_R-2},\qquad q_R=1-q_L,
\]

where `epsilon = 1e-6`, with a one-half fallback when the denominator is too small. The same split-probability bounds apply to the final tail split and the overlap-count estimator (`allocated = 0`). Child masses are parent mass times these branch proportions; the bounds constrain conditional split probabilities, not the resulting node masses. This is a mode-style update with safeguards, not the usual posterior-mean formula `(alpha + n_L) / (2 * alpha + n_L + n_R)`.

## Conditional base distribution inside terminal intervals

Let `A = (a,b]` be a terminal interval, using `[0,b]` when `a = 0`, `p_A` its fitted mass, and `F_0`, `f_0`, `S_0` the base CDF, density, and survival function. For exponential and Weibull bases,

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

The fitted density and survival function mix these conditional distributions using `p_A`. Finite uniform leaves use their interval's uniform distribution. The uniform tail `(T,Inf)` instead holds a point mass `p_inf` at infinity: its density is zero at every finite time, and its contribution to survival is `p_inf` at every finite time. Consequently `S(t) = p_inf` for finite `t >= T`, and the continuous density integrates to `1-p_inf`. At `Inf`, the reported survival and density are zero. Right-censored intervals include this point mass, while finite observation intervals do not.

EM integrates these same terminal distributions when allocating an interval-censored observation to child nodes:

\[
w_{i,B}=\frac{P(T_i\in B\cap I_i)}{P(T_i\in I_i)}.
\]

Exact observations are assigned directly by their event times. At a positive split point, including `T`, both allocation and density use the interval on the left; zero is included in intervals starting at zero. Density is zero at infinite evaluation times. Exponential and Weibull conditional probabilities use cumulative-hazard differences and `expm1`, so the new tail retains its conditional mass even when unconditioned survival at `T` underflows. Splits use `log1p` and `expm1` to preserve small truncated probabilities; ordinary splitting stops if no interior median is representable.

For `allocated = 1`, the same conditional-allocation rule is used initially and in every E-step. If leaf `A_j` has current mass `p_j` and conditional observation probability `h_ij`, its weight for observation `i` is `p_j * h_ij / sum_k(p_k * h_ik)`. Initial masses are products of one-half branch probabilities; subsequent E-steps use the preceding M-step's masses. Base probabilities are evaluated only inside terminal intervals. An internal node's weight sums those of its descendant leaves. These completed-tree initial counts replace the provisional construction weights, without changing the partition or its equal initial branch probabilities.

For the uniform tail `(T,Inf)`, `h_ij` is one for every right-censored observation and zero for every finite observation. All finite observation endpoints are at most `T`, so this rule requires no uniform distribution truncated to an infinite interval. A right-censored observation at `L` contributes `p_tail / S(L)` to the tail; it need not be allocated entirely there. Finite uniform leaves use overlap length divided by leaf width. Exponential and Weibull leaves use their conditional base probabilities. The separate `allocated = 0` option retains overlap-count initialization, its mode update and no EM; uniform right-censored observations are entirely assigned to its tail.

During EM, the partition, observation intervals and base parameters stay fixed. The implementation caches observation/leaf conditional probabilities and child links, and computes each observation's denominator once per iteration. A node's numerator sums its descendant leaves in the original `terminal_nodes` order. Observation accumulation, branch updates and the convergence calculation retain their original order; no child-total regrouping or prefix-sum subtraction is used. The additional dense probability cache is limited to 64 MiB per fit; larger problems and degenerate numerical partitions use the uncached allocation calculation.

## Truncation and tail conventions

1. **Tree truncation:** stopping ordinary recursive partitioning at a finite depth or low initial count, followed by the mandatory tail split.
2. **Split distribution truncation:** restricting the base to the node's intersection with `[0,T]` to select its ordinary median split.
3. **Conditional terminal distribution:** restricting and renormalizing the base inside a finite leaf or an exponential/Weibull tail; the uniform tail is instead an atom at infinity.

There is no finite tail cutoff. The exponential tail has distribution `largest_obs + Exp(lambda)`, while the Weibull tail follows the original Weibull law conditional on exceeding `largest_obs`. None of these operations implements delayed entry or a likelihood for sampling-truncated data.

## Implementation map

| Source | Responsibility |
| --- | --- |
| `R/PT.R` | Input wrapper, defaults, result names, conversion of terminal indices to R indexing |
| `src/PT.cpp` | PT tree, initialization, EM, shared conditional-base probabilities, output |
| `R/MLE_estimate.R` | Parametric likelihoods, gradients, and optimization |
| `R/data_example.R` | Simulation generators |
| `tests/pt-conditional-base.R` | Conditional allocation, output, mass, and tail regression checks |
| `tests/pt-tree-structure.R` | Truncated-median partitions, final tail layer, Beta parameters, and numerical tail checks |

`survPT()` provides point estimates. It does not return posterior samples, credible intervals, automatic tuning selection, or proof of convergence.
