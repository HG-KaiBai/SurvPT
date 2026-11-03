# Examples and research scripts

Start with [`quickstart.R`](quickstart.R). It uses synthetic data, fits `survPT()` with uniform, exponential, and Weibull bases, checks basic probability properties, and does not write files.

[`posterior.R`](posterior.R) compares PT with uniform and exponential centers,
fiducial inference, mixed Polya trees (MPT), and transformed Bernstein polynomials
(TBP). It uses the current `test_function.R` data example: `data_rc(seed=12,
n=200,K=1,Lambda=0.12,q=0.3)`. All methods use the same 201 target times in the
observed range. PT retains MAP-EM curves alongside posterior means and pointwise
95% credible intervals. Predictions above a PT fit's largest finite endpoint
remain unspecified and return `NA`.

The comparison requires `devtools`, `survival`, `spBayesSurv` **1.1.9**, `pracma`
and `quadprog`. MPT/TBP share the fitting helpers in `bayes_helpers.R` with
`test_function.R`, including its baseline-only MPT initializer adapter. Their
package-default priors and MCMC controls are preserved (3000 burn-in sweeps and
2000 retained draws). PT uses 2000 sweeps with 500 burn-in; fiducial uses 1000
retained draws with 100 burn-in. Each method has an explicit seed.

Run `Rscript Examples/posterior.R` from the repository root. It writes these
files under `output/posterior-comparison/`:

- `survival-comparison.png`: an overlay and individual plots with 95% intervals;
  PT panels also show MAP-EM curves.
- `comparison.csv`: survival SSE, RMSE, integrated squared error and average
  interval width. PT MAP-EM and posterior-mean results are separate rows.
- `survival-curves.csv`: point estimates, interval limits, PT MAP curves and truth.
- `pt-convergence.csv`, `pt-posterior-diagnostics.csv`: EM convergence, posterior
  summaries, ESS and MCSE. Reaching the iteration cap is recorded as unconverged.
- `timings.csv`, `session-info.txt`: elapsed fit/prediction times and R versions.

The fiducial point estimate is its median; its survival limits are
`1 - fid$upper` and `1 - fid$lower`. MPT/TBP use the package's `Shat`, `Shatlow`
and `Shatup`. Fiducial intervals and Bayesian credible intervals have different
interpretations. These are pointwise intervals; errors on this one dataset do
not estimate repeated-sampling coverage or establish a ranking of methods.

## PT calls in these examples

`survPT()` only fits the model. Supply the numeric vector of target times to
`predict()` and read curves from its returned object:

```r
dat <- data.frame(L=c(0,0.4,1,1.5,2,3,4), R=c(0.6,0.4,2,Inf,3,4,Inf))
fit <- survPT(dat, prior=2, verbose=FALSE)
x <- seq(0,fit$Tau,length.out=101)

pred <- predict(fit,x,posterior=FALSE)
pred$Survival_Probs  # deterministic MAP-EM survival
pred$Density
fit$Posterior_Probs[fit$Terminal_Nodes]  # fitted tree masses

set.seed(2026)
post <- predict(fit,x,0.95)  # fresh sampling at prediction time
post$Posterior$Survival[,c("time","mean","lower","upper")]

draws <- samplePT(fit)  # standalone sampling, without a time grid
predict(fit,c(3,1,2),0.90,samples=draws)$Posterior$Survival
```

Bootstrap and deterministic simulation curves use `posterior=FALSE`.
`posterior.R` runs posterior prediction for both PT centers. Target times, interval levels and sampler
controls belong to `predict()` or `samplePT()`; `survPT()` no longer accepts
`test_seq`, `support`, or sampler arguments. Tree diagnostics belong to `fit`;
`Test_seq`, `Survival_Probs`, `Density` and `Posterior` belong to `pred` or `post`.

From a checkout, load the current source with `devtools::load_all(".")` before
running examples interactively. With an installed package, load `library(SurvPT)`
and use version 0.3.2 or later.

## Research scripts

The research scripts below support simulation and data-analysis work. Some execute experiments when sourced, use local interactive state, or write results into the working directory. Review their settings and dependencies before running them.

| File | Purpose | Requirements and execution notes |
| --- | --- | --- |
| `compare_tail_split.R` | Light-group comparison with and without the final split at `largest_obs` | `Rcpp`, `RcppArmadillo`, `ggplot2`, C++17 toolchain; compiles a temporary comparison implementation, reads only embedded data from `Data Analysis.R`, and writes plots and fitted results under `output/tail-split-comparison` |
| `compare_partition_versions.R` | Light-group comparison of archived 0.2.0 and both 0.2.1 tail-split variants, for exponential and uniform bases | Also requires `SurvPT_0.2.0.tar.gz` in the repository root; compiles temporary classes from the archived and current sources, checks probability conservation, and writes curves, partition paths and fitted results under `output/partition-version-comparison` |
| `diagnose_uniform_tail_em.R` | Audit uniform tail allocations and the effect of the Beta terms during EM | Run `compare_partition_versions.R` first; independently replays both saved fits with leaf-mixture algebra, runs likelihood-only EM on the same fixed partitions, and profiles the current tail mass under likelihood and penalized objectives; writes results under `output/uniform-tail-em-diagnostic` |
| `compare_tail_regularization.R` | Controlled tail-prior experiments on the revised light-group partitions | Requires `output/probability-allocation/light_check.rds`; holds data and terminal distributions fixed, compares the current update with weakened tail-path priors and a separate tail probability, verifies objective ascent and two starting fits, and writes `output/tail-regularization-comparison` without altering package code |
| `inspect_legacy_tail_partitions.R` | Inspect every archived 0.2.0 leaf entirely beyond the light-group maximum endpoint | Reads the saved version comparison, verifies the 0.2.0 archive checksum, separates full-tail and boundary-crossing contributions, and writes leaf probabilities, disjoint subtree totals and plots under `output/legacy-tail-partitions` |
| `test_function.R` | Historical estimator comparison and plotting | `devtools`, `interval`, `spBayesSurv` 1.1.9; PT models are `fit_unif`/`fit_exp` and predictions are `result_unif`/`result_exp`; optional fiducial plots require constructing `fid` |
| `bayes_helpers.R` | Shared baseline MPT/TBP fitting helpers | `survival`, `spBayesSurv` 1.1.9; sourced by `test_function.R` and `posterior.R`; defining the helpers does not run MCMC |
| `Simulation.R` | Repeated comparisons with parametric fits and interval-censoring NPMLE | `survival`, `interval`, `parallel`; launches loops immediately, requests 300 workers, writes CSV files under `n*` |
| `Simulation_CI.R` | Bootstrap survival-coverage experiments at selected quantile times | `devtools`, `interval`, `parallel`, `boot`, `bootstrap.R`; launches large nested experiments and writes `CI_n*` output |
| `bootstrap.R` | Helpers used by bootstrap comparisons | Assumes SurvPT, `interval`, and `boot` functions are attached; helper is tailored to three evaluation times and three estimators |
| `Fiducial.R` | Experimental fiducial comparison sampler | Sources its interpolation helper separately; returns a list including `point_li` and `density`; not an exported package estimator |
| `LinInterpolation.R` | Quadratic-programming interpolation for the fiducial routine | `pracma`; helper for `Fiducial.R` |
| `Fiducial_Simu.R` | Large fiducial/PT/NPMLE simulation study | `interval`, `parallel`, `pracma`, both fiducial helpers; requests 300 workers, writes `fiducial/`; contains calls expecting a numeric return from a now list-returning `fiducial()` |
| `Data Analysis.R` | Historical treatment-group analysis and bootstrap plots | `devtools`, `interval`, `parallel`, `ggplot2`, `pracma`; embeds data, defines both groups sequentially, writes RData and PNG files; review group selection and data-use permission before execution |

## Running research scripts deliberately

Use the repository root as the working directory because helper paths are relative to it. Install optional packages only for the comparison being studied. Reduce replication and bootstrap counts for a pilot run. Replace the fixed worker count with a setting appropriate to your system; Windows does not support the forked multiworker `mclapply()` workflow used in the historical simulations.

Use a dedicated output directory and set seeds explicitly. Resolve the documented interface mismatches before interpreting a script's results. The package's [regression tests](../tests/scpt.R), [quickstart](quickstart.R) and [posterior example](posterior.R) are the current validation entry points. The version-comparison and tail-diagnostic scripts in the historical inventory above are not included in this checkout and do not validate 0.3.0.

## Validation

After installing SurvPT, run these commands from the repository root:

```sh
Rscript Examples/quickstart.R
Rscript Examples/posterior.R
Rscript tests/scpt.R
Rscript tests/interfaces.R
```

The regression suite checks reference masses and asymmetric Beta parameters,
reference-count pruning, probability conservation, both MAP objective histories,
final allocation consistency, exact Beta posterior sampling, an ambiguous-censoring
posterior against numerical integration, unrestricted draws, reproducibility,
retention controls and deterministic-output preservation. It also checks automatic
exponential rates and numerical tail behavior. It does not run the full simulation studies.
The interface suite also checks target-time ordering, automatic fresh sampling
at prediction time, explicit sample reuse, serialization, and unchanged fitted state.

## Attribution

The fiducial comparison code and scenario labels relate to Cui, Hannig, and Kosorok (2024), *A Unified Nonparametric Fiducial Approach to Interval-Censored Data*, JASA 119(547), 2230–2241, [doi:10.1080/01621459.2023.2252143](https://doi.org/10.1080/01621459.2023.2252143). Preserve this attribution when adapting the routines. The local versions are modified research code, not a verified reproduction of the published experiments. See [Citations](../docs/citations.md).
