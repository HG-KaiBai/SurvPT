# Examples and research scripts

Start with [`quickstart.R`](quickstart.R). It uses synthetic data, fits `survPT()` with uniform, exponential, and Weibull bases, checks basic probability properties, and does not write files.

The research scripts below support simulation and data-analysis work. Some execute experiments when sourced, use local interactive state, or write results into the working directory. Review their settings and dependencies before running them.

| File | Purpose | Requirements and execution notes |
| --- | --- | --- |
| `compare_tail_split.R` | Light-group comparison with and without the final split at `largest_obs` | `Rcpp`, `RcppArmadillo`, `ggplot2`, C++17 toolchain; compiles a temporary comparison implementation, reads only embedded data from `Data Analysis.R`, and writes plots and fitted results under `output/tail-split-comparison` |
| `compare_partition_versions.R` | Light-group comparison of archived 0.2.0 and both 0.2.1 tail-split variants, for exponential and uniform bases | Also requires `SurvPT_0.2.0.tar.gz` in the repository root; compiles temporary classes from the archived and current sources, checks probability conservation, and writes curves, partition paths and fitted results under `output/partition-version-comparison` |
| `diagnose_uniform_tail_em.R` | Audit uniform tail allocations and the effect of the Beta terms during EM | Run `compare_partition_versions.R` first; independently replays both saved fits with leaf-mixture algebra, runs likelihood-only EM on the same fixed partitions, and profiles the current tail mass under likelihood and penalized objectives; writes results under `output/uniform-tail-em-diagnostic` |
| `compare_tail_regularization.R` | Controlled tail-prior experiments on the revised light-group partitions | Requires `output/probability-allocation/light_check.rds`; holds data and terminal distributions fixed, compares the current update with weakened tail-path priors and a separate tail probability, verifies objective ascent and two starting fits, and writes `output/tail-regularization-comparison` without altering package code |
| `inspect_legacy_tail_partitions.R` | Inspect every archived 0.2.0 leaf entirely beyond the light-group maximum endpoint | Reads the saved version comparison, verifies the 0.2.0 archive checksum, separates full-tail and boundary-crossing contributions, and writes leaf probabilities, disjoint subtree totals and plots under `output/legacy-tail-partitions` |
| `test_function.R` | Historical estimator comparison and plotting | `devtools`, `interval`; calls `View()` and refers to `fid` although its construction is commented out |
| `Simulation.R` | Repeated comparisons with parametric fits and interval-censoring NPMLE | `survival`, `interval`, `parallel`; launches loops immediately, requests 300 workers, writes CSV files under `n*` |
| `Simulation_CI.R` | Bootstrap survival-coverage experiments at selected quantile times | `devtools`, `interval`, `parallel`, `boot`, `bootstrap.R`; launches large nested experiments and writes `CI_n*` output |
| `bootstrap.R` | Helpers used by bootstrap comparisons | Assumes SurvPT, `interval`, and `boot` functions are attached; helper is tailored to three evaluation times and three estimators |
| `Fiducial.R` | Experimental fiducial comparison sampler | Sources its interpolation helper separately; returns a list including `point_li` and `density`; not an exported package estimator |
| `LinInterpolation.R` | Quadratic-programming interpolation for the fiducial routine | `pracma`; helper for `Fiducial.R` |
| `Fiducial_Simu.R` | Large fiducial/PT/NPMLE simulation study | `interval`, `parallel`, `pracma`, both fiducial helpers; requests 300 workers, writes `fiducial/`; contains calls expecting a numeric return from a now list-returning `fiducial()` |
| `Data Analysis.R` | Historical treatment-group analysis and bootstrap plots | `devtools`, `interval`, `parallel`, `ggplot2`, `pracma`; embeds data, defines both groups sequentially, writes RData and PNG files; review group selection and data-use permission before execution |

## Running research scripts deliberately

Use the repository root as the working directory because helper paths are relative to it. Install optional packages only for the comparison being studied. Reduce replication and bootstrap counts for a pilot run. Replace the fixed worker count with a setting appropriate to your system; Windows does not support the forked multiworker `mclapply()` workflow used in the historical simulations.

Use a dedicated output directory and set seeds explicitly. Resolve the documented interface mismatches before interpreting a script's results. The package's [regression tests](../tests/pt-conditional-base.R) and [quickstart](quickstart.R) are the current validation entry points.

## Validation

After installing SurvPT, run these commands from the repository root:

```sh
Rscript Examples/quickstart.R
Rscript tests/pt-conditional-base.R
Rscript tests/pt-tree-structure.R
```

The regression tests check initial and subsequent conditional allocations, survival and density output, probability conservation, and tail behavior across the three PT bases. They cover 12 mixed-censoring configurations, three trees with only the final tail split, three exact-data configurations, endpoint conventions, and default grids. Tree tests check truncated medians, the additional tail layer and its Beta parameters, equal conditional branch probabilities at initialization, early stopping, numerical tail stability, uniform tail factors for finite and right-censored observations, and the separate overlap-count option. They do not run the full simulation studies.

## Attribution

The fiducial comparison code and scenario labels relate to Cui, Hannig, and Kosorok (2024), *A Unified Nonparametric Fiducial Approach to Interval-Censored Data*, JASA 119(547), 2230–2241, [doi:10.1080/01621459.2023.2252143](https://doi.org/10.1080/01621459.2023.2252143). Preserve this attribution when adapting the routines. The local versions are modified research code, not a verified reproduction of the published experiments. See [Citations](../docs/citations.md).
