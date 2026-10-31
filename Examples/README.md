# Examples and research scripts

Start with [`quickstart.R`](quickstart.R). It uses synthetic data, fits `survPT()` with uniform, exponential, and Weibull bases, checks basic probability properties, and does not write files.

The research scripts below support simulation and data-analysis work. Some execute experiments when sourced, use local interactive state, or write results into the working directory. Review their settings and dependencies before running them.

| File | Purpose | Requirements and execution notes |
| --- | --- | --- |
| `test_function.R` | Uniform/exponential PT, OPT, NPMLE, MPT, and TBP comparison and plotting | Run from the repository root; requires `devtools`, `interval`, `survival`, and `spBayesSurv` 1.1.9. Defines data, both PT fits, and a common grid through the largest finite endpoint; MPT/TBP use package-default MCMC. Fiducial calls are commented out. |
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
Rscript tests/pt-tail-root.R
```

The regression tests check conditional allocations, survival and density output, probability conservation, and tail behavior across the three PT bases. They cover 12 mixed-censoring configurations, three root-only trees, and three exact-data configurations. They do not run the full simulation studies.

## Attribution

The fiducial comparison code and scenario labels relate to Cui, Hannig, and Kosorok (2024), *A Unified Nonparametric Fiducial Approach to Interval-Censored Data*, JASA 119(547), 2230–2241, [doi:10.1080/01621459.2023.2252143](https://doi.org/10.1080/01621459.2023.2252143). Preserve this attribution when adapting the routines. The local versions are modified research code, not a verified reproduction of the published experiments. See [Citations](../docs/citations.md).
