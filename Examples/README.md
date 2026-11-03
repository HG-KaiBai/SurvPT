# Examples and research scripts

Start with [`quickstart.R`](quickstart.R). It uses synthetic data, fits `survPT()` with uniform, exponential, Weibull, and spliced uniform-exponential bases, checks basic probability properties, and does not write files.

For posterior intervals, run `source("Examples/posterior.R")` in an interactive R session after installing SurvPT 0.2.3, `spBayesSurv` 1.1.9, `survival`, and `ggplot2`. The example compares three PT bases with the MPT and TBP models from `spBayesSurv` on the same simulated data and time grid. PT chains retain 5,000 draws after 3,000 burn-in sweeps. Both this example and `test_function.R` use the package-default MCMC settings for spBayesSurv: in version 1.1.9, 2,000 retained draws after 3,000 burn-in sweeps, `nskip = 0`, and progress every 500 retained draws. Its default parametric initialization chain is also retained. The posterior example prints the main-chain settings, selected-time intervals and width summaries, and displays a faceted ggplot with credible ribbons and simulated truth. Full results remain in `posterior_intervals`, with the plot in `posterior_plot`; no result or image files are saved. Noninteractive runs leave the plot in memory without opening a graphics device.

The research scripts below support simulation and data-analysis work. Some execute experiments when sourced, use local interactive state, or write results into the working directory. Review their settings and dependencies before running them.

| File | Purpose | Requirements and execution notes |
| --- | --- | --- |
| `posterior.R` | PT, MPT, and TBP posterior survival interval comparison | SurvPT >= 0.2.3, `spBayesSurv` 1.1.9, `survival`, `ggplot2`; reuses the initialization adaptation from `test_function.R`; prints tables and an interactive ggplot, writes no result files |
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
Rscript tests/pt-uniform-tail-factor.R
Rscript tests/pt-uniform-exp.R
Rscript tests/pt-predict.R
```

The regression tests check conditional allocations, survival and density output, probability conservation, and tail behavior. They cover 20 mixed-censoring configurations for priors 1-3, uniform cutoff controls, and 24 configurations of the spliced uniform-exponential base, along with endpoint and interface checks. Prediction checks compare exact and censored posterior draws with analytic distributions and verify reproducibility and fit/predict separation. They do not run the full simulation studies.

## Attribution

The fiducial comparison code and scenario labels relate to Cui, Hannig, and Kosorok (2024), *A Unified Nonparametric Fiducial Approach to Interval-Censored Data*, JASA 119(547), 2230–2241, [doi:10.1080/01621459.2023.2252143](https://doi.org/10.1080/01621459.2023.2252143). Preserve this attribution when adapting the routines. The local versions are modified research code, not a verified reproduction of the published experiments. See [Citations](../docs/citations.md).
