# Example and research-script inventory

Start with [`quickstart.R`](quickstart.R). It uses synthetic data, runs all three PT bases and the univariate right-censored OPT interface, checks basic probability properties, and does not write files. OPT methodology: Seok, Tian, and Wong (2014), [doi:10.1093/biostatistics/kxt025](https://doi.org/10.1093/biostatistics/kxt025).

The other files are **historical research scripts**, not installation tests or turnkey reproductions of published results. Read them before sourcing: some execute their experiments immediately, use local interactive state, or write results into the working directory. Their statistical calculations have not been revalidated by this documentation update.

| File | Purpose | Requirements and execution notes |
| --- | --- | --- |
| `test_function.R` | Interactive PT/OPT comparison and plotting | `devtools`, `interval`; calls `View()` and refers to `fid` although its construction is commented out |
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

## Attribution

The fiducial comparison code and scenario labels relate to Cui, Hannig, and Kosorok (2024), *A Unified Nonparametric Fiducial Approach to Interval-Censored Data*, JASA 119(547), 2230–2241, [doi:10.1080/01621459.2023.2252143](https://doi.org/10.1080/01621459.2023.2252143). Preserve this attribution when adapting the routines. The local versions are modified research code, not a verified reproduction of the published experiments. See [Citations](../docs/citations.md).
