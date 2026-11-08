# Examples and research scripts

Run examples from the repository root in R. They use the 0.2.3 fit/predict
workflow with the 0.4.0 reference tree, alpha parameters and unspecified tail.
Results stay as R objects. Tables print in the console; plots display only
in interactive sessions. No result, workspace or image files are saved.

```r
# Load this checkout, or load an installed build with library(SurvPT).
devtools::load_all(".")
source("Examples/quickstart.R")
```

## Fit, estimate and sample

```r
dat <- data.frame(L=c(0,0.4,1,1.5,2,3,4), R=c(0.6,0.4,2,Inf,3,4,Inf))
fit <- survPT(dat, prior=2, flat=FALSE, verbose=FALSE)
x <- seq(0,fit$Tau,length.out=101)
pred <- predict(fit,x,sample=FALSE)
pred$estimate
attr(pred,"tail")

set.seed(2026)
post <- predict(fit,x,0.95)  # defaults: n_samples=2000, burn_in=3000
post[,c("x","estimate","mean","lower","upper")]
attr(post,"tail")
```

`flat=TRUE` selects Beta(1,1) throughout fitting and posterior sampling.
There is one EM fit starting at the centering distribution. Above `Tau`,
survival and density predictions and all posterior summaries are `NA`, including
at `Inf`. Tail probability `S(Tau)` is returned separately. `samplePT()` is removed.
Prediction returns summaries without raw draws. At child level `m`, the default
MAP probability bound is `1e-6 * 2^(-m)`; posterior draws are unrestricted.

## Small examples and comparisons

- `quickstart.R` fits all three body families, demonstrates flat-prior fitting,
  checks probability properties, and prints tail estimates and posterior intervals.
  Objects: `fits`, `predictions`, `flat_fit`, `checks`, `tail_probabilities`,
  `tail_predictions` (including `NA` beyond `Tau` and at `Inf`).
- `posterior.R` follows the 0.2.3 comparison of PT, spBayesSurv MPT and TBP.
  Requires `survival`, `spBayesSurv` 1.1.9 and `ggplot2`. The PT families are
  uniform, exponential and Weibull; all use the unspecified 0.4.0 tail. Editable
  PT controls default to 2000 retained draws after 3000 burn-in sweeps. MPT/TBP
  retain their package-default priors and MCMC settings. The shared initialization
  adapter is in `bayes_helpers.R`. Tables print in R, and `posterior_plot` displays
  interactively. Objects: `pt_models`, `pt_predictions`, `posterior_intervals`,
  `tail_probabilities`, `sampling_settings`, `interval_summary`, `posterior_plot`.
- `test_function.R` compares point estimates with interval NPMLE, OPT, MPT and
  TBP. It preserves the locally selected simulation and truth settings. Requires
  `devtools`, `interval` and `spBayesSurv` 1.1.9. PT curves and density predictions
  use separate `predict()` calls. `plot_comparison()` displays the plots in R.
- `Data Analysis.R` contains the treatment-group analysis and bootstrap comparison.
  The second group assignment selects the light group. Requires `devtools`,
  `interval`, `parallel`, `ggplot2`, `pracma` and `quadprog`. Results remain in
  `analysis_results` and `plot`; the existing bootstrap settings are substantial.

Posterior intervals condition on each fitted model's settings. PT, MPT and TBP
use different priors. Interval width on one dataset is not repeated-sampling
coverage, and finite chains require mixing assessment.

## Simulation scripts

These historical studies execute large loops when sourced. Reduce replications
for a pilot run. Windows uses one worker; other platforms default to two.

| File | Purpose | Results in R |
| --- | --- | --- |
| `Simulation.R` | PT/parametric/NPMLE error comparisons | `simulation_results` |
| `Simulation_CI.R` | Bootstrap survival coverage | `coverage_results` |
| `Fiducial_Simu.R` | Fiducial/PT/NPMLE comparisons | `fiducial_results`, `final` |
| `bootstrap.R` | Three-time bootstrap helper | Returned data frame |
| `Fiducial.R`, `LinInterpolation.R` | Fiducial sampler and interpolation | Returned list |

Simulation scripts require `interval` and `parallel`; bootstrap studies also
require `boot`, and fiducial comparisons require `pracma` and `quadprog`.
Bootstrap resamples can have a smaller `Tau`; predictions beyond that resample's
observed range are `NA`. Coverage at such times is unavailable for that resample.

## Validation

After loading the current package:

```r
source("tests/pt-tail-na.R")
source("Examples/quickstart.R")
source("Examples/posterior.R")
```

Regression checks cover centering initialization, asymmetric alpha, flat priors,
one-stage EM, analytic Beta and numerical posterior targets, tail output,
thinning, ordering, serialization and unchanged model state.

## Attribution

The fiducial routines relate to Cui, Hannig, and Kosorok (2024), *A Unified
Nonparametric Fiducial Approach to Interval-Censored Data*, JASA 119(547),
2230–2241, [doi:10.1080/01621459.2023.2252143](https://doi.org/10.1080/01621459.2023.2252143).
The local routines are modified research code. See [Citations](../docs/citations.md).
