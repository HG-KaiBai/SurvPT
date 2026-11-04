# SurvPT

[![Install from GitHub](https://img.shields.io/badge/Install-GitHub-181717?logo=github&logoColor=white)](#install-from-github)
[![Downloads](https://img.shields.io/github/downloads/USERNAME/REPO/total?label=Downloads&logo=github&color=181717)](https://github.com/HG-KaiBai/SurvPT/releases)
[![Views](https://view-counter.tobyhagan.com/?user=USERNAME/REPO&base=555555&accent=181717&text=white&flat=true)](https://github.com/HG-KaiBai/SurvPT)
  

SurvPT is an R package developed for an accompanying research manuscript currently in preparation. It estimates univariate survival functions from censored observations using self-consistent finite Polya trees. Version **0.4.1** uses `survPT()` for fitting and `predict()` for estimation and posterior sampling.

The source code is available on [GitHub](https://github.com/HG-KaiBai/SurvPT) under the [MIT license](LICENSE.md).

## Fit with SurvPT

| Function | Input | Model and scope |
| --- | --- | --- |
| `survPT()` | Numeric left and right interval endpoints | Finite Polya tree with uniform, exponential, or Weibull conditional base distributions; supports exact, interval-censored, right-censored, and current-status observations |
| `predict(fit, times, 0.95)` | Fitted model and a numeric vector of target times | MAP survival predictions plus posterior means and pointwise credible intervals; runs sampling during prediction by default |
| `estimate()`, `estimate_exp()` | Left and right endpoints | Parametric Weibull and exponential maximum-likelihood estimates, useful for choosing PT base parameters |

Start with `survPT()` and encode each observation as a pair of interval endpoints. Choose a uniform or exponential body center; Weibull remains available as an extension. Every family estimates only the total probability beyond the largest finite endpoint, without assuming a tail shape. Predictions beyond that endpoint return `NA`. See [Methods](docs/methods.md) for the reference prior and the two fitting algorithms.

## Install from GitHub

Run the following commands in R to install SurvPT directly from GitHub using [`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html):

```r
install.packages(c("remotes", "Rcpp", "RcppArmadillo"))
remotes::install_github("HG-KaiBai/SurvPT", upgrade = "never")
library(SurvPT)
packageVersion("SurvPT")
```

Installation compiles the package from source and requires a C++17 toolchain. On Windows, install the [Rtools version matching your R installation](https://cran.r-project.org/bin/windows/Rtools/). Installation has been tested with R 4.5.1 on Windows x64.

For installation from a local source archive or checkout, reproducible version selection, and platform notes, see [Installation](docs/installation.md). Source archives are build artifacts and may not be tracked by Git.

## First fit

```r
library(SurvPT)

# Each row represents one event time: exact, interval, or right censored.
dat <- data.frame(
  L = c(0.3, 0.7, 1.2, 1.8, 2.4, 3.0),
  R = c(0.3, 1.1, 1.2, 2.5, Inf, Inf)
)
x <- seq(0, 3, length.out = 101)  # target times, not observation intervals

fit <- survPT(
  dat[, c("L", "R")],
  prior = 2, prior_params = 0.5,
  max_M = 6, verbose = FALSE
)

pred <- predict(fit, x, sample = FALSE)
plot(pred$x, pred$estimate, type = "l",
     xlab = "Time", ylab = "Survival probability", ylim = c(0, 1))
sum(fit$Posterior_Probs[fit$Terminal_Nodes])  # approximately 1
fit$Tail_Prob                              # S(tau)
fit$Convergence$map$converged
```

`prior_params = 0.5` is an exponential **rate**, not a mean. `survPT()` reads the first two columns of `X` as interval endpoints; it does not read a status column. Use `as.numeric()` when a returned Armadillo vector appears as a one-column matrix.

## Add posterior means and uncertainty

```r
set.seed(2026)
pred <- predict(fit, x, 0.95) # samples now, then evaluates at exactly x
pred                        # x, estimate, lower, upper, mean, median
attr(pred, "tail")           # fitted and posterior summaries of S(tau)

# Defaults: 2000 posterior draws after 3000 burn-in sweeps.
# Only curve and tail summaries are returned; raw draws are not retained.
```

Sampling alternates discrete leaf labels and unrestricted conjugate Beta draws.
Intervals are pointwise and conditional on the retained tree and fixed base
parameters. `survPT()` only fits the model. `predict()` runs a fresh chain
unless `sample=FALSE` or `n_samples=0`. Predictions preserve
the order of repeated or unsorted target times and do not refit the model.
Posterior draws are used internally to compute summaries. Tail reference
mass `s0`, prior strength `gamma`, calibration `kappa`, and pruning `min_count`
can be configured. See [the guide](docs/user-guide.md) and
[posterior example](Examples/posterior.R).

Allocation starts directly from the centering reference distribution. There is
one MAP-EM fit with the selected prior. Set `flat=TRUE` in `survPT()` to use
Beta(1,1) at every split for both fitting and posterior sampling; the default
`flat=FALSE` uses the depth- and reference-dependent 0.4.0 alpha parameters.
At child level `m`, fitted branch probabilities are bounded by
`split_epsilon * 2^(-m)` and `1 - split_epsilon * 2^(-m)`; the base scale
defaults to `1e-6`. Posterior Beta draws remain unrestricted.
Examples print results and keep objects and plots in the R session without
automatically saving result files.

A runnable example covering all three PT bases is in [Examples/quickstart.R](Examples/quickstart.R). Simulation and data-analysis scripts, their dependencies, and validation commands are described in [Examples/README.md](Examples/README.md).

## Documentation

- [Installation from GitHub or source](docs/installation.md)
- [Data encoding, fitting, and interpreting output](docs/user-guide.md)
- [Finite Polya tree implementation details](docs/methods.md)
- [Known limitations and troubleshooting](docs/limitations.md)
- [Citations and BibTeX](docs/citations.md)

Installed help is available with `help(package = "SurvPT")`, `?survPT`, and `?predict.survPT`.

## Cite SurvPT

Until the manuscript citation is available, use `citation("SurvPT")` for the software reference. Record the package version and Git commit or source-archive checksum used in your analysis. See [Citations and BibTeX](docs/citations.md) for citation metadata and attribution for the comparison routines.

## Authors and contact

Authors: Yixing Zhao and Liqun Diao. Developer and maintainer: Yixing Zhao.

For questions about the package, contact Yixing Zhao at the email address in [DESCRIPTION](DESCRIPTION).
