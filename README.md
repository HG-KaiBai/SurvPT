# SurvPT

[![Install from GitHub](https://img.shields.io/badge/Install-GitHub-181717?logo=github&logoColor=white)](#install-from-github)

SurvPT is an R package developed for an accompanying research manuscript currently in preparation. It estimates a univariate survival function and density from censored observations using finite Polya trees. The current package version is **0.2.0**.

The source code is available on [GitHub](https://github.com/HG-KaiBai/SurvPT) under the [MIT license](LICENSE.md).

## Fit with SurvPT

| Function | Input | Model and scope |
| --- | --- | --- |
| `survPT()` | Numeric left and right interval endpoints | Finite Polya tree with uniform, exponential, or Weibull conditional base distributions; supports exact, interval-censored, right-censored, and current-status observations |
| `estimate()`, `estimate_exp()` | Left and right endpoints | Parametric Weibull and exponential maximum-likelihood estimates, useful for choosing PT base parameters |

Start with `survPT()` and encode each observation as a pair of interval endpoints. Choose a uniform, exponential, or Weibull base distribution, fit the tree, and evaluate the survival function on a time grid. The estimator returns point estimates; it does not accept covariates or return posterior draws, confidence intervals, or a multivariate fitted distribution. Read the [method and tail definitions](docs/methods.md) before interpreting extrapolation.

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
x <- data.frame(
  L = c(0.3, 0.7, 1.2, 1.8, 2.4, 3.0),
  R = c(0.3, 1.1, 1.2, 2.5, Inf, Inf)
)
grid <- seq(0, 6, length.out = 101)

fit <- survPT(
  x[, c("L", "R")], test_seq = grid,
  prior = 2, prior_params = 0.5,
  max_M = 6, max_EM_iter = 100, verbose = FALSE
)

plot(fit$Test_seq, as.numeric(fit$Survival_Probs), type = "l",
     xlab = "Time", ylab = "Survival probability", ylim = c(0, 1))
sum(fit$Posterior_Probs[fit$Terminal_Nodes])  # approximately 1
```

`prior_params = 0.5` is an exponential **rate**, not a mean. `survPT()` reads the first two columns of `X` as interval endpoints; it does not read a status column. Use `as.numeric()` when a returned Armadillo vector appears as a one-column matrix.

A runnable example covering all three PT bases is in [Examples/quickstart.R](Examples/quickstart.R). Simulation and data-analysis scripts, their dependencies, and validation commands are described in [Examples/README.md](Examples/README.md).

## Documentation

- [Installation from GitHub or source](docs/installation.md)
- [Data encoding, fitting, and interpreting output](docs/user-guide.md)
- [Finite Polya tree implementation details](docs/methods.md)
- [Known limitations and troubleshooting](docs/limitations.md)
- [Citations and BibTeX](docs/citations.md)

Installed help is available with `help(package = "SurvPT")`, `?survPT`, and `?SurvPT` (the low-level native interface).

## Cite SurvPT

Until the manuscript citation is available, use `citation("SurvPT")` for the software reference. Record the package version and Git commit or source-archive checksum used in your analysis. See [Citations and BibTeX](docs/citations.md) for citation metadata and attribution for the comparison routines.

## Authors and contact

Authors: Yixing Zhao and Liqun Diao. Developer and maintainer: Yixing Zhao.

For questions about the package, contact Yixing Zhao at the email address in [DESCRIPTION](DESCRIPTION).
