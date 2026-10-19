# SurvPT

SurvPT is an R package for estimating a univariate survival function and density from censored observations using finite Polya trees. The current package version is **0.2.0**.

This is a **private research repository**. Downloads and GitHub installation require access from the maintainer. No software license has been selected; obtain the maintainer's permission before using or redistributing the code. See [LICENSE](LICENSE). These documents do not change the repository's visibility or grant access.

## Choose an estimator

| Function | Input | Model and scope |
| --- | --- | --- |
| `survPT()` | Numeric left and right interval endpoints | Finite Polya tree with uniform, exponential, or Weibull conditional base distributions; supports exact, interval-censored, right-censored, and current-status observations |
| `survOPT()` | One observed time and one event indicator per subject | Optional Polya tree for **univariate right-censored data**, with uniform finite terminal intervals |
| `estimate()`, `estimate_exp()` | Left and right endpoints | Parametric Weibull and exponential maximum-likelihood estimates, useful for choosing PT base parameters |

Neither tree interface accepts covariates or returns posterior draws, confidence intervals, or a multivariate fitted distribution. Read the [method and tail definitions](docs/methods.md) before interpreting extrapolation.

## Install

Install R, a compatible C++17 compilation toolchain, and the package dependencies. On Windows, use the [Rtools version matching your R installation](https://cran.r-project.org/bin/windows/Rtools/). Installation has been tested with R 4.5.1 on Windows x64.

```r
install.packages(c("Rcpp", "RcppArmadillo"))
install.packages("SurvPT_0.2.0.tar.gz", repos = NULL, type = "source")
library(SurvPT)
packageVersion("SurvPT")
```

Run this after obtaining the archive through an authorized channel, from the directory containing the archive. For installation from a private GitHub checkout, authentication, and platform notes, see [Installation](docs/installation.md). Source archives are build artifacts and may not be tracked by Git.

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

A complete runnable example, including all three PT bases and OPT, is in [Examples/quickstart.R](Examples/quickstart.R). The other historical scripts have additional dependencies and may launch large simulations; consult [Examples/README.md](Examples/README.md) before running them.

## Documentation

- [Installation and private GitHub access](docs/installation.md)
- [Data encoding, fitting, and interpreting output](docs/user-guide.md)
- [PT and OPT implementation details](docs/methods.md)
- [Known limitations and troubleshooting](docs/limitations.md)
- [Development, checks, and source-package builds](docs/development.md)
- [Citations and BibTeX](docs/citations.md)
- [Contribution workflow](CONTRIBUTING.md)
- [Changes in version 0.2.0](NEWS)

Installed help is available with `help(package = "SurvPT")`, `?survPT`, `?survOPT`, and `?SurvPT` (the low-level native interface).

## Cite the OPT methodology

When using or discussing the OPT implementation, cite:

Seok, J., Tian, L., and Wong, W. H. (2014). **Density estimation on multivariate censored data with optional Pólya tree.** *Biostatistics*, 15(1), 182–195. [doi:10.1093/biostatistics/kxt025](https://doi.org/10.1093/biostatistics/kxt025).

The article studies a multivariate method. This repository exposes a univariate right-censoring implementation related to that methodology; it is not a claim of a complete reproduction of the article. Use `citation("SurvPT")` for the software and OPT references, and record your package version and Git commit.

## Support

Authorized collaborators can report a reproducible example in the repository's private issue tracker. For access, permission, or an issue that cannot be posted there, contact Yixing Zhao at the maintainer address in [DESCRIPTION](DESCRIPTION). Include `sessionInfo()`, the function call, and synthetic data that reproduce the problem.
