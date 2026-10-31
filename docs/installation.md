# Installation

[Back to README](../README.md)

## Requirements

SurvPT is a source R package containing C++17 code. It requires `Rcpp` (at least 1.1.0), `RcppArmadillo`, and a working compiler toolchain. R's standard `stats` package supplies simulation and optimization functions. A local validation used R 4.5.1, Windows x64, and Rtools45; other platform and version combinations have not been validated in this repository.

- **Windows:** install the [Rtools release appropriate for your R version](https://cran.r-project.org/bin/windows/Rtools/).
- **macOS:** install the compiler tools required by your R distribution; consult [R for macOS tools](https://mac.r-project.org/tools/).
- **Linux:** install R development headers, `make`, and a C++17 compiler through your distribution. BLAS, LAPACK, and the Fortran runtime must be available to R's linker.

The package uses R's BLAS/LAPACK and OpenMP make variables. Its current makefiles also add debug compilation flags. Compiler performance and compatibility can therefore differ across platforms; see [package-check notes](limitations.md#package-check-notes).

## Install from GitHub

Install the package directly from the repository:

```r
install.packages(c("remotes", "Rcpp", "RcppArmadillo"))
remotes::install_github("HG-KaiBai/SurvPT", upgrade = "never")
library(SurvPT)
```

For reproducibility, pass an existing commit SHA with `ref = "YOUR_COMMIT_SHA"`. SurvPT is distributed under the [MIT license](../LICENSE.md).

## Install a local source archive

If you have a source archive such as `SurvPT_0.3.2.tar.gz`, set the R working directory to its location or supply an absolute path:

```r
install.packages(c("Rcpp", "RcppArmadillo"))
install.packages("SurvPT_0.3.2.tar.gz", repos = NULL, type = "source")
library(SurvPT)
packageVersion("SurvPT")
```

`repos = NULL` does not fetch missing dependencies for this local archive, so install those first. Restart R before replacing a loaded package, particularly on Windows where a loaded DLL can be locked.

## Install from a checkout

Clone the repository, then install from its root directory:

```sh
git clone https://github.com/HG-KaiBai/SurvPT.git
cd SurvPT
R CMD INSTALL --preclean .
```

Run the dependency-installation command above first. A GitHub ZIP download is a checkout snapshot, not an R source-package archive: extract it and run `R CMD INSTALL` on the directory containing `DESCRIPTION`.

## Confirm the installation

```r
library(SurvPT)
packageVersion("SurvPT")
find.package("SurvPT")
help(package = "SurvPT")
citation("SurvPT")
```

Run [the quickstart](../Examples/quickstart.R) from a checkout using `Rscript Examples/quickstart.R`. It needs only the installed package and standard R packages and prints a compact validation summary.

## Reproducibility

Version 0.2.0 changes censored-data allocations for exponential and Weibull PT bases and for the uniform infinite terminal interval. Refit analyses made with 0.1.0 instead of assuming identical estimates. Record `packageVersion("SurvPT")`, `sessionInfo()`, the chosen support, base parameters, stopping controls, and source commit or archive checksum.

Version 0.3.0 adds an internal Beta-prior body/tail split at the largest finite observed endpoint for unbounded support. It replaces the uniform cutoff and leaves tail shape unspecified for all bases. Cached leaf probabilities accelerate EM while preserving its update rules.

Version 0.3.1 carries uniform survival forward as `S(tau)` at finite times beyond `tau` and uses a conditional exponential tail for the exponential base. Weibull tail predictions remain unchanged. Refit earlier analyses to use the updated tail predictions; see [Methods](methods.md).

Version 0.3.2 bounds all PT split probabilities to `[1e-6, 1 - 1e-6]`, uses a `Beta(1, 1)` body/tail prior, and sets finite-tree prior parameters to `alpha_m = max(1, 0.5 * m^2)`. These changes affect fitted masses for both bounded and unbounded support; refit analyses and record the source commit or archive checksum.
