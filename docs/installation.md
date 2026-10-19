# Installation

[Back to README](../README.md)

## Requirements

SurvPT is a source R package containing C++17 code. It requires `Rcpp` (at least 1.1.0), `RcppArmadillo`, and a working compiler toolchain. R's standard `stats` package supplies simulation and optimization functions. A local validation used R 4.5.1, Windows x64, and Rtools45; other platform and version combinations have not been validated in this repository.

- **Windows:** install the [Rtools release appropriate for your R version](https://cran.r-project.org/bin/windows/Rtools/).
- **macOS:** install the compiler tools required by your R distribution; consult [R for macOS tools](https://mac.r-project.org/tools/).
- **Linux:** install R development headers, `make`, and a C++17 compiler through your distribution. BLAS, LAPACK, and the Fortran runtime must be available to R's linker.

The package uses R's BLAS/LAPACK and OpenMP make variables. Its current makefiles also add debug compilation flags. Compiler performance and compatibility can therefore differ across platforms; see [development notes](development.md).

## Install an authorized source archive

Download the source archive supplied by the maintainer. In R, set the working directory to its location or supply an absolute path:

```r
install.packages(c("Rcpp", "RcppArmadillo"))
install.packages("SurvPT_0.2.0.tar.gz", repos = NULL, type = "source")
library(SurvPT)
packageVersion("SurvPT")
```

`repos = NULL` does not fetch missing dependencies for this local archive, so install those first. Restart R before replacing a loaded package, particularly on Windows where a loaded DLL can be locked.

## Install from a checkout

After the maintainer grants repository access, authenticate with your normal Git credential manager or SSH configuration, then clone:

```sh
git clone https://github.com/HG-KaiBai/SurvPT.git
cd SurvPT
R CMD INSTALL --preclean .
```

Run the dependency-installation command above first. A GitHub ZIP download is a checkout snapshot, not an R source-package archive: extract it and run `R CMD INSTALL` on the directory containing `DESCRIPTION`.

## Optional installation with remotes

The [`remotes` documentation](https://remotes.r-lib.org/reference/install_github.html) describes private-repository authentication via `GITHUB_PAT`. Configure a credential with access to this repository outside your script; do not paste tokens into tracked files, notebooks, or issue reports.

```r
install.packages("remotes")
remotes::install_github("HG-KaiBai/SurvPT", upgrade = "never")
```

For reproducibility, pass an existing commit SHA with `ref = "YOUR_COMMIT_SHA"`. The command does not make the repository public. A 404 may mean that your account or token lacks access; it is not a reason to change repository visibility.

## Confirm the installation

```r
library(SurvPT)
packageVersion("SurvPT")
find.package("SurvPT")
help(package = "SurvPT")
citation("SurvPT")
```

Run [the quickstart](../Examples/quickstart.R) from a checkout using `Rscript Examples/quickstart.R`. It needs only the installed package and standard R packages and prints a compact validation summary.

## Updating an existing installation

Version 0.2.0 changes censored-data allocations for exponential and Weibull PT bases and for the uniform infinite terminal interval. Refit analyses made with 0.1.0 instead of assuming identical estimates. Record `packageVersion("SurvPT")`, `sessionInfo()`, the chosen support, base parameters, stopping controls, and source commit or archive checksum.
