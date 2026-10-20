# Citations

[Back to README](../README.md)

## Software

SurvPT was developed for an accompanying research manuscript currently in preparation. Until its bibliographic details are available, use `citation("SurvPT")` for the software reference. Record the package version and Git commit or source-archive checksum used in your analysis. Software authorship is listed in `DESCRIPTION`.

The repository also provides [CITATION.cff](../CITATION.cff) for GitHub and [inst/CITATION](../inst/CITATION) for R. SurvPT is distributed under the [MIT license](../LICENSE.md).

To export the installed software citation as BibTeX:

```r
toBibtex(citation("SurvPT"))
```

## Fiducial comparison code

The historical comparison scripts and simulation-scenario names refer to:

**Cui, Yifan; Hannig, Jan; and Kosorok, Michael R. (2024).** A Unified Nonparametric Fiducial Approach to Interval-Censored Data. *Journal of the American Statistical Association*, **119**(547), 2230–2241. [doi:10.1080/01621459.2023.2252143](https://doi.org/10.1080/01621459.2023.2252143).

The [publisher record](https://www.tandfonline.com/doi/full/10.1080/01621459.2023.2252143) gives online publication in December 2023 and a 2024 journal issue. Cite this work when using or adapting those comparison routines. The local routines include modifications and have not been certified as reproductions of the authors' supplementary software or reported experiments. The exported `fiducial_1()` through `fiducial_4()` functions generate data; they do not run fiducial inference.
