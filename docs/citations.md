# Citations

[Back to README](../README.md)

## Software

Use `citation("SurvPT")` after installing the package. Record version 0.2.0 and the Git commit or source-archive checksum used in your analysis. The package's software authorship is listed in `DESCRIPTION`; the authors of the methodological references are not being identified as authors or maintainers of this repository.

The repository also provides [CITATION.cff](../CITATION.cff) for GitHub and [inst/CITATION](../inst/CITATION) for R. A citation is not permission to redistribute a private repository; see [LICENSE](../LICENSE).

## OPT methodology

**Seok, Junhee; Tian, Lu; and Wong, Wing H. (2014).** Density estimation on multivariate censored data with optional Pólya tree. *Biostatistics*, **15**(1), 182–195. DOI: [10.1093/biostatistics/kxt025](https://doi.org/10.1093/biostatistics/kxt025).

Verified against the [publisher record](https://academic.oup.com/biostatistics/article-abstract/15/1/182/244923) and [PubMed record, PMID 23902636](https://pubmed.ncbi.nlm.nih.gov/23902636/). The journal issue year is **2014**; online publication was July 31, 2013. Use 2014 with the volume, issue, and page citation above.

```bibtex
@article{SeokTianWong2014,
  author  = {Seok, Junhee and Tian, Lu and Wong, Wing H.},
  title   = {Density estimation on multivariate censored data with optional {P\'olya} tree},
  journal = {Biostatistics},
  year    = {2014},
  volume  = {15},
  number  = {1},
  pages   = {182--195},
  doi     = {10.1093/biostatistics/kxt025}
}
```

Cite this paper when reporting results from `survOPT()` or discussing the associated censored-data OPT methodology. The implementation in this package is univariate and right-censored; the paper's multivariate scope should not be presented as a feature of the package. Implementation differences and the infinite-tail convention are described in [Methods](methods.md).

## Fiducial comparison code

The historical comparison scripts and simulation-scenario names refer to:

**Cui, Yifan; Hannig, Jan; and Kosorok, Michael R. (2024).** A Unified Nonparametric Fiducial Approach to Interval-Censored Data. *Journal of the American Statistical Association*, **119**(547), 2230–2241. [doi:10.1080/01621459.2023.2252143](https://doi.org/10.1080/01621459.2023.2252143).

The [publisher record](https://www.tandfonline.com/doi/full/10.1080/01621459.2023.2252143) gives online publication in December 2023 and a 2024 journal issue. Cite this work when using or adapting those comparison routines. The local routines include modifications and have not been certified as reproductions of the authors' supplementary software or reported experiments. The exported `fiducial_1()` through `fiducial_4()` functions generate data; they do not run fiducial inference.
