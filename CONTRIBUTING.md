# Contributing to SurvPT

This repository remains private. Contributions are for authorized collaborators and do not change its visibility or grant a software license. Follow the permission statement in [LICENSE](LICENSE).

## Discuss and reproduce

Report a concise problem statement, expected and observed behavior, package version, `sessionInfo()`, and a minimal synthetic dataset. Explain the censoring encoding and the chosen support and base parameters. Do not attach identifiable research data or credentials.

## Prepare a change

1. Work on a focused branch in an authorized checkout.
2. For a behavior change, add a regression that fails on the previous implementation and checks a statistical property or an independently computed result.
3. Edit help text in the roxygen comments under `R/`, then regenerate `man/` and `NAMESPACE`. Do not hand-edit Rcpp-generated interfaces.
4. Update the relevant guide, runnable example, and `NEWS` when input, output, assumptions, or numerical results change.
5. Follow [Development](docs/development.md) to build and check an isolated source archive.

Write documentation and comments in English. Distinguish verified behavior from intended methods. For OPT, retain the Seok, Tian, and Wong (2014) citation and the univariate scope qualification. Preserve attribution when adapting the fiducial comparison code.

## Review checklist

- Inputs, returned fields, index conventions, and parameter units are documented.
- Exact observations, censoring intervals, boundary values, and tails are covered where relevant.
- Node probabilities and allocation counts conserve mass under the tested model.
- Examples run without private data, interactive state, or optional large simulations.
- Citation metadata and version numbers agree.
- No source archive, compiled object, private dataset, token, or local IDE state is accidentally staged.

Submit the change through the private repository's usual review process. Keep release archives and any published documentation within channels approved by the maintainer; do not enable public hosting or change repository visibility as part of a documentation or packaging change.
