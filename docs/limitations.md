# Known limitations and troubleshooting

[Back to README](../README.md)

## Statistical and interface limitations

| Limitation | Practical consequence |
| --- | --- |
| PT and OPT return point estimates only | No built-in confidence bands, posterior draws, or convergence flag |
| Inputs are univariate with no covariates | No multivariate joint survival or regression interface |
| No sampling-truncation or delayed-entry arguments | Do not treat the `support` setting as a delayed-entry correction |
| Uniform PT uses a data-dependent upper cutoff | Its estimated survival becomes zero beyond that cutoff |
| OPT can retain mass at infinity | Infinite-support output can have a positive survival plateau and a finite density integral below one |
| Partition is fixed after initialization | EM does not revisit early stopping or split locations |
| Direct exponentials and CDF differences are used | Extreme rates, shapes, or tail locations can lose precision |
| Input checks are incomplete | Validate missing values, numeric types, interval ordering, support, and iteration controls yourself |
| Exact events at the root's lower endpoint are not assigned by the strict left-bound rule | Use positive exact event times with the default zero lower support; assess another model for structural zero times |
| Finite support and observation bounds can disagree during initialization | Use a justified support and encode intervals consistently; investigate any allocation-sum diagnostic |

The field name `Posterior_Probs` is historical. Read [Methods](methods.md) for the actual PT branch-mass update and OPT mass convention.

## Simulation and legacy examples

- `data_cs()` overwrites its event-by-inspection `status` with an exact-event test. Fit current-status data using `L` and `R`; the returned `status` must not be passed to `survOPT()`.
- `data_ic()` uses a legacy inspection loop `1:Ri`. If `Ri` is zero, R evaluates `1:0` as two indices rather than an empty sequence. Do not claim exact reproduction of an intended zero-inspection design without correcting and validating that generator.
- `estimate()` and `estimate_exp()` return only `optim()$par`. Inspect a separately run optimization if diagnostics are needed; small samples or heavily censored samples may not identify stable parameters.
- Historical scripts in `Examples/` include large simulations, hard-coded worker counts, interactive viewers, and experimental comparison code. The validated starting point is `Examples/quickstart.R`; see the [example inventory](../Examples/README.md).

These limitations are documented to distinguish observed software behavior from the intended statistical model. A documentation update does not change those algorithms.

## Common problems

**GitHub returns 404 or authentication fails.** Ask the maintainer to confirm collaborator access and check your credential's repository permissions. Do not change the repository to public.

**Compiler, `make`, or headers are missing.** Install the toolchain matching your R distribution, then reinstall dependencies. See [Installation](installation.md).

**The loaded version is not the installed version.** Restart R, run `find.package("SurvPT")`, inspect `.libPaths()`, and reinstall into the intended library. `devtools::load_all()` uses a checkout, whereas `library(SurvPT)` loads an installed package.

**An evaluation grid cannot be generated.** Supply a finite sorted `test_seq`. The PT default grid relies on rows with a finite right endpoint; all-right-censored data need an explicit grid and may also contain insufficient information for the intended fit.

**“Sum of allocations exceeds parent sample size.”** This diagnostic is emitted when child and parent counts differ, even if the difference is a deficit. Check support/interval compatibility and model settings. Do not interpret the result as validated merely because the fit returned.

**A uniform curve reaches zero unexpectedly early.** Inspect the largest finite endpoint and the factor-1.1 cutoff. Use a scientifically appropriate exponential or Weibull base if an unbounded PT tail is required.

**An OPT curve does not reach zero.** Inspect the mass of its infinite terminal interval. The plateau is documented implementation behavior, not a plotting error.

**A historical simulation fails on Windows.** Scripts using `parallel::mclapply()` with multiple forked workers require adaptation. Start serially with one worker, or implement and validate a socket-cluster workflow.

## Package-check notes

The development guide records the commands used for validation. Remaining issues can include debug flags overriding user settings, GNU make extensions, and native `std::cout` output. Those are packaging and portability concerns and are distinct from passing the PT regression tests. This repository is not presented as a CRAN-ready submission.
