# Known limitations and troubleshooting

[Back to README](../README.md)

## Statistical and interface limitations

| Limitation | Practical consequence |
| --- | --- |
| `survPT()` returns point estimates | No built-in confidence bands or posterior draws; its convergence flag only checks terminal-mass changes |
| Inputs are univariate with no covariates | No multivariate joint survival or regression interface |
| No sampling-truncation or delayed-entry arguments | Do not treat the `support` setting as a delayed-entry correction |
| Tail shape beyond `tau` depends on the base | Uniform survival stays at `S(tau)` with density `NA`; exponential predictions use the conditional base with a fixed supplied rate; Weibull survival and density remain `NA` |
| Partition is fixed after initialization | EM does not revisit early stopping or split locations |
| Direct exponentials and CDF differences are used | Extreme rates, shapes, or tail locations can lose precision |
| Input checks are incomplete | Validate missing values, numeric types, interval ordering, support, and iteration controls yourself |
| Exact events at the root's lower endpoint are not assigned by the strict left-bound rule | Use positive exact event times with the default zero lower support; assess another model for structural zero times |
| Finite support and observation bounds can disagree during initialization | Use a justified support and encode intervals consistently; investigate any allocation-sum diagnostic |

The field name `Posterior_Probs` is historical. Read [Methods](methods.md) for the actual PT branch-mass update.

## Simulation and legacy examples

- `data_cs()` overwrites its event-by-inspection `status` with an exact-event test. Fit current-status data using `L` and `R`; the returned `status` is not an event-by-inspection indicator.
- `data_ic()` uses a legacy inspection loop `1:Ri`. If `Ri` is zero, R evaluates `1:0` as two indices rather than an empty sequence. Do not claim exact reproduction of an intended zero-inspection design without correcting and validating that generator.
- `estimate()` and `estimate_exp()` return only `optim()$par`. Inspect a separately run optimization if diagnostics are needed; small samples or heavily censored samples may not identify stable parameters.
- Historical scripts in `Examples/` include large simulations, hard-coded worker counts, interactive viewers, and experimental comparison code. The validated starting point is `Examples/quickstart.R`; see the [example inventory](../Examples/README.md).

These limitations are documented to distinguish observed software behavior from the intended statistical model. A documentation update does not change those algorithms.

## Common problems

**GitHub installation fails.** Check the repository URL, the requested commit or tag, and the error message. For compiler or dependency failures, follow [Installation](installation.md). Report a reproducible installation problem in the [issue tracker](https://github.com/HG-KaiBai/SurvPT/issues).

**Compiler, `make`, or headers are missing.** Install the toolchain matching your R distribution, then reinstall dependencies. See [Installation](installation.md).

**The loaded version is not the installed version.** Restart R, run `find.package("SurvPT")`, inspect `.libPaths()`, and reinstall into the intended library. `devtools::load_all()` uses a checkout, whereas `library(SurvPT)` loads an installed package.

**An evaluation grid cannot be generated.** Check for valid finite left endpoints. The default grid includes right-censoring times. Infinite-support fits require the largest finite endpoint to exceed the lower support bound. All-right-censored data do not identify a finite-body distribution; its fitted mass depends on the Beta prior and mode-style updates.

**“Sum of allocations exceeds parent sample size.”** This diagnostic is emitted when child and parent counts differ, even if the difference is a deficit. Check support/interval compatibility and model settings. Do not interpret the result as validated merely because the fit returned.

**A curve stops or flattens at `tau`.** Uniform survival stays at `S(tau)` beyond `tau`, while its density is `NA`. Exponential survival and density use the conditional exponential tail with the supplied rate. Weibull survival and density remain `NA` beyond `tau`; this does not mean the tail mass is zero.

**A historical simulation fails on Windows.** Scripts using `parallel::mclapply()` with multiple forked workers require adaptation. Start serially with one worker, or implement and validate a socket-cluster workflow.

## Package-check notes

Validation commands are listed in the [example guide](../Examples/README.md#validation). Package checks can report debug flags overriding user settings, GNU make extensions, nonportable example filenames, citation-file placement, undocumented exports, and native `std::cout` output. These findings are separate from the numerical regression checks.
