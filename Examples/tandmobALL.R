# Teeth 14, 24, 34 and 44: compare SurvPT, Fiducial, MPT and TBP in boys and girls.
# Run from the repository root:
#   devtools::load_all(".")
#   source("Examples/tandmobALL.R")
# Reuse the functions with L/R data frames (NA denotes a censored endpoint):
#   fit_tandmob_group(boy, group = "boy")
#   compare_tandmob(boy, girl, tooth = 14)
#   analyze_tandmob(tandmobAll, teeth = c(14, 24, 34, 44))
# Optional reuse of completed teeth before sourcing:
#   options(SurvPT.tandmob_existing = list("44" = completed_tooth_44))
# Optional export before sourcing:
#   options(SurvPT.tandmob_output_dir = "output/tandmobALL")

source("Examples/bayes_helpers.R", local = TRUE)
source("Examples/LinInterpolation.R", local = TRUE)
source("Examples/Fiducial.R", local = TRUE)

fit_tandmob_group <- function(data, group, x = NULL, seed = 20261109L,
                              pt_samples = 2000L, pt_burn_in = 3000L, pt_thin = 1L,
                              fid_samples = 1000L, fid_burn_in = 100L) {
  # The existing fiducial helper returns fixed 95% equal-tailed intervals.
  level <- 0.95
  uninformative <- is.na(data$L) & is.na(data$R)
  excluded <- sum(uninformative)
  data <- data[!uninformative, c("L", "R"), drop = FALSE]
  data$L[is.na(data$L)] <- 0
  data$R[is.na(data$R)] <- Inf
  rownames(data) <- NULL
  n <- nrow(data)
  censoring_summary <- data.frame(
    group = group, n = n, left = sum(data$L == 0 & is.finite(data$R)),
    interval = sum(data$L > 0 & is.finite(data$R) & data$L < data$R),
    right = sum(is.infinite(data$R)), exact = sum(data$L == data$R),
    excluded_both_missing = excluded)

  endpoints <- unlist(data, use.names = FALSE)
  tau <- max(endpoints[is.finite(endpoints)])
  if (is.null(x)) x <- seq(0, tau, length.out = 201L)
  grid <- seq(0, tau, length.out = 100L)
  method_order <- c("SurvPT", "Fiducial", "MPT", "TBP")
  seeds <- setNames(seed + 0:3, method_order)
  run_times <- setNames(numeric(4L), method_order)

  # Elapsed time includes fitting, sampling and interval summaries.
  # Shared data preparation, printing and graphics are outside these timers.
  cat("\n", group, ": fitting and sampling SurvPT ...\n", sep = "")
  set.seed(seeds[["SurvPT"]])
  started <- proc.time()[["elapsed"]]
  exp_par <- SurvPT::estimate_exp(data$L, data$R)
  pt_model <- SurvPT::survPT(
    data, prior = 2, prior_params = exp_par, max_M = 8L,
    max_EM_iter = 2000L, EM_tol = 1e-6, verbose = FALSE)
  pt_prediction <- predict(pt_model, x, level = level, sample = TRUE,
                           n_samples = pt_samples, burn_in = pt_burn_in,
                           thin = pt_thin)
  run_times["SurvPT"] <- proc.time()[["elapsed"]] - started

  cat("\n", group, ": sampling Fiducial ...\n", sep = "")
  set.seed(seeds[["Fiducial"]])
  started <- proc.time()[["elapsed"]]
  fid <- fiducial(data$L, data$R, nsample = n, grid = grid, testgrid = x,
                  mfid = fid_samples, mburn = fid_burn_in)
  run_times["Fiducial"] <- proc.time()[["elapsed"]] - started

  # Surv(interval2) uses NA for the unobserved endpoint.
  bayes_data <- data
  bayes_data$L[bayes_data$L == 0 & bayes_data$R > 0] <- NA_real_
  bayes_data$R[is.infinite(bayes_data$R)] <- NA_real_
  positive <- x > 0

  # MPT/TBP retain package-default priors and MCMC settings.
  # Their timers include the package's parametric initialization.
  cat("\n", group, ": fitting and sampling MPT ...\n", sep = "")
  set.seed(seeds[["MPT"]])
  started <- proc.time()[["elapsed"]]
  mpt_model <- fit_mpt(bayes_data)
  mpt_curves <- spBayesSurv::plot.survregbayes2(
    mpt_model, xpred = matrix(0, 1L, 1L), tgrid = x[positive],
    CI = level, PLOT = FALSE)
  run_times["MPT"] <- proc.time()[["elapsed"]] - started

  cat("\n", group, ": fitting and sampling TBP ...\n", sep = "")
  set.seed(seeds[["TBP"]])
  started <- proc.time()[["elapsed"]]
  tbp_model <- fit_tbp(bayes_data)
  tbp_curves <- spBayesSurv::plot.survregbayes(
    tbp_model, tgrid = x[positive], CI = level, PLOT = FALSE)
  run_times["TBP"] <- proc.time()[["elapsed"]] - started

  bayes_table <- function(curves, method) {
    result <- data.frame(method = method, x = x, center = 1, lower = 1, upper = 1)
    result$center[positive] <- as.numeric(curves$Shat[, 1L])
    result$lower[positive] <- as.numeric(curves$Shatlow[, 1L])
    result$upper[positive] <- as.numeric(curves$Shatup[, 1L])
    result
  }
  survival_intervals <- rbind(
    data.frame(method = "SurvPT", x = x, center = pt_prediction$mean,
               lower = pt_prediction$lower, upper = pt_prediction$upper),
    # S = 1 - F reverses the CDF bounds.
    data.frame(method = "Fiducial", x = x, center = 1 - fid$point_li,
               lower = 1 - fid$upper, upper = 1 - fid$lower),
    bayes_table(mpt_curves, "MPT"),
    bayes_table(tbp_curves, "TBP"))
  survival_intervals$group <- group
  survival_intervals$method <- factor(survival_intervals$method, levels = method_order)
  survival_intervals$interval_type <- ifelse(
    survival_intervals$method == "Fiducial", "Fiducial", "Posterior credible")
  survival_intervals$level <- level
  survival_intervals$width <- survival_intervals$upper - survival_intervals$lower
  sampling_settings <- data.frame(
    group = group, method = method_order, seed = unname(seeds),
    retained = c(pt_samples, fid_samples, mpt_model$mcmc$nsave, tbp_model$mcmc$nsave),
    burn_in = c(pt_burn_in, fid_burn_in, mpt_model$mcmc$nburn, tbp_model$mcmc$nburn),
    thin = c(pt_thin, 1L, mpt_model$mcmc$nskip + 1L, tbp_model$mcmc$nskip + 1L),
    elapsed_seconds = unname(run_times))
  cat("\n", group, ": completed all four methods.\n", sep = "")
  print(sampling_settings, row.names = FALSE)

  list(
    data = data, censoring_summary = censoring_summary,
    models = list(SurvPT = pt_model, Fiducial = fid, MPT = mpt_model, TBP = tbp_model),
    pt_prediction = pt_prediction, survival_intervals = survival_intervals,
    tail_probability = data.frame(group = group, attr(pt_prediction, "tail")),
    sampling_settings = sampling_settings)
}

plot_tandmob <- function(survival_intervals, censoring_summary) {
  palette <- c(SurvPT = "#0072B2", Fiducial = "#D55E00",
               MPT = "#009E73", TBP = "#CC79A7")
  point_times <- sort(unique(survival_intervals$x))
  point_times <- point_times[seq(1L, length(point_times), by = 16L)]
  point_data <- survival_intervals[
    survival_intervals$x %in% point_times & survival_intervals$x >= 6, ]
  tooth_labels <- vapply(split(censoring_summary, censoring_summary$tooth), function(d) {
    paste0("Tooth ", d$tooth[1L], "\n",
           paste(d$group, "n =", d$n, collapse = "; "))
  }, character(1L))

  ggplot2::ggplot(survival_intervals,
                  ggplot2::aes(x = x, colour = method,
                               group = interaction(tooth, method, group))) +
    ggplot2::geom_line(ggplot2::aes(y = lower), linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = upper), linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = center), linetype = "solid", linewidth = 0.8) +
    ggplot2::geom_point(data = point_data, ggplot2::aes(y = center, shape = group),
                        size = 2.2) +
    ggplot2::scale_colour_manual(values = palette, name = "Method") +
    ggplot2::scale_shape_manual(values = c(boy = 16, girl = 17), name = "Group") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, override.aes = list(shape = NA)),
      shape = ggplot2::guide_legend(order = 2, override.aes = list(colour = "black"))) +
    ggplot2::coord_cartesian(xlim = range(survival_intervals$x), ylim = c(0, 1)) +
    ggplot2::facet_wrap(~ tooth, ncol = 2L,
                        labeller = ggplot2::as_labeller(tooth_labels)) +
    ggplot2::labs(
      title = "Tooth emergence in boys and girls",
      subtitle = "Solid: estimates; dashed: pointwise 95% intervals",
      x = "Age (years)", y = "Probability tooth has not emerged",
      caption = paste("SurvPT: posterior mean; Fiducial: interpolated median; MPT/TBP: package Shat.",
                      "Intervals: posterior credible intervals, except Fiducial (fiducial intervals).",
                      sep = "\n")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   plot.caption = ggplot2::element_text(hjust = 0))
}

compare_tandmob <- function(boy, girl, seed = 123,
                            pt_samples = 2000L, pt_burn_in = 3000L, pt_thin = 1L,
                            fid_samples = 1000L, fid_burn_in = 100L, tooth = 44L) {
  groups <- list(boy = boy, girl = girl)
  # Both groups use the same grid, within the observable range of each PT fit.
  group_tau <- vapply(groups, function(data) {
    endpoints <- unlist(data[, c("L", "R")], use.names = FALSE)
    max(endpoints[is.finite(endpoints)])
  }, numeric(1L))
  x <- seq(0, min(group_tau), length.out = 201L)

  # Sequential execution keeps the runtime measurements free of competing fits.
  group_results <- lapply(seq_along(groups), function(i) {
    fit_tandmob_group(
      groups[[i]], group = names(groups)[i], x = x, seed = seed + 4L * (i - 1L),
      pt_samples = pt_samples, pt_burn_in = pt_burn_in, pt_thin = pt_thin,
      fid_samples = fid_samples, fid_burn_in = fid_burn_in)
  })
  names(group_results) <- names(groups)
  combine <- function(name) {
    out <- do.call(rbind, lapply(group_results, function(result) result[[name]]))
    rownames(out) <- NULL
    out$tooth <- tooth
    out
  }
  survival_intervals <- combine("survival_intervals")
  survival_intervals$group <- factor(survival_intervals$group, levels = names(groups))
  censoring_summary <- combine("censoring_summary")
  sampling_settings <- combine("sampling_settings")
  runtime_summary <- data.frame(
    tooth = tooth,
    method = group_results$boy$sampling_settings$method,
    boy_seconds = group_results$boy$sampling_settings$elapsed_seconds,
    girl_seconds = group_results$girl$sampling_settings$elapsed_seconds)
  # Arithmetic mean across the two groups, one complete run per group.
  runtime_summary$mean_seconds <- rowMeans(runtime_summary[, c("boy_seconds", "girl_seconds")])

  list(
    tooth = tooth, groups = group_results, censoring_summary = censoring_summary,
    survival_intervals = survival_intervals, sampling_settings = sampling_settings,
    runtime_summary = runtime_summary, tail_probabilities = combine("tail_probability"),
    comparison_plot = plot_tandmob(survival_intervals, censoring_summary),
    session_info = utils::sessionInfo())
}

analyze_tandmob <- function(data, teeth = c(14L, 24L, 34L, 44L), seed = 123L,
                            existing_results = list(), output_dir = NULL, ...) {
  tooth_results <- existing_results
  for (i in seq_along(teeth)) {
    tooth <- as.character(teeth[i])
    if (tooth %in% names(existing_results)) next
    columns <- paste0(c("EBEG.", "EEND."), tooth)
    boy <- data[which(data$GENDERNum == 0), columns]
    girl <- data[which(data$GENDERNum == 1), columns]
    names(boy) <- names(girl) <- c("L", "R")
    cat("\n========== Tooth ", tooth, " ==========\n", sep = "")
    tooth_results[[tooth]] <- compare_tandmob(
      boy, girl, tooth = teeth[i], seed = seed + 8L * (i - 1L), ...)
    # Save each completed tooth so it can be explicitly reused in a later call.
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(tooth_results[[tooth]], file.path(output_dir, paste0("tooth_", tooth, ".rds")))
    }
  }
  tooth_results <- tooth_results[as.character(teeth)]
  combine <- function(name) {
    out <- do.call(rbind, lapply(tooth_results, function(result) result[[name]]))
    rownames(out) <- NULL
    out
  }
  survival_intervals <- combine("survival_intervals")
  censoring_summary <- combine("censoring_summary")
  sampling_settings <- combine("sampling_settings")
  runtime_by_tooth <- combine("runtime_summary")
  method_order <- c("SurvPT", "Fiducial", "MPT", "TBP")
  runtime_summary <- do.call(rbind, lapply(method_order, function(method) {
    d <- sampling_settings[sampling_settings$method == method, ]
    data.frame(method = method, runs = nrow(d),
               boy_mean_seconds = mean(d$elapsed_seconds[d$group == "boy"]),
               girl_mean_seconds = mean(d$elapsed_seconds[d$group == "girl"]),
               mean_seconds = mean(d$elapsed_seconds))
  }))

  list(
    teeth = tooth_results, censoring_summary = censoring_summary,
    survival_intervals = survival_intervals, sampling_settings = sampling_settings,
    runtime_by_tooth = runtime_by_tooth, runtime_summary = runtime_summary,
    tail_probabilities = combine("tail_probabilities"),
    comparison_plot = plot_tandmob(survival_intervals, censoring_summary),
    session_info = utils::sessionInfo())
}

utils::data("tandmobAll", package = "icensBKL")
# GENDERNum: 0 = boy, 1 = girl; EBEG/EEND: ages in years.
# NA in EBEG denotes left censoring; NA in EEND denotes right censoring.
# https://cloud.r-project.org/web/packages/icensBKL/icensBKL.pdf
teeth <- c(14L, 24L, 34L, 44L)
output_dir <- getOption("SurvPT.tandmob_output_dir", NULL)
tandmob_results <- analyze_tandmob(
  tandmobAll, teeth = teeth, seed = 123L,
  existing_results = getOption("SurvPT.tandmob_existing", list()),
  output_dir = output_dir)
survival_intervals <- tandmob_results$survival_intervals
sampling_settings <- tandmob_results$sampling_settings
runtime_summary <- tandmob_results$runtime_summary
runtime_by_tooth <- tandmob_results$runtime_by_tooth
comparison_plot <- tandmob_results$comparison_plot
cat("\nCensoring summary:\n")
print(tandmob_results$censoring_summary, row.names = FALSE)
cat("\nRuntime in seconds: fitting + sampling + interval calculation; excludes plotting.\n")
cat("Per-tooth mean: one boy run and one girl run per method.\n")
print(runtime_by_tooth, row.names = FALSE, digits = 5)
cat("\nOverall mean: all tooth/group runs per method, with equal weight.\n")
print(runtime_summary, row.names = FALSE, digits = 5)

if (!is.null(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(survival_intervals, file.path(output_dir, "survival_intervals.csv"),
                   row.names = FALSE)
  utils::write.csv(sampling_settings, file.path(output_dir, "sampling_settings.csv"),
                   row.names = FALSE)
  utils::write.csv(runtime_summary, file.path(output_dir, "runtime_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(runtime_by_tooth, file.path(output_dir, "runtime_by_tooth.csv"),
                   row.names = FALSE)
  utils::write.csv(tandmob_results$censoring_summary,
                   file.path(output_dir, "censoring_summary.csv"), row.names = FALSE)
  ggplot2::ggsave(file.path(output_dir, "survival_comparison.png"), comparison_plot,
                  width = 12, height = 10, dpi = 180)
  saveRDS(tandmob_results, file.path(output_dir, "tandmob_results.rds"))
  message("Combined plot, runtime table and sampled results saved in ", normalizePath(output_dir))
}
if (interactive()) {
  print(comparison_plot)
} else {
  message("Combined plot available as comparison_plot (also in tandmob_results).")
}
