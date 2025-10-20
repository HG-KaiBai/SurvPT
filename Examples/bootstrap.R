fit_turnbull <- function(times, mod) {
  last_partition <- mod$intmap[, ncol(mod$intmap)]
  if (is.finite(last_partition[2])) {
    # single-stratum: getsurv() returns a list with element S
    return(getsurv(times, mod, "interpolation")[[1]]$S)
  } else {
    res <- numeric(length(times))

    idx <- which(times < last_partition[1])
    if (length(idx) > 0) {
      res[idx] <- getsurv(times[idx], mod, "interpolation")[[1]]$S
    }

    idx <- which(times >= last_partition[1])
    if (length(idx) > 0) {
      len <- max(0, 2 * last_partition[1] - times[idx])
      res[idx] <- mod$pf[length(mod$pf)] * len / (last_partition[1])
    }
    return(res)
  }
}




# --- bootstrap wrapper given raw interval-censored data ---
# df must contain: L, R ; optional logical Lin, Rin (see ?icfit)
icfit_boot_quantile <- function(df, times, B = 1000) {
  stat_fn <- function(dat, indices) {
    d <- dat[indices, , drop = FALSE]
    fit <- icfit(d$L, d$R)
    ep = estimate_exp(d$L, d$R)
    ans = c()
    for (prior in c(1, 2)){
      tree = survPT(d, test_seq = times, prior = prior, prior_params = ep, EM_tol = 1e-4, verbose = FALSE)
      ans = c(ans, tree$Survival_Probs)
    }
    np = fit_turnbull(times, fit)
    ans = c(ans, np)
    ans
  }
  bt <- boot(data = df, statistic = stat_fn, R = B)
  fit0 <- icfit(df$L, df$R)
  ep0 = estimate_exp(df$L, df$R)
  t_hat = c()
  for (prior in c(1, 2)){
    tree0 = survPT(df, test_seq = times, prior = prior, prior_params = ep0, EM_tol = 1e-4, verbose = FALSE)
    t_hat = c(t_hat, tree0$Survival_Probs)
  }
  np_hat = fit_turnbull(times, fit0)
  t_hat <- c(t_hat, np_hat)

  # Percentile CI for each quantile
  ci_percentile <- t(sapply(seq_along(1:9), function(i) {
    quant <- boot.ci(bt, type = "perc", index = i)
    if (inherits(quant, "try-error") || is.null(quant$percent)) c(NA_real_, NA_real_) else quant$percent[4:5]
  }))
  colnames(ci_percentile) <- c("lower", "upper")
  rownames(ci_percentile) <- paste0("S", rep(c(0.75,0.5,0.25), 3))
  data.frame(estimate = t_hat, lower = ci_percentile[,1], upper = ci_percentile[,2])
}

