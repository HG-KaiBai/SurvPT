# Historical large simulation: executes experiment loops when sourced.
# See Examples/README.md; reduce replication/worker counts before a pilot run.
# Requires SurvPT, survival, interval, and parallel; results remain in R.
library(SurvPT)
library(survival)
library(interval)
library(parallel)

simu <- function(t, n, K, Lambda, q, type){
  if (type == "right") {
    data = data_rc(t, n, K, Lambda, q)
  } else if (type == "interval"){
    data = data_ic(t, n, K, Lambda, q)
  } else if (type == "current"){
    data = data_cs(t, n, K, Lambda, q)
  } else {
    message("Type must be one of right, interval and current")
    return(1)
  }

  exp_par = estimate_exp(data[,1], data[,2])
  weibull_par = estimate(data[,1], data[,2])

  obs.mat = c(data$L[!is.infinite(data$L)], data$R[!is.infinite(data$R)])
  x = seq(0, max(obs.mat), by = max(obs.mat)/200)
  true.surv = 1 - pweibull(x, K, 1/Lambda)

  time = c()
  sse = c()

  for (prior in c(1, 2)){
    start_time <- Sys.time()
    fit = survPT(data, prior = prior, prior_params = exp_par, EM_tol = 1e-4, verbose = T)
    prediction = predict(fit, x, sample = FALSE)
    pt_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    time = c(time, pt_time)
    err = sum((prediction$estimate - true.surv)^2)
    sse = c(sse, err)
  }

  start_time <- Sys.time()
  kmf = icfit(data$L, data$R)
  turnbull_surv = getsurv(x, kmf, "interpolation")[[1]]$S
  turnbull_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  time = c(time, turnbull_time)

  sse = c(sse, NP = sum((turnbull_surv - true.surv)^2))

  weibull_surv = 1 - pweibull(x, weibull_par[2], 1/weibull_par[1])
  sse = c(sse, Weibull = sum((weibull_surv - true.surv)^2))

  exp_surv = 1 - pexp(x, exp_par)
  sse = c(sse, Exp = sum((exp_surv - true.surv)^2))

  return(c(sse, time))
}


output <- function(R, n, K, Lambda, q, type){
  num_cores <- if (.Platform$OS.type == "windows") 1L else 2L
  result = mclapply(1:R, function(t) simu(t, n, K, Lambda, q, type), mc.cores = num_cores)
  result = do.call(rbind, result)
  colnames(result) =c("Uniform_PT", "Exp_PT", "NP", "Weibull", "Exp", "Uniform_PT_time", "Exp_PT_time",  "NP_time")
  result = rbind(result, colMeans(result))
  rownames(result) = c(seq(1:R), "Ave")
  print(tail(result, 1))
  return(result)
}

# for (type in c("right", "current", "interval")){
#   for (q in c(0.3, 0.5)){
#     message(type, q)
#     output(1000, 200, 1, 0.12, q, type)
#     output(1000, 200, 2, 0.12, q, type)
#   }
# }

simulation_results <- list()
for (type in c("right", "current", "interval")){
  for (q in c(0.3, 0.5)){
    message(type, q)
    for (K in c(1, 2)) {
      key <- paste(type, "n500", paste0("k", K), paste0("q", q), sep="_")
      simulation_results[[key]] <- output(1000, 500, K, 0.12, q, type)
    }
  }
}
