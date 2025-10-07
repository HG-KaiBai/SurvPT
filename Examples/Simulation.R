library(SurvPT)
library(interval)
library(parallel)

simu <- function(t, n, K, Lambda, q, type){
    message(t)
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

    kmf = icfit(data$L, data$R)
    obs.mat = c(data$L[!is.infinite(data$L)], data$R[!is.infinite(data$R)])
    x = seq(0, max(obs.mat), by = max(obs.mat)/200)
    true.surv = 1 - pweibull(x, K, 1/Lambda)

    sse = c()
    for (prior in c(1, 2)){
        fit = survPT(data, test_seq = x, prior = prior, prior_params = exp_par)
        err = sum((fit$Survival_Probs - true.surv)^2)
        sse = c(sse, err)
    }

    turnbull_surv = getsurv(x, kmf, "interpolation")[[1]]$S
    sse = c(sse, NP = sum((turnbull_surv - true.surv)^2))

    weibull_surv = 1 - pweibull(x, weibull_par[2], 1/weibull_par[1])
    sse = c(sse, Weibull = sum((weibull_surv - true.surv)^2))

    exp_surv = 1 - pexp(x, exp_par)
    sse = c(sse, Exp = sum((exp_surv - true.surv)^2))

    return(sse)
}


output <- function(R, n, K, Lambda, q, type){
  num_cores <- 1
  result = mclapply(1:R, function(t) simu(t, n, K, Lambda, q, type), mc.cores = num_cores)
  result = do.call(rbind, result)
  colnames(result) =c("Uniform_PT", "Exp_PT", "NP", "Weibull", "Exp")
  result = rbind(result, colMeans(result))
  rownames(result) = c(seq(1:R), "Ave")
  # filename = paste0(type, n, "_k", K,"_lambda", Lambda, "_rate", q,".csv")
  # write.csv(result, file = filename, row.names = TRUE)
  return(result)
}


for (type in c("right", "current", "interval")){
  for (q in c(0.3, 0.5)){
    message(type, q)
    output(1000, 200, 1, 0.12, q, type)
    output(1000, 200, 2, 0.12, q, type)
  }
}
