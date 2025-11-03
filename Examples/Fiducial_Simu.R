library(SurvPT)
library(interval)
library(parallel)
source('./Examples/Fiducial.R')
source('./Examples/LinInterpolation.R')

simu <- function(t, n, type){
  if (type == 1) {
    data = fiducial_1(t, n)
    grid = seq(0, 5, length.out = 100)
  } else if (type == 2){
    data = fiducial_2(t, n)
    grid = seq(0, 8, length.out = 100)
  } else if (type == 3){
    data = fiducial_3(t, n)
    grid = seq(0, 5, length.out = 100)
  } else if (type == 4){
    data = fiducial_4(t, n)
    grid = seq(0, 3, length.out = 100)
  } else {
    message("Type must be one of 1, 2, 3 and 4")
    return(1)
  }

  exp_par = estimate_exp(data[,1], data[,2])
  weibull_par = estimate(data[,1], data[,2])

  time = c()
  sse = c()

  obs.mat = c(data$L[!is.infinite(data$L)], data$R[!is.infinite(data$R)])
  x = seq(0, max(obs.mat), by = max(obs.mat)/200)
  grid = seq(0, max(obs.mat), length.out = 100)
  if (type == 1){
    true.surv = 1 - pexp(x, 1)
  } else if (type == 2){
    true.surv = 1 - pgamma(x, 3, 1)
  } else if (type ==3){
    true.surv = 1 - pgamma(x, 2, 1)
  } else if (type ==4){
    true.surv = 1 - pexp(x, 1)
  }

  for (prior in c(1, 2)){
    start_time <- Sys.time()
    fit = survPT(data, test_seq = x, prior = prior, prior_params = exp_par, EM_tol = 1e-6, verbose = T)
    pt_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    time = c(time, pt_time)
    err = sum((fit$Survival_Probs - true.surv)^2)
    sse = c(sse, err)
  }

  start_time <- Sys.time()
  fiducial_surv = 1 - fiducial(data$L, data$R, n, grid, x)
  fid_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  time = c(time, fid_time)
  err = sum((fiducial_surv - true.surv)^2)
  sse = c(sse, err)

  kmf = icfit(data$L, data$R)
  turnbull_surv = getsurv(x, kmf, "interpolation")[[1]]$S
  sse = c(sse, NP = sum((turnbull_surv - true.surv)^2))

  weibull_surv = 1 - pweibull(x, weibull_par[2], 1/weibull_par[1])
  sse = c(sse, Weibull = sum((weibull_surv - true.surv)^2))

  exp_surv = 1 - pexp(x, exp_par)
  sse = c(sse, Exp = sum((exp_surv - true.surv)^2))
  return(c(sse, time))
}


output <- function(R, n, type){
  num_cores <- 300
  result = mclapply(1:R, function(t) simu(t, n, type), mc.cores = num_cores)
  result = do.call(rbind, result)
  # colnames(result) =c("Uniform_PT", "Exp_PT", "Fiducial", "NP", "Weibull", "Exp", "Uniform_PT_time", "Exp_PT_time",  "Fiducial_time")
  result = rbind(result, colMeans(result))
  # rownames(result) = c(seq(1:R), "Ave")
  filename = paste0("Scenario_", type, "_n", n, ".csv")
  dirname = paste0("Scenario_", type)
  dir.create(file.path("fiducial", dirname), showWarnings = FALSE)
  write.csv(result, file = file.path("fiducial", dirname, filename), row.names = TRUE)

  sd = apply(result, 2, sd)
  out = c(result[nrow(result), 1:6], sd[1:6])
  out = out[c(1,7,2,8,3,9,4,10,5,11,6,12)]
  return(out)
}


 final = matrix(nrow=0, ncol= 6 * 2 * 4)
 for (n in c(50, 75, 100, 200, 500, 800, 1000)){
    type_ans = c()
    for (type in c(1,2,3,4)){
        message("S_", type, "_", n)
        type_ans = c(type_ans, output(1000, n, type))
    }
    final = rbind(final, type_ans)
}
rownames(final) = c("50", "75", "100", "200", "500", "800", "1000")
colnames(final) = c("S1_Unif_PT_SSE", "S1_Unif_PT_SD", "S1_Exp_PT_SSE", "S1_Exp_PT_SD", "S1_Fiducial_SSE", "S1_Fiducial_SD",
                        "S1_NP_SSE", "S1_NP_SD", "S1_Weibull_SSE", "S1_Weibull_SD", "S1_Exp_SSE", "S1_Exp_SD",
                    "S2_Unif_PT_SSE", "S2_Unif_PT_SD", "S2_Exp_PT_SSE", "S2_Exp_PT_SD", "S2_Fiducial_SSE", "S2_Fiducial_SD",
                        "S2_NP_SSE", "S2_NP_SD", "S2_Weibull_SSE", "S2_Weibull_SD", "S2_Exp_SSE", "S2_Exp_SD",
                    "S3_Unif_PT_SSE", "S3_Unif_PT_SD", "S3_Exp_PT_SSE", "S3_Exp_PT_SD", "S3_Fiducial_SSE", "S3_Fiducial_SD",
                        "S3_NP_SSE", "S3_NP_SD", "S3_Weibull_SSE", "S3_Weibull_SD", "S3_Exp_SSE", "S3_Exp_SD",
                    "S4_Unif_PT_SSE", "S4_Unif_PT_SD", "S4_Exp_PT_SSE", "S4_Exp_PT_SD", "S4_Fiducial_SSE", "S4_Fiducial_SD" ,
                        "S4_NP_SSE", "S4_NP_SD", "S4_Weibull_SSE", "S4_Weibull_SD", "S4_Exp_SSE", "S4_Exp_SD")

write.csv(final, file = file.path("fiducial","final_results.csv"), row.names = TRUE)
