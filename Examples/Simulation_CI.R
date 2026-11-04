# Historical nested bootstrap study: executes experiment loops when sourced.
# See Examples/README.md for dependencies. Results remain in R.
# These are bootstrap survival intervals, not native PT posterior intervals.
devtools::load_all()
library(interval)
library(parallel)
source("./Examples/bootstrap.R")
library(boot)

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

  true_quant = qweibull(c(0.25, 0.5, 0.75), K, 1/Lambda)
  ans = c()
  boot_quant = icfit_boot_quantile(data, true_quant, B = 1000)
  ans_names = c("Unif-PT", "Exp-PT", "NPMLE")
  for (i in 1:(nrow(boot_quant)/3)){
    ans = cbind(ans, c(boot_quant[(3*(i-1)+1):(3*i),1] - c(0.75,0.5,0.25),
                       boot_quant[(3*(i-1)+1):(3*i),3] - boot_quant[(3*(i-1)+1):(3*i),2],
                       as.numeric(boot_quant[(3*(i-1)+1):(3*i),2] <= c(0.75,0.5,0.25) & boot_quant[(3*(i-1)+1):(3*i),3] >= c(0.75,0.5,0.25))))
  }
  colnames(ans) = ans_names

  return(ans)
}


output <- function(R, n, K, Lambda, q, type){
  num_cores <- if (.Platform$OS.type == "windows") 1L else 2L
  result = mclapply(1:R, function(t) simu(t, n, K, Lambda, q, type), mc.cores = num_cores)
  result = do.call(rbind, result)
  colnames(result) =c("Unif-PT", "Exp-PT", "NPMLE")
  print(head(result))
  return(result)
}


coverage_results <- list()
for (type in c("right", "current", "interval")){
  for (n in c(200, 500)){
    for (q in c(0.3, 0.5)){
      message(type, n, "_", q)
      for (K in c(1, 2)) {
        key <- paste(type, paste0("n", n), paste0("k", K), paste0("q", q), sep="_")
        coverage_results[[key]] <- output(1000, n, K, 0.12, q, type)
      }
    }
  }
}
