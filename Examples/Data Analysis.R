# Heavily treatment
L = c(15, 15, 10, 12, 9, 15, 5, 1, 15, 1, 9, 11, 12, 7, 11, 10, 12, 10, 1, 5, 10, 8,
           16, 16, 16, 1, 14, 14, 15, 15, 14, 9, 14, 14, 10, 13, 1, 11, 12, 9, 14, 3, 9, 13, 9,
           17, 17, 17, 12, 13, 1, 1, 14, 10, 10, 3, 15, 1, 10, 9, 11, 7, 13, 7, 10, 1, 1, 10,
           13, 14, 12, 10, 12, 10, 1, 12, 1, 9, 13, 7, 1, 9, 12, 13, 9, 10, 10,
           14, 13, 11, 5, 12, 8, 1, 14, 1, 1, 1, 14, 8, 3, 13, 10, 5, 10)

R = c(Inf, Inf, 11, 14, 12, 16, 7, 14, 16, 13, 10, 13, 13, 10, 13, 11, 13, 11, 7, 7, 11, 10,
      Inf, Inf, Inf, 16, 15, 15, 16, 15, 15, 10, 15, 15, 11, 14, 15, 13, 13, 13, 15, 7, 11, 15, 12,
      Inf, Inf, Inf, 13, 16, 11, 13, 15, 11, 12, 15, 16, 7, 15, 12, 13, 9, 14, 9, 11, 10, 7, 12,
      15, 15, 14, 11, 13, 11, 14, 13, 15, 12, 15, 15, 6, 12, 13, 15, 13, 12, 14,
      16, 15, 12, 7, 13, 10, 15, 15, 14, 11, 12, 16, 15, 14, 14, 12, 8, 12)
data <- data.frame(L = L, R = R)

# Lightly treatment
L = c(1, 10, 10, 8, 1, 1, 1, 12, 1, 6, 12, 1, 11, 11, 1, 1, 1, 15, 11, 1, 7,
      rep(15, 19), 12, 1, 1, 9, 1, 13, 9, 14, 13, 1, 13, 13, 13, 12, 12, 1, 13, 1, 12,
      rep(16, 31), 1, 15, 13, 8, 11, 9, 13, 1, 15, 12, 8, 3, 12, 1, 9, 12, 13, 11, 6,
      rep(17, 10), 1, 3, 5, 11, 8, 1, 4, 14, 7, 13, 10, 6, 7, 10, 12, 14, 10, 1, 11,
      18, 1, 8, 14, 1, 9, 1, 1, 1, 12, 15, 14, 15, 1, 11, 11, 1, 6, 10)

R = c(Inf, 15, 12, 13, 11, 18, 15, 15, 15, 13, 14, 15, 15, 14, 15, 14, 11, 16, 12, 12, 13,
      rep(Inf, 19), 14, 16, 12, 14, 15, 14, 12, 15, 14, 13, 15, 15, 14, 13, 13, 16, 13, 3, 13,
      rep(Inf, 31), 15, 16, 14, 16, 16, 14, 15, 12, 16, 13, 16, 16, 14, 15, 14, 13, 14, 14, 13,
      rep(Inf, 10), 15, 10, 11, 12, 12, 5, 11, 15, 12, 15, 12, 8, 10, 16, 14, 15, 12, 5, 14,
      Inf, 15, 15, 16, 17, 13, 16, 16, 14, 14, 16, 15, 16, 12, 14 ,12, 15, 12, 11)
data_light <- data.frame(L = L, R = R)

devtools::load_all()
library(interval)
source('./Examples/Fiducial.R')
source('./Examples/LinInterpolation.R')

n = nrow(data)
exp_par = estimate_exp(data[,1], data[,2])
weibull_par = estimate(data[,1], data[,2])

time = c()

obs.mat = c(data$L[!is.infinite(data$L)], data$R[!is.infinite(data$R)])
x = seq(0, max(obs.mat), by = max(obs.mat)/200)
grid = seq(0, max(obs.mat), length.out = 100)


for (prior in c(1, 2)){
  start_time <- Sys.time()
  fit = survPT(data, test_seq = x, prior = prior, prior_params = exp_par, EM_tol = 1e-6, verbose = T)
  pt_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  time = c(time, pt_time)
}

start_time <- Sys.time()
fiducial_surv = 1 - fiducial(data$L, data$R, n, grid, x)
fid_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
time = c(time, fid_time)

kmf = icfit(data$L, data$R)
turnbull_surv = getsurv(x, kmf, "interpolation")[[1]]$S


plot(kmf)
lines(x, fit$Survival_Probs, col = "green", lwd = 2)
lines(x, fiducial_surv, col = "red", lwd = 2)






