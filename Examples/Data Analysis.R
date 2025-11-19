devtools::load_all()
library(interval)
library(parallel)
library(ggplot2)
source('./Examples/Fiducial.R')
source('./Examples/LinInterpolation.R')

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
data <- data.frame(L = L, R = R)



pt_boot <- function(data, B = 1000){
  set.seed(123)
  n = nrow(data)
  obs.mat = c(data$L[!is.infinite(data$L)], data$R[!is.infinite(data$R)])
  x = seq(0, max(obs.mat), by = max(obs.mat)/200)
  exp_par = estimate_exp(data[,1], data[,2])

  fit = survPT(data, test_seq = x, prior = 1, prior_params = exp_par, EM_tol = 1e-6, verbose = F)
  surv_probs = fit$Survival_Probs


  n_cores <- parallel::detectCores()
  boot_surv_list <- mclapply(1:B, function(b) {
    indices = sample(1:n, n, replace = TRUE)
    boot_data = data[indices, ]
    exp_par = estimate_exp(data[,1], data[,2])
    boot_fit = survPT(boot_data, test_seq = x, prior = 1, prior_params = exp_par, EM_tol = 1e-6, verbose = F)
    return(c(boot_fit$Survival_Probs))
  }, mc.cores = 1)

  boot_surv = do.call(rbind, boot_surv_list)

  lower_ci = apply(boot_surv, 2, quantile, probs = 0.025)
  upper_ci = apply(boot_surv, 2, quantile, probs = 0.975)

  return(list(surv_probs = surv_probs, lower_ci = lower_ci, upper_ci = upper_ci))
}


n = nrow(data)
exp_par = estimate_exp(data[,1], data[,2])

time = c()

obs.mat = c(data$L[!is.infinite(data$L)], data$R[!is.infinite(data$R)])
x = seq(0, max(obs.mat), by = max(obs.mat)/200)
grid = seq(0, max(obs.mat), length.out = 100)

tryCatch({
  start_time <- Sys.time()
  fit = pt_boot(data, B = 100)
  pt_surv = fit$surv_probs
  pt_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  time = c(time, pt_time)

  pt_lower = fit$lower_ci
  pt_upper = fit$upper_ci
}, error = function(e) {
  cat("Error in PT Bootstrap:\n")
  print(e)
})



start_time <- Sys.time()
fid = fiducial(data$L, data$R, n, grid, x)
fiducial_surv = 1 - fid[[1]]
fid_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
time = c(time, fid_time)


fid_lower = 1 - fid[[2]]
fid_upper = 1 - fid[[3]]

plot = ggplot()+
  geom_step(aes(x = x, y = pt_surv, color = "PT Estimate"), size = 1)+
  #geom_ribbon(aes(x = x, ymin = pt_lower, ymax = pt_upper, fill = "PT Estimate"), alpha = 0.2)+
  geom_step(aes(x = x, y = fiducial_surv, color = "Fiducial Estimate"), size = 1)+
  geom_ribbon(aes(x = x, ymin = fid_lower, ymax = fid_upper, fill = "Fiducial Estimate"), alpha = 0.2)+
  labs(x = "Time",
       y = "Survival Probability")+
  scale_color_manual(name = "Method", values = c("PT Estimate" = "blue", "Fiducial Estimate" = "red")) +
  scale_fill_manual(name = "Method", values = c("PT Estimate" = "blue", "Fiducial Estimate" = "red")) +
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,max(x)))+
  theme_minimal()+
  theme(text = element_text(size = 16))+
  theme(legend.position = "right")


save(fit, fid, file = "heavy_treat_estimates.RData")
ggsave("Heavy_treat.png", plot, width = 3.4, height = 3.0, units = "in")

save(fit, fid, file = "light_treat_estimates.RData")
ggsave("Light_treat.png", plot, width = 3.4, height = 3.0, units = "in")





