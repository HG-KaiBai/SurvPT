# Historical interactive comparison; see Examples/README.md before running.
# Requires devtools and interval. View() calls need an interactive session;
# later fiducial plots also require constructing the currently commented-out fid.
devtools::load_all()
library(interval)
K = 1
Lambda = 0.12
n = 200
test_data = data_rc(seed = 354, n = n, K = K, Lambda = Lambda, q = 0.3)
exp_par = estimate_exp(test_data[,1], test_data[,2])
weib_par = estimate(test_data[,1], test_data[,2])

obs.mat = c(test_data$L[!is.infinite(test_data$L)], test_data$R[!is.infinite(test_data$R)])
x = seq(0, max(obs.mat), by = max(obs.mat)/200)
true.surv = 1 - pweibull(x, K, 1/Lambda)

# Try calling the function with basic parameters
tryCatch({
  cat("Calling SurvPT with uniform prior...\n")
  result_unif <- survPT(test_data, test_seq = x, prior = 1, prior_params = exp_par, max_M = 8, max_EM_iter = 100, EM_tol = 1e-4)
  cat("Success! Result structure:\n")
  print(str(result_unif))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

for (l in 0:9){
  print(sum(result_unif$Posterior_Probs[result_unif$Level == l]))
}
sum(result_unif$Posterior_Probs[result_unif$Terminal_Nodes])

message("SSE Error For Unif Prior")
sum((result_unif$Survival_Probs - true.surv)^2)
message("RMSE Error For Unif Prior")
sqrt(mean((abs(diff(c(result_unif$Survival_Probs, 0))) - abs(diff(c(true.surv,0))))^2))

tryCatch({
  cat("Calling SurvPT with uniform prior...\n")
  result_exp <- survPT(test_data, test_seq = x, prior = 2, prior_params = exp_par, max_M = 8, max_EM_iter = 100, EM_tol = 1e-4)

  cat("Success! Result structure:\n")
  print(str(result_exp))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

for (l in 0:9){
  print(sum(result_exp$Posterior_Probs[result_exp$Level == l]))
}
sum(result_exp$Posterior_Probs[result_exp$Terminal_Nodes])

message("SSE Error For Exp Prior")
sum((result_exp$Survival_Probs - true.surv)^2)
message("RMSE Error For Exp Prior")
sqrt(mean((abs(diff(c(result_exp$Survival_Probs, 0))) - abs(diff(c(true.surv,0))))^2))

# SSE error view
matrix = cbind(x, true.surv, result_exp$Survival_Probs, result_unif$Survival_Probs,
               (result_exp$Survival_Probs - true.surv), (result_unif$Survival_Probs - true.surv))
colnames(matrix) = c("Time","True_S", "Exp_S", "Unif_S", "Exp_Diff", "Unif_Diff")
View(matrix)

# RMSE error View
matrix = cbind(x[-1],-diff(true.surv), -diff(result_exp$Survival_Probs), -diff(result_unif$Survival_Probs),
               (-diff(result_exp$Survival_Probs) - -diff(true.surv)), (-diff(result_unif$Survival_Probs) - -diff(true.surv)))
colnames(matrix) = c("Time","True_S", "Exp_S", "Unif_S", "Exp_Diff", "Unif_Diff")
View(matrix)










tryCatch({
  cat("Calling SurvOPT with uniform prior...\n")
  result_opt <- survOPT(test_data$L, test_data$status, test_seq = x, max_M = 8)

  cat("Success! Result structure:\n")
  print(str(result_opt))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

for (l in 0:9){
  print(sum(result_opt$Posterior_Probs[result_opt$Level == l]))
}
print(sum(result_opt$Posterior_Probs[result_opt$Terminal_Nodes]))

grid = result_unif$Test_seq
#fid = fiducial(test_data$L, test_data$R, n = n, grid, grid)


kmf = icfit(test_data$L, test_data$R)
plot(kmf)
lines(result_unif$Test_seq, result_unif$Survival_Probs, col = "blue", lwd = 2)
lines(result_exp$Test_seq, result_exp$Survival_Probs, col = "orange", lwd = 2)
lines(result_opt$Test_seq, result_opt$Survival_Probs, col = "green", lwd = 2)
lines(result_unif$Test_seq, 1-pweibull(result_unif$Test_seq, K, 1/Lambda), col = "red", lwd = 2)
#lines(grid, 1 - fid$point_li, col = "purple", lwd = 2)
legend("topright", legend = c("Kaplan-Meier", "SurvPT Unif", "SurvPT Exp", "OPT", "True Survival", "Fiducial"),
       col = c("black", "blue","orange", "green" ,"red", "purple"), lwd = 2)


plot(result_unif$Test_seq, dweibull(result_unif$Test_seq, K, 1/Lambda), type = "l", col = "red", lwd = 2)
lines(result_unif$Test_seq, result_unif$Density, col = "blue", lwd = 2)
lines(result_exp$Test_seq, result_exp$Density, col = "orange", lwd = 2)
lines(result_opt$Test_seq, result_opt$Density, col = "green", lwd = 2)
lines(grid, fid$density, col = "purple", lwd = 2)







