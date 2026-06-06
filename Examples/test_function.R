devtools::load_all()
library(interval)
K = 1
Lambda = 0.12

test_data = data_rc(seed = 123, n = 200, K = K, Lambda = Lambda, q = 0.3)
exp_par = estimate_exp(test_data[,1], test_data[,2])
weib_par = estimate(test_data[,1], test_data[,2])
# Try calling the function with basic parameters
tryCatch({
  cat("Calling SurvPT with uniform prior...\n")
  result <- survPT(test_data, prior = 1, prior_params = exp_par, max_M = 8, max_EM_iter = 100, EM_tol = 1e-4)

  cat("Success! Result structure:\n")
  print(str(result))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

for (l in 0:9){
  print(sum(result$Posterior_Probs[result$Level == l]))
}
sum(result$Posterior_Probs[result$Terminal_Nodes])

tryCatch({
  cat("Calling SurvOPT with uniform prior...\n")
  result1 <- survOPT(test_data$L, test_data$status, max_M = 8)

  cat("Success! Result structure:\n")
  print(str(result1))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

for (l in 0:9){
  print(sum(result1$Posterior_Probs[result1$Level == l]))
}
print(sum(result1$Posterior_Probs[result1$Terminal_Nodes]))



kmf = icfit(test_data$L, test_data$R)
plot(kmf)
lines(result$Test_seq, result$Survival_Probs, col = "blue", lwd = 2)
lines(result1$Test_seq, result1$Survival_Probs, col = "green", lwd = 2)
lines(result$Test_seq, 1-pweibull(result$Test_seq, K, 1/Lambda), col = "red", lwd = 2)
legend("topright", legend = c("Kaplan-Meier", "SurvPT Estimate", "OPT Estimate", "True Survival"), col = c("black", "blue", "green" ,"red"), lwd = 2)


plot(result$Test_seq, dweibull(result$Test_seq, K, 1/Lambda), type = "l", col = "red", lwd = 2)
lines(result$Test_seq, result$Density, col = "blue", lwd = 2)
lines(result1$Test_seq, result1$Density, col = "green", lwd = 2)







