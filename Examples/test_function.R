devtools::load_all()
library(interval)
K = 1
Lambda = 0.12
test_data = data_ic(seed = 123, n = 1000, K = K, Lambda = Lambda, q = 0.3)
exp_par = estimate_exp(test_data[,1], test_data[,2])
weib_par = estimate(test_data[,1], test_data[,2])
# Try calling the function with basic parameters
tryCatch({
  cat("Calling SurvPT with uniform prior...\n")
  result <- survPT(test_data, max_EM_iter = 100, EM_tol = 1e-4)

  cat("Success! Result structure:\n")
  print(str(result))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

for (l in 0:8){
  print(sum(result$Posterior_Probs[result$Level == l]))
}

kmf = icfit(test_data$L, test_data$R)
plot(kmf)
lines(result$Test_seq, result$Survival_Probs, col = "blue", lwd = 2)
lines(result$Test_seq, 1-pweibull(result$Test_seq, K, 1/Lambda), col = "red", lwd = 2)
legend("topright", legend = c("Kaplan-Meier", "SurvPT Estimate", "True Survival"), col = c("black", "blue", "red"), lwd = 2)

