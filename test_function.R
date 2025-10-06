devtools::load_all()
test_data = data_rc(seed = 123, n = 200, K = 1, Lambda = 0.12, q = 0.3)
exp_par = estimate_exp(test_data[,1], test_data[,2])
weib_par = estimate(test_data[,1], test_data[,2])
# Try calling the function with basic parameters
tryCatch({
  cat("Calling SurvPT with uniform prior...\n")
  result <- survPT(test_data)

  cat("Success! Result structure:\n")
  print(str(result))

}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
})

for (l in 0:8){
  print(sum(result$Posterior_Probs[result$Level == l]))
}

