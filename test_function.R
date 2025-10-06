test_data = data_cs(seed = 123, n = 200, K = 1, Lambda = 0.12, q = 0.3)

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

