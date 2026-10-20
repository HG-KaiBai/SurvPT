# SurvPT 0.2.0: small, reproducible examples using only installed dependencies.
# Run from a checkout: Rscript Examples/quickstart.R
# No files are written; the fitted objects remain available when sourced.
library(SurvPT)

d <- data_rc(seed = 123, n = 40, K = 1.5, Lambda = 0.4, q = 0.3)
x <- d[, c("L", "R")]
grid <- seq(0, 12, length.out = 121)

fits <- list(
  uniform = survPT(x, test_seq = grid, prior = 1,
                   max_M = 5, max_EM_iter = 50, verbose = FALSE),
  exponential = survPT(x, test_seq = grid, prior = 2, prior_params = 0.4,
                       max_M = 5, max_EM_iter = 50, verbose = FALSE),
  weibull = survPT(x, test_seq = grid, prior = 3, prior_params = c(0.4, 1.5),
                   max_M = 5, max_EM_iter = 50, verbose = FALSE)
)

checks <- do.call(rbind, lapply(names(fits), function(name) {
  fit <- fits[[name]]
  s <- as.numeric(fit$Survival_Probs)
  density <- as.numeric(fit$Density)
  mass <- sum(fit$Posterior_Probs[fit$Terminal_Nodes])
  stopifnot(abs(mass-1) < 1e-8, all(is.finite(s)),
            all(diff(s) <= 1e-10), all(s >= 0), all(s <= 1+1e-10),
            all(is.finite(density)), all(density >= 0))
  data.frame(method = name, terminal_nodes = length(fit$Terminal_Nodes),
             total_mass = mass, survival_at_last_grid_point = tail(s, 1))
}))

print(checks, row.names = FALSE)
cat("SurvPT version:", as.character(packageVersion("SurvPT")), "\n")
cat("Uniform PT has a finite cutoff; exponential and Weibull bases retain unbounded tails.\n")

# Optional plotting in an interactive R session:
if (interactive()) {
  plot(grid, as.numeric(fits$exponential$Survival_Probs), type = "l",
       xlab = "Time", ylab = "Survival probability", ylim = c(0, 1))
  lines(grid, pweibull(grid, shape = 1.5, scale = 1/0.4, lower.tail = FALSE),
        col = "red", lty = 2)
  legend("topright", c("PT exponential base", "Simulation truth"),
         col = c("black", "red"), lty = c(1, 2), bty = "n")
}
