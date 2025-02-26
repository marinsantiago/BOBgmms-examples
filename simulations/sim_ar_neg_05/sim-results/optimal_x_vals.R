# ------------------------------------------------------------------------------
# Sim setting 1-9: Optimal x-values derived from BOB
# ------------------------------------------------------------------------------

folder <- "simulations/sim_ar_neg_05"
for (s in seq_len(9)) {
  file <- paste(
    folder,
    paste("/setting", s, sep = ""),
    paste("/sim", s, "_results.Rdata", sep = ""),
    sep = ""
  )
  load(file)
}

results <- list(
  results_sim1, results_sim2, results_sim3,
  results_sim4, results_sim5, results_sim6,
  results_sim7, results_sim8, results_sim9
)

optimal_x_median <- sapply(results, \(r) median(r$optimal.x.values))
optimal_x_iqr <- sapply(results, \(r) IQR(r$optimal.x.values))

optimal_x.df <- data.frame(
  optimal_x_median,
  optimal_x_iqr
) |> t()

