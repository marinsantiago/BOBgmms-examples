# ------------------------------------------------------------------------------
# Sim setting 4: Results
# ------------------------------------------------------------------------------
folder_path <- "simulations/sim_ar_neg_05/setting4"
load(paste(folder_path, "/wbb1_sim4.Rdata", sep = ""))
load(paste(folder_path, "/wbb2_sim4.Rdata", sep = ""))
load(paste(folder_path, "/bob_sim4.Rdata", sep = ""))
load(paste(folder_path, "/nuts_sim4.Rdata", sep = ""))
load(paste(folder_path, "/advi_sim4.Rdata", sep = ""))

SW2.distance <- list(
  WBB1 = sapply(wbb1_setting4, \(m) m$sw2.dist),
  WBB2 = sapply(wbb2_setting4, \(m) m$sw2.dist),
  BOB = sapply(bob_setting4, \(m) m$sw2.dist),
  NUTS = sapply(nuts_setting4, \(m) m$sw2.dist),
  ADVI = sapply(advi_setting4, \(m) m$sw2.dist)
)

run.times <- list(
  WBB1 = sapply(wbb1_setting4, \(m) m$run.time),
  WBB2 = sapply(wbb2_setting4, \(m) m$run.time),
  BOB = sapply(bob_setting4, \(m) m$run.time),
  NUTS = sapply(nuts_setting4, \(m) m$run.time),
  ADVI = sapply(advi_setting4, \(m) m$run.time)
)

optimal.x.values <- sapply(bob_setting4, \(m) m$x.optim)

results_sim4 <- list(
  SW2.distance = SW2.distance, 
  run.times = run.times,
  optimal.x.values = optimal.x.values
)

# Save the results
save(results_sim4, file = paste(folder_path, "/sim4_results.Rdata", sep = ""))
