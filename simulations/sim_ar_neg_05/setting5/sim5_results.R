# ------------------------------------------------------------------------------
# Sim setting 1: Results
# ------------------------------------------------------------------------------
folder_path <- "simulations/sim_ar_neg_05/setting5"
load(paste(folder_path, "/wbb1_sim5.Rdata", sep = ""))
load(paste(folder_path, "/wbb2_sim5.Rdata", sep = ""))
load(paste(folder_path, "/bob_sim5.Rdata", sep = ""))
load(paste(folder_path, "/nuts_sim5.Rdata", sep = ""))
load(paste(folder_path, "/advi_sim5.Rdata", sep = ""))

SW2.distance <- list(
  WBB1 = sapply(wbb1_setting5, \(m) m$sw2.dist),
  WBB2 = sapply(wbb2_setting5, \(m) m$sw2.dist),
  BOB = sapply(bob_setting5, \(m) m$sw2.dist),
  NUTS = sapply(nuts_setting5, \(m) m$sw2.dist),
  ADVI = sapply(advi_setting5, \(m) m$sw2.dist)
)

run.times <- list(
  WBB1 = sapply(wbb1_setting5, \(m) m$run.time),
  WBB2 = sapply(wbb2_setting5, \(m) m$run.time),
  BOB = sapply(bob_setting5, \(m) m$run.time),
  NUTS = sapply(nuts_setting5, \(m) m$run.time),
  ADVI = sapply(advi_setting5, \(m) m$run.time)
)

optimal.x.values <- sapply(bob_setting5, \(m) m$x.optim)

results_sim5 <- list(
  SW2.distance = SW2.distance, 
  run.times = run.times,
  optimal.x.values = optimal.x.values
)

# Save the results
save(results_sim5, file = paste(folder_path, "/sim5_results.Rdata", sep = ""))
