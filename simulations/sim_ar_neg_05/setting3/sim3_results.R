# ------------------------------------------------------------------------------
# Sim setting 3: Results
# ------------------------------------------------------------------------------
folder_path <- "simulations/sim_ar_neg_05/setting3"
load(paste(folder_path, "/wbb1_sim3.Rdata", sep = ""))
load(paste(folder_path, "/wbb2_sim3.Rdata", sep = ""))
load(paste(folder_path, "/bob_sim3.Rdata", sep = ""))
load(paste(folder_path, "/nuts_sim3.Rdata", sep = ""))
load(paste(folder_path, "/advi_sim3.Rdata", sep = ""))

SW2.distance <- list(
  WBB1 = sapply(wbb1_setting3, \(m) m$sw2.dist),
  WBB2 = sapply(wbb2_setting3, \(m) m$sw2.dist),
  BOB = sapply(bob_setting3, \(m) m$sw2.dist),
  NUTS = sapply(nuts_setting3, \(m) m$sw2.dist),
  ADVI = sapply(advi_setting3, \(m) m$sw2.dist)
)

run.times <- list(
  WBB1 = sapply(wbb1_setting3, \(m) m$run.time),
  WBB2 = sapply(wbb2_setting3, \(m) m$run.time),
  BOB = sapply(bob_setting3, \(m) m$run.time),
  NUTS = sapply(nuts_setting3, \(m) m$run.time),
  ADVI = sapply(advi_setting3, \(m) m$run.time)
)

optimal.x.values <- sapply(bob_setting3, \(m) m$x.optim)

results_sim3 <- list(
  SW2.distance = SW2.distance, 
  run.times = run.times,
  optimal.x.values = optimal.x.values
)

# Save the results
save(results_sim3, file = paste(folder_path, "/sim3_results.Rdata", sep = ""))
