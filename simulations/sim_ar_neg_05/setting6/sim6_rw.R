# ------------------------------------------------------------------------------
# Sim setting 6: Random weighting
# ------------------------------------------------------------------------------

library(BOBgmms)

source("R/generate_data.R")
source("R/init.R")
source("R/tune.R")
source("R/bayes_post.R")
source("R/predictive.R")
source("R/dens_plots.R")
#source("R/advi.R")
#source("R/nuts.R")

load("simulations/sim_ar_neg_05/setting6/data_setting6.Rdata")

# Simulation settings
total_sims <- length(data_setting6)
K <- length(data_setting6$rep1$params.prior$alphas)
data_dims <- dim(data_setting6$rep1$y)
n <- data_dims[1]
d <- data_dims[2]


# Initialize data structures to store the data
wbb1_setting6 <- vector("list", total_sims)
wbb2_setting6 <- vector("list", total_sims)
bob_setting6 <- vector("list", total_sims)

# Run the algorithms -----------------------------------------------------------
for (sim in seq_len(total_sims)) {
  
  cat("------------------------------------- \n")
  cat("Starting simulation:", sim, "\n")
  cat("------------------------------------- \n")
  
  # Extract data and hyper-parameters from current replication
  current_rep <- data_setting6[[sim]]
  y <- current_rep$y
  alphas <- current_rep$params.prior$alphas
  betas <- current_rep$params.prior$betas
  psis <- current_rep$params.prior$psis
  lambdas <- current_rep$params.prior$lambdas
  nus <- current_rep$params.prior$nus
  
  a <- current_rep$params.temp$a
  b <- current_rep$params.temp$b
  c <- current_rep$params.temp$c
  r <- current_rep$params.temp$r
  
  # Extract initial parameters
  means.init <- current_rep$params.init$means.init
  covs.init <- current_rep$params.init$covs.init
  probs.init <- current_rep$params.init$probs.init
  
  # Weighted Bayesian Bootstrap 1 ----------------------------------------------
  
  set.seed(1)
  start.wbb1 <- Sys.time()
  wbb1.out <- BOBgmms::wbb.gmm(
    y = y, means.init = means.init, covs.init = covs.init,
    probs.init = probs.init, betas = betas, lambdas = lambdas,
    nus = nus, psis = psis, alphas = alphas, wbb.scheme = "wbb1",
    max.iters = 20000, a = a, b = b, c = c, r = r
  )
  end.wbb1 <- Sys.time()
  (wbb1.time <- difftime(end.wbb1, start.wbb1, units = "secs"))
  
  # Weighted Bayesian Bootstrap 2 ----------------------------------------------
  
  set.seed(1)
  start.wbb2 <- Sys.time()
  wbb2.out <- BOBgmms::wbb.gmm(
    y = y, means.init = means.init, covs.init = covs.init,
    probs.init = probs.init, betas = betas, lambdas = lambdas,
    nus = nus, psis = psis, alphas = alphas, wbb.scheme = "wbb2",
    max.iters = 20000, a = a, b = b, c = c, r = r
  )
  end.wbb2 <- Sys.time()
  (wbb2.time <- difftime(end.wbb2, start.wbb2, units = "secs"))
  
  # Bayesian Optimized Bootstrap -----------------------------------------------
  
  lower_bound <- 1e-05
  upper_bound <- 1.5
  
  set.seed(1)
  start.bob <- Sys.time()
  bob.out <- BOBgmms::bob.gmm(
    y = y, means.init = means.init, covs.init = covs.init,
    probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
    psis = psis, alphas = alphas, lower_bound = lower_bound,
    upper_bound = upper_bound, max.iters = 20000, size.batch = 1000,
    bo.iters = 100, a = a, b = b, c = c, r = r
  )
  end.bob <- Sys.time()
  (bob.time <- difftime(end.bob, start.bob, units = "secs"))
  bob_x.optim <- bob.out$x.optim # optimal x value - 1e-05  
  
  # Posterior predictive distributions -----------------------------------------
  wbb1.pred <- post.predictive(wbb1.out, d = d, K = K)
  wbb2.pred <- post.predictive(wbb2.out, d = d, K = K)
  bob.pred <- post.predictive(bob.out$post.draws, d = d, K = K)
  
  # Distances  -----------------------------------------------------------------
  set.seed(1)
  bayes.pred <- current_rep$bayes.pred
  dist.wbb1 <- tryCatch(
    {
      T4transport::swdist(bayes.pred, wbb1.pred)$distance
    }, error = function(e) {
      NA
    }
  )
  dist.wbb2 <- tryCatch(
    {
      T4transport::swdist(bayes.pred, wbb2.pred)$distance
    }, error = function(e) {
      NA
    }
  )
  dist.bob <- tryCatch(
    {
      T4transport::swdist(bayes.pred, bob.pred)$distance
    }, error = function(e) {
      NA
    }
  )
  
  # Save the results
  results.wbb1 <- list(
    sw2.dist = dist.wbb1,
    run.time = wbb1.time
  )
  
  results.wbb2 <- list(
    sw2.dist = dist.wbb2,
    run.time = wbb2.time
  )
  
  results.bob <- list(
    sw2.dist = dist.bob,
    run.time = bob.time,
    x.optim = bob_x.optim
  )
  
  wbb1_setting6[[sim]] <- results.wbb1
  wbb2_setting6[[sim]] <- results.wbb2
  bob_setting6[[sim]] <- results.bob
  
  rm(
    wbb1.out, wbb2.out, bob.out,
    wbb1.pred, wbb2.pred, bob.pred
  )
  gc()
}

names(wbb1_setting6) <- paste0("rep", seq_len(total_sims))
names(wbb2_setting6) <- paste0("rep", seq_len(total_sims))
names(bob_setting6) <- paste0("rep", seq_len(total_sims))

# Export the RW results
folder_path <- "simulations/sim_ar_neg_05/setting6"
save(wbb1_setting6, file = paste(folder_path, "/wbb1_sim6.Rdata", sep = ""))
save(wbb2_setting6, file = paste(folder_path, "/wbb2_sim6.Rdata", sep = ""))
save(bob_setting6, file = paste(folder_path, "/bob_sim6.Rdata", sep = ""))
