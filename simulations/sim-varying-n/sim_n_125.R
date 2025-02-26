# ------------------------------------------------------------------------------
# Simulations for a varying sample size (n = 50, 125, 250, 375, 500)
# ------------------------------------------------------------------------------

library(BOBgmms)

source("R/generate_data.R")
source("R/init.R")
source("R/tune.R")
source("R/bayes_post.R")
source("R/predictive.R")
source("R/dens_plots.R")

# n = 125 ----------------------------------------------------------------------

total_sims <- 30

# Initialize data structures to store the output
dists_wbb1_n125 <- dists_wbb2_n125 <- vector("list", total_sims)
dists_bob_n125 <- optimals_x_bob_n125 <- vector("list", total_sims)

# Simulation settings
n <- 125
d <- 15
K <- 2

for (sim in seq_len(total_sims)) {
  
  cat("------------------------------------------- \n")
  cat("Replication number: ", sim, "\n")
  cat("------------------------------------------- \n")
  
  # Generate the data
  set.seed(sim)
  data_sim <- generate.data(
    sample.size = n, data.dim = d, 
    n.clust = K, rho = -0.5
  )
  y <- data_sim$y
  Z <- data_sim$Z
  
  # Set prior and tuning hyper-parameters
  alphas <- rep(1.1, K)
  betas <- lapply(1:K, function(x) rep(0, d))
  psis <- lapply(1:K, function(x) diag(d))
  
  set.seed(1)
  shrink.params <- shrinkage.cv.tuning(y, betas, psis, alphas)
  
  (lambdas <- rep(shrink.params$best.lambda, K))
  (nus <- rep(shrink.params$best.nu, K))
  
  init.params <- initial.values(y, K)
  means.init <- init.params$values$means
  covs.init <- init.params$values$covs
  probs.init <- init.params$values$probs
  
  range.a <- seq(0, 0)
  range.b <- seq(0, 4, by = 0.6)
  range.c <- seq(1, 4, by = 0.6)
  range.r <- seq(1, 4, by = 0.6)
  
  temper <- temperature.tuning(
    y = y, means.init = means.init, covs.init = covs.init,
    probs.init = probs.init, range.a = range.a, range.b = range.b,
    range.c = range.c, range.r = range.r, betas = betas, lambdas = lambdas,
    nus = nus, psis = psis, alphas = alphas
  )
  
  (a <- temper$best.a)
  (b <- temper$best.b)
  (c <- temper$best.c)
  (r <- temper$best.r)
  
  # Weighted Bayesian Bootstrap 1
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
  
  # Weighted Bayesian Bootstrap 2
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
  
  
  # Bayesian Optimized Bootstrap
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
  
  # BAYES
  bayes.true <- bayesian_posterior(y, Z, betas, lambdas, nus, psis, alphas)
  mean(bayes.true)
  
  # Posterior predictive distributions
  wbb1.pred <- post.predictive(wbb1.out, d = d, K = K)
  wbb2.pred <- post.predictive(wbb2.out, d = d, K = K)
  bob.pred <- post.predictive(bob.out$post.draws, d = d, K = K)
  bayes.pred <- post.predictive(bayes.true, d = d, K = K)
  
  # Distances
  sw2wbb1 <- T4transport::swdist(bayes.pred, wbb1.pred)$distance
  sw2wbb2 <- T4transport::swdist(bayes.pred, wbb2.pred)$distance
  sw2bob <- T4transport::swdist(bayes.pred, bob.pred)$distance
  
  # Store the results
  dists_wbb1_n125[[sim]] <- sw2wbb1
  dists_wbb2_n125[[sim]] <- sw2wbb2
  dists_bob_n125[[sim]] <- sw2bob
  
  optimals_x_bob_n125[[sim]] <- bob.out$x.optim
}

results_n125 <- list(
  dist_wbb1 = dists_wbb1_n125,
  dist_wbb2 = dists_wbb2_n125,
  dist_bob = dists_bob_n125,
  opts_x_bob = optimals_x_bob_n125
)

# Export the results -----------------------------------------------------------

save(results_n125, file = "simulations/sim-varying-n/results_n125.Rdata")
