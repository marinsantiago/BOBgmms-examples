#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/bayesPost.R"); source("./R/tuning.R"); source("./R/generateData.R")
rm(shhh)

ncores <- parallel::detectCores()

totalSims <- 10

# Initialize data structures to store the output
wbb1Sim.n2 <- vector("list", totalSims)
wbb2Sim.n2 <- vector("list", totalSims)
bobSim.n2  <- vector("list", totalSims)

# Simulation settings
n <- 125
d <- 15
K <- 2

# Set prior hyper-parameters
alphas <- rep(1.1, K)
betas  <- lapply(1:K, function(x) rep(0, d))
psis   <- lapply(1:K, function(x) diag(d))

# Grid tempering hyper-parameters
range.a <- seq(0, 0)
range.b <- seq(0, 4, by = 0.2)
range.c <- seq(1, 4, by = 0.2)
range.r <- seq(1, 4, by = 0.2)

coresWBB   <- as.integer(ncores - 1)
coresBO    <- as.integer(ncores - 1)
coresSamp  <- as.integer(ncores - 1)

# Posterior draws
S <- 20000; Sbatch <- 4000

# BO bounds 
lowerBound <- c(1, rep(1e-5, K), rep(1e-5, K), 1e-5)
upperBound <- c(1.5, rep(1.5, K),  rep(1.5,  K), 1.5)

# Run simulations
for (i in 1:totalSims) {
  
  set.seed(i) # For replication purposes
  
  # Generate data
  Data <- generate.data(d, n, K)
  y    <- Data$y
  Z    <- Data$Z
  
  # The regularization path for the HTK-means shrinakge parameter
  clusterHD::diagPlot(clusterHD::HTKmeans(y, K)) # 0.75
  
  # Select lambda and nu via cross-validation
  shrink.params <- cv.tuning(y, alphas, betas, psis, 0.75, init = "Pool", mult = 2)
  shrink.params$lambda
  shrink.params$nu
  
  lambdas <- rep(shrink.params$lambda, K)
  nus     <- rep(shrink.params$nu,     K)
  
  # Initialize model parameters
  init.params <- initialVals(y, alphas, betas, lambdas, nus, psis, 0.75)
  init.params$algorithm
  
  mu    <- init.params$means
  sigma <- init.params$covs
  prob  <- init.params$probs
  
  # Grid search tempering hyper-parameters
  temp.params <- temp.tuning(y, mu, sigma, prob,
                             range.a, range.b, range.c, range.r,
                             alphas, betas, lambdas, nus, psis, ncores - 1)
  (a.temp <- temp.params$a)
  (b.temp <- temp.params$b)
  (c.temp <- temp.params$c)
  (r.temp <- temp.params$r)
  
  # -----------------------------------------
  # Weighted Bayesian Bootstrap
  # -----------------------------------------
  
  # WBB1 - Random prior weights
  start.wbb1 <- Sys.time()
  wbb1.draws <- wbb.Sampler(y, mu, sigma, prob, S, a.temp, b.temp, c.temp,
                            r.temp, alphas, betas, lambdas, nus, psis,
                            randomPrior = T, cores = coresWBB)
  end.wbb1   <- Sys.time()
  (wbb1.time <- difftime(end.wbb1, start.wbb1, units = "mins"))
  
  # WBB2 - Fixed prior weights
  start.wbb2 <- Sys.time()
  wbb2.draws <- wbb.Sampler(y, mu, sigma, prob, S, a.temp, b.temp, c.temp,
                            r.temp, alphas, betas, lambdas, nus, psis,
                            randomPrior = F, cores = coresWBB)
  end.wbb2   <- Sys.time()
  (wbb2.time <- difftime(end.wbb2, start.wbb2, units = "mins"))
  
  # -----------------------------------------
  # Bayesian Optimized Bootstrap
  # -----------------------------------------
  
  start.bob <- Sys.time()
  bob.fit   <- bob.Sampler(y, mu, sigma, prob, S, Sbatch, a.temp, b.temp, c.temp,
                           r.temp, alphas, betas, lambdas, nus, psis,
                           lowerBound, upperBound, coresBO, coresSamp)
  end.bob    <- Sys.time()
  (bob.time  <- difftime(end.bob, start.bob, units = "mins"))
  
  bob.draws <- bob.fit$Draws
  bob.xvals <- bob.fit$x.optim
  
  rm(bob.fit)
  
  # -----------------------------------------
  # True Bayesian Posterior
  # -----------------------------------------
  
  bayes.draws <- bayesianPosterior(y, Z, 2 * S, alphas, betas, lambdas, nus, psis)
  
  # -----------------------------------------
  # Posterior predictive distributions
  # -----------------------------------------
  
  y.new.wbb1  <- post.predictive(wbb1.draws,  d, K)
  y.new.wbb2  <- post.predictive(wbb2.draws,  d, K)
  y.new.bob   <- post.predictive(bob.draws,   d, K)
  y.new.bayes <- post.predictive(bayes.draws, d, K)
  
  # -----------------------------------------
  # Distances
  # -----------------------------------------
  
  # TV distances
  (wbb1.tv <- distance.tv(y.new.bayes, y.new.wbb1))
  (wbb2.tv <- distance.tv(y.new.bayes, y.new.wbb2))
  (bob.tv  <- distance.tv(y.new.bayes, y.new.bob)) 
  
  # KS distances
  (wbb1.ks <- distance.ks(y.new.bayes, y.new.wbb1))
  (wbb2.ks <- distance.ks(y.new.bayes, y.new.wbb2))
  (bob.ks  <- distance.ks(y.new.bayes, y.new.bob))
  
  # -----------------------------------------
  # Distances
  # -----------------------------------------
  
  results.wbb1 <- list(Time.wbb1 = wbb1.time,
                       TV.wbb1   = wbb1.tv,
                       KS.wbb1   = wbb1.ks)
  
  results.wbb2 <- list(Time.wbb2 = wbb2.time,
                       TV.wbb2   = wbb2.tv,
                       KS.wbb2   = wbb2.ks)
  
  results.bob  <- list(Time.bob  = bob.time,
                       TV.bob    = bob.tv,
                       KS.bob    = bob.ks)
  
  wbb1Sim.n2[[i]] <- results.wbb1
  wbb2Sim.n2[[i]] <- results.wbb2
  bobSim.n2[[i]]  <- results.bob
  
  rm(wbb1.draws, wbb2.draws, bob.draws, bayes.draws)
  rm(y.new.wbb1, y.new.wbb2, y.new.bob, y.new.bayes)
}

names(wbb1Sim.n2) <- paste0("iteration", 1:totalSims)
names(wbb2Sim.n2) <- paste0("iteration", 1:totalSims)
names(bobSim.n2)  <- paste0("iteration", 1:totalSims)

# Export wbb results
save(wbb1Sim.n2, file = "simulations/sim-varying-n/sim_n2/wbb1Sim_n2.Rdata")
save(wbb2Sim.n2, file = "simulations/sim-varying-n/sim_n2/wbb2Sim_n2.Rdata")
save(bobSim.n2,  file = "simulations/sim-varying-n/sim_n2/bobSim_n2.Rdata")
