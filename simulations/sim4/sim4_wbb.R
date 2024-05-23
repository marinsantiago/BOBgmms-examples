#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)

ncores <- parallel::detectCores()

# ------------------------------------------
# Sim4 - Weighted Bayesian Bootstrap
# ------------------------------------------

set.seed(1)

load("simulations/sim4/dataSim4.Rdata")

totalSims <- length(dataSim4)

K <- length(dataSim4$iteration1$params.prior$alphas)

# Initialize data structures to store the WBB output
wbb1Sim4 <- vector("list", totalSims)
wbb2Sim4 <- vector("list", totalSims)

# Total number of posterior draws
S <- 20000

# Cores
coresWBB   <- as.integer(ncores - 1)

# ----------------------------------
# Run WBB algorithms
# ----------------------------------

for (i in 1:totalSims) {
  
  cat("Starting simulation:", i, "\n")
  
  # Extract data for current iteration
  iteration_i <- dataSim4[[i]]
  y           <- iteration_i$y
  n           <- nrow(y)
  d           <- ncol(y)
  
  # Extract prior hyper-parameters
  alphas  <- iteration_i$params.prior$alphas
  betas   <- iteration_i$params.prior$betas
  psis    <- iteration_i$params.prior$psis
  lambdas <- iteration_i$params.prior$lambdas
  nus     <- iteration_i$params.prior$nus
  
  # Extract initial parameters
  mu    <- iteration_i$params.init$mu
  sigma <- iteration_i$params.init$sigma
  prob  <- iteration_i$params.init$prob
  
  init.params <- list(mu, sigma, prob)
  
  # Extract tempering parameters
  a.temp <- iteration_i$params.temp$a.temp
  b.temp <- iteration_i$params.temp$b.temp
  c.temp <- iteration_i$params.temp$c.temp
  r.temp <- iteration_i$params.temp$r.temp
  
  # ----------------------------------
  # Run algorithms
  # ----------------------------------
  
  # WBB1 - Random prior weights
  start.wbb1 <- Sys.time()
  wbb1.draws <- wbb.Sampler(y, mu, sigma, prob, S, a.temp, b.temp, c.temp,
                            r.temp, alphas, betas, lambdas, nus, psis,
                            randomPrior = T, cores = coresWBB)
  end.wbb1   <- Sys.time()
  Time.wbb1  <- difftime(end.wbb1, start.wbb1, units = "mins")
  
  # WBB2 - Fixed prior weights
  start.wbb2 <- Sys.time()
  wbb2.draws <- wbb.Sampler(y, mu, sigma, prob, S, a.temp, b.temp, c.temp,
                            r.temp, alphas, betas, lambdas, nus, psis,
                            randomPrior = F, cores = coresWBB)
  end.wbb2   <- Sys.time()
  Time.wbb2  <- difftime(end.wbb2, start.wbb2, units = "mins")
  
  cat("------------------------------------- \n")
  cat("Wall-clock time:                      \n")
  cat("Time WBB1:", Time.wbb1, "mins         \n")
  cat("Time WBB2:", Time.wbb2, "mins         \n")
  cat("------------------------------------- \n")
  
  # ----------------------------------------------
  # Posterior predictive distribution
  # ----------------------------------------------
  
  y.new.wbb1 <- post.predictive(wbb1.draws, d, K)
  y.new.wbb2 <- post.predictive(wbb2.draws, d, K)
  
  # ----------------------------------------------
  # Distances w.r.t Bayes Post. Pred.
  # ----------------------------------------------
  
  y.new.bayes <- iteration_i$bayes.pred
  
  TV.wbb1 <- distance.tv(y.new.bayes, y.new.wbb1)
  TV.wbb2 <- distance.tv(y.new.bayes, y.new.wbb2)
  
  cat("------------------------------------- \n")
  cat("TV Distances:                         \n")
  cat("TV-WBB1:", TV.wbb1,                  "\n")
  cat("TV-WBB2:", TV.wbb2,                  "\n")
  cat("------------------------------------- \n")
  
  KS.wbb1 <- distance.ks(y.new.bayes, y.new.wbb1)
  KS.wbb2 <- distance.ks(y.new.bayes, y.new.wbb2)
  
  cat("------------------------------------- \n")
  cat("KS Distances:                         \n")
  cat("KS-WBB1:", KS.wbb1,                  "\n")
  cat("KS-WBB2:", KS.wbb2,                  "\n")
  cat("------------------------------------- \n")
  
  # Save results
  results.wbb1 <- list(Time.wbb1 = Time.wbb1,
                       TV.wbb1   = TV.wbb1,
                       KS.wbb1   = KS.wbb1)
  
  results.wbb2 <- list(Time.wbb2 = Time.wbb2,
                       TV.wbb2   = TV.wbb2,
                       KS.wbb2   = KS.wbb2)
  
  wbb1Sim4[[i]] <- results.wbb1
  wbb2Sim4[[i]] <- results.wbb2
  
  rm(wbb1.draws, wbb2.draws)
  rm(y.new.wbb1, y.new.wbb2)
}

names(wbb1Sim4) <- paste0("iteration", 1:totalSims)
names(wbb2Sim4) <- paste0("iteration", 1:totalSims)

# Export wbb results
save(wbb1Sim4, file = "simulations/sim4/wbb1Sim4.Rdata")
save(wbb2Sim4, file = "simulations/sim4/wbb2Sim4.Rdata")
