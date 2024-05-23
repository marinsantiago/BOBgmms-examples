#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/nuts.R"); source("./R/advi.R")
rm(shhh)

# ------------------------------------------
# Sim1 - NUTS and ADVI
# ------------------------------------------

set.seed(1)

load("simulations/sim1/dataSim1.Rdata")
load("simulations/sim1/wbb1Sim1.Rdata")
load("simulations/sim1/wbb2Sim1.Rdata")
load("simulations/sim1/bobSim1.Rdata")

totalSims <- length(dataSim1)

# Initialize data structures to store the nuts output
nutsSim1 <- vector("list", totalSims)
adviSim1 <- vector("list", totalSims)

S <- 20000

# Run nuts algorithm
for (i in 1:totalSims) {
  
  cat("Starting simulation:", i, "\n")
  
  # Extract data for current iteration
  iteration_i <- dataSim1[[i]]
  y           <- iteration_i$y
  n           <- nrow(y)
  d           <- ncol(y)
  K           <- length(iteration_i$params.prior$alphas)
  
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
  
  # ----------------------------------
  # Run NUTS - Stan
  # ----------------------------------

  start.nuts <- Sys.time()
  nuts.draws <- nuts.Sampler(y, init.params, S * 2,
                             alphas, betas, lambdas, nus, psis)
  end.nuts   <- Sys.time()
  Time.nuts  <- difftime(end.nuts, start.nuts, units = "mins")

  cat("------------------------------------- \n")
  cat("Wall-clock time:                      \n")
  cat("Time NUTS:", Time.nuts,         "mins \n")
  cat("------------------------------------- \n")
  
  # ----------------------------------
  # Run ADVI - Stan
  # ----------------------------------
  
  start.advi <- Sys.time()
  advi.draws <- advi(y, init.params, S, alphas, betas, lambdas, nus, psis)
  end.advi   <- Sys.time()
  Time.advi  <- difftime(end.advi, start.advi, units = "mins")
  
  cat("------------------------------------- \n")
  cat("Wall-clock time:                      \n")
  cat("Time ADVI:", Time.advi,         "mins \n")
  cat("------------------------------------- \n")
  
  # ----------------------------------------------
  # Posterior predictive distribution
  # ----------------------------------------------

  y.new.nuts <- post.predictive(nuts.draws, d, K)
  y.new.advi <- post.predictive(advi.draws, d, K)
  
  # ----------------------------------------------
  # Distances w.r.t Bayes Post. Pred.
  # ----------------------------------------------
  y.new.bayes <- iteration_i$bayes.pred
  
  # TV distances  
  nuts.tv <- distance.tv(y.new.bayes, y.new.nuts)
  advi.tv <- distance.tv(y.new.bayes, y.new.advi)
  
  # KS distances
  nuts.ks <- distance.ks(y.new.bayes, y.new.nuts)
  advi.ks <- distance.ks(y.new.bayes, y.new.advi)

  # Save results
  results.nuts <- list(Time.nuts = Time.nuts,
                       TV.nuts   = nuts.tv,
                       KS.nuts   = nuts.ks)
  
  results.advi <- list(Time.advi = Time.advi,
                       TV.advi   = advi.tv,
                       KS.advi   = advi.ks)
  
  nutsSim1[[i]] <- results.nuts
  adviSim1[[i]] <- results.advi
  
  rm(nuts.draws, advi.draws)
  rm(y.new.nuts, y.new.advi)
}

names(nutsSim1) <- paste0("iteration", 1:totalSims)
names(adviSim1) <- paste0("iteration", 1:totalSims)

# Export NUTS and ADVI results
save(nutsSim1, file = "simulations/sim1/nutsSim1.Rdata")
save(adviSim1, file = "simulations/sim1/adviSim1.Rdata")
