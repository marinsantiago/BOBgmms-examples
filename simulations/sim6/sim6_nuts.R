#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/nuts.R"); source("./R/advi.R")
rm(shhh)

# ------------------------------------------
# Sim6 - NUTS and ADVI
# ------------------------------------------

set.seed(1)

load("simulations/sim6/dataSim6.Rdata")
load("simulations/sim6/wbb1Sim6.Rdata")
load("simulations/sim6/wbb2Sim6.Rdata")
load("simulations/sim6/bobSim6.Rdata")

totalSims <- length(dataSim6)

# Initialize data structures to store the nuts output
nutsSim6 <- vector("list", totalSims)
adviSim6 <- vector("list", totalSims)

S <- 20000

# Run nuts algorithm
for (i in 1:totalSims) {
  
  cat("Starting simulation:", i, "\n")
  
  # Extract data for current iteration
  iteration_i <- dataSim6[[i]]
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
  nuts.draws <- tryCatch(
    {
      nuts.Sampler(y, init.params, S * 2, alphas, betas, lambdas, nus, psis)
    }, error = function(e){
      NA
    }
  )
  end.nuts   <- Sys.time()
  if (sum(is.na(nuts.draws)) > 0) {
    Time.nuts <- NA
  } else {
    Time.nuts  <- difftime(end.nuts, start.nuts, units = "mins")
  }
  
  cat("------------------------------------- \n")
  cat("Wall-clock time:                      \n")
  cat("Time NUTS:", Time.nuts,         "mins \n")
  cat("------------------------------------- \n")
  
  # ----------------------------------
  # Run ADVI - Stan
  # ----------------------------------
  
  start.advi <- Sys.time()
  advi.draws <- tryCatch(
    {
      advi(y, init.params, S, alphas, betas, lambdas, nus, psis)
    }, error = function(e){
      NA
    }
  )
  end.advi   <- Sys.time()
  if (sum(is.na(advi.draws)) > 0) {
    Time.advi <- NA
  } else {
    Time.advi  <- difftime(end.advi, start.advi, units = "mins")
  }
  
  cat("------------------------------------- \n")
  cat("Wall-clock time:                      \n")
  cat("Time ADVI:", Time.advi,         "mins \n")
  cat("------------------------------------- \n")
  
  # ----------------------------------------------
  # Posterior predictive distribution
  # ----------------------------------------------
  
  if (sum(is.na(nuts.draws)) > 0) {
    y.new.nuts <- NA
  } else {
    y.new.nuts <- post.predictive(nuts.draws, d, K)
  }
  
  if (sum(is.na(advi.draws)) > 0) {
    y.new.advi <- NA
  } else {
    y.new.advi <- post.predictive(advi.draws, d, K)
  }
  
  # ----------------------------------------------
  # Distances w.r.t Bayes Post. Pred.
  # ----------------------------------------------
  y.new.bayes <- iteration_i$bayes.pred
  
  # TV distances  
  if (sum(is.na(nuts.draws)) > 0) {
    nuts.tv <- NA
  } else {
    nuts.tv <- distance.tv(y.new.bayes, y.new.nuts)
  }
  
  if (sum(is.na(advi.draws)) > 0) {
    advi.tv <- NA
  } else {
    advi.tv <- distance.tv(y.new.bayes, y.new.advi)
  }
  
  # KS distances
  if (sum(is.na(nuts.draws)) > 0) {
    nuts.ks <- NA
  } else {
    nuts.ks <- distance.ks(y.new.bayes, y.new.nuts)
  }
  
  if (sum(is.na(advi.draws)) > 0) {
    advi.ks <- NA
  } else {
    advi.ks <- distance.ks(y.new.bayes, y.new.advi)
  }
  
  # Save results
  results.nuts <- list(Time.nuts = Time.nuts,
                       TV.nuts   = nuts.tv,
                       KS.nuts   = nuts.ks)
  
  results.advi <- list(Time.advi = Time.advi,
                       TV.advi   = advi.tv,
                       KS.advi   = advi.ks)
  
  nutsSim6[[i]] <- results.nuts
  adviSim6[[i]] <- results.advi
  
  rm(nuts.draws, advi.draws)
  rm(y.new.nuts, y.new.advi)
}

names(nutsSim6) <- paste0("iteration", 1:totalSims)
names(adviSim6) <- paste0("iteration", 1:totalSims)

# Export NUTS and ADVI results
save(nutsSim6, file = "simulations/sim6/nutsSim6.Rdata")
save(adviSim6, file = "simulations/sim6/adviSim6.Rdata")
