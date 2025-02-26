# ------------------------------------------------------------------------------
# Sim setting 2: Stan
# ------------------------------------------------------------------------------

library(BOBgmms)

source("R/generate_data.R")
source("R/init.R")
source("R/tune.R")
source("R/bayes_post.R")
source("R/predictive.R")
source("R/dens_plots.R")
source("R/advi.R")
source("R/nuts.R")

load("simulations/sim_ar_neg_05/setting2/data_setting2.Rdata")

# Simulation settings
total_sims <- length(data_setting2)
K <- length(data_setting2$rep1$params.prior$alphas)
data_dims <- dim(data_setting2$rep1$y)
n <- data_dims[1]
d <- data_dims[2]

# Initialize data structures to store the data
nuts_setting2 <- vector("list", total_sims)
advi_setting2 <- vector("list", total_sims)

for (sim in seq_len(total_sims)) {
  
  cat("------------------------------------- \n")
  cat("Starting simulation:", sim, "\n")
  cat("------------------------------------- \n")
  
  # Extract data and hyper-parameters from current replication
  current_rep <- data_setting2[[sim]]
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
  
  # NUTS -----------------------------------------------------------------------
  set.seed(1)
  start.nuts <- Sys.time()
  nuts.out <- tryCatch(
    {
      nuts(
        y = y, means.init = means.init, covs.init = covs.init,
        probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
        psis = psis, alphas = alphas, max.iters = 40000, warm_init = TRUE
      )
    }, error = function(e) {
      NA
    }
  )
  end.nuts <- Sys.time()
  if (sum(is.na(nuts.out)) > 0) {
    nuts.time <- NA
  } else {
    nuts.time <- difftime(end.nuts, start.nuts, units = "secs")
  }
  
  # ADVI -----------------------------------------------------------------------
  set.seed(1)
  start.advi <- Sys.time()
  advi.out <- tryCatch(
    {
      advi(
        y = y, means.init = means.init, covs.init = covs.init,
        probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
        psis = psis, alphas = alphas, max.iters = 20000, vb.model = stan.model
      )
    }, error = function(e) {
      NA
    }
  )
  end.advi <- Sys.time()
  if (sum(is.na(advi.out)) > 0) {
    advi.time <- NA
  } else {
    advi.time <- difftime(end.advi, start.advi, units = "secs")
  }
  
  # Posterior predictive distributions -----------------------------------------
  if (sum(is.na(nuts.out)) > 0) {
    nuts.pred <- NA
  } else {
    nuts.pred <- post.predictive(nuts.out, d = d, K = K)
  }
  
  if (sum(is.na(advi.out)) > 0) {
    advi.pred <- NA
  } else {
    advi.pred <- post.predictive(advi.out, d = d, K = K)
  }
  
  # Distances  -----------------------------------------------------------------
  set.seed(1)
  bayes.pred <- current_rep$bayes.pred
  if (sum(is.na(nuts.out)) > 0) {
    dist.nuts <- NA
  } else {
    dist.nuts <- T4transport::swdist(bayes.pred, nuts.pred)$distance
  }
  
  if (sum(is.na(advi.out)) > 0) {
    dist.advi <- NA
  } else {
    dist.advi <- T4transport::swdist(bayes.pred, advi.pred)$distance
  }
  
  # Save the results
  results.nuts <- list(
    sw2.dist = dist.nuts,
    run.time = nuts.time
  )
  
  results.advi <- list(
    sw2.dist = dist.advi,
    run.time = advi.time
  )
  
  nuts_setting2[[sim]] <- results.nuts
  advi_setting2[[sim]] <- results.advi
  
  rm(
    nuts.out, advi.out, nuts.pred, advi.pred
  )
  gc()
}

names(nuts_setting2) <- paste0("rep", seq_len(total_sims))
names(advi_setting2) <- paste0("rep", seq_len(total_sims))

# Export the Stan results
folder_path <- "simulations/sim_ar_neg_05/setting2"
save(nuts_setting2, file = paste(folder_path, "/nuts_sim2.Rdata", sep = ""))
save(advi_setting2, file = paste(folder_path, "/advi_sim2.Rdata", sep = ""))
