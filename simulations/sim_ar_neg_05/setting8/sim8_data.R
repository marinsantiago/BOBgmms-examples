# ------------------------------------------------------------------------------
# Sim setting 8: Generate the data
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

# Setting 8 --------------------------------------------------------------------
n <- 150
d <- 10
K <- 4
rho <- -0.5

total_sims <- 30

# Initialize data structures to store the data
data_setting8 <- vector("list", total_sims)

# Generate the data ------------------------------------------------------------
for (sim in seq_len(total_sims)) {
  
  cat("------------------------------------- \n")
  cat("Starting simulation:", sim, "\n")
  cat("------------------------------------- \n")
  
  set.seed(sim)
  data_sim <- generate.data(
    sample.size = n, data.dim = d, 
    n.clust = K, rho = rho
  )
  y <- data_sim$y
  Z <- data_sim$Z
  
  # Set prior and tuning hyper-parameters
  betas <- lapply(1:K, \(k) rep(0, d))
  psis <- lapply(1:K, \(k) diag(d))
  alphas <- rep(1.1, K)
  
  set.seed(1)
  shrink.params <- tryCatch(
    {
      shrinkage.cv.tuning(y, betas, psis, alphas)
    }, error = function(e){
      list(
        best.lambda = 4.5 * d/n,
        best.nu = 4.5 * d/n + (d + 1)
      )
    }
  )
  
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
  
  # Target Bayesian Posterior
  bayes.true <- bayesian_posterior(y, Z, betas, lambdas, nus, psis, alphas)
  bayes.pred <- post.predictive(bayes.true, d = d, K = K)
  
  # Store the results
  params.prior <- list(
    alphas = alphas, betas = betas, psis = psis,
    lambdas = lambdas, nus = nus
  )
  
  params.init <- list(
    means.init = means.init, covs.init = covs.init, probs.init = probs.init
  )
  
  params.temp <- list(a = a, b = b, c = c, r = r)
  
  iteration_i <- list(
    params.prior = params.prior, params.init = params.init,
    params.temp = params.temp, y = y, bayes.pred = bayes.pred
  )
  
  data_setting8[[sim]] <- iteration_i
}

names(data_setting8) <- paste0("rep", seq_len(total_sims))

# Export generated data --------------------------------------------------------
folder_path <- "simulations/sim_ar_neg_05/setting8/data_setting8.Rdata"
save(data_setting8, file = folder_path)
