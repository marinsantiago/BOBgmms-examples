#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/bayesPost.R"); source("./R/tuning.R")
source("./R/generateData.R")
rm(shhh)

ncores <- parallel::detectCores()

set.seed(1)

# ----------------------------------------
# Sim8 - Data generation
# ----------------------------------------

n   <- 150
d   <- 10
K   <- 4

# Set prior hyper-parameters
alphas <- rep(1.1, K)
betas  <- lapply(1:K, function(x) rep(0, d))
psis   <- lapply(1:K, function(x) diag(d))

# Grid tempering hyper-parameters
range.a <- seq(0, 0)
range.b <- seq(0, 4, by = 0.2)
range.c <- seq(1, 4, by = 0.2)
range.r <- seq(1, 4, by = 0.2)

S <- 20000

totalSims <- 10

# Initialize data structures to store the data
dataSim8 <- vector("list", totalSims)

# ----------------------------------
# Generate the data
# ----------------------------------

for (i in 1:totalSims) {
  
  cat("Starting simulation:", i, "\n")
  
  # Generate the data
  Data <- generate.data(d, n, K)
  y    <- Data$y
  Z    <- Data$Z
  
  # The regularization path for the HTK-means shrinakge parameter
  clusterHD::diagPlot(clusterHD::HTKmeans(y, K)) # 0.75
  abline(v = 1 - 0.75, col = 2)
  
  # Select lambda and nu via cross-validation
  shrink.params <- cv.tuning(y, alphas, betas, psis, 0.75, init = "HT", mult = 2)
  lambdas       <- rep(shrink.params$lambda, K)
  nus           <- rep(shrink.params$nu,     K)
  
  cat("------------------------------------- \n")
  cat("Shrinkage hyper-params:               \n")
  cat("lambda:", shrink.params$lambda,      "\n")
  cat("nu:",     shrink.params$nu,          "\n")
  cat("------------------------------------- \n")
  
  # Initialize model parameters
  init.params <- initialVals.HT(y, alphas, betas, lambdas, nus, psis, 0.75)
  
  mu    <- init.params$means
  sigma <- init.params$covs
  prob  <- init.params$probs
  
  # Grid search tempering hyper-parameters
  temp.params <- temp.tuning(y, mu, sigma, prob,
                             range.a, range.b, range.c, range.r,
                             alphas, betas, lambdas, nus, psis, ncores - 1)
  a.temp <- temp.params$a 
  b.temp <- temp.params$b
  c.temp <- temp.params$c
  r.temp <- temp.params$r
  
  cat("------------------------------------- \n")
  cat("Tempering hyper-params:               \n")
  cat("a:", a.temp,                         "\n")
  cat("b:", b.temp,                         "\n")
  cat("c:", c.temp,                         "\n")
  cat("r:", r.temp,                         "\n")
  cat("------------------------------------- \n")
  
  # -----------------------------------------
  # True Bayesian Posterior
  # -----------------------------------------
  
  bayes.draws <- bayesianPosterior(y, Z, 2 * S, alphas, betas, lambdas, nus, psis)
  y.new.bayes <- post.predictive(bayes.draws, d, K)
  
  # Store the results
  
  params.prior <- list(alphas  = alphas, 
                       betas   = betas,
                       psis    = psis,
                       lambdas = lambdas,
                       nus     = nus)
  
  params.init  <- list(sigma   = sigma,
                       mu      = mu,
                       prob    = prob)
  
  params.temp  <- list(a.temp  = a.temp,
                       b.temp  = b.temp,
                       c.temp  = c.temp,
                       r.temp  = r.temp)
  
  iteration_i  <- list(params.prior = params.prior,
                       params.init  = params.init,
                       params.temp  = params.temp,
                       y            = y,
                       bayes.pred   = y.new.bayes)
  
  dataSim8[[i]] <- iteration_i
  
}

names(dataSim8) <- paste0("iteration", 1:totalSims)

# Export generated data
save(dataSim8, file = "simulations/sim8/dataSim8.Rdata")
