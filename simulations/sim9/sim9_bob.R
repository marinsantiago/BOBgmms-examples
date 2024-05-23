#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)

ncores <- parallel::detectCores()

# ------------------------------------------
# Sim9 - Bayesian Optimized Bootstrap
# ------------------------------------------

set.seed(1)

load("simulations/sim9/dataSim9.Rdata")
load("simulations/sim9/wbb1Sim9.Rdata")
load("simulations/sim9/wbb2Sim9.Rdata")

totalSims <- length(dataSim9)

K <- length(dataSim9$iteration1$params.prior$alphas)

# Initialize data structures to store the BOB output
bobSim9  <- vector("list", totalSims)

# Total number of posterior draws
S <- 20000

# Mini-batch size
Sbatch <- 4000

# BO bounds
lowerBound <- c(1, rep(1e-5, K), rep(1e-5, K), 1e-5)
upperBound <- c(1.5, rep(1.5, K),  rep(1.5,  K), 1.5)

# Cores
coresBO    <- as.integer(ncores - 1)
coresSamp  <- as.integer(ncores - 1)

# ----------------------------------
# Run RW algorithms
# ----------------------------------

for (i in 1:totalSims) {
  
  cat("Starting simulation:", i, "\n")
  
  # Extract data for current iteration
  iteration_i <- dataSim9[[i]]
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
  
  # BOB - Bayesian Optimized Bootstrap
  start.bob <- Sys.time()
  bob.fit   <- bob.Sampler(y, mu, sigma, prob, S, Sbatch, a.temp, b.temp, c.temp,
                           r.temp, alphas, betas, lambdas, nus, psis,
                           lowerBound, upperBound, coresBO, coresSamp)
  end.bob    <- Sys.time()
  (bob.time  <- difftime(end.bob, start.bob, units = "mins"))
  
  bob.draws <- bob.fit$Draws
  bob.xvals <- bob.fit$x.optim
  
  rm(bob.fit)
  
  cat("------------------------------------- \n")
  cat("Wall-clock time:                      \n")
  cat("Time BOB:",  bob.time,          "mins \n")
  cat("------------------------------------- \n")
  
  # ----------------------------------------------
  # Posterior predictive distribution
  # ----------------------------------------------
  
  y.new.bob <- post.predictive(bob.draws, d, K)
  
  # ----------------------------------------------
  # Distances w.r.t Bayes Post. Pred.
  # ----------------------------------------------
  
  y.new.bayes <- iteration_i$bayes.pred
  
  TV.bob <- distance.tv(y.new.bayes, y.new.bob)  
  
  cat("------------------------------------- \n")
  cat("TV Distances:                         \n")
  cat("TV-BOB:",  TV.bob,                   "\n")
  cat("------------------------------------- \n")
  
  KS.bob <- distance.ks(y.new.bayes, y.new.bob)
  
  cat("------------------------------------- \n")
  cat("KS Distances:                         \n")
  cat("KS-BOB:",  KS.bob,                   "\n")
  cat("------------------------------------- \n")
  
  # Save results
  results.bob <- list(Time.bob = bob.time,
                      TV.bob   = TV.bob,
                      KS.bob   = KS.bob)
  
  bobSim9[[i]]  <- results.bob
  
  rm(bob.draws, y.new.bob)
}

names(bobSim9) <- paste0("iteration", 1:totalSims)

# Export bob results
save(bobSim9, file = "simulations/sim9/bobSim9.Rdata")
