#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/bayesPost.R"); source("./R/tuning.R"); source("./R/generateData.R")
rm(shhh)

ncores <- parallel::detectCores()

set.seed(101)

# --------------------------------------------
# Generate the data and set hyper-parameters
# --------------------------------------------

n   <- 50
d   <- 10
K   <- 2

# Set prior hyper-parameters
alphas <- rep(1.1, K)
betas  <- lapply(1:K, function(x) rep(0, d))
psis   <- lapply(1:K, function(x) diag(d))

# Grid tempering hyper-parameters
range.a <- seq(0, 0)
range.b <- seq(0, 4, by = 0.2)
range.c <- seq(1, 4, by = 0.2)
range.r <- seq(1, 4, by = 0.2)

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

coresWBB <- as.integer(ncores - 1)

S <- 20000

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

Sbatch <- 4000

# BO bounds 
lowerBound <- c(1, rep(1e-5, K), rep(1e-5, K), 1e-5)
upperBound <- c(1.5, rep(1.5, K),  rep(1.5,  K), 1.5)

coresBO    <- as.integer(ncores - 1)
coresSamp  <- as.integer(ncores - 1)

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

# Export results

priors <- list(alphas = alphas, betas = betas, lambdas = lambdas,
               nus = nus, psis = psis)

save(y, file = "simulations/sim-illustrative/y_data.Rdata")
save(priors,      file = "simulations/sim-illustrative/priors.Rdata")
save(init.params, file = "simulations/sim-illustrative/init_params.Rdata")

save(bob.draws,   file = "simulations/sim-illustrative/bob_draws.Rdata")
save(wbb1.draws,  file = "simulations/sim-illustrative/wbb1_draws.Rdata")
save(wbb2.draws,  file = "simulations/sim-illustrative/wbb2_draws.Rdata")
save(bayes.draws, file = "simulations/sim-illustrative/bayes_draws.Rdata")

save(y.new.bob,   file = "simulations/sim-illustrative/y_new_bob.Rdata")
save(y.new.wbb1,  file = "simulations/sim-illustrative/y_new_wbb1.Rdata")
save(y.new.wbb2,  file = "simulations/sim-illustrative/y_new_wbb2.Rdata")
save(y.new.bayes, file = "simulations/sim-illustrative/y_new_bayes.Rdata")

save(bob.time,   file = "simulations/sim-illustrative/bob_time.Rdata")
save(wbb1.time,  file = "simulations/sim-illustrative/wbb1_time.Rdata")
save(wbb2.time,  file = "simulations/sim-illustrative/wbb2_time.Rdata")
