#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/scaleFeatures.R"); source("./R/bayesPost.R"); source("./R/tuning.R")
shhh(library(coda, include.only = "effectiveSize"))
rm(shhh)

ncores <- parallel::detectCores()

set.seed(100)

# -----------------------------------------
# Wine data-set
# -----------------------------------------

load("wines/wine.Rdata")

sum(wines[,1] == "Barolo")           # 59 Barolo type
sum(wines[,1] == "Grignolino")       # 71 Grignolino type
sum(wines[,1] == "Barbera")          # 48 Barbera type

# Split the data into train and held-out
n.test <- 78
n <- nrow(wines) - n.test 

shuflle.indxs <- sample(1:nrow(wines), nrow(wines), replace = F)
train.indxs   <- shuflle.indxs[1:n]
test.indxs    <- shuflle.indxs[(n + 1):nrow(wines)]

# Prepare and standardize the data
Data.scaled  <- train.testScale(wines[train.indxs, - 1], wines[test.indxs, - 1])
labels.train <- wines[train.indxs, 1]
labels.test  <- wines[test.indxs,  1]

y <- Data.scaled$y.train.scaled; rownames(y) <- NULL
K <- length(levels(labels.train))
d <- ncol(y)
Z <- do.call(cbind, lapply(1:K, \(k) ifelse(as.numeric(labels.train) == k, 1, 0)))
sum(Z[,1])                            # 34 Barolo type - train
sum(Z[,2])                            # 40 Grignolino type - train
sum(Z[,3])                            # 26 Barbera type - train

y.test <- Data.scaled$y.test.scaled; row.names(y.test) <- NULL
Z.test <- do.call(cbind, lapply(1:K, \(k) ifelse(as.numeric(labels.test) == k, 1, 0)))
sum(Z.test[,1])                       # 25 Barolo type - test
sum(Z.test[,2])                       # 31 Grignolino type - test
sum(Z.test[,3])                       # 22 Barbera type - test

rm(Data.scaled); rm(wines); rm(shuflle.indxs)
rm(test.indxs); rm(train.indxs)

# Set prior hyper-parameters
alphas <- rep(1.1, K)
betas  <- lapply(1:K, function(x) rep(0, d))
psis   <- lapply(1:K, function(x) diag(d))

# Grid tempering hyper-parameters
range.a <- seq(0, 0)
range.b <- seq(0, 4, by = 0.2)
range.c <- seq(1, 4, by = 0.2)
range.r <- seq(1, 4, by = 0.2)

# Total number of posterior draws
S <- 20000

# Mini-batch size
Sbatch <- 4000

# Select lambda and nu via cross-validation
shrink.params <- cv.tuning(y, alphas, betas, psis)
(lambdas <- rep(shrink.params$lambda, K))
(nus     <- rep(shrink.params$nu,     K))

# Initialize model parameters
init.params <- initialVals(y, alphas, betas, lambdas, nus, psis)
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

# WBB1 - Random prior weights
start.wbb1 <- Sys.time()
wbb1.draws <- wbb.Sampler(y, mu, sigma, prob, S, a.temp, b.temp, c.temp,
                          r.temp, alphas, betas, lambdas, nus, psis,
                          randomPrior = T, cores = coresWBB)
end.wbb1   <- Sys.time()
(wbb1.time <- difftime(end.wbb1, start.wbb1, units = "mins"))

#mean(wbb1.draws[,1]) # 0.8824324

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

# BO bounds 
lowerBound <- c(1.0, rep(1e-5, K), rep(1e-5, K), 1e-5)
upperBound <- c(1.5, rep(1.5, K),  rep(1.5,  K), 1.5)

coresBO    <- as.integer(ncores - 1)
coresSamp  <- as.integer(ncores - 1)

# BOB - Bayesian Optimized Bootstrap
start.bob  <- Sys.time()
bob.fit    <- bob.Sampler(y, mu, sigma, prob, S, Sbatch, a.temp, b.temp, c.temp,
                          r.temp, alphas, betas, lambdas, nus, psis,
                          lowerBound, upperBound, coresBO, coresSamp, seed = 10)
end.bob    <- Sys.time()
(bob.time  <- difftime(end.bob, start.bob, units = "mins"))

bob.draws <- bob.fit$Draws
bob.xvals <- bob.fit$x.optim

rm(bob.fit)

#mean(bob.draws[,1]) #0.8822172

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

save(y, file = "wines/y_data.Rdata")
save(priors,      file = "wines/priors.Rdata")
save(init.params, file = "wines/init_params.Rdata")
save(y.test,      file = "wines/y_data_test.Rdata")
save(labels.test, file = "wines/labels_test.Rdata")

save(bob.draws,   file = "wines/bob_draws.Rdata")
save(wbb1.draws,  file = "wines/wbb1_draws.Rdata")
save(wbb2.draws,  file = "wines/wbb2_draws.Rdata")
save(bayes.draws, file = "wines/bayes_draws.Rdata")

save(y.new.bob,   file = "wines/y_new_bob.Rdata")
save(y.new.wbb1,  file = "wines/y_new_wbb1.Rdata")
save(y.new.wbb2,  file = "wines/y_new_wbb2.Rdata")
save(y.new.bayes, file = "wines/y_new_bayes.Rdata")

save(bob.time,   file = "wines/bob_time.Rdata")
save(wbb1.time,  file = "wines/wbb1_time.Rdata")
save(wbb2.time,  file = "wines/wbb2_time.Rdata")
