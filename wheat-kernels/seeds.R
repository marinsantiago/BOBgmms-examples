#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/nuts.R"); source("./R/lppd.R"); source("./R/scaleFeatures.R")
source("./R/densPlots.R"); source("./R/bayesPost.R"); source("./R/tuning.R")
#source("./R/advi.R")
shhh(library(coda, include.only = "effectiveSize"))
rm(shhh)

ncores <- parallel::detectCores()

set.seed(100)

# -----------------------------------------
# Seeds data-set
# -----------------------------------------

load("wheat-kernels/seeds.Rdata")

sum(seeds[,8] == "Kama")             # 70 Kama variety
sum(seeds[,8] == "Rosa")             # 70 Rosa variety
sum(seeds[,8] == "Canadian")         # 70 Canadian 

# Split the data into train and held-out
n.test <- 100
n <- nrow(seeds) - n.test 

shuflle.indxs <- sample(1:nrow(seeds), nrow(seeds), replace = F)
train.indxs   <- shuflle.indxs[1:n]
test.indxs    <- shuflle.indxs[(n + 1):nrow(seeds)]

# Prepare and standardize the data
Data.scaled  <- train.testScale(seeds[train.indxs, - 8], seeds[test.indxs, - 8])
labels.train <- seeds[train.indxs, 8]
labels.test  <- seeds[test.indxs,  8]

y <- Data.scaled$y.train.scaled; rownames(y) <- NULL
K <- length(levels(labels.train))
d <- ncol(y)
Z <- do.call(cbind, lapply(1:K, \(k) ifelse(as.numeric(labels.train) == k, 1, 0)))
sum(Z[,1])                            # 33 Kama type - train
sum(Z[,2])                            # 38 Rosa type - train
sum(Z[,3])                            # 39 Canadian type - train

y.test <- Data.scaled$y.test.scaled; row.names(y.test) <- NULL
Z.test <- do.call(cbind, lapply(1:K, \(k) ifelse(as.numeric(labels.test) == k, 1, 0)))
sum(Z.test[,1])                       # 37 Kama type - test
sum(Z.test[,2])                       # 32 Rosa type - test
sum(Z.test[,3])                       # 31 Canadian type - test

rm(Data.scaled); rm(seeds); rm(shuflle.indxs)
rm(test.indxs); rm(train.indxs)

pairs(y)
pairs(y.test)

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

#mean(wbb1.draws[,1]) # -0.1168917

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
#upperBound <- c(d/3.25, rep(d/3.25, K),  rep(d/3.25,  K), d/3.25)
upperBound <- c(1.5, rep(1.5, K),  rep(1.5,  K), 1.5)

coresBO    <- as.integer(ncores - 1)
coresSamp  <- as.integer(ncores - 1)

# BOB - Bayesian Optimized Bootstrap
start.bob  <- Sys.time()
bob.fit    <- bob.Sampler(y, mu, sigma, prob, S, Sbatch, a.temp, b.temp, c.temp,
                          r.temp, alphas, betas, lambdas, nus, psis,
                          lowerBound, upperBound, coresBO, coresSamp, seed = 20)
end.bob    <- Sys.time()
(bob.time  <- difftime(end.bob, start.bob, units = "mins"))

bob.draws <- bob.fit$Draws
bob.xvals <- bob.fit$x.optim

# Loss: 2.891837
#mean(bob.draws[,1]) # -0.1152939

rm(bob.fit)

# -----------------------------------------
# No-U-Turn-Sampler (NUTS)
# -----------------------------------------
set.seed(1)

# NUTS - Stan
start.nuts <- Sys.time()
nuts.draws <- nuts.Sampler(y, init.params, S * 2,
                           alphas, betas, lambdas, nus, psis)
end.nuts   <- Sys.time()
(nuts.time <- difftime(end.nuts, start.nuts, units = "mins"))


# -----------------------------------------
# ADVI - Mean Field Variational Bayes
# -----------------------------------------

# ADVI - Stan
start.advi <- Sys.time()
advi.draws <- advi(y, init.params, S, alphas, betas, lambdas, nus, psis)
end.advi   <- Sys.time()
(advi.time <- difftime(end.advi, start.advi, units = "mins"))

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
y.new.nuts  <- post.predictive(nuts.draws,  d, K)
y.new.advi  <- post.predictive(advi.draws,  d, K)
y.new.bayes <- post.predictive(bayes.draws, d, K)

# -----------------------------------------
# Distances
# -----------------------------------------

# TV distances
(wbb1.tv <- distance.tv(y.new.bayes, y.new.wbb1))
(wbb2.tv <- distance.tv(y.new.bayes, y.new.wbb2))
(bob.tv  <- distance.tv(y.new.bayes, y.new.bob)) 
(nuts.tv <- distance.tv(y.new.bayes, y.new.nuts)) 
(advi.tv <- distance.tv(y.new.bayes, y.new.advi)) 

# KS distances
(wbb1.ks <- distance.ks(y.new.bayes, y.new.wbb1))
(wbb2.ks <- distance.ks(y.new.bayes, y.new.wbb2))
(bob.ks  <- distance.ks(y.new.bayes, y.new.bob))
(nuts.ks <- distance.ks(y.new.bayes, y.new.nuts))
(advi.ks <- distance.ks(y.new.bayes, y.new.advi))

# -----------------------------------------
# Effective Sample size
# -----------------------------------------

(wbb1.ess <- median(apply(wbb1.draws, 2, coda::effectiveSize)))
(wbb2.ess <- median(apply(wbb2.draws, 2, coda::effectiveSize)))
(bob.ess  <- median(apply(bob.draws,  2, coda::effectiveSize)))
(nuts.ess <- median(apply(nuts.draws, 2, coda::effectiveSize)))
(advi.ess <- median(apply(advi.draws, 2, coda::effectiveSize)))

# -----------------------------------------
# Density plots
# -----------------------------------------

densPlot(as.vector(nuts.draws[,seq(from = 5, by = d, length.out = K)]),
         as.vector(nuts.draws[,seq(from = 7, by = d, length.out = K)]),
         "NUTS", "Width", "Area",
         colPlot = "Blues", xlim = c(-1.3, 1.4),  ylim = c(-1.3, 1.4))
addContours(as.vector(bayes.draws[,seq(from = 5, by = d, length.out = K)]),
            as.vector(bayes.draws[,seq(from = 7, by = d, length.out = K)]))
legend("bottomleft", "True posterior", col = 2, lwd = 2, bty = "n")


densPlot(as.vector(wbb1.draws[,seq(from = 5, by = d, length.out = K)]),
         as.vector(wbb1.draws[,seq(from = 7, by = d, length.out = K)]),
         "WBB1", "Width", "Area",
         colPlot = "Blues", xlim = c(-1.3, 1.4),  ylim = c(-1.3, 1.4))
addContours(as.vector(bayes.draws[,seq(from = 5, by = d, length.out = K)]),
            as.vector(bayes.draws[,seq(from = 7, by = d, length.out = K)]))
legend("bottomleft", "True posterior", col = 2, lwd = 2, bty = "n")


densPlot(as.vector(wbb2.draws[,seq(from = 5, by = d, length.out = K)]),
         as.vector(wbb2.draws[,seq(from = 7, by = d, length.out = K)]),
         "WBB2", "Width", "Area",
         colPlot = "Blues", xlim = c(-1.3, 1.4),  ylim = c(-1.3, 1.4))
addContours(as.vector(bayes.draws[,seq(from = 5, by = d, length.out = K)]),
            as.vector(bayes.draws[,seq(from = 7, by = d, length.out = K)]))
legend("bottomleft", "True posterior", col = 2, lwd = 2, bty = "n")


densPlot(as.vector(bob.draws[,seq(from = 5, by = d, length.out = K)]),
         as.vector(bob.draws[,seq(from = 7, by = d, length.out = K)]),
         "BOB", "Width", "Area",
         colPlot = "Blues", xlim = c(-1.3, 1.4),  ylim = c(-1.3, 1.4))
addContours(as.vector(bayes.draws[,seq(from = 5, by = d, length.out = K)]),
            as.vector(bayes.draws[,seq(from = 7, by = d, length.out = K)]))
legend("bottomleft", "True posterior", col = 2, lwd = 2, bty = "n")


densPlot(as.vector(advi.draws[,seq(from = 5, by = d, length.out = K)]),
         as.vector(advi.draws[,seq(from = 7, by = d, length.out = K)]),
         "ADVI", "Bottom", "Diagonal",
         colPlot = "Blues", xlim = c(-1.3, 1.4),  ylim = c(-1.3, 1.4))
addContours(as.vector(bayes.draws[,seq(from = 5, by = d, length.out = K)]),
            as.vector(bayes.draws[,seq(from = 7, by = d, length.out = K)]))
legend("bottomleft", "True posterior", col = 2, lwd = 2, bty = "n")


# Posterior Predictives

densPlot(y.new.wbb1[,5], y.new.wbb1[,7], "WBB1", "Width", "Area",
         colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
         x.points = y.test[,5], y.points = y.test[,7], z.clust = labels.test)
addContours(y.new.bayes[,5], y.new.bayes[,7])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot(y.new.wbb2[,5], y.new.wbb2[,7], "WBB2", "Width", "Area",
         colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
         x.points = y.test[,5], y.points = y.test[,7], z.clust = labels.test)
addContours(y.new.bayes[,5], y.new.bayes[,7])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot(y.new.bob[,5], y.new.bob[,7], "BOB", "Width", "Area",
         colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
         x.points = y.test[,5], y.points = y.test[,7], z.clust = labels.test)
addContours(y.new.bayes[,5], y.new.bayes[,7])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot(y.new.nuts[,5], y.new.nuts[,7], "NUTS", "Width", "Area",
         colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
         x.points = y.test[,5], y.points = y.test[,7], z.clust = labels.test)
addContours(y.new.bayes[,5], y.new.bayes[,7])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot(y.new.advi[,5], y.new.advi[,7], "ADVI", "Width", "Area",
         colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
         x.points = y.test[,5], y.points = y.test[,7], z.clust = labels.test)
addContours(y.new.bayes[,5], y.new.bayes[,7])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

bayesDens(y.new.bayes[,5], y.new.bayes[,7], 
          xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1),
          main = "Bayesian Posterior Predictive", 
          xlab = "Width", ylab = "Area")




pairs(y)


densPlot(y.new.wbb1[,5], y.new.wbb1[,4], "WBB1", "Width", "Area",
         colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
         x.points = y.test[,5], y.points = y.test[,4], z.clust = labels.test)
addContours(y.new.bayes[,5], y.new.bayes[,4])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot(y.new.bob[,5], y.new.bob[,4], "BOB", "Width", "Area",
         colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
         x.points = y.test[,5], y.points = y.test[,4], z.clust = labels.test)
addContours(y.new.bayes[,5], y.new.bayes[,4])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

