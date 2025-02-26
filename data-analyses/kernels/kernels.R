# ------------------------------------------------------------------------------
# Kernels data-set
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

load("data-analyses/kernels/seeds.Rdata")

sum(seeds[,8] == "Kama")             # 70 K variety
sum(seeds[,8] == "Rosa")             # 70 R variety
sum(seeds[,8] == "Canadian")         # 70 C variety

# Prepare and standardize the data ---------------------------------------------
labels <- seeds[,8]
y <- base::scale(seeds[,-8], center = TRUE, scale = FALSE)
attr(y, "scaled:center") <- NULL; attr(y, "scaled:scale")  <- NULL
K <- length(levels(labels))
d <- ncol(y)
Z <- lapply(1:K, \(k) ifelse(as.numeric(labels) == k, 1, 0)) |> 
  unlist(x = _) |> matrix(data = _, ncol = K, byrow = FALSE)

sum(Z[,1])                            # 59 K type - train
sum(Z[,2])                            # 71 R type - train
sum(Z[,3])                            # 48 C type - train

# Set prior and tuning hyper-parameters ----------------------------------------
alphas <- rep(1.1, K)
betas <- lapply(1:K, function(x) rep(0, d))
psis <- lapply(1:K, function(x) diag(d))

set.seed(1)
shrink.params <- shrinkage.cv.tuning(y, betas, psis, alphas)

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

# Weighted Bayesian Bootstrap 1 ------------------------------------------------

set.seed(1)
start.wbb1 <- Sys.time()
wbb1.out <- BOBgmms::wbb.gmm(
  y = y, means.init = means.init, covs.init = covs.init,
  probs.init = probs.init, betas = betas, lambdas = lambdas,
  nus = nus, psis = psis, alphas = alphas, wbb.scheme = "wbb1",
  max.iters = 20000, a = a, b = b, c = c, r = r
)
end.wbb1 <- Sys.time()
(wbb1.time <- difftime(end.wbb1, start.wbb1, units = "secs"))
mean(wbb1.out) # 0.0988211

# Weighted Bayesian Bootstrap 2 ------------------------------------------------

set.seed(1)
start.wbb2 <- Sys.time()
wbb2.out <- BOBgmms::wbb.gmm(
  y = y, means.init = means.init, covs.init = covs.init,
  probs.init = probs.init, betas = betas, lambdas = lambdas,
  nus = nus, psis = psis, alphas = alphas, wbb.scheme = "wbb2",
  max.iters = 20000, a = a, b = b, c = c, r = r
)
end.wbb2 <- Sys.time()
(wbb2.time <- difftime(end.wbb2, start.wbb2, units = "secs"))
mean(wbb2.out) # 0.09651176

# Bayesian Optimized Bootstrap -------------------------------------------------

lower_bound <- 1e-05
upper_bound <- 1.5

set.seed(1)
start.bob <- Sys.time()
bob.out <- BOBgmms::bob.gmm(
  y = y, means.init = means.init, covs.init = covs.init,
  probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
  psis = psis, alphas = alphas, lower_bound = lower_bound,
  upper_bound = upper_bound, max.iters = 20000, size.batch = 1000,
  bo.iters = 100, a = a, b = b, c = c, r = r
)
end.bob <- Sys.time()
(bob.time <- difftime(end.bob, start.bob, units = "secs"))
mean((bob.out$post.draws)) # 0.09774964
bob.out$x.optim # optimal x value - 0.9051748

# ADVI -------------------------------------------------------------------------

set.seed(1)
start.advi <- Sys.time()
advi.out <- advi(
  y = y, means.init = means.init, covs.init = covs.init,
  probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
  psis = psis, alphas = alphas, max.iters = 20000, vb.model = stan.model
)
end.advi <- Sys.time()
(advi.time <- difftime(end.advi, start.advi, units = "secs"))
mean((advi.out)) # 0.2017045

# NUTS -------------------------------------------------------------------------

set.seed(1)
start.nuts <- Sys.time()
nuts.out <- nuts(
  y = y, means.init = means.init, covs.init = covs.init,
  probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
  psis = psis, alphas = alphas, max.iters = 40000, warm_init = TRUE
)
end.nuts <- Sys.time()
(nuts.time <- difftime(end.nuts, start.nuts, units = "secs"))
mean((nuts.out)) # 0.3223127

# BAYES ------------------------------------------------------------------------

bayes.true <- bayesian_posterior(y, Z, betas, lambdas, nus, psis, alphas)
mean(bayes.true) # 0.1158226

# Posterior predictive distributions -------------------------------------------

wbb1.pred <- post.predictive(wbb1.out, d = d, K = K)
wbb2.pred <- post.predictive(wbb2.out, d = d, K = K)
bob.pred <- post.predictive(bob.out$post.draws, d = d, K = K)
advi.pred <- post.predictive(advi.out, d = d, K = K)
nuts.pred <- post.predictive(nuts.out, d = d, K = K)
bayes.pred <- post.predictive(bayes.true, d = d, K = K)

# Export the results -----------------------------------------------------------

save(wbb1.out, file = "data-analyses/kernels/wbb1_out.Rdata")
save(wbb2.out, file = "data-analyses/kernels/wbb2_out.Rdata")
save(bob.out, file = "data-analyses/kernels/bob_out.Rdata")
save(advi.out, file = "data-analyses/kernels/advi_out.Rdata")
save(nuts.out, file = "data-analyses/kernels/nuts_out.Rdata")
save(bayes.true, file = "data-analyses/kernels/bayes_true.Rdata")

save(wbb1.pred, file = "data-analyses/kernels/wbb1_pred.Rdata")
save(wbb2.pred, file = "data-analyses/kernels/wbb2_pred.Rdata")
save(bob.pred, file = "data-analyses/kernels/bob_pred.Rdata")
save(advi.pred, file = "data-analyses/kernels/advi_pred.Rdata")
save(nuts.pred, file = "data-analyses/kernels/nuts_pred.Rdata")
save(bayes.pred, file = "data-analyses/kernels/bayes_pred.Rdata")

# Upload (back) the results ----------------------------------------------------

load("data-analyses/kernels/wbb1_out.Rdata"); mean(wbb1.out) # 0.0988211
load("data-analyses/kernels/wbb2_out.Rdata"); mean(wbb2.out) # 0.09651176
load("data-analyses/kernels/bob_out.Rdata"); mean((bob.out$post.draws)) # 0.0977
load("data-analyses/kernels/advi_out.Rdata"); mean((advi.out)) # 0.2017045
load("data-analyses/kernels/nuts_out.Rdata"); mean((nuts.out)) # 0.3223127
load("data-analyses/kernels/bayes_true.Rdata"); mean(bayes.true) # 0.1158226

load("data-analyses/kernels/wbb1_pred.Rdata")
load("data-analyses/kernels/wbb2_pred.Rdata")
load("data-analyses/kernels/bob_pred.Rdata")
load("data-analyses/kernels/advi_pred.Rdata")
load("data-analyses/kernels/nuts_pred.Rdata")
load("data-analyses/kernels/bayes_pred.Rdata")


# Distances --------------------------------------------------------------------

set.seed(101)
round(T4transport::swdist(bayes.pred, wbb1.pred)$distance, 3)     # 0.071
round(T4transport::swdist(bayes.pred, wbb2.pred)$distance, 3)     # 0.071
round(T4transport::swdist(bayes.pred, bob.pred)$distance, 3)      # 0.07
round(T4transport::swdist(bayes.pred, advi.pred)$distance, 3)     # 0.228
round(T4transport::swdist(bayes.pred, nuts.pred)$distance, 3)     # 0.357

# Running times ----------------------------------------------------------------

round(bob.time, 3)     # 24.681
round(wbb1.time, 3)    # 0.377
round(wbb2.time, 3)    # 0.505
round(nuts.time, 3)    # 2224.107
round(advi.time, 3)    # 19.204

# Post. predictive plots -------------------------------------------------------
# Size PDF plots 6 x 6

# Main plots -------------------------------------------------------------------

var1 <- 2
var2 <- 6
colnames(y)[var1]
colnames(y)[var2]
xlab <- colnames(y)[var1]
ylab <- colnames(y)[var2]
xlim <- range(bob.pred[,var1])
ylim <- range(bob.pred[,var2])


# Bayes
bayes_dens(
  x = bayes.pred[,var1], y = bayes.pred[,var2], col = 2,
  xlim = xlim,
  ylim = ylim,
  main = "Bayes",
  xlab = xlab,
  ylab = ylab
)


# BOB
dens_plot(
  x = bob.pred[,var1], y = bob.pred[,var2], main = "BOB",
  xlab = xlab,
  ylab = ylab,
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# WBB1
dens_plot(
  x = wbb1.pred[,var1], y = wbb1.pred[,var2], main = "WBB1",
  xlab = xlab,
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# WBB2
dens_plot(
  x = wbb2.pred[,var1], y = wbb2.pred[,var2], main = "WBB2",
  xlab = xlab, 
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# ADVI
dens_plot(
  x = advi.pred[,var1], y = advi.pred[,var2], main = "ADVI",
  xlab = xlab, 
  ylab = ylab,
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# NUTS
dens_plot(
  x = nuts.pred[,var1], y = nuts.pred[,var2], main = "NUTS",
  xlab = xlab,
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# Supp plots 1 -----------------------------------------------------------------

var1 <- 5
var2 <- 7
colnames(y)[var1]
colnames(y)[var2]
xlab <- colnames(y)[var1]
ylab <- colnames(y)[var2]
xlim <- range(bob.pred[,var1])
ylim <- range(bob.pred[,var2])


# Bayes
bayes_dens(
  x = bayes.pred[,var1], y = bayes.pred[,var2], col = 2,
  xlim = xlim,
  ylim = ylim,
  main = "Bayes",
  xlab = xlab,
  ylab = ylab
)


# BOB
dens_plot(
  x = bob.pred[,var1], y = bob.pred[,var2], main = "BOB",
  xlab = xlab,
  ylab = ylab,
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# WBB1
dens_plot(
  x = wbb1.pred[,var1], y = wbb1.pred[,var2], main = "WBB1",
  xlab = xlab,
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# WBB2
dens_plot(
  x = wbb2.pred[,var1], y = wbb2.pred[,var2], main = "WBB2",
  xlab = xlab, 
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# ADVI
dens_plot(
  x = advi.pred[,var1], y = advi.pred[,var2], main = "ADVI",
  xlab = xlab, 
  ylab = ylab,
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# NUTS
dens_plot(
  x = nuts.pred[,var1], y = nuts.pred[,var2], main = "NUTS",
  xlab = xlab,
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# Supp plots 2 -----------------------------------------------------------------

var1 <- 1
var2 <- 4
colnames(y)[var1]
colnames(y)[var2]
xlab <- colnames(y)[var1]
ylab <- colnames(y)[var2]
xlim <- range(bob.pred[,var1])
ylim <- range(bob.pred[,var2])


# Bayes
bayes_dens(
  x = bayes.pred[,var1], y = bayes.pred[,var2], col = 2,
  xlim = xlim,
  ylim = ylim,
  main = "Bayes",
  xlab = xlab,
  ylab = ylab
)


# BOB
dens_plot(
  x = bob.pred[,var1], y = bob.pred[,var2], main = "BOB",
  xlab = xlab,
  ylab = ylab,
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# WBB1
dens_plot(
  x = wbb1.pred[,var1], y = wbb1.pred[,var2], main = "WBB1",
  xlab = xlab,
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# WBB2
dens_plot(
  x = wbb2.pred[,var1], y = wbb2.pred[,var2], main = "WBB2",
  xlab = xlab, 
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# ADVI
dens_plot(
  x = advi.pred[,var1], y = advi.pred[,var2], main = "ADVI",
  xlab = xlab, 
  ylab = ylab,
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


# NUTS
dens_plot(
  x = nuts.pred[,var1], y = nuts.pred[,var2], main = "NUTS",
  xlab = xlab,
  ylab = ylab, 
  col.palette = "Blues",
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

