# ------------------------------------------------------------------------------
# Wine data-set
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

load("data-analyses/wine/wine.Rdata")

sum(wines[,1] == "Barolo")           # 59 Barolo type
sum(wines[,1] == "Grignolino")       # 71 Grignolino type
sum(wines[,1] == "Barbera")          # 48 Barbera type

# Prepare and standardize the data ---------------------------------------------
labels <- wines[,1]
y <- base::scale(wines[,-1], center = TRUE, scale = FALSE)
attr(y, "scaled:center") <- NULL; attr(y, "scaled:scale")  <- NULL
K <- length(levels(labels))
d <- ncol(y)
Z <- lapply(1:K, \(k) ifelse(as.numeric(labels) == k, 1, 0)) |> 
  unlist(x = _) |> matrix(data = _, ncol = K, byrow = FALSE) 

sum(Z[,1])                            # 59 Barolo type - train
sum(Z[,2])                            # 71 Grignolino type - train
sum(Z[,3])                            # 48 Barbera type - train

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
mean(wbb1.out) # 123.7781

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
mean(wbb2.out) # 114.2872

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
mean((bob.out$post.draws)) # 141.3791
bob.out$x.optim # optimal x value - 0.42001

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
mean((advi.out)) # 628.9826

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
mean((nuts.out)) # 187.7145

# BAYES ------------------------------------------------------------------------

bayes.true <- bayesian_posterior(y, Z, betas, lambdas, nus, psis, alphas)
mean(bayes.true) # 198.8007

# Posterior predictive distributions -------------------------------------------

wbb1.pred <- post.predictive(wbb1.out, d = d, K = K)
wbb2.pred <- post.predictive(wbb2.out, d = d, K = K)
bob.pred <- post.predictive(bob.out$post.draws, d = d, K = K)
advi.pred <- post.predictive(advi.out, d = d, K = K)
nuts.pred <- post.predictive(nuts.out, d = d, K = K)
bayes.pred <- post.predictive(bayes.true, d = d, K = K)

# Export the results -----------------------------------------------------------

save(wbb1.out, file = "data-analyses/wine/wbb1_out.Rdata")
save(wbb2.out, file = "data-analyses/wine/wbb2_out.Rdata")
save(bob.out, file = "data-analyses/wine/bob_out.Rdata")
save(advi.out, file = "data-analyses/wine/advi_out.Rdata")
save(nuts.out, file = "data-analyses/wine/nuts_out.Rdata")
save(bayes.true, file = "data-analyses/wine/bayes_true.Rdata")

save(wbb1.pred, file = "data-analyses/wine/wbb1_pred.Rdata")
save(wbb2.pred, file = "data-analyses/wine/wbb2_pred.Rdata")
save(bob.pred, file = "data-analyses/wine/bob_pred.Rdata")
save(advi.pred, file = "data-analyses/wine/advi_pred.Rdata")
save(nuts.pred, file = "data-analyses/wine/nuts_pred.Rdata")
save(bayes.pred, file = "data-analyses/wine/bayes_pred.Rdata")

# Upload (back) the results ----------------------------------------------------

load("data-analyses/wine/wbb1_out.Rdata"); mean(wbb1.out)            # 123.7781
load("data-analyses/wine/wbb2_out.Rdata"); mean(wbb2.out)            # 114.2872
load("data-analyses/wine/bob_out.Rdata"); mean((bob.out$post.draws)) # 141.3791
load("data-analyses/wine/advi_out.Rdata"); mean((advi.out))          # 628.9826
load("data-analyses/wine/nuts_out.Rdata"); mean((nuts.out))          # 187.7145
load("data-analyses/wine/bayes_true.Rdata"); mean(bayes.true)        # 166.0071

load("data-analyses/wine/wbb1_pred.Rdata")
load("data-analyses/wine/wbb2_pred.Rdata")
load("data-analyses/wine/bob_pred.Rdata")
load("data-analyses/wine/advi_pred.Rdata")
load("data-analyses/wine/nuts_pred.Rdata")
load("data-analyses/wine/bayes_pred.Rdata")

# Distances --------------------------------------------------------------------

set.seed(1)
round(T4transport::swdist(bayes.pred, wbb1.pred)$distance, 3)     # 4.694
round(T4transport::swdist(bayes.pred, wbb2.pred)$distance, 3)     # 5.957
round(T4transport::swdist(bayes.pred, bob.pred)$distance, 3)      # 3.201
round(T4transport::swdist(bayes.pred, advi.pred)$distance, 3)     # 18.906
round(T4transport::swdist(bayes.pred, nuts.pred)$distance, 3)     # 16.65

# Running times ----------------------------------------------------------------

round(bob.time, 3)     # 31.509
round(wbb1.time, 3)    # 1.283
round(wbb2.time, 3)    # 0.571
round(nuts.time, 3)    # 7501.166
round(advi.time, 3)    # 21.307


# Post. predictive plots -------------------------------------------------------
# Size PDF plots 6 x 6

# Main plots -------------------------------------------------------------------
var1 <- 1
var2 <- 12
colnames(y)[var1]
colnames(y)[var2]
xlab = paste(colnames(y)[var1], "hol", sep = "")
ylab = tolower(colnames(y)[var2])
xlim <- range(y[,var1]) + c(-0.5, 0.5)
ylim <- range(y[,var2]) + c(-0.5, 0.5)


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
var1 <- 6
var2 <- 10
colnames(y)[var1]
colnames(y)[var2]
xlab = paste(colnames(y)[var1], "ols", sep = "")
ylab = paste(colnames(y)[var2], "or", sep = "")
xlim <- range(y[,var1]) + c(-0.5, 0.5)
ylim <- range(y[,var2]) + c(-0.5, 0.5)


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
var1 <- 7
var2 <- 13
colnames(y)[7]
colnames(y)[13]
xlim <- range(y[,var1]) + c(-1.0, 0.0)
ylim <- range(y[,var2]) + c(-1.0, 1.0)
ylim <- c(-550, ylim[2])
xlab <- "flavanoids"
ylab <- "proline"

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
  #xlim = c(-2.5, 2.5),
  #ylim = c(-4.5, 7.5),
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
  #xlim = c(-2.5, 2.5),
  #ylim = c(-4.5, 7.5),
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
  #xlim = c(-2.5, 2.5),
  #ylim = c(-4.5, 7.5),
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
  #xlim = c(-2.5, 2.5),
  #ylim = c(-4.5, 7.5),
  xlim = xlim,
  ylim = ylim,
  x.points = y[,var1], y.points = y[,var2], z.clust = labels
)
add_contours(bayes.pred[,var1], bayes.pred[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")
