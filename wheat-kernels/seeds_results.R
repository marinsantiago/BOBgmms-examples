library(BOBgmms)
source("./R/densPlots.R")
shhh(library(coda, include.only = "effectiveSize"))
rm(shhh)

load("wheat-kernels/y_new_bayes.Rdata")
load("wheat-kernels/y_new_wbb1.Rdata")
load("wheat-kernels/y_new_wbb2.Rdata")
load("wheat-kernels/y_new_advi.Rdata")
load("wheat-kernels/y_new_nuts.Rdata")
load("wheat-kernels/y_new_bob.Rdata")

#load("wheat-kernels/bayes_draws.Rdata")
#load("wheat-kernels/wbb1_draws.Rdata")
#load("wheat-kernels/wbb2_draws.Rdata")
#load("wheat-kernels/advi_draws.Rdata")
#load("wheat-kernels/nuts_draws.Rdata")
#load("wheat-kernels/bob_draws.Rdata")

load("wheat-kernels/advi_time.Rdata")
load("wheat-kernels/nuts_time.Rdata")
load("wheat-kernels/wbb1_time.Rdata")
load("wheat-kernels/wbb2_time.Rdata")
load("wheat-kernels/bob_time.Rdata")

load("wheat-kernels/y_data_test.Rdata")
load("wheat-kernels/labels_test.Rdata")

# TV distances
(wbb1.tv <- distance.tv(y.new.bayes, y.new.wbb1))
(wbb2.tv <- distance.tv(y.new.bayes, y.new.wbb2))
(bob.tv  <- distance.tv(y.new.bayes, y.new.bob)) 
(nuts.tv <- distance.tv(y.new.bayes, y.new.nuts)) 
(advi.tv <- distance.tv(y.new.bayes, y.new.advi))

round(bob.tv, 3)
round(wbb1.tv, 3)
round(wbb2.tv, 3)
round(nuts.tv, 3)
round(advi.tv, 3)

# KS distances
(wbb1.ks <- distance.ks(y.new.bayes, y.new.wbb1))
(wbb2.ks <- distance.ks(y.new.bayes, y.new.wbb2))
(bob.ks  <- distance.ks(y.new.bayes, y.new.bob))
(nuts.ks <- distance.ks(y.new.bayes, y.new.nuts))
(advi.ks <- distance.ks(y.new.bayes, y.new.advi))

round(bob.ks, 3)
round(wbb1.ks, 3)
round(wbb2.ks, 3)
round(nuts.ks, 3)
round(advi.ks, 3)

# Median Effective Sample size
#(wbb1.ess <- median(apply(wbb1.draws, 2, coda::effectiveSize)))
#(wbb2.ess <- median(apply(wbb2.draws, 2, coda::effectiveSize)))
#(bob.ess  <- median(apply(bob.draws,  2, coda::effectiveSize)))
#(nuts.ess <- median(apply(nuts.draws, 2, coda::effectiveSize)))
#(advi.ess <- median(apply(advi.draws, 2, coda::effectiveSize)))

# Elapsed
round(bob.time, 3)
round(wbb1.time, 3)
round(wbb2.time, 3)
round(nuts.time, 3)
round(advi.time, 3)

# -----------------------------------------
# Posterior Predictive plots
# -----------------------------------------

var1 <- 5
var2 <- 7

colnames(y.test)[var1]
colnames(y.test)[var2]

densPlot.reduced(y.new.wbb1[,var1], y.new.wbb1[,var2], "WBB1", "width", "k.groove.len",
    colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.3, 2.3), 
      x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.wbb2[,var1], y.new.wbb2[,var2], "WBB2", "width", "k.groove.len",
    colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.3, 2.3), 
      x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.bob[,var1], y.new.bob[,var2], "BOB", "width", "k.groove.len",
    colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.3, 2.3),
      x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.nuts[,var1], y.new.nuts[,var2], "NUTS", "width", "k.groove.len",
    colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.3, 2.3),
      x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.advi[,var1], y.new.advi[,var2], "ADVI", "width", "k.groove.len",
    colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.3, 2.3),
      x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


bayesDens(y.new.bayes[,var1], y.new.bayes[,var2], 
          xlim = c(-2.1, 2.1),  ylim = c(-2.3, 2.3),
          main = "Bayes", xlab = "width", ylab = "k.groove.len")


# ----------------------------------------------
# Supplementary Posterior predictive plots
# ----------------------------------------------

var1 <- 2
var2 <- 3

colnames(y.test)[var1]
colnames(y.test)[var2]


densPlot.reduced(y.new.wbb1[,var1], y.new.wbb1[,var2], "WBB1", "perimeter", "compactness",
                 colPlot = "Blues", xlim = c(-2.1, 2.3),  ylim = c(-2.8, 2.3), 
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.wbb2[,var1], y.new.wbb2[,var2], "WBB2", "perimeter", "compactness",
                 colPlot = "Blues", xlim = c(-2.1, 2.3),  ylim = c(-2.8, 2.3),
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.bob[,var1], y.new.bob[,var2], "BOB", "perimeter", "compactness",
                 colPlot = "Blues", xlim = c(-2.1, 2.3),  ylim = c(-2.8, 2.3),
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.nuts[,var1], y.new.nuts[,var2], "NUTS", "perimeter", "compactness",
                 colPlot = "Blues", xlim = c(-2.1, 2.3),  ylim = c(-2.8, 2.3),
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.advi[,var1], y.new.advi[,var2], "ADVI", "perimeter", "compactness",
                 colPlot = "Blues", xlim = c(-2.1, 2.3),  ylim = c(-2.8, 2.3),
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


bayesDens(y.new.bayes[,var1], y.new.bayes[,var2], 
          xlim = c(-2.1, 2.3),  ylim = c(-2.8, 2.3),
          main = "Bayes", xlab = "perimeter", ylab = "compactness")
