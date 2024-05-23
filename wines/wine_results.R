library(BOBgmms)
source("./R/densPlots.R")
shhh(library(coda, include.only = "effectiveSize"))
rm(shhh)

load("wines/y_new_bayes.Rdata")
load("wines/y_new_wbb1.Rdata")
load("wines/y_new_wbb2.Rdata")
load("wines/y_new_advi.Rdata")
load("wines/y_new_nuts.Rdata")
load("wines/y_new_bob.Rdata")

#load("wines/bayes_draws.Rdata")
#load("wines/wbb1_draws.Rdata")
#load("wines/wbb2_draws.Rdata")
#load("wines/advi_draws.Rdata")
#load("wines/nuts_draws.Rdata")
#load("wines/bob_draws.Rdata")

load("wines/advi_time.Rdata")
load("wines/nuts_time.Rdata")
load("wines/wbb1_time.Rdata")
load("wines/wbb2_time.Rdata")
load("wines/bob_time.Rdata")

load("wines/y_data_test.Rdata")
load("wines/labels_test.Rdata")

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
# Posterior predictive plots
# -----------------------------------------

var1 <- 1
var2 <- 12

colnames(y.test)[var1]
colnames(y.test)[var2]

densPlot.reduced(y.new.wbb1[,var1], y.new.wbb1[,var2], "WBB1", "alcohol", "oddw",
    colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.7, 2.2), 
        x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.wbb2[,var1], y.new.wbb2[,var2], "WBB2", "alcohol", "oddw",
    colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.7, 2.2), 
      x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.bob[,var1], y.new.bob[,var2], "BOB", "alcohol", "oddw",
  colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.7, 2.2), 
    x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.nuts[,var1], y.new.nuts[,var2], "NUTS", "alcohol", "oddw",
    colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.7, 2.2), 
      x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.advi[,var1], y.new.advi[,var2], "ADVI", "alcohol", "oddw",
    colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.7, 2.2), 
      x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


bayesDens(y.new.bayes[,var1], y.new.bayes[,var2], 
          xlim = c(-2.5, 2.5),  ylim = c(-2.7, 2.2),
          main = "Bayes", xlab = "alcohol", ylab = "oddw")


# ----------------------------------------------
# Supplementary Posterior predictive plots
# ----------------------------------------------

var1 <- 6
var2 <- 10

colnames(y.test)[var1]
colnames(y.test)[var2]

densPlot.reduced(y.new.wbb1[,var1], y.new.wbb1[,var2], "WBB1", "phenols", "colour",
                 colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.0, 2.7), 
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.wbb2[,var1], y.new.wbb2[,var2], "WBB2", "phenols", "colour",
                 colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.0, 2.7), 
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.bob[,var1], y.new.bob[,var2], "BOB", "phenols", "colour",
                 colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.0, 2.7),  
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.nuts[,var1], y.new.nuts[,var2], "NUTS", "phenols", "colour",
                 colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.0, 2.7), 
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


densPlot.reduced(y.new.advi[,var1], y.new.advi[,var2], "ADVI", "phenols", "colour",
                 colPlot = "Blues", xlim = c(-2.5, 2.5),  ylim = c(-2.0, 2.7), 
                 x.points = y.test[,var1], y.points = y.test[,var2], z.clust = labels.test)
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("bottomleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")


bayesDens(y.new.bayes[,var1], y.new.bayes[,var2], 
          xlim = c(-2.5, 2.5),  ylim = c(-2.0, 2.7), 
          main = "Bayes", xlab = "phenols", ylab = "colour")

