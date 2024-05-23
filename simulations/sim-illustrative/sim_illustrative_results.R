library(BOBgmms)
source("./R/densPlots.R")

load("simulations/sim-illustrative/y_new_bayes.Rdata")
load("simulations/sim-illustrative/y_new_wbb1.Rdata")
load("simulations/sim-illustrative/y_new_wbb2.Rdata")
load("simulations/sim-illustrative/y_new_advi.Rdata")
load("simulations/sim-illustrative/y_new_nuts.Rdata")
load("simulations/sim-illustrative/y_new_bob.Rdata")

load("simulations/sim-illustrative/advi_time.Rdata")
load("simulations/sim-illustrative/nuts_time.Rdata")
load("simulations/sim-illustrative/wbb1_time.Rdata")
load("simulations/sim-illustrative/wbb2_time.Rdata")
load("simulations/sim-illustrative/bob_time.Rdata")

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
var2 <- 2

label.x <- as.expression(bquote(y [.(var1)]))
label.y <- as.expression(bquote(y [.(var2)]))

densPlot.reduced(y.new.wbb1[,1], y.new.wbb1[,2], "WBB1", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,1], y.new.bayes[,2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

densPlot.reduced(y.new.wbb2[,1], y.new.wbb2[,2], "WBB2", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,1], y.new.bayes[,2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

densPlot.reduced(y.new.bob[,1], y.new.bob[,2], "BOB", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,1], y.new.bayes[,2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

densPlot.reduced(y.new.nuts[,1], y.new.nuts[,2], "NUTS", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,1], y.new.bayes[,2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

densPlot.reduced(y.new.advi[,1], y.new.advi[,2], "ADVI", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,1], y.new.bayes[,2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

bayesDens(y.new.bayes[,1], y.new.bayes[,2], 
          xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
          main = "Bayes", xlab = label.x, ylab = label.y)

# Size PDF plots 6 x 6. 

# ----------------------------------------------
# Supplementary Posterior predictive plots
# ----------------------------------------------

var1 <- 3
var2 <- 4


label.x <- as.expression(bquote(y [.(var1)]))
label.y <- as.expression(bquote(y [.(var2)]))

densPlot.reduced(y.new.wbb1[,var1], y.new.wbb1[,var2], "WBB1", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

densPlot.reduced(y.new.wbb2[,var1], y.new.wbb2[,var2], "WBB2", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

densPlot.reduced(y.new.bob[,var1], y.new.bob[,var2], "BOB", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

densPlot.reduced(y.new.nuts[,var1], y.new.nuts[,var2], "NUTS", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

densPlot.reduced(y.new.advi[,var1], y.new.advi[,var2], "ADVI", label.x, label.y,
                 colPlot = "Blues", xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1))
addContours(y.new.bayes[,var1], y.new.bayes[,var2])
legend("topleft", "True posterior pred.", col = 2, lwd = 2, bty = "n")

bayesDens(y.new.bayes[,var1], y.new.bayes[,var2], 
          xlim = c(-2.1, 2.1),  ylim = c(-2.1, 2.1), 
          main = "Bayes", xlab = label.x, ylab = label.y)
