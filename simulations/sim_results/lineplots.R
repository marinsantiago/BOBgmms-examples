# ------------------------------------------#
#    Simulation plots varying sample size   #
# ------------------------------------------#

load("simulations/sim-varying-n/sim_n1/bobSim_n1.Rdata")
load("simulations/sim-varying-n/sim_n1/wbb1Sim_n1.Rdata")
load("simulations/sim-varying-n/sim_n1/wbb2Sim_n1.Rdata")

load("simulations/sim-varying-n/sim_n2/bobSim_n2.Rdata")
load("simulations/sim-varying-n/sim_n2/wbb1Sim_n2.Rdata")
load("simulations/sim-varying-n/sim_n2/wbb2Sim_n2.Rdata")

load("simulations/sim-varying-n/sim_n3/bobSim_n3.Rdata")
load("simulations/sim-varying-n/sim_n3/wbb1Sim_n3.Rdata")
load("simulations/sim-varying-n/sim_n3/wbb2Sim_n3.Rdata")

load("simulations/sim-varying-n/sim_n4/bobSim_n4.Rdata")
load("simulations/sim-varying-n/sim_n4/wbb1Sim_n4.Rdata")
load("simulations/sim-varying-n/sim_n4/wbb2Sim_n4.Rdata")

load("simulations/sim-varying-n/sim_n5/bobSim_n5.Rdata")
load("simulations/sim-varying-n/sim_n5/wbb1Sim_n5.Rdata")
load("simulations/sim-varying-n/sim_n5/wbb2Sim_n5.Rdata")

bob.results  <- list(bobSim.n1, bobSim.n2, bobSim.n3, bobSim.n4, bobSim.n5)
wbb1.results <- list(wbb1Sim.n1, wbb1Sim.n2, wbb1Sim.n3, wbb1Sim.n4, wbb1Sim.n5)
wbb2.results <- list(wbb2Sim.n1, wbb2Sim.n2, wbb2Sim.n3, wbb2Sim.n4, wbb2Sim.n5)

TV.bob  <- lapply(bob.results,  \(r) sapply(r, \(b) b$TV.bob))
TV.wbb1 <- lapply(wbb1.results, \(r) sapply(r, \(b) b$TV.wbb1))
TV.wbb2 <- lapply(wbb2.results, \(r) sapply(r, \(b) b$TV.wbb2))

KS.bob  <- lapply(bob.results,  \(r) sapply(r, \(b) b$KS.bob))
KS.wbb1 <- lapply(wbb1.results, \(r) sapply(r, \(b) b$KS.wbb1))
KS.wbb2 <- lapply(wbb2.results, \(r) sapply(r, \(b) b$KS.wbb2))

# Compute the medians
TV.bob.medians  <- sapply(TV.bob, median)
TV.wbb1.medians <- sapply(TV.wbb1, median)
TV.wbb2.medians <- sapply(TV.wbb2, median)

KS.bob.medians  <- sapply(KS.bob, median)
KS.wbb1.medians <- sapply(KS.wbb1, median)
KS.wbb2.medians <- sapply(KS.wbb2, median)

# Compute the upper quantiles
TV.bob.uq  <- sapply(TV.bob,  \(d) quantile(d, 0.75))
TV.wbb1.uq <- sapply(TV.wbb1, \(d) quantile(d, 0.75))
TV.wbb2.uq <- sapply(TV.wbb2, \(d) quantile(d, 0.75))

KS.bob.uq  <- sapply(KS.bob,  \(d) quantile(d, 0.75))
KS.wbb1.uq <- sapply(KS.wbb1, \(d) quantile(d, 0.75))
KS.wbb2.uq <- sapply(KS.wbb2, \(d) quantile(d, 0.75))

# Compute the lower quantiles
TV.bob.lq  <- sapply(TV.bob,  \(d) quantile(d, 0.25))
TV.wbb1.lq <- sapply(TV.wbb1, \(d) quantile(d, 0.25))
TV.wbb2.lq <- sapply(TV.wbb2, \(d) quantile(d, 0.25))

KS.bob.lq  <- sapply(KS.bob,  \(d) quantile(d, 0.25))
KS.wbb1.lq <- sapply(KS.wbb1, \(d) quantile(d, 0.25))
KS.wbb2.lq <- sapply(KS.wbb2, \(d) quantile(d, 0.25))

n_sizes <- c(50, 125, 250, 375, 500)

lwd <- 1.7
par(mfrow = c(1, 2), mar = c(4, 4.6, 4, 2) + 0.1)

# TV plot
plot(n_sizes, TV.bob.medians, type = "l", ylim = c(0, 0.1), col = "tan2",
     lwd = lwd, ylab = as.expression(bquote("" ~ hat(TV))), lty = 1,
     xlab = "n", main = as.expression(bquote("" ~ hat(TV) ~ "Distance")),
     cex = 1, cex.main = 1.5, cex.axis = 1.3, cex.lab = 1.3)
lines(n_sizes, TV.wbb1.medians, type = "l", col = "#3D5A80", lwd = lwd, lty = 2)
lines(n_sizes, TV.wbb2.medians, type = "l", col = "#9B1D20", lwd = lwd, lty = 3)
segments(x0 = n_sizes,
         y0 = TV.wbb2.lq,
         x1 = n_sizes, 
         y1 = TV.wbb2.uq,
         col = "#9B1D20", lwd = 2)
segments(x0 = n_sizes,
         y0 = TV.wbb1.lq,
         x1 = n_sizes, 
         y1 = TV.wbb1.uq,
         col = "#3D5A80", lwd = 2)
segments(x0 = n_sizes,
         y0 = TV.bob.lq,
         x1 = n_sizes, 
         y1 = TV.bob.uq,
         col = "tan2", lwd = 2)
legend(xy.coords(355, 0.105), c("BOB", "WBB1", "WBB2"), lty = c(1, 2, 3),
       col = c("tan2", "#3D5A80", "#9B1D20"), lwd = 1.7, bty = "n",
       cex = 1, pt.cex = 0.6, x.intersp = 0.5)


# KS plot
plot(n_sizes, KS.bob.medians, type = "l", ylim = c(0, 0.1), col = "tan2",
     lwd = lwd, ylab = as.expression(bquote("" ~ hat(KS))),
     xlab = "n", main = as.expression(bquote("" ~ hat(KS) ~ "Distance")),
     cex = 1, cex.main = 1.5, cex.axis = 1.3, cex.lab = 1.3)
lines(n_sizes, KS.wbb1.medians, type = "l", col = "#3D5A80", lwd = lwd, lty = 2)
lines(n_sizes, KS.wbb2.medians, type = "l", col = "#9B1D20", lwd = lwd, lty = 3)
segments(x0 = n_sizes,
         y0 = KS.wbb2.lq,
         x1 = n_sizes, 
         y1 = KS.wbb2.uq,
         col = "#9B1D20", lwd = lwd)
segments(x0 = n_sizes,
         y0 = KS.wbb1.lq,
         x1 = n_sizes, 
         y1 = KS.wbb1.uq,
         col = "#3D5A80", lwd = lwd)
segments(x0 = n_sizes,
         y0 = KS.bob.lq,
         x1 = n_sizes, 
         y1 = KS.bob.uq,
         col = "tan2", lwd = lwd)
legend(xy.coords(350, 0.105), c("BOB", "WBB1", "WBB2"), lty = c(1, 2, 3),
       col = c("tan2", "#3D5A80", "#9B1D20"), lwd = 1.7, bty = "n",
       cex = 1, pt.cex = 0.6, x.intersp = 0.5)

# Size of plot:
# 5.55 X 11.8 inches - Landscape mode. 
