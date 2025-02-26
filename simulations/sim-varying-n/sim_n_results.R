# ------------------------------------------------------------------------------
# Simulations for a varying sample size (n = 50, 125, 250, 375, 500)
# ------------------------------------------------------------------------------

# Load the results -------------------------------------------------------------

load("simulations/sim-varying-n/results_n50.Rdata")
load("simulations/sim-varying-n/results_n125.Rdata")
load("simulations/sim-varying-n/results_n250.Rdata")
load("simulations/sim-varying-n/results_n375.Rdata")
load("simulations/sim-varying-n/results_n500.Rdata")

results <- list(
  results_n50, results_n125, results_n250, results_n375, results_n500
)

# WBB1 -------------------------------------------------------------------------
results.wbb1 <- lapply(results, \(r) unlist(r$dist_wbb1))
# Medians
SW2.wbb1.median <- sapply(results.wbb1, median)
# Upper quantiles
SW2.wbb1.uq <- sapply(results.wbb1, \(r) quantile(r, 0.75))
# Lower quantiles
SW2.wbb1.lq <- sapply(results.wbb1, \(r) quantile(r, 0.25))

# WBB2 -------------------------------------------------------------------------
results.wbb2 <- lapply(results, \(r) unlist(r$dist_wbb2))
# Medians
SW2.wbb2.median <- sapply(results.wbb2, median)
# Upper quantiles
SW2.wbb2.uq <- sapply(results.wbb2, \(r) quantile(r, 0.75))
# Lower quantiles
SW2.wbb2.lq <- sapply(results.wbb2, \(r) quantile(r, 0.25))

# BOB --------------------------------------------------------------------------
results.bob <- lapply(results, \(r) unlist(r$dist_bob))
# Medians
SW2.bob.median <- sapply(results.bob, median)
# Upper quantiles
SW2.bob.uq <- sapply(results.bob, \(r) quantile(r, 0.75))
# Lower quantiles
SW2.bob.lq <- sapply(results.bob, \(r) quantile(r, 0.25))

# Line-plots -------------------------------------------------------------------
par(mar = c(4, 4.6, 4, 2) + 0.1)
n_sizes <- c(50, 125, 250, 375, 500)
lwd <- 1.7
SW2.dist <- as.expression(bquote("" ~ SW[2] ~ "  Distance"))
plot(
  x = n_sizes, y = SW2.bob.median, type = "l", ylim = c(0, 0.35), col = "tan2",
  lwd = lwd, xlab = "n", ylab = SW2.dist, main = SW2.dist, lty = 1,
  cex = 1, cex.main = 1.5, cex.axis = 1.3, cex.lab = 1.3
)
lines(
  x = n_sizes, y = SW2.wbb1.median, 
  type = "l", col = "#3D5A80", lwd = lwd, lty = 2
)
lines(
  x = n_sizes, y = SW2.wbb2.median,
  type = "l", col = "#9B1D20", lwd = lwd, lty = 3
)
segments(
  x0 = n_sizes, y0 = SW2.wbb2.lq,
  x1 = n_sizes, y1 = SW2.wbb2.uq, col = "#9B1D20", lwd = 2
)
segments(
  x0 = n_sizes, y0 = SW2.wbb1.lq,
  x1 = n_sizes, y1 = SW2.wbb1.uq, col = "#3D5A80", lwd = 2
)
segments(
  x0 = n_sizes, y0 = SW2.bob.lq,
  x1 = n_sizes, y1 = SW2.bob.uq, col = "tan2", lwd = 2
)
legend(
  x = "topright", legend = c("BOB", "WBB1", "WBB2"), lty = c(1, 2, 3),
  col = c("tan2", "#3D5A80", "#9B1D20"), lwd = 1.7, bty = "n",
  cex = 1, pt.cex = 0.6, x.intersp = 0.5
)

# Size of plot:
# 9.14 X 6.26 inches - Landscape mode. 