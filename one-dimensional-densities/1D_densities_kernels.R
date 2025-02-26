# ------------------------------------------------------------------------------
# One-dimensional density plots: Kernels data-set
# ------------------------------------------------------------------------------

folder <- "data-analyses/kernels"
load(paste(folder, "/seeds.Rdata", sep = ""))
load(paste(folder, "/bayes_pred.Rdata", sep = ""))
load(paste(folder, "/bob_pred.Rdata", sep = ""))
load(paste(folder, "/wbb1_pred.Rdata", sep = ""))
load(paste(folder, "/wbb2_pred.Rdata", sep = ""))
load(paste(folder, "/nuts_pred.Rdata", sep = ""))
load(paste(folder, "/advi_pred.Rdata", sep = ""))

var.names <- colnames(seeds[,-8])
d <- length(var.names)


ylims <- list(
  c(0, 0.20), c(0, 0.45), c(0, 5.50), c(0, 1.10),
  c(0, 1.10), c(0, 0.35), c(0, 1.10)
)

par(mfrow = c(3, 3))
for (j in seq_len(d)) {
  variable <- j
  plot(
    x = density(bayes.pred[,variable]),
    ylim = ylims[[variable]], type = "n", ylab = "", 
    main = var.names[variable], xlab = var.names[variable]
  )
  lines(density(advi.pred[,variable]), col = "grey", lwd = 2, lty = 6)
  lines(density(nuts.pred[,variable]), col = "aquamarine3", lwd = 2, lty = 5)
  lines(density(wbb2.pred[,variable]), col = "magenta", lwd = 2, lty = 4)
  lines(density(wbb1.pred[,variable]), col = "goldenrod3", lwd = 2, lty = 3)
  lines(density(bayes.pred[,variable]), lwd = 2, lty = 1, col = 2)
  lines(density(bob.pred[,variable]), col = 4, lwd = 2, lty = 2)
}

plot.new()
legend(
  x = "center", legend = c("Bayes", "BOB", "WBB1", "WBB2", "NUTS", "ADVI"), 
  xpd = TRUE, horiz = FALSE, bty = "n", lty = c(1,2,3, 4, 5, 6), 
  col = c(2, 4, "goldenrod3", "magenta", "aquamarine3", "grey"), 
  lwd = 2, cex = 1.5, 
  x.intersp = 1, y.intersp = 0.5
)

# Size of plot:
# 13 X 11.8 inches - Portrait mode. 
