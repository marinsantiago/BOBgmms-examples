# ---------------------------------------------------------#
#    Boxplots for TV and KS distances across 9 settings    #
# ---------------------------------------------------------#

load("simulations/sim1/resultsSim1.Rdata")
load("simulations/sim2/resultsSim2.Rdata")
load("simulations/sim3/resultsSim3.Rdata")
load("simulations/sim4/resultsSim4.Rdata")
load("simulations/sim5/resultsSim5.Rdata")
load("simulations/sim6/resultsSim6.Rdata")
load("simulations/sim7/resultsSim7.Rdata")
load("simulations/sim8/resultsSim8.Rdata")
load("simulations/sim9/resultsSim9.Rdata")

results <- list(
  resultsSim1, resultsSim2, resultsSim3, resultsSim4, resultsSim5,
  resultsSim6, resultsSim7, resultsSim8, resultsSim9
)

# Create data frames containing the results

create.df <- function(TV_BOB, KS_BOB, TV_WBB1, KS_WBB1, TV_WBB2, KS_WBB2){
  Total_BOB  <- length(TV_BOB) + length(KS_BOB)
  Total_WBB1 <- length(TV_WBB1) + length(KS_WBB1)
  Total_WBB2 <- length(TV_WBB2) + length(KS_WBB2)
  DF <- data.frame(
    x = c(c(TV_BOB,  KS_BOB),
          c(TV_WBB1, KS_WBB1),
          c(TV_WBB2, KS_WBB2)),
    y = rep(c("BOB", "WBB1", "WBB2"), c(Total_BOB, Total_WBB1, Total_WBB2)),
    z = rep(rep(1:2, each = length(TV_BOB)),  3),
    stringsAsFactors = FALSE
  )
  return(DF)
}

DFs <- lapply(results, \(r)
              create.df(r$TVdistances$BOB,  r$KSdistances$BOB,
                        r$TVdistances$WBB1, r$KSdistances$WBB1,
                        r$TVdistances$WBB2, r$KSdistances$WBB2))

# Generate the boxplots

cols <- c("#9B1D20", "#3D5A80")
par(oma = c(7,1,1,1), mfrow = c(3, 3), mar = c(4, 4.25, 2.5, 1))
for (i in 1:length(DFs)){
  boxplot(x ~ z + y, data = DFs[[i]],
          col = cols, at = c(1:2, 5:6, 9:10), xaxt = "n", xlab = "",
          ylab = "", main = paste("Setting", i), pch = 16,
          cex = 1, cex.main = 1.7, cex.axis = 1.3, cex.lab = 1.5)
  axis(side = 1, at = c(1.5, 5.5, 9.5), labels = c("BOB", "WBB1", "WBB2"), cex.axis = 1.8)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c(as.expression(bquote("" ~ hat(TV) ~ "Dist.")),
                           as.expression(bquote("" ~ hat(KS) ~ "Dist."))),
       #fill = cols, 
       xpd = T, horiz = T, seg.len = 1, bty = 'n',
       cex = 1.5, pt.cex = 3.5, x.intersp = 1, pch = 15, col = cols)
# xpd = TRUE makes the legend plot to the figure
mtext("        Distance between the Bayesian and the approximate posterior predictive distributions",
      side = 2, line = -1.8, outer = TRUE, cex = 1.4, las = 0)

# Size of plot:
# 13 X 11.8 inches - Portrait mode. 