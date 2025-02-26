# ------------------------------------------------------------------------------
# Sim setting 1-9: Boxplots
# ------------------------------------------------------------------------------

folder <- "simulations/sim_ar_neg_05"
for (s in seq_len(9)) {
  file <- paste(
    folder,
    paste("/setting", s, sep = ""),
    paste("/sim", s, "_results.Rdata", sep = ""),
    sep = ""
  )
  load(file)
}

results <- list(
  results_sim1, results_sim2, results_sim3,
  results_sim4, results_sim5, results_sim6,
  results_sim7, results_sim8, results_sim9
)

#col.palette <- c("#228B22", "#9B1D20", "#3D5A80")
col.palette <- c("#3D5A80", "#9B1D20", "tan2")
par(mfrow = c(3, 3), oma = c(0, 2, 0, 0))
for (s in seq_len(9)) {
  boxplot(
    x = results[[s]]$SW2.distance[1:3], 
    col = col.palette, pch = 16,
    main = paste("Setting", s),
    cex = 1, cex.main = 1.7, cex.axis = 1.3, cex.lab = 1.5,
    outline = TRUE
  )
}

mtext(
  text = as.expression(
    bquote(
      "" ~ SW[2] ~ " distance between the Bayesian and the approximate posterior predictive distributions"
    )
  ),
  side = 2, line = -0.5, outer = TRUE, cex = 1.2, las = 0
)

# Size of plot:
# 13 X 11.8 inches - Portrait mode. 