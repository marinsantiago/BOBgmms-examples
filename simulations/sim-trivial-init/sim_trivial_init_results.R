# ------------------------------------------------------------------------------
# Simulations under a trivial initialization - Results
# ------------------------------------------------------------------------------

# Load the results -------------------------------------------------------------

load("simulations/sim-trivial-init/results_d5.Rdata")
load("simulations/sim-trivial-init/results_d10.Rdata")
load("simulations/sim-trivial-init/results_d15.Rdata")

results_d5 <- lapply(results_d5, \(r) unlist(r))
results_d10 <- lapply(results_d10, \(r) unlist(r))
results_d15 <- lapply(results_d15, \(r) unlist(r))
methods_names <- c("WBB1", "WBB2", "BOB")
names(results_d5) <- names(results_d10) <- names(results_d15) <- methods_names

results <- list(results_d5, results_d10, results_d15)

# Boxplots ---------------------------------------------------------------------

col.palette <- c("#3D5A80", "#9B1D20", "tan2")
par(mfrow = c(1, 3), oma = c(0, 2, 0, 0))
for (s in seq_len(length(results))) {
  boxplot(
    x = results[[s]], 
    col = col.palette, pch = 16,
    main = paste("Setting TI", s),
    cex = 1, cex.main = 1.7, cex.axis = 1.3, cex.lab = 1.5,
    outline = TRUE
  )
}

mtext(
  text = as.expression(bquote("" ~ SW[2] ~ " Distance")),
  side = 2, line = -0.5, outer = TRUE, cex = 1.2, las = 0
)

# Size of plot:
# 10 X 4 inches - Landscape mode. 


# table of medians -------------------------------------------------------------

bob.median <- sapply(results, \(r) median(r$BOB))
bob.iqr <- sapply(results, \(r) IQR(r$BOB))

wbb1.median <- sapply(results, \(r) median(r$WBB1))
wbb1.iqr <- sapply(results, \(r) IQR(r$WBB1))

wbb2.median <- sapply(results, \(r) median(r$WBB2))
wbb2.iqr <- sapply(results, \(r) IQR(r$WBB2))

medians.df <- data.frame(
  bob.median, bob.iqr,
  wbb1.median, wbb1.iqr,
  wbb2.median, wbb2.iqr
) |> t()

xtable::xtable(medians.df, digits = 3)
