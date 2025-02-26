# ------------------------------------------------------------------------------
# Sim setting 1-9: Table Medians
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
rr <- results

# SW2 Distances ----------------------------------------------------------------

# BOB
SW2.bob.median <- sapply(rr, \(r) median(r$SW2.distance$BOB, na.rm = TRUE))
SW2.bob.iqr <- sapply(rr, \(r) IQR(r$SW2.distance$BOB, na.rm = TRUE))

# WBB1
SW2.wbb1.median <- sapply(rr, \(r) median(r$SW2.distance$WBB1, na.rm = TRUE))
SW2.wbb1.iqr <- sapply(rr, \(r) IQR(r$SW2.distance$WBB1, na.rm = TRUE))

# WBB2
SW2.wbb2.median <- sapply(rr, \(r) median(r$SW2.distance$WBB2, na.rm = TRUE))
SW2.wbb2.iqr <- sapply(rr, \(r) IQR(r$SW2.distance$WBB2, na.rm = TRUE))

# NUTS
SW2.nuts.median <- sapply(rr, \(r) quantile(r$SW2.distance$NUTS, 0.5))
SW2.nuts.iqr <- sapply(rr, \(r) IQR(r$SW2.distance$NUTS, na.rm = TRUE))

# ADVI
SW2.advi.median <- sapply(rr, \(r) median(r$SW2.distance$ADVI, na.rm = TRUE))
SW2.advi.iqr <- sapply(rr, \(r) IQR(r$SW2.distance$ADVI, na.rm = TRUE))

# SW2 output
SW2.df <- data.frame(
  SW2.bob.median, SW2.bob.iqr,
  SW2.wbb1.median, SW2.wbb1.iqr,
  SW2.wbb2.median, SW2.wbb2.iqr,
  SW2.nuts.median, SW2.nuts.iqr,
  SW2.advi.median, SW2.advi.iqr
) |> t()

xtable::xtable(SW2.df, digits = 3)


# Elapsed Times ----------------------------------------------------------------

# BOB
times.bob.median <- sapply(rr, \(r) median(r$run.times$BOB, na.rm = TRUE))
times.bob.iqr <- sapply(rr, \(r) IQR(r$run.times$BOB, na.rm = TRUE))

# WBB1
times.wbb1.median <- sapply(rr, \(r) median(r$run.times$WBB1, na.rm = TRUE))
times.wbb1.iqr <- sapply(rr, \(r) IQR(r$run.times$WBB1, na.rm = TRUE))

# WBB2
times.wbb2.median <- sapply(rr, \(r) median(r$run.times$WBB2, na.rm = TRUE))
times.wbb2.iqr <- sapply(rr, \(r) IQR(r$run.times$WBB2, na.rm = TRUE))

# NUTS
times.nuts.median <- sapply(rr, \(r) quantile(r$run.times$NUTS, 0.5, na.rm = T))
times.nuts.iqr <- sapply(rr, \(r) IQR(r$run.times$NUTS, na.rm = TRUE))

# ADVI
times.advi.median <- sapply(rr, \(r) median(r$run.times$ADVI, na.rm = TRUE))
times.advi.iqr <- sapply(rr, \(r) IQR(r$run.times$ADVI, na.rm = TRUE))

# Elapsed output
times.df <- data.frame(
  times.bob.median, times.bob.iqr,
  times.wbb1.median, times.wbb1.iqr,
  times.wbb2.median, times.wbb2.iqr,
  times.nuts.median, times.nuts.iqr,
  times.advi.median, times.advi.iqr
) |> t()

xtable::xtable(times.df, digits = 3)
