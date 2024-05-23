# ---------------------------------------------------------#
#    Medians for TV and KS distances across 9 settings     #
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

# -----------------------------------------
# TV Distances
# -----------------------------------------

# BOB
(TV.bob.median <- sapply(results, \(r) median(r$TVdistances$BOB)))
(TV.bob.iqr    <- sapply(results, \(r) IQR(r$TVdistances$BOB)))

# WBB1
(TV.wbb1.median <- sapply(results, \(r) median(r$TVdistances$WBB1)))
(TV.wbb1.iqr    <- sapply(results, \(r) IQR(r$TVdistances$WBB1)))

# WBB2
(TV.wbb2.median <- sapply(results, \(r) median(r$TVdistances$WBB2)))
(TV.wbb2.iqr    <- sapply(results, \(r) IQR(r$TVdistances$WBB2)))

# NUTS
(TV.nuts.median <- sapply(results, \(r) quantile(r$TVdistances$NUTS, 0.5)))
(TV.nuts.iqr    <- sapply(results, \(r) IQR(r$TVdistances$NUTS)))

# ADVI
(TV.advi.median <- sapply(results, \(r) median(r$TVdistances$ADVI, na.rm = T)))
(TV.advi.iqr    <- sapply(results, \(r) IQR(r$TVdistances$ADVI, na.rm = T)))

TV.df <- t(data.frame(TV.bob.median,  TV.bob.iqr,
                      TV.wbb1.median, TV.wbb1.iqr, 
                      TV.wbb2.median, TV.wbb2.iqr, 
                      TV.nuts.median, TV.nuts.iqr,
                      TV.advi.median, TV.advi.iqr))

xtable::xtable(TV.df, digits = 3)


# -----------------------------------------
# KS Distances
# -----------------------------------------

# BOB
(KS.bob.median <- sapply(results, \(r) median(r$KSdistances$BOB)))
(KS.bob.iqr    <- sapply(results, \(r) IQR(r$KSdistances$BOB)))

# WBB1
(KS.wbb1.median <- sapply(results, \(r) median(r$KSdistances$WBB1)))
(KS.wbb1.iqr    <- sapply(results, \(r) IQR(r$KSdistances$WBB1)))

# WBB2
(KS.wbb2.median <- sapply(results, \(r) median(r$KSdistances$WBB2)))
(KS.wbb2.iqr    <- sapply(results, \(r) IQR(r$KSdistances$WBB2)))

# NUTS
(KS.nuts.median <- sapply(results, \(r) quantile(r$KSdistances$NUTS, 0.5)))
(KS.nuts.iqr    <- sapply(results, \(r) IQR(r$KSdistances$NUTS)))

# ADVI
(KS.advi.median <- sapply(results, \(r) median(r$KSdistances$ADVI, na.rm = T)))
(KS.advi.iqr    <- sapply(results, \(r) IQR(r$KSdistances$ADVI, na.rm = T)))

KS.df <- t(data.frame(KS.bob.median,  KS.bob.iqr,
                      KS.wbb1.median, KS.wbb1.iqr, 
                      KS.wbb2.median, KS.wbb2.iqr, 
                      KS.nuts.median, KS.nuts.iqr,
                      KS.advi.median, KS.advi.iqr))

xtable::xtable(KS.df, digits = 3)


# -----------------------------------------
# Elapsed
# -----------------------------------------

# BOB
(Time.bob.median <- sapply(results, \(r) median(r$runTimes$BOB)))
(Time.bob.iqr    <- sapply(results, \(r) IQR(r$runTimes$BOB)))

# WBB1
(Time.wbb1.median <- sapply(results, \(r) median(r$runTimes$WBB1)))
(Time.wbb1.iqr    <- sapply(results, \(r) IQR(r$runTimes$WBB1)))

# WBB2
(Time.wbb2.median <- sapply(results, \(r) median(r$runTimes$WBB2)))
(Time.wbb2.iqr    <- sapply(results, \(r) IQR(r$runTimes$WBB2)))

# NUTS
(Time.nuts.median <- sapply(results, \(r) quantile(r$runTimes$NUTS, 0.5)))
(Time.nuts.iqr    <- sapply(results, \(r) IQR(r$runTimes$NUTS)))

# ADVI
(Time.advi.median <- sapply(results, \(r) median(r$runTimes$ADVI, na.rm = T)))
(Time.advi.iqr    <- sapply(results, \(r) IQR(r$runTimes$ADVI, na.rm = T)))

Time.df <- t(data.frame(Time.bob.median,  Time.bob.iqr,
                        Time.wbb1.median, Time.wbb1.iqr, 
                        Time.wbb2.median, Time.wbb2.iqr, 
                        Time.nuts.median, Time.nuts.iqr,
                        Time.advi.median, Time.advi.iqr))

xtable::xtable(Time.df, digits = 3)
