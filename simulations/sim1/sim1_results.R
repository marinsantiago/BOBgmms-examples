load("simulations/sim1/wbb1Sim1.Rdata")
load("simulations/sim1/wbb2Sim1.Rdata")
load("simulations/sim1/nutsSim1.Rdata")
load("simulations/sim1/adviSim1.Rdata")
load("simulations/sim1/bobSim1.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim1, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim1, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim1,  \(m) m$TV.bob ),
                    NUTS = sapply(nutsSim1, \(m) m$TV.nuts),
                    ADVI = sapply(adviSim1, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim1, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim1, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim1,  \(m) m$KS.bob),
                    NUTS = sapply(nutsSim1, \(m) m$KS.nuts),
                    ADVI = sapply(adviSim1, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim1, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim1, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim1,  \(m) m$Time.bob),
                    NUTS = sapply(nutsSim1, \(m) m$Time.nuts),
                    ADVI = sapply(adviSim1, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim1 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim1, file = "simulations/sim1/resultsSim1.Rdata")
