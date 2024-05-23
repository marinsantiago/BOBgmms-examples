load("simulations/sim5/wbb1Sim5.Rdata")
load("simulations/sim5/wbb2Sim5.Rdata")
load("simulations/sim5/nutsSim5.Rdata")
load("simulations/sim5/adviSim5.Rdata")
load("simulations/sim5/bobSim5.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim5, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim5, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim5,  \(m) m$TV.bob ),
                    NUTS = sapply(nutsSim5, \(m) m$TV.nuts),
                    ADVI = sapply(adviSim5, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim5, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim5, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim5,  \(m) m$KS.bob),
                    NUTS = sapply(nutsSim5, \(m) m$KS.nuts),
                    ADVI = sapply(adviSim5, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim5, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim5, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim5,  \(m) m$Time.bob),
                    NUTS = sapply(nutsSim5, \(m) m$Time.nuts),
                    ADVI = sapply(adviSim5, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim5 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim5, file = "simulations/sim5/resultsSim5.Rdata")
