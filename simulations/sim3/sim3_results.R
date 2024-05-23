load("simulations/sim3/wbb1Sim3.Rdata")
load("simulations/sim3/wbb2Sim3.Rdata")
load("simulations/sim3/nutsSim3.Rdata")
load("simulations/sim3/adviSim3.Rdata")
load("simulations/sim3/bobSim3.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim3, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim3, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim3,  \(m) m$TV.bob ),
                    NUTS = sapply(nutsSim3, \(m) m$TV.nuts),
                    ADVI = sapply(adviSim3, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim3, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim3, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim3,  \(m) m$KS.bob),
                    NUTS = sapply(nutsSim3, \(m) m$KS.nuts),
                    ADVI = sapply(adviSim3, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim3, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim3, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim3,  \(m) m$Time.bob),
                    NUTS = sapply(nutsSim3, \(m) m$Time.nuts),
                    ADVI = sapply(adviSim3, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim3 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim3, file = "simulations/sim3/resultsSim3.Rdata")
