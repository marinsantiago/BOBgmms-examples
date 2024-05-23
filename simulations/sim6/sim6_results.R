load("simulations/sim6/wbb1Sim6.Rdata")
load("simulations/sim6/wbb2Sim6.Rdata")
load("simulations/sim6/nutsSim6.Rdata")
load("simulations/sim6/adviSim6.Rdata")
load("simulations/sim6/bobSim6.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim6, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim6, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim6,  \(m) m$TV.bob ),
                    NUTS = sapply(nutsSim6, \(m) m$TV.nuts),
                    ADVI = sapply(adviSim6, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim6, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim6, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim6,  \(m) m$KS.bob),
                    NUTS = sapply(nutsSim6, \(m) m$KS.nuts),
                    ADVI = sapply(adviSim6, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim6, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim6, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim6,  \(m) m$Time.bob),
                    NUTS = sapply(nutsSim6, \(m) m$Time.nuts),
                    ADVI = sapply(adviSim6, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim6 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim6, file = "simulations/sim6/resultsSim6.Rdata")
