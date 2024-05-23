load("simulations/sim8/wbb1Sim8.Rdata")
load("simulations/sim8/wbb2Sim8.Rdata")
load("simulations/sim8/nutsSim8.Rdata")
load("simulations/sim8/adviSim8.Rdata")
load("simulations/sim8/bobSim8.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim8, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim8, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim8,  \(m) m$TV.bob ),
                    NUTS = sapply(nutsSim8, \(m) m$TV.nuts),
                    ADVI = sapply(adviSim8, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim8, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim8, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim8,  \(m) m$KS.bob),
                    NUTS = sapply(nutsSim8, \(m) m$KS.nuts),
                    ADVI = sapply(adviSim8, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim8, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim8, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim8,  \(m) m$Time.bob),
                    NUTS = sapply(nutsSim8, \(m) m$Time.nuts),
                    ADVI = sapply(adviSim8, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim8 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim8, file = "simulations/sim8/resultsSim8.Rdata")

