load("simulations/sim7/wbb1Sim7.Rdata")
load("simulations/sim7/wbb2Sim7.Rdata")
load("simulations/sim7/nutsSim7.Rdata")
load("simulations/sim7/adviSim7.Rdata")
load("simulations/sim7/bobSim7.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim7, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim7, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim7,  \(m) m$TV.bob ),
                    NUTS = sapply(nutsSim7, \(m) m$TV.nuts),
                    ADVI = sapply(adviSim7, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim7, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim7, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim7,  \(m) m$KS.bob),
                    NUTS = sapply(nutsSim7, \(m) m$KS.nuts),
                    ADVI = sapply(adviSim7, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim7, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim7, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim7,  \(m) m$Time.bob),
                    NUTS = sapply(nutsSim7, \(m) m$Time.nuts),
                    ADVI = sapply(adviSim7, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim7 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim7, file = "simulations/sim7/resultsSim7.Rdata")
