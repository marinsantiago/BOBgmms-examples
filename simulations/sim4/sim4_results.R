load("simulations/sim4/wbb1Sim4.Rdata")
load("simulations/sim4/wbb2Sim4.Rdata")
load("simulations/sim4/nutsSim4.Rdata")
load("simulations/sim4/adviSim4.Rdata")
load("simulations/sim4/bobSim4.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim4, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim4, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim4,  \(m) m$TV.bob ),
                    NUTS = sapply(nutsSim4, \(m) m$TV.nuts),
                    ADVI = sapply(adviSim4, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim4, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim4, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim4,  \(m) m$KS.bob),
                    NUTS = sapply(nutsSim4, \(m) m$KS.nuts),
                    ADVI = sapply(adviSim4, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim4, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim4, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim4,  \(m) m$Time.bob),
                    NUTS = sapply(nutsSim4, \(m) m$Time.nuts),
                    ADVI = sapply(adviSim4, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim4 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim4, file = "simulations/sim4/resultsSim4.Rdata")
