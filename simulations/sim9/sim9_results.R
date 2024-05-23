load("simulations/sim9/wbb1Sim9.Rdata")
load("simulations/sim9/wbb2Sim9.Rdata")
load("simulations/sim9/adviSim9.Rdata")
load("simulations/sim9/bobSim9.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim9, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim9, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim9,  \(m) m$TV.bob ),
                    ADVI = sapply(adviSim9, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim9, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim9, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim9,  \(m) m$KS.bob),
                    ADVI = sapply(adviSim9, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim9, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim9, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim9,  \(m) m$Time.bob),
                    ADVI = sapply(adviSim9, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim9 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim9, file = "simulations/sim9/resultsSim9.Rdata")
