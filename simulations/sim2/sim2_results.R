load("simulations/sim2/wbb1Sim2.Rdata")
load("simulations/sim2/wbb2Sim2.Rdata")
load("simulations/sim2/nutsSim2.Rdata")
load("simulations/sim2/adviSim2.Rdata")
load("simulations/sim2/bobSim2.Rdata")

TVdistances <- list(WBB1 = sapply(wbb1Sim2, \(m) m$TV.wbb1),
                    WBB2 = sapply(wbb2Sim2, \(m) m$TV.wbb2),
                    BOB  = sapply(bobSim2,  \(m) m$TV.bob ),
                    NUTS = sapply(nutsSim2, \(m) m$TV.nuts),
                    ADVI = sapply(adviSim2, \(m) m$TV.advi))

KSdistances <- list(WBB1 = sapply(wbb1Sim2, \(m) m$KS.wbb1),
                    WBB2 = sapply(wbb2Sim2, \(m) m$KS.wbb2),
                    BOB  = sapply(bobSim2,  \(m) m$KS.bob),
                    NUTS = sapply(nutsSim2, \(m) m$KS.nuts),
                    ADVI = sapply(adviSim2, \(m) m$KS.advi))

runTimes    <- list(WBB1 = sapply(wbb1Sim2, \(m) m$Time.wbb1),
                    WBB2 = sapply(wbb2Sim2, \(m) m$Time.wbb2),
                    BOB  = sapply(bobSim2,  \(m) m$Time.bob),
                    NUTS = sapply(nutsSim2, \(m) m$Time.nuts),
                    ADVI = sapply(adviSim2, \(m) m$Time.advi)) 

boxplot(TVdistances)
boxplot(KSdistances)
boxplot(runTimes)

resultsSim2 <- list(TVdistances = TVdistances, 
                    KSdistances = KSdistances,
                    runTimes = runTimes)

save(resultsSim2, file = "simulations/sim2/resultsSim2.Rdata")
