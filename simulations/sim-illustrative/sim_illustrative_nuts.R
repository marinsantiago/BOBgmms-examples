#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/nuts.R"); source("./R/advi.R")
rm(shhh)

load("simulations/sim-illustrative/priors.Rdata")
load("simulations/sim-illustrative/y_data.Rdata")
load("simulations/sim-illustrative/init_params.Rdata")

set.seed(101)

alphas  <- priors$alphas
betas   <- priors$betas
lambdas <- priors$lambdas
nus     <- priors$nus
psis    <- priors$psis 

# -----------------------------------------
# No-U-Turn-Sampler (NUTS)
# -----------------------------------------

S <- 20000

# NUTS - Stan
start.nuts <- Sys.time()
nuts.draws <- nuts.Sampler(y, init.params, S * 2,
                           alphas, betas, lambdas, nus, psis)
end.nuts   <- Sys.time()
(nuts.time <- difftime(end.nuts, start.nuts, units = "mins"))


# -----------------------------------------
# ADVI - Mean Field Variational Bayes
# -----------------------------------------

# ADVI - Stan
start.advi <- Sys.time()
advi.draws <- advi(y, init.params, S, alphas, betas, lambdas, nus, psis)
end.advi   <- Sys.time()
(advi.time <- difftime(end.advi, start.advi, units = "mins"))

# -----------------------------------------
# Posterior predictive distributions
# -----------------------------------------

d <- 10
K <- 2

y.new.nuts <- post.predictive(nuts.draws,  d, K)
y.new.advi <- post.predictive(advi.draws,  d, K)

# Export the results
save(nuts.draws, file = "simulations/sim-illustrative/nuts_draws.Rdata")
save(advi.draws, file = "simulations/sim-illustrative/advi_draws.Rdata")

save(y.new.nuts, file = "simulations/sim-illustrative/y_new_nuts.Rdata")
save(y.new.advi, file = "simulations/sim-illustrative/y_new_advi.Rdata")

save(nuts.time, file = "simulations/sim-illustrative/nuts_time.Rdata")
save(advi.time, file = "simulations/sim-illustrative/advi_time.Rdata")
