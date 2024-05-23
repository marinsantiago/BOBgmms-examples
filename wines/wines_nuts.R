#devtools::install_github("marinsantiago/BOBgmms")
library(BOBgmms)
source("./R/nuts.R"); source("./R/advi.R")
rm(shhh)

load("wines/priors.Rdata")
load("wines/y_data.Rdata")
load("wines/init_params.Rdata")

set.seed(100)

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

d <- 13
K <- 3

y.new.nuts <- post.predictive(nuts.draws,  d, K)
y.new.advi <- post.predictive(advi.draws,  d, K)

# Export the results
save(nuts.draws, file = "wines/nuts_draws.Rdata")
save(advi.draws, file = "wines/advi_draws.Rdata")

save(y.new.nuts, file = "wines/y_new_nuts.Rdata")
save(y.new.advi, file = "wines/y_new_advi.Rdata")

save(nuts.time, file = "wines/nuts_time.Rdata")
save(advi.time, file = "wines/advi_time.Rdata")
