#install.packages("rstan")
shhh <- suppressPackageStartupMessages
shhh(library(rstan, include.only = c("stan", "extract")))

nuts.Sampler <- function(y, init.params, Snuts,
                         alphas, betas, lambdas, nus, psis, warm_init = TRUE){
  
  dims.y <- dim(y)
  K      <- length(betas)
  n      <- dims.y[1]
  p      <- dims.y[2]
  
  stan.data <- list(K = K, n = n, p = p, Y = y, beta = betas, lambda = lambdas,
                    nu = nus, Psi = psis, alpha = alphas)
  
  if(warm_init){
    init        <- init.params
    names(init) <- c("mu", "Sigma", "theta") 
    init        <- list(init)
  } else {
    init        <- "random"
  }
  
  seed <- sample.int(.Machine$integer.max, 1)
  
  # Fit the model
  suppressWarnings(
    stan.fit <- rstan::stan("./stan/GMM.stan", data = stan.data, chains = 1,
                            iter = Snuts, init = init, seed = seed)
  )
  
  # Extract model parameters
  stan.draws <- rstan::extract(stan.fit)
  pi.stan    <- stan.draws$theta
  mu.stan    <- do.call(cbind, lapply(1:K, \(k) stan.draws$mu[,k,]))
  sigma.stan <- do.call(cbind, lapply(1:K, \(k) 
                do.call(cbind, lapply(1:p, \(j) stan.draws$Sigma[,k,,j]))))
  
  cbind(mu.stan, sigma.stan, pi.stan)
  
}
