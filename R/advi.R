#install.packages("rstan")
shhh <- suppressPackageStartupMessages
shhh(library(rstan, include.only = c("stan_model", "vb")))

stan.model <- rstan::stan_model("./stan/GMM.stan")

advi <- function(y, init.params, Sadvi, alphas, betas, lambdas, nus, psis, 
                 vb.model = stan.model, warm_init = TRUE){
  
  dims.y <- dim(y)
  K      <- length(betas)
  n      <- dims.y[1]
  p      <- dims.y[2]
  
  stan.data <- list(K = K, n = n, p = p, Y = y, beta = betas, lambda = lambdas,
                    nu = nus, Psi = psis, alpha = alphas)
  
  if(warm_init){
    init        <- lapply(1:3, \(i) init.params[[i]])
    names(init) <- c("mu", "Sigma", "theta") 
    #init        <- list(init)
  } else {
    init        <- "random"
  }
  
  seed <- sample.int(.Machine$integer.max, 1)
  
  # Fit the model
  suppressWarnings(
    vb.fit <- rstan::vb(vb.model, data = stan.data, init = init, 
                        output_samples = Sadvi, seed = seed)
    # The argument `output_samples` wil be deprecated in the future. In that case,
    # use `draws` insetad.
  )
  
  # Extract model parameters
  advi.draws <- extract(vb.fit)
  pi.advi    <- advi.draws$theta
  mu.advi    <- do.call(cbind, lapply(1:K, \(k) advi.draws$mu[,k,]))
  sigma.advi <- do.call(cbind, lapply(1:K, \(k) 
                do.call(cbind, lapply(1:p, \(j) advi.draws$Sigma[,k,,j]))))
  
  cbind(mu.advi, sigma.advi, pi.advi)

}
