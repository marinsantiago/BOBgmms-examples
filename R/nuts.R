# ------------------------------------------------------------------------------
# No U-Turn-Sampler (NUTS)
# ------------------------------------------------------------------------------

shhh <- \(o) o |> suppressPackageStartupMessages() |> suppressWarnings()
shhh(library(rstan, include.only = c("stan_model", "stan", "extract")))


nuts <- function(y, means.init, covs.init, probs.init,
                 betas, lambdas, nus, psis, alphas,
                 max.iters = 40000, warm_init = TRUE,
                 seed = sample.int(.Machine$integer.max, 1)) {
  
  dims.y <- dim(y)
  n <- dims.y[1]
  p <- dims.y[2]
  K <- length(betas)
  
  # Input data for stan model --------------------------------------------------
  stan.data <- list(
    K = K, n = n, p = p, Y = y, 
    beta = betas, lambda = lambdas, 
    nu = nus, Psi = psis, alpha = alphas
  )
  
  # Initialization -------------------------------------------------------------
  if (warm_init) {
    init <- list(means.init, covs.init, probs.init)
    names(init) <- c("mu", "Sigma", "theta")
    init <- list(init)
  } else {
    init <- "random"
  }
  
  # Fit the model --------------------------------------------------------------
  stan.fit <- rstan::stan(
    file = "./stan/GMM.stan", data = stan.data, chains = 1,
    iter = max.iters, init = init, seed = seed
  ) |> suppressWarnings() |> suppressMessages()
  
  # Extract model parameters ---------------------------------------------------
  mat_list <- \(x) do.call(cbind, x)
  out_nuts <- rstan::extract(stan.fit)
  means.nuts <- lapply(
    1:K, \(k) out_nuts$mu[,k,]
  ) |> mat_list(x = _)
  covs.nuts <- lapply(
    1:K, \(k) {
      lapply(1:p, \(j) out_nuts$Sigma[,k,,j]) |> mat_list(x = _)
    }
  ) |> mat_list(x = _)
  probs.nuts <- out_nuts$theta
  rm(stan.fit, out_nuts)
  gc()
  
  cbind(means.nuts, covs.nuts, probs.nuts)
}
