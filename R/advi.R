# ------------------------------------------------------------------------------
# Automatic Differentiation Variational Inference (ADVI)
# ------------------------------------------------------------------------------

shhh <- \(o) o |> suppressPackageStartupMessages() |> suppressWarnings()
shhh(library(rstan, include.only = c("stan_model", "vb", "extract")))

# Set up the stan model
stan.model <- rstan::stan_model("./stan/GMM.stan")

advi <- function(y, means.init, covs.init, probs.init,
                 betas, lambdas, nus, psis, alphas,
                 max.iters = 20000, vb.model = stan.model,
                 warm_init = TRUE,
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
  } else {
    init <- "random"
  }
  
  # Fit the model --------------------------------------------------------------
  vb.fit <- rstan::vb(
    object = vb.model, data = stan.data, init = init,
    output_samples = max.iters, seed = seed
  ) |> suppressWarnings() |> suppressMessages()
  # The argument `output_samples` wil be deprecated in the future. In that case,
  # use `draws` instead.
  
  # Extract model parameters ---------------------------------------------------
  mat_list <- \(x) do.call(cbind, x)
  out_advi <- rstan::extract(vb.fit)
  means.advi <- lapply(
    1:K, \(k) out_advi$mu[,k,]
  ) |> mat_list(x = _)
  covs.advi <- lapply(
    1:K, \(k) {
      lapply(1:p, \(j) out_advi$Sigma[,k,,j]) |> mat_list(x = _)
    }
  ) |> mat_list(x = _)
  probs.advi <- out_advi$theta
  rm(vb.fit, out_advi)
  gc()
  
  cbind(means.advi, covs.advi, probs.advi)
}
