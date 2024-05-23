shhh <- suppressPackageStartupMessages
shhh(library(mvnfast, include.only = "rmvn"))

generate.data <- function(p, n, K, probs = rep(1/K, K), 
                          gam = 5, sparsity = 0.6, scale = TRUE){
  
  # Number of important parameters
  imp.params <- ceiling(p * sparsity)
  
  # Location of the clusters
  gammas <- seq(from = 1, by = gam, length.out = K)
  
  # Generate mean vectors
  mus <- lapply(gammas, \(g) c(rep(g, imp.params), rep(0, p - imp.params)))
  
  # Generate covariance matrices
  sigmas <- lapply(1:K, \(k) diag(p))
  
  # Generate latent clustering variables
  z <- sample(1:K, n, replace = T, prob = probs)
  
  # Generate a single observation
  genOne <- function(mus, sigmas, z) mvnfast::rmvn(1, mus[[z]], sigmas[[z]])
  
  y <- do.call(rbind, lapply(1:n, \(i) genOne(mus, sigmas, z[i])))
  Z <- do.call(rbind, lapply(z,   \(x) as.numeric(1:K == x)))
  
  if(scale) {
    
    y <- scale(y)
    attr(y, "scaled:center") <- NULL 
    attr(y, "scaled:scale")  <- NULL
    
  }
  
  list("y" = y, "means" = mus, "covs" = sigmas, "probs" = probs, "Z" = Z)
  
}
