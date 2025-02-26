# ------------------------------------------------------------------------------
# Data generating mechanism
# ------------------------------------------------------------------------------

shhh <- \(o) o |> suppressPackageStartupMessages() |> suppressWarnings()
shhh(library(mvnfast, include.only = "rmvn"))

library(BOBgmms)

ar.matrix <- \(rho, d) rho^abs(matrix(1:d - 1, d, d, byrow = T) - (1:d - 1))

generate.data <- function(sample.size, data.dim, n.clust,
                          probs = rep(1/n.clust, n.clust), gamm = 5, 
                          sparsity = 0.6, rho = 0, 
                          center = TRUE, scale = FALSE,
                          seed = sample.int(.Machine$integer.max, 1)) {
  
  #if (!(n.clust %in% c(2, 3, 4))) stop("n.clust should be 2, 3, or 4.")
  set.seed(seed)
  # Pre-compute constants
  d <- data.dim
  n <- sample.size
  
  # Locations of the clusters --------------------------------------------------
  #ciel_1 <- ceiling(d * sparsity/1)
  #ciel_2 <- ceiling(d * sparsity/2)
  #means <- vector("list", n.clust)
  #if (n.clust == 2) {
  #  means[[1]] <- c(rep(gamm, ciel_1), rep(0, d - ciel_1))
  #  means[[2]] <- c(rep(-gamm, ciel_1), rep(0, d - ciel_1))
  #} else if (n.clust == 3) {
  #  means[[1]] <- c(rep(gamm, ciel_1), rep(0, d - ciel_1))
  #  means[[2]] <- c(rep(-gamm, ciel_2), rep(gamm, ciel_2), rep(0, d - 2*ciel_2))
  #  means[[3]] <- c(rep(-gamm, ciel_1), rep(0, d - ciel_1))
  #} else if (n.clust == 4) {
  #  means[[1]] <- c(rep(gamm, ciel_1), rep(0, d - ciel_1))
  #  means[[2]] <- c(rep(-gamm, ciel_2), rep(gamm, ciel_2), rep(0, d - 2*ciel_2))
  #  means[[3]] <- c(rep(gamm, ciel_2), rep(-gamm, ciel_2), rep(0, d - 2*ciel_2))
  #  means[[4]] <- c(rep(-gamm, ciel_1), rep(0, d - ciel_1))
  #}
  ciel_d <- ceiling(d * sparsity)
  gammas <- seq(from = 1, by = gamm, length.out = n.clust)
  means <- lapply(gammas, \(g) c(rep(g, ciel_d), rep(0, d - ciel_d)))
  
  # Covariance matrices --------------------------------------------------------
  covs <- lapply(1:n.clust, \(k) ar.matrix(rho, d))
  
  # Generate latent clustering variables ---------------------------------------
  z <- sample(1:n.clust, n, replace = T, prob = probs)
  Z <- lapply(z, \(x) as.numeric(1:n.clust == x)) |> do.call(rbind, args = _)
  
  # Generate the observations --------------------------------------------------
  y <- lapply(
    seq_len(n), \(i) mvnfast::rmvn(1, means[[z[i]]], covs[[z[i]]])
  ) |> do.call(rbind, args = _)
  # Center and scale the data (if needed)
  y <- base::scale(y, center = center, scale = scale)
  attr(y, "scaled:center") <- NULL 
  attr(y, "scaled:scale")  <- NULL
  
  list(y = y, means = means, covs = covs, probs = probs, Z = Z)
}
