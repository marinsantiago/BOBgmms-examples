# ------------------------------------------------------------------------------
# Posterior predictive distribution of a GMM
# ------------------------------------------------------------------------------

post.predictive <- function(post.draws, d, K,
                            cores = parallel::detectCores() - 1,
                            seed = sample.int(.Machine$integer.max, 1)) {
  
  # Pre-compute constants
  dK = d * K
  m.idxs <- seq_len(dK)                # Mean indices
  c.idxs <- seq_len(dK * d) + dK       # Covariance indices
  p.idxs <- seq_len(K) + dK * (1 + d)  # Probs. indices
  max.iters <- nrow(post.draws)
  
  # Posterior predictive sampler (in parallel)
  base::RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  # Iterate
  post.pred.out <- pbmcapply::pbmclapply(
    seq_len(max.iters), \(iter) {
      theta <- post.draws[iter,]
      # Extract current parameters
      means.mat <- theta[m.idxs] |> matrix(data = _, nrow = K, byrow = T)
      covs.raw <- theta[c.idxs] |> matrix(data = _, nrow = K, byrow = T)
      covs.list <- lapply(seq_len(K), \(k) matrix(covs.raw[k,], ncol = d))
      probs.vec <- theta[p.idxs]
      # Sample from the model
      mvnfast::rmixn(1, means.mat, covs.list, probs.vec)
    },
    mc.set.seed = TRUE,
    mc.cores = cores
  ) |> unlist(x = _) |> matrix(data = _, nrow = max.iters, byrow = TRUE)
  # Re-set random generator
  base::RNGkind("default", "default", "default")
  
  post.pred.out
}
