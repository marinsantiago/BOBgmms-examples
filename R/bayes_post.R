# ------------------------------------------------------------------------------
# Bayesian posterior distribution of a GMM
# ------------------------------------------------------------------------------

shhh <- \(o) o |> suppressPackageStartupMessages() |> suppressWarnings()
shhh(library(LaplacesDemon, include.only = c("rinvwishart", "rdirichlet")))
shhh(library(mvnfast, include.only = "dmvn"))
shhh(library(pbmcapply)); shhh(library(parallel))

library(BOBgmms)

bayesian_posterior <- function(y, true.Z.mat,
                               betas, lambdas, nus, psis, alphas,
                               max.iters = 20000,
                               cores = parallel::detectCores(),
                               seed = sample.int(.Machine$integer.max, 1)) {
  
  # Pre-compute constants
  K <- ncol(true.Z.mat)
  seq_K <- seq_len(K)
  
  # Summary statistics ---------------------------------------------------------
  n_ks <- apply(true.Z.mat, 2, sum)
  # The "== 1" is correct, as we want to extract the Z equal to one for each "k"
  clusters <- lapply(seq_K, \(k) y[true.Z.mat[,k] == 1, ]) 
  y.bar_ks <- lapply(clusters, colMeans)
  SSR_ks <- lapply(seq_K, \(k) var(clusters[[k]]) * (n_ks[k] - 1))
  
  # Posterior hyper-parameters -------------------------------------------------
  alphas.hat <- alphas + n_ks
  nus.hat <- nus + n_ks
  lambdas.hat <- lambdas + n_ks
  betas.hat <- lapply(
    seq_K, \(k) {
      (lambdas[k] / lambdas.hat[k] * betas[[k]]) +
        (n_ks[k] / lambdas.hat[k] * y.bar_ks[[k]])
    }
  )
  psis.hat <- lapply(
    seq_K, \(k) {
      (psis[[k]] + SSR_ks[[k]]) + ((lambdas[k] * n_ks[k]) / lambdas.hat[k]) *
        tcrossprod(y.bar_ks[[k]] - betas[[k]])
    }
  )
  
  # Posterior sampler (in parallel) --------------------------------------------
  base::RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  # Iterate
  bayes.out <- pbmcapply::pbmclapply(
    seq_len(max.iters), \(iter) {
      Sigmas_k <- lapply(seq_K, \(k) rinvwishart(nus.hat[[k]], psis.hat[[k]]))
      means_k <- lapply(
        seq_K, \(k) {
          c(rmvn(1, betas.hat[[k]], Sigmas_k[[k]] / lambdas.hat[k]))
        }
      )
      probs_k <- c(rdirichlet(1, alphas.hat))
      # Theta vector
      list(means_k, Sigmas_k, probs_k) |> unlist()
    },
    mc.set.seed = TRUE,
    mc.cores = cores
  ) |> unlist(x = _) |> matrix(data = _, nrow = max.iters, byrow = TRUE)
  # Re-set random generator
  base::RNGkind("default", "default", "default")
  
  bayes.out
}
