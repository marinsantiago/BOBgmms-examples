# ------------------------------------------------------------------------------
# Warm initialization of model parameters
# ------------------------------------------------------------------------------

shhh <- \(o) o |> suppressPackageStartupMessages() |> suppressWarnings()
shhh(library(mclust)); shhh(library(sparcl)); shhh(library(clusterHD))

library(BOBgmms)

initial.values <- function(y, K) {
  # Pre-compute constants
  dims_y <- dim(y)
  n <- dims_y[1]
  d <- dims_y[2]
  # Pool of candidate initial values -------------------------------------------
  pool.vals <- initial_pool(y = y, K = K, n = n, d = d)
  # Maximize the log-likelihood ------------------------------------------------
  loglik.vals <- sapply(
    pool.vals, \(vals) {
      BOBgmms::loglik(
        means = vals$means, covs = vals$covs, probs = vals$probs, y = y
      )
    }
  )
  best.init <- which.max(loglik.vals)
  initial.vals <- pool.vals[[best.init]]
  list(values = initial.vals, algorithm = names(best.init))
}


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------


# Pool of candidate initial values
initial_pool <- function(y, K, n, d, HTKshrinkage = 0.75) {
  
  # Pre-compute constants
  adj <- diag(d) * 1e-10
  seq_K <- seq_len(K)
  
  # 1. HT K-means -- regularization path plot ----------------------------------
  htk.model <- clusterHD::HTKmeans(y, K)
  #clusterHD::diagPlot(clusterHD::HTKmeans(y, K)) # 0.75
  #abline(v = 1 - 0.75, col = 2)
  htk.lambda.idx.plot <- which.min(abs(htk.model$lambdas - HTKshrinkage))
  # Clustering
  htk.plot.clust <- c(htk.model$HTKmeans.out[[htk.lambda.idx.plot]]$cluster)
  # GMM parameters
  htk.plot.means <- lapply(seq_K, \(k) colMeans(y[htk.plot.clust == k, ]))
  htk.plot.covs <- lapply(seq_K, \(k) var(y[htk.plot.clust == k, ]) + adj)
  htk.plot.probs <- sapply(seq_K, \(k) sum(htk.plot.clust == k)) / n
  htk.plot.probs <- htk.plot.probs / sum(htk.plot.probs)
  htk.plot.out <- list(
    means = htk.plot.means, covs = htk.plot.covs, probs = htk.plot.probs
  )
  
  # 2. HT K-means -- with AIC --------------------------------------------------
  htk.model <- clusterHD::HTKmeans(y, K)
  htk.lambda.idx_aic <- which(
    htk.model$lambdas == clusterHD::getLambda(htk.model, type = "AIC")
  )
  # Clustering
  htk.aic.clust <- c(htk.model$HTKmeans.out[[htk.lambda.idx_aic]]$cluster)
  # GMM parameters
  htk.aic.means <- lapply(seq_K, \(k) colMeans(y[htk.aic.clust == k, ]))
  htk.aic.covs <- lapply(seq_K, \(k) var(y[htk.aic.clust == k, ]) + adj)
  htk.aic.probs <- sapply(seq_K, \(k) sum(htk.aic.clust == k)) / n
  htk.aic.probs <- htk.aic.probs / sum(htk.aic.probs)
  htk.aic.out <- list(
    means = htk.aic.means, covs = htk.aic.covs, probs = htk.aic.probs
  )
  
  # 3. HT K-means -- with BIC --------------------------------------------------
  htk.lambda.idx_bic <- which(
    htk.model$lambdas == clusterHD::getLambda(htk.model, type = "BIC")
  )
  # Clustering
  htk.bic.clust <- c(htk.model$HTKmeans.out[[htk.lambda.idx_bic]]$cluster)
  # GMM parameters
  htk.bic.means <- lapply(seq_K, \(k) colMeans(y[htk.bic.clust == k, ]))
  htk.bic.covs <- lapply(seq_K, \(k) var(y[htk.bic.clust == k, ]) + adj)
  htk.bic.probs <- sapply(seq_K, \(k) sum(htk.bic.clust == k)) / n
  htk.bic.probs <- htk.bic.probs / sum(htk.bic.probs)
  htk.bic.out <- list(
    means = htk.bic.means, covs = htk.bic.covs, probs = htk.bic.probs
  )
  
  # 4. GMM via mclust ----------------------------------------------------------
  mclust.bic <- mclust::mclustBIC(y, verbose = F, warn = F)
  mclust.mod <- mclust::Mclust(y, G = K, x = mclust.bic, verbose = F, warn = F)
  # Clustering
  mclust.clust <- mclust.mod$classification
  # GMM parameters
  mclust.means <- lapply(seq_K, \(k) colMeans(y[mclust.clust == k, ]))
  mclust.covs <- lapply(seq_K, \(k) var(y[mclust.clust == k, ]) + adj)
  mclust.probs <- sapply(seq_K, \(k)sum(mclust.clust == k)) / n
  mclust.probs <- mclust.probs / sum(mclust.probs)
  mclust.out <- list(
    means = mclust.means, covs = mclust.covs, probs = mclust.probs
  )
  
  # 5. Sparse K-means ----------------------------------------------------------
  sparsekm.tunparam <- tryCatch(
    {
      sparcl::KMeansSparseCluster.permute(
        y, K, nperms = 30, silent = TRUE
      )$bestw
    }, error = function(e){
      d/n * 50
    }
  )
  # Clustering
  sparsekm.clust <- sparcl::KMeansSparseCluster(
    y, K, wbounds = sparsekm.tunparam, nstart = 100, silent = TRUE
  )[[1]]$Cs
  # GMM parameters
  sparsekm.means <- lapply(seq_K, \(k) colMeans(y[sparsekm.clust == k, ]))
  sparsekm.covs <- lapply(seq_K, \(k) var(y[sparsekm.clust == k, ]) + adj)
  sparsekm.probs <- sapply(seq_K, \(k) sum(sparsekm.clust == k)) / n
  sparsekm.probs <- sparsekm.probs / sum(sparsekm.probs)
  sparsekm.out <- list(
    means = sparsekm.means, covs = sparsekm.covs, probs = sparsekm.probs
  )
  
  # 6. K-means -----------------------------------------------------------------
  km.clust <- stats::kmeans(y, K, iter.max = 100, nstart = 100)$cluster
  # GMM parameters
  km.means <- lapply(seq_K, \(k) colMeans(y[km.clust == k, ]))
  km.covs <- lapply(seq_K, \(k) var(y[km.clust == k, ]) + adj)
  km.probs <- sapply(seq_K, \(k) sum(km.clust == k)) / n
  km.probs <- km.probs / sum(km.probs)
  km.out <- list(means = km.means, covs = km.covs, probs = km.probs)
  gc() # Collect garbage
  
  list(
    "HTKmeans-plot" = htk.plot.out, "HTKmeans-AIC" = htk.aic.out, 
    "HTKmeans-BIC" = htk.bic.out, "GMM-mclust" = mclust.out,
    "SparseKmeans" = sparsekm.out, "Kmeans" = km.out
  )
}
