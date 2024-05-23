library(BOBgmms)
shhh <- suppressPackageStartupMessages
shhh(library(clusterHD)); shhh(library(parallel))
shhh(library(mclust, include.only = c("mclustBIC", "Mclust")))
shhh(library(sparcl, include.only = c("KMeansSparseCluster.permute", 
                                      "KMeansSparseCluster")))

# --------------------------------------------
# Tune shrinkage hyper-parameters via CV
# --------------------------------------------

cv.tuning <- function(y, alphas, betas, psis,
                      HTKshrinkage = 0.75, init = "Pool", 
                      mult = min(length(alphas), 3)){

  n <- nrow(y)
  d <- ncol(y)
  K <- length(alphas)
  
  n.train <- floor(n * 0.9)  
  n.test  <- n - n.train
  
  # Randomly split the data between train and validation sets
  shuffle.indices <- sample(n, n)
  train.indices   <- shuffle.indices[1:n.train]
  test.indices    <- shuffle.indices[(n.train + 1):n]
  y.train         <- y[train.indices,]
  y.test          <- y[test.indices,]
  
  rm(shuffle.indices); rm(train.indices); rm(test.indices)
  
  # Initial values cv procedure
  if (init == "Pool") {
    
    pool.train <- initialPool(y.train, K, n.train, d, HTKshrinkage)
    best.init  <- which.max(sapply(pool.train, \(v) 
                                   BOBgmms::log_lik(v, y.train, d, K)))
    optimVals  <- pool.train[[best.init]]
    
    # Extract the mean vectors
    mu.train <- split(optimVals[1:(d * K)], unlist(lapply(1:K, \(k) rep(k, d))))
    names(mu.train) <- NULL
    
    # Extract the covariance matrices
    sigma.train <- optimVals[(K * d + ((1 - 1) * d^(2)) + 1):(K * d + (K * d^(2)))]
    sigma.train <- split(sigma.train, unlist(lapply(1:K, \(k) rep(k, d * d))))
    sigma.train <- lapply(sigma.train, \(c) matrix(c, ncol = d, nrow = d))
    names(sigma.train) <- NULL
    
    # Extract the probabilities
    prob.train <- tail(optimVals, K)
    
    rm(pool.train); rm(optimVals); rm(best.init)
    
  } else if (init == "HT") {
    
    adj  <- diag(d) * 1e-10
    mod  <- clusterHD::HTKmeans(y.train, K)
    best <- which(abs(mod$lambdas - HTKshrinkage) < 1e-5)
    
    # Extract model parameters
    best.mod    <- mod$HTKmeans.out[[best]]
    mu.train    <- lapply(1:K, \(k) best.mod$centers[k,])
    prob.train  <- sapply(1:K, \(k) sum(best.mod$cluster == k)) / n.train
    sigma.train <- lapply(1:K, \(k) var(y.train[best.mod$cluster == k,]) + adj)
    
    rm(adj); rm(mod); rm(best); rm(best.mod)
    
  }
  
  # Grid for lambda and nu
  Lambdas <- sort(c((d * mult) / n, seq(0.01, 10, length.out = 9) * (d * mult) / n))
  Nus     <- sort(c((d^(2) * mult) / n, 
                    seq(0.01, 10, length.out = 9) * (d^(2) * mult) / n)) + d + 1
  
  log.lik.best <- -Inf
  best.lambda  <- best.nu <- NULL
  
  for (lambda in Lambdas) {
    for (nu in Nus) {
      
      # Run EM algorithm with current values of lambda and nu
      theta <- BOBgmms::em.optim(y.train, mu.train, sigma.train, prob.train,
        0, 0, 0, 1, rep(1, n), rep(1, K), rep(1, K), 1,
          alphas, betas, rep(lambda, K), rep(nu, K), psis)
      
      # Evaluate the log-likelihood over the validation set
      log.lik.current <- BOBgmms::log_lik(theta, y.test, d, K)
      
      # Update best log-likelihood
      if (log.lik.current > log.lik.best) {
        
        log.lik.best <- log.lik.current
        
        # Combination of lambda and nu that yields the best likelihood
        best.lambda <- lambda
        best.nu     <- nu
        
      }
    }
  }

  if (is.null(c(best.lambda, best.nu))) {
    out <- list("lambda" = (d * mult) / n, "nu" = (d^(2) * mult) / n + d + 1)
  } else {
    out <- list("lambda" = best.lambda, "nu" = best.nu)
  }
  
  out
}

# --------------------------------------------
# Pool of initial values
# --------------------------------------------

initialPool <- function(y, K, n, d, HTKshrinkage){
  
  adj <- diag(d) * 1e-10
  
  # 0. HT K-means -- Tuned via regularization path plot
  mod0   <- clusterHD::HTKmeans(y, K)
  best0  <- which(abs(mod0$lambdas - HTKshrinkage) < 1e-5)
  loc0   <- lapply(1:K, \(k) colMeans(y[mod0$HTKmeans.out[[best0]]$cluster == k,]))
  prob0  <- sapply(1:K, \(k) sum(mod0$HTKmeans.out[[best0]]$cluster == k)) / n
  covs0  <- lapply(1:K, \(k) var(y[mod0$HTKmeans.out[[best0]]$cluster == k,]) + adj)
  theta0 <- c(unlist(loc0), unlist(covs0), prob0)
  
  # 1. HT K-means -- AIC
  best1  <- which(mod0$lambdas == clusterHD::getLambda(mod0, type = "AIC"))
  loc1   <- lapply(1:K, \(k) colMeans(y[mod0$HTKmeans.out[[best1]]$cluster == k,]))
  prob1  <- sapply(1:K, \(k) sum(mod0$HTKmeans.out[[best1]]$cluster == k)) / n
  covs1  <- lapply(1:K, \(k) var(y[mod0$HTKmeans.out[[best1]]$cluster == k,]) + adj)
  theta1 <- c(unlist(loc1), unlist(covs1), prob1)
  
  # 2. HT K-means -- BIC
  best2  <- which(mod0$lambdas == clusterHD::getLambda(mod0, type = "BIC"))
  loc2   <- lapply(1:K, \(k) colMeans(y[mod0$HTKmeans.out[[best2]]$cluster == k,]))
  prob2  <- sapply(1:K, \(k) sum(mod0$HTKmeans.out[[best2]]$cluster == k)) / n
  covs2  <- lapply(1:K, \(k) var(y[mod0$HTKmeans.out[[best2]]$cluster == k,]) + adj)
  theta2 <- c(unlist(loc2), unlist(covs2), prob2)
  
  # 3. GMM -- mclust
  mcBIC  <- mclust::mclustBIC(y, verbose = F, warn = F)
  mod3   <- mclust::Mclust(y, G = K, x = mcBIC, verbose = F, warn = F)
  loc3   <- lapply(1:K, \(k) mod3$parameters$mean[,k])
  prob3  <- sapply(1:K, \(k) sum(mod3$classification == k)) / n
  covs3  <- lapply(1:K, \(k) var(y[mod3$classification == k,]) + adj)
  theta3 <- c(unlist(loc3), unlist(covs3), prob3)
  
  # 4. Sparse K-means
  best3  <- sparcl::KMeansSparseCluster.permute(y, K, nperms = 30, silent = T)$bestw
  mod4   <- sparcl::KMeansSparseCluster(y, K, wbounds = best3, nstart = 100, silent = T)
  loc4   <- lapply(1:K, \(k) colMeans(y[mod4[[1]]$Cs == k, ]))
  prob4  <- sapply(1:K, \(k) sum(mod4[[1]]$Cs == k)) / n
  covs4  <- lapply(1:K, \(k) var(y[mod4[[1]]$Cs == k,]) + adj)
  theta4 <- c(unlist(loc4), unlist(covs4), prob4)
  
  # 5. K-means
  mod5   <- stats::kmeans(y, K, iter.max = 100, nstart = 100)
  loc5   <- lapply(1:K, \(k) mod5$centers[k,])
  prob5  <- mod5$size / n
  covs5  <- lapply(1:K, \(k) var(y[mod5$cluster == k,]) + adj)
  theta5 <- c(unlist(loc5), unlist(covs5), prob5)
  
  # 6. Partitioning Around Medoids
  mod6   <- cluster::pam(y, K, nstart = 100)

  list("HTKmeans-Regularization-Path" = theta0,
       "HTKmeans-AIC" = theta1, "HTKmeans-BIC" = theta2,
       "GMM-mclust"   = theta3, "SparKmeans"   = theta4,
       "Kmeans"       = theta5)
}

# -----------------------------------------
# Initialization of model parameters
# -----------------------------------------

initialVals <- function(y, alphas, betas, lambdas, nus, psis,
                        HTKshrinkage = 0.75){
  
  n <- nrow(y)
  d <- ncol(y)
  K <- length(alphas)
  
  # Pool of candidate initial values
  vals <- initialPool(y, K, n, d, HTKshrinkage)
  
  # Maximize the log-posterior
  best.init <- which.max(
    sapply(vals, \(v) BOBgmms::log_post(v, y, alphas, betas, lambdas, psis, nus))
  )
  
  optimalVals <- vals[[best.init]]
  
  # Extract the mean vectors
  means <- split(optimalVals[1:(d * K)], unlist(lapply(1:K, \(k) rep(k, d))))
  names(means) <- NULL
  
  # Extract the covariance matrices
  covs <- optimalVals[(K * d + ((1 - 1) * d^(2)) + 1):(K * d + (K * d^(2)))]
  covs <- split(covs, unlist(lapply(1:K, \(k) rep(k, d * d))))
  covs <- lapply(covs, \(c) matrix(c, ncol = d, nrow = d))
  names(covs) <- NULL
  
  # Extract the probabilities
  probs <- tail(optimalVals, K)
  
  list(means = means, covs = covs, probs = probs, algorithm = names(best.init))
}

# -----------------------------------------
# Tune tempering hyper-parameters 
# -----------------------------------------

temp.tuning <- function(y, means, covs, probs, 
                        range.a, range.b, range.c, range.r, 
                        alphas, betas, lambdas, nus, psis, temp.cores){
  
  temper.grid <- expand.grid(range.a, range.b, range.c, range.r)
  colnames(temper.grid) <- c("a", "b", "c", "r")
  
  log.postTemp <- parallel::mclapply(asplit(temper.grid, 1), function(temp){
    
    theta.temp <- BOBgmms::em.optim(y, means, covs, probs, temp["a"], temp["b"], 
      temp["c"], temp["r"], rep(1, n), rep(1, K), rep(1, K), 1, 
        alphas, betas, lambdas, nus, psis)
    
    BOBgmms::log_post(theta.temp, y, alphas, betas, lambdas, psis, nus)
    
  }, mc.cores = temp.cores)
  
  as.list(asplit(temper.grid, 1)[[which.max(log.postTemp)]])
}


# -----------------------------------------
# HT Initialization of model parameters
# -----------------------------------------

initialVals.HT <- function(y, alphas, betas, lambdas, nus, psis, 
                           HTKshrinkage = 0.75){
  
  n <- nrow(y)
  d <- ncol(y)
  K <- length(alphas)
  
  # HT K-means -- Tuned via regularization path plot
  adj  <- diag(d) * 1e-10
  mod  <- clusterHD::HTKmeans(y, K)
  best <- which(abs(mod$lambdas - HTKshrinkage) < 1e-5)
  locs <- lapply(1:K, \(k) mod$HTKmeans.out[[best]]$centers[k,])
  prob <- sapply(1:K, \(k)sum(mod$HTKmeans.out[[best]]$cluster == k)) / n
  covs <- lapply(1:K, \(k) var(y[mod$HTKmeans.out[[best]]$cluster == k,]) + adj)
    
  list(means = locs, covs = covs, probs = prob, algorithm = "HT K-means")
}
