# ------------------------------------------------------------------------------
# Tune hyper-parameters
# ------------------------------------------------------------------------------

shhh <- \(o) o |> suppressPackageStartupMessages() |> suppressWarnings()
shhh(library(pbmcapply)); shhh(library(parallel))

library(BOBgmms)

source("~/Desktop/BOBgmms-examples-Eigen/R/init.R")

# 10-fold cross-validation to select shrinkage parameters (in parallel)
shrinkage.cv.tuning <- function(y, betas, psis, alphas,
                                cv.folds = 10L,
                                cores = parallel::detectCores(),
                                seed = sample.int(.Machine$integer.max, 1)) {
  
  set.seed(seed)
  # Pre-compute constants
  dims_y <- dim(y)
  n <- dims_y[1]
  d <- dims_y[2]
  K <- length(alphas)
  
  # Grid for lambda and nu
  lambdas.seq <- seq(0.1, 10, length.out = 10) * (d / n)
  nus.seq <- lambdas.seq + (d + 1) # Make sure the posterior mean exists
  grid. <- expand.grid(lambdas.seq, nus.seq) |> asplit(x = _, MARGIN = 1)
  
  # Data structure to store output. Each row is a cross validation fold
  test.logliks <- matrix(NA, nrow = cv.folds, ncol = 10 * 10)
  
  # Cross-validation 
  folds <- sample(rep(seq_len(cv.folds), each = n / cv.folds))
  for (fold. in seq_len(cv.folds)) {
    cat("\rFold number: ", fold.)
    
    # Split the data into train and held-out sets-------------------------------
    train.idx <- folds != fold.
    test.idx <- folds == fold.
    y.train <- y[train.idx, ]
    y.test <- y[test.idx, ]
    n.train <- sum(train.idx)
    n.test <- sum(test.idx)
    
    # Initialize the model parameters ------------------------------------------
    init.vals <- initial.values(y.train, K)
    means.init <- init.vals$values$means
    covs.init <- init.vals$values$covs
    probs.init <- init.vals$values$probs
    
    # EM optimization over each (lambda, nu) pair and evaluation of the 
    # held-out log-likelihood --------------------------------------------------
    test.logliks[fold., ] <- parallel::mclapply(
      grid., \(g) {
        # EM algorithm
        current.lambdas <- rep(g[1], K) 
        current.nus <- rep(g[2], K)
        current.params <- BOBgmms::em.optim(
          y = y.train, means.init = means.init, covs.init = covs.init, 
          probs.init = probs.init, betas = betas, 
          lambdas = current.lambdas, nus = current.nus,
          psis = psis, alphas = alphas
        )
        current.means <- current.params$means
        current.covs <- lapply(current.params$covs, \(c) round(c, 5))
        current.probs <- current.params$probs
        # Held-out log-likelihood
        BOBgmms::loglik(
          means = current.means, covs = current.covs, 
          probs = current.probs, y = y.test
        )
      },
      mc.cores = cores
    ) |> unlist()
  }
  
  # Optimal hyper-parameters
  best.idx <- which.max(colMeans(test.logliks))
  best.lambda <- grid.[[best.idx]][1]; names(best.lambda) <- NULL
  best.nu <- grid.[[best.idx]][2]; names(best.nu) <- NULL
  gc() # Collect garbage
  
  list(best.lambda = best.lambda, best.nu = best.nu)
}


# ------------------------------------------------------------------------------
# Tune tempering hyper-parameters 
# ------------------------------------------------------------------------------


temperature.tuning <- function(y, means.init, covs.init, probs.init,
                               range.a, range.b, range.c, range.r,
                               betas, lambdas, nus, psis, alphas,
                               cores = parallel::detectCores()) {
  
  # Define the grid of possible values
  temp.grid <- expand.grid(
    a = range.a, b = range.b, c = range.c, r = range.r
  ) |> asplit(x =, MARGIN = 1)
  # Iterate over the grid
  log.posts <- pbmcapply::pbmclapply(
    temp.grid, \(temp) {
      current_a <- temp["a"]
      current_b <- temp["b"]
      current_c <- temp["c"]
      current_r <- temp["r"]
      # EM optimization
      current.params <- BOBgmms::em.optim(
        y = y, means.init = means.init, covs.init = covs.init,
        probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
        psis = psis, alphas = alphas, lik.weights = rep(1, nrow(y)),
        means.weights = rep(1, length(means.init)), 
        covs.weights = rep(1, length(means.init)),
        probs.weights = rep(1, length(means.init)), 
        a = current_a, b = current_b, c = current_c, r = current_r
      )
      # Log-posterior
      BOBgmms::lpost(
        means = means.init, covs = covs.init, probs = probs.init,
        y = y, betas = betas, lambdas = lambdas, nus = nus,
        psis = psis, alphas = alphas
      )
    },
    mc.cores = cores
  ) |> unlist()
  
  best.idx <- which.max(log.posts)
  best.a <- temp.grid[[best.idx]]["a"]; names(best.a) <- NULL
  best.b <- temp.grid[[best.idx]]["b"]; names(best.b) <- NULL
  best.c <- temp.grid[[best.idx]]["c"]; names(best.c) <- NULL
  best.r <- temp.grid[[best.idx]]["r"]; names(best.r) <- NULL
  
  list(best.a = best.a, best.b = best.b, best.c = best.c, best.r = best.r)
}
