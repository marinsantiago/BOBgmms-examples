shhh <- suppressPackageStartupMessages
shhh(library(LaplacesDemon, include.only = c("rinvwishart", "rdirichlet")))
shhh(library(mvnfast, include.only = "rmvt"))

bayesianPosterior <- function(y, true.Z, S, alphas, betas, lambdas, nus, psis){
  
  nks      <- apply(true.Z, 2, sum)
  clusters <- lapply(1:K, \(k) y[true.Z[,k] == 1,])
  yBars    <- lapply(clusters, colMeans)
  SSRs     <- lapply(1:K, \(k) var(clusters[[k]]) * (nks[k] - 1))
  
  # Posterior hyper-parameters
  aHats    <- alphas + nks
  nuHats   <- nus + nks
  lambHats <- lambdas + nks
  betaHats <- lapply(1:K, \(k) 
                     lambdas[k]/lambHats[k] * betas[[k]] + 
                     nks[k]/lambHats[k] * yBars[[k]])
  psiHats  <- lapply(1:K, \(k) 
                     psis[[k]] + SSRs[[k]] + ((lambdas[k] * nks[k])/lambHats[k] *
                     (yBars[[k]] - betas[[k]]) %*% t((yBars[[k]] - betas[[k]]))))
  
  # Posterior parameters
  t.dfs    <- nuHats - d + 1
  t.scales <- lapply(1:K, \(k) psiHats[[k]]/(lambHats[k] * (t.dfs[k])))
  
  postOneShot <- function(t.scales, t.dfs, betaHats, nuHats, psiHats, aHats){
    
    mus.post    <- lapply(1:K, \(k) c(
      mvnfast::rmvt(1, betaHats[[k]], t.scales[[k]], t.dfs[k])))
    
    sigmas.post <- lapply(1:K, \(k) as.vector(
      LaplacesDemon::rinvwishart(nuHats[[k]], psiHats[[k]])))
    
    probs.post  <- c(LaplacesDemon::rdirichlet(1, aHats))
    
    c(do.call(c, mus.post), do.call(c, sigmas.post), probs.post)
    
  }
  
  out <- lapply(1:S, \(s) 
                postOneShot(t.scales, t.dfs, betaHats, nuHats, psiHats, aHats))
  
  matrix(unlist(out), nrow = S, byrow = T)
  
}
