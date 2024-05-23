
// Gaussian Mixture Model - Stan

data {
  int<lower=1> K;              // number of mixture components
  int<lower=1> n;              // number of data points
  int<lower=1> p;              // number of variables
  matrix[n, p] Y;              // data matrix
  vector[p] beta[K];           // prior means
  vector[K] lambda;            // lambda hyper-parameters
  vector[K] nu;                // nu hyper-parameters
  cov_matrix[p] Psi[K];        // Psi matrices hyper-parameters
  vector[K] alpha;             // alpha hyper-parameters
}

parameters {
  simplex[K] theta;            // mixing proportions
  vector[p] mu[K];             // mean vectors of mixture components
  cov_matrix[p] Sigma[K];      // covariance matrices of mixture components
}

model {
  vector[K] log_theta = log(theta);  // cache log calculation
  for (k in 1:K){
    
    // Normal prior for mean vectors
    mu[k] ~ multi_normal(beta[k], Sigma[k]/lambda[k]);
    
    // Inverse Wishart prior for cov matrices
    Sigma[k] ~ inv_wishart(nu[k], Psi[k]);
    
  }
  // Dirichlet prior for mixing probabilities
  theta ~ dirichlet(alpha);
  
  // Sampling distribution
  for (i in 1:n) {
    
    vector[K] lps = log_theta;
    for (k in 1:K) {
      
      lps[k] += multi_normal_lpdf(Y[i] | mu[k], Sigma[k]);
      
    }
    target += log_sum_exp(lps);
    
  }
}
