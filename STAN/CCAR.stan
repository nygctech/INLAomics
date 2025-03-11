functions {
  real sparse_car_lpdf(vector phi, vector mu, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
      vector[n] loc = (phi-mu);
    
      phit_D = (loc .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + loc[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + loc[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m((alpha-1)/alpha * lambda[i]);

      return 0.5 * (n * log(tau) + n * log(alpha)
                    + sum(ldet_terms)
                    - tau * ((1-alpha)*(phit_D * loc - (phit_W * loc)) 
                    + alpha*dot_product(loc,loc)));   
  }
}
data {
  int<lower = 1> n;
  int<lower = 1> N;
  int<lower = 0> y[N];
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
  int W_n;                // number of adjacent region pairs
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[n] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[n] lambda;       // eigenvalues of D-W
}

parameters {
  vector[n] phi1;
  vector[n] phi2;
  real theta1;
  real theta2;
  real eta0;
  real eta1;
  real nu1;
  real nu2;
}


transformed parameters{
  real<lower = 0> tau1;
  real<lower = 0> tau2;
  real<lower = 0, upper = 1> alpha1;
  real<lower = 0, upper = 1> alpha2;
  tau1 = exp(theta1);
  tau2 = exp(theta2);
  alpha1 = inv_logit(nu1);
  alpha2 = inv_logit(nu2);
}

model {
  phi1 ~ sparse_car(rep_vector(0, n), tau1, alpha1, W_sparse, D_sparse, lambda, n, W_n);
  phi2 ~ sparse_car((eta0 * identity_matrix(n) + eta1 * W)*phi1, tau2, alpha2, W_sparse, D_sparse, lambda, n, W_n);
  tau1 ~ chi_square(1);
  tau2 ~ chi_square(1);
  alpha1 ~ uniform(0,1);
  alpha2 ~ uniform(0,1);
  eta0 ~ normal(0, 10);
  eta1 ~ normal(0, 10);
  target += nu1 - 2 * log(1 + exp(nu1)) + nu2 - 2 * log(1 + exp(nu2)) + theta1 + theta2;
  for (i in 1:N){
    target += poisson_log_lpmf(y[i]|i <= n ? phi1[i] : phi2[i-n]);
  }
}
