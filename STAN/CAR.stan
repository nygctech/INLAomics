functions {
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m((alpha-1)/alpha * lambda[i]);
      return 0.5 * (n * log(tau) + n * log(alpha)
                    + sum(ldet_terms)
                    - tau * ((1-alpha)*(phit_D * phi - (phit_W * phi)) 
                    + alpha*dot_product(phi,phi)));
  }
}

data {
  int<lower = 1> n;
  int<lower = 0> y[n];
  matrix<lower = 0, upper = 1>[n, n] W; 
  int W_n; 
  int W_sparse[W_n, 2];
  vector[n] D_sparse;
  vector[n] lambda;
}

parameters {
  vector[n] phi;
  real theta;
  real nu;
  //real zipt;
}

transformed parameters{
  real<lower = 0> tau;
  real<lower = 0, upper = 1> alpha;
  //real<lower = 0, upper = 1> zip;
  tau = exp(theta);
  alpha = inv_logit(nu);
  //zip = inv_logit(zipt);
}

model {
  phi ~ sparse_car(tau, alpha, W_sparse, D_sparse, lambda, n, W_n);
  tau ~ chi_square(1);
  //zip ~ beta(1,1);
  alpha ~ uniform(0,1);
  //b ~ normal(0, 10);
  target += nu - 2 * log(1 + exp(nu)) + theta; //+ zipt - 2 * log(1 + exp(zipt));
  for (i in 1:n){
    /*
    if (y[i] == 0)
      target += log_sum_exp(bernoulli_lpmf(1|zip), bernoulli_lpmf(0|zip) + poisson_log_lpmf(y[i]|phi[i]+ dot_product(X[i,],b)));
    else{
      target += bernoulli_lpmf(0|zip) + poisson_log_lpmf(y[i]|phi[i]+ dot_product(X[i,],b));
    }
    */
    target += poisson_log_lpmf(y[i]|phi[i]);
  }
}
