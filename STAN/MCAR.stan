functions {
  real logit2(real x) {
    return logit((x+1)/2);
  }

  real inv_logit2(real x) {
    return (exp(x)-1)/(exp(x)+1);
  }

  real sparse_car_lpdf(vector phi, real tau1, real tau2, real rho, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int N, int W_n) {
      row_vector[n] phit_D_1; // phi_1' * D
      row_vector[n] phit_W_1; // phi_1' * W
      row_vector[n] phit_D_2; // phi_2' * D
      row_vector[n] phit_W_2; // phi_2' * W
      vector[n] ldet_terms;
      vector[n] phi_1;
      vector[n] phi_2;
      vector[n] ldet_terms;
      real ldet_lambda;

      phi_1 = phi[1:n];
      phi_2 = phi[(n+1):N];
      phit_D_1 = (phi_1 .* D_sparse)';
      phit_D_2 = (phi_2 .* D_sparse)';
      phit_W_1 = rep_row_vector(0, n);
      phit_W_2 = rep_row_vector(0, n);

      for (i in 1:W_n) {
        phit_W_1[W_sparse[i, 1]] = phit_W_1[W_sparse[i, 1]] + phi_1[W_sparse[i, 2]];
        phit_W_1[W_sparse[i, 2]] = phit_W_1[W_sparse[i, 2]] + phi_1[W_sparse[i, 1]];
        phit_W_2[W_sparse[i, 1]] = phit_W_2[W_sparse[i, 1]] + phi_2[W_sparse[i, 2]];
        phit_W_2[W_sparse[i, 2]] = phit_W_2[W_sparse[i, 2]] + phi_2[W_sparse[i, 1]];
      }
      for (i in 1:n) ldet_terms[i] = log1m((alpha-1)/alpha * lambda[i]);

      ldet_lambda = tau1*tau2 / (1-square(rho));
      return 0.5 * (2*sum(ldet_terms) + n*log(ldet_lambda)
                    - 1/(1-square(rho)) * (tau1*(phit_D_1 * phi_1 - alpha * phit_W_1 * phi_1)
                    + tau2*(phit_D_2 * phi_2 - alpha * phit_W_2 * phi_2)
                    -2*rho*sqrt(tau1*tau2)*(phit_D_1 * phi_2 - alpha * phit_W_1 * phi_2)));
  }
}

data {
  int<lower = 1> n;
  int<lower = 0> y[n];
  int<lower=1> p;
  matrix[n, p] X;
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
  int W_n;                // number of adjacent region pairs
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[n] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[n] lambda;       // eigenvalues of D-W          // number of adjacent region pairs
}

parameters {
  vector[N] phi;
  real theta1;
  real theta2;
  real nu;
  real zipt1;
  real zipt2;
}

transformed parameters{
  real<lower = 0> tau;
  real<lower = 0, upper = 1> alpha;
  real<lower = 0, upper = 1> zip1;
  real<lower = 0, upper = 1> zip2;
  tau = exp(theta);
  alpha = inv_logit2(nu);
  zip1 = inv_logit(zipt1);
  zip2 = inv_logit(zipt2);
}

model {
  phi ~ sparse_car(tau1, tau2, rho, alpha, W_sparse, D_sparse, lambda, n, N, W_n);
  tau1 ~ gamma(2, 2);
  tau2 ~ gamma(2, 2);
  if (zip) {
    theta1 ~ beta(1,2);
    theta2 ~ beta(1,2);
  }
  if (zip){
    for (i in 1:N){
      if (y[i] == 0)
        target += log_sum_exp(bernoulli_lpmf(1|i <= n ? theta1 : theta2),
                              bernoulli_lpmf(0|i <= n ? theta1 : theta2)
                              +poisson_log_lpmf(y[i]|phi[i]));
      else{
        target += bernoulli_lpmf(0|i <= n ? theta1 : theta2)
                    + poisson_log_lpmf(y[i]|phi[i]);
      }
    }
  } else {
    y ~ poisson_log(phi);
  }
}

