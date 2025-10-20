// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec sample_from_inverse_kron(const arma::mat& Lambda, const arma::mat& Q) {
  arma::mat KronProd = arma::kron(Lambda, Q);
  arma::mat Sigma = arma::inv(KronProd);
  int dim = Sigma.n_rows;
  arma::vec z = arma::randn<arma::vec>(dim);
  arma::mat L = arma::chol(Sigma, "lower");
  return L * z;
}

// [[Rcpp::export]]
arma::vec sample_from_inverse( const arma::mat& Q, const arma::vec& mu) {
  arma::mat Sigma = arma::inv(Q);
  int dim = Sigma.n_rows;
  arma::vec z = arma::randn<arma::vec>(dim);
  arma::mat L = arma::chol(Sigma, "lower");
  return L * z + mu;
}