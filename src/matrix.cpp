// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP matMult(const Eigen::Map<Eigen::MatrixXd>& A,
             const Eigen::Map<Eigen::MatrixXd>& B) {
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP matInv(const Eigen::Map<Eigen::MatrixXd>& A) {
  Eigen::MatrixXd A_inv = A.inverse();
  return Rcpp::wrap(A_inv);
}

// [[Rcpp::export]]
SEXP matSolve(const Eigen::Map<Eigen::MatrixXd>& A,
              const Eigen::Map<Eigen::MatrixXd>& B) {
  Eigen::MatrixXd X = A.colPivHouseholderQr().solve(B);
  return Rcpp::wrap(X);
}
