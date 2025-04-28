#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat boot_stat(const arma::mat& mat,
                    const arma::umat& boot_indices,
                    std::string method) {
  int B       = boot_indices.n_cols;  //
  int n_cols  = mat.n_cols;          //
  arma::mat boot_results(B, n_cols); //

  //
  int n_threads = omp_get_max_threads();

#pragma omp parallel for num_threads(n_threads)
  for (int b = 0; b < B; b++) {
    arma::uvec idx = boot_indices.col(b);
    arma::mat mat_sample = mat.rows(idx);
    arma::rowvec stats(n_cols);

    if (method == "median") {
      for (int j = 0; j < n_cols; j++) {
        arma::vec col = mat_sample.col(j);
        arma::vec good = col.elem(find_finite(col));
        stats(j) = (good.n_elem>0 ? arma::median(good) : arma::datum::nan);
      }
    }
    else if (method == "mean") {
      for (int j = 0; j < n_cols; j++) {
        arma::vec col = mat_sample.col(j);
        arma::vec good = col.elem(find_finite(col));
        stats(j) = (good.n_elem>0 ? arma::mean(good) : arma::datum::nan);
      }
    }
    else if (method == "sum") {
      for (int j = 0; j < n_cols; j++) {
        arma::vec col = mat_sample.col(j);
        arma::vec good = col.elem(find_finite(col));
        stats(j) = (good.n_elem>0 ? arma::sum(good) : arma::datum::nan);
      }
    }
    else {
      Rcpp::stop("Unknown method ' %s '. Use 'mean','median' or 'sum'", method);
    }

    boot_results.row(b) = stats;
  }

  return boot_results;
}
