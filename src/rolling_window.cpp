#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rolling_window_sum(NumericVector positions, NumericVector expressions,
                                 double window, String type) {
  int n = positions.size();
  NumericVector result(n);

  for(int i = 0; i < n; i++) {
    double pos_i = positions[i];
    double sum = 0.0;

    if(type == "sum1") {
      // window 1: [pos_i, pos_i + window]
      for(int j = 0; j < n; j++) {
        if(positions[j] >= pos_i && positions[j] <= pos_i + window) {
          sum += expressions[j];
        }
      }
    } else {
      // window 2: [pos_i + window/2, pos_i + 1.5*window]
      for(int j = 0; j < n; j++) {
        if(positions[j] >= pos_i + window/2 && positions[j] <= pos_i + 1.5*window) {
          sum += expressions[j];
        }
      }
    }

    result[i] = sum;
  }

  return result;
}
