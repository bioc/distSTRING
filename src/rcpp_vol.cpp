#include <Rcpp.h>
#include <string.h>
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
using namespace Rcpp;

//' @useDynLib distSTRING, .registration = TRUE
//' @import Rcpp
//' @import RcppThread
//' @export rcpp_vol
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_vol( Rcpp::NumericMatrix inputMatrix, int ncores = 1 ) {
  int nRows = inputMatrix.nrow();
  Rcpp::NumericVector v (nRows);
  RcppThread::parallelFor(0, nRows, [&] (int i) {
    double vol = 0;
    double x = inputMatrix(i,0);
    double y = inputMatrix(i,1);
    double z = inputMatrix(i,2);
    vol = x * y * z;
    v(i) = vol;
  }, ncores);
  return Rcpp::NumericVector(v);
}
