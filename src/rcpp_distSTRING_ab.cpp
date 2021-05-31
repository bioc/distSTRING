#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

//' @useDynLib distSTRING, .registration = TRUE
//' @import Rcpp
//' @export rcpp_distSTRING_ab
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_distSTRING_ab( std::string a, std::string b, Rcpp::NumericMatrix scoreMatrix, int nsites ) {
  std::unordered_map<std::string, double> dist_mat;
  int nCols = scoreMatrix.ncol();
  int nRows = scoreMatrix.nrow();
  CharacterVector Colnames = colnames(scoreMatrix);
  CharacterVector Rownames = rownames(scoreMatrix);
  for( int is = 0; is < nCols; is++ ){
    for( int js = 0; js < nRows; js++){
      std::string isName = "";
      isName = Colnames(is);
      std::string jsName = "";
      jsName = Rownames(js);
      dist_mat[isName+jsName] = scoreMatrix(is,js);
    }
  }
  double eqnum = 0;
  int ab_n = nsites;
  for( int s=0; s < nsites; s++){
    std::string as;
    std::string bs;
    as = a[s];
    bs = b[s];
    double ab_dist;
    ab_dist = dist_mat[as+bs];
    if(ab_dist >= 0.0){
      eqnum = eqnum + ab_dist;
    } else {
      ab_n = ab_n -1;
    };
  }
  return Rcpp::NumericVector::create(Rcpp::Named("distSTRING") = eqnum, Rcpp::Named("sitesUsed") = ab_n);
}
