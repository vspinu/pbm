// #include <tr1/random>
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
// using std::vector;

//[[Rcpp::export]]
IntegerVector rdunif(const IntegerVector min, const IntegerVector max){
  // No need for GetRNGstate() and PutRNGstate(); rcpp does it for us in higher
  // order functions.

  int N = min.length();
  if( N != max.length() ) Rf_error("lenghts of min and max vectors must be the same");

  IntegerVector out(N);
  for( int i = 0; i < N; i++){
    int m = min[i];
    int n = max[i] - m + 1;
    if(n < 0) Rf_error("max[i] - min[i] is less than 1 for i=%d", i + 1);
    out[i] = int(unif_rand() * n) + m;
  }
  return(out);
}

//[[Rcpp::export]]
NumericVector ddunif(const IntegerVector min, const IntegerVector max){
  // No need for GetRNGstate() and PutRNGstate(); rcpp does it for us in higher
  // order functions.

  int N = min.length();
  if( N != max.length() ) Rf_error("lenghts of min and max vectors must be the same");

  NumericVector out(N);
  for( int i = 0; i < N; i++){
    out[i] = 1.0/(max[i] - min[i] + 1);
  }
  return(out);
}
