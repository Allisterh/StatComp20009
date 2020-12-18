#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk chain generator using Rcpp
//' @description A random walk chain generator using Rcpp
//' @param sigma The value of variance
//' @param x0 the start point of the random walk chain
//' @param N repeat times
//' @return A list containing the random walk chain and the number
//' of accepeted times
//' @examples
//' \dontrun{
//' rw <- rwRcpp(2,15,10000)
//' rw$k
//' plot(rw$x,type='l')
//' }
//' @export
// [[Rcpp::export]]
List rwRcpp(double sigma, double x0, int N) {
     NumericVector x(N);
     x[0]=x0;
     auto u = runif(N);
     int k = 0;
     for (int i = 1;i < N; i++) {
         auto y = rnorm(1, x[i - 1], sigma)[0];
         if (u[i] <= (exp(-abs(y))/exp(-abs(x[i-1])))) {
             x(i) = y;
         } else {
             x[i] = x[i-1];
             k++;
         }
     }
     return List::create(Named("x") = x, Named("k") = k);
 }