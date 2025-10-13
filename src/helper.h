#ifndef HELPER_H
#define HELPER_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

double newton_raphson(const double& s_null, const int& level, const double& theta);

#endif // HELPER_H

