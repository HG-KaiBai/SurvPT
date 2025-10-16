// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

double newton_raphson(const double& s_null, const int& level, const double& theta) {
    if (theta == 0 || theta == 1) {
        return theta; // boundary case
    }
    double target =  R::digamma(pow(level, 2.0) * s_null) - R::digamma(pow(level, 2.0) * (1-s_null)) - log(theta) + log(1-theta);
    double s_new = s_null;
    int count = 0;
    while (count < 10) { // Limit iterations to 10, usually sufficient
        double f_old = R::digamma(pow(level, 2.0) * s_new) - R::digamma(pow(level, 2.0) * (1-s_new)) - log(theta) + log(1-theta);
        double f_deriv = pow(level, 2.0) * (R::trigamma(pow(level, 2.0) * s_new) + R::trigamma(pow(level, 2.0) * (1-s_new)));
        s_new -= f_old / f_deriv;
        double f_new = R::digamma(pow(level, 2.0) * s_new) - R::digamma(pow(level, 2.0) * (1-s_new)) - log(theta) + log(1-theta);
        if (abs(f_new) <= abs(target)) {
            return min(1.0, max(s_new, 0.0)); // ensure within [0, 1]
        }
        count++;
    }
    return s_new;
}
