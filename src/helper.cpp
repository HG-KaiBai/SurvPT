// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

double newton_raphson(const double& s_hat, const double& s_null, const int& level, const double& theta) {
    if (s_hat == 0 || s_hat == 1) {
        return s_hat; // boundary case
    }
    double s_old = s_null;
    double s_new = s_hat;
    double f_old = R::digamma(pow(level, 2.0) * s_old) - R::digamma(pow(level, 2.0) * (1-s_old)) - log(theta) + log(1-theta);
    int count = 0;
    while (count < 10) { // Limit iterations to 10, usually sufficient
        if (count == 9) {
            std::cout << "Final Count: " << count << ", s_old: " << s_old << ", s_new: " << s_new << std::endl; // Test if convergence is achieved
        }
        double f_new = R::digamma(pow(level, 2.0) * s_new) - R::digamma(pow(level, 2.0) * (1-s_new)) - log(theta) + log(1-theta);
        if (abs(f_new) <= abs(f_old) || abs(s_new - s_old) < 1e-6) {
            return s_new;
        }
        double f_deriv = pow(level, 2.0) * (R::trigamma(pow(level, 2.0) * s_new) + R::trigamma(pow(level, 2.0) * (1-s_new)));
        s_new -= f_new / f_deriv;
        count++;
    }
    return s_new;
}