// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "PT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
List SurvPT(const arma::mat X, //Data (n by d matrix)
            const int prior,
            const arma::vec support, //Support domain (length 2 vector)
            const int max_M, //Maximum level of tree
            const arma::vec prior_params, //Prior parameters (only for Exponential and Weibull
            const int allocated //0: not allocated (MW) method, 1: allocated (PT) method
){
    //Initialize the class object
    PT my_tree(X, prior, support, max_M, prior_params, allocated);

    my_tree.re_allocation();

    List out = my_tree.output();

    my_tree.clear();

    return out;
}
