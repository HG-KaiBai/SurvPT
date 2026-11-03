// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "PT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export(rng = false)]]
List SurvPT(const arma::mat X, //Data (n by d matrix)
            const int prior,
            const arma::vec support, //Support domain (length 2 vector)
            const int max_M, //Maximum level of tree
            const arma::vec prior_params, //Base parameters: rate for priors 2/4, rate and shape for 3
            const int allocated, //0: not allocated (MW) method, 1: allocated (PT) method
            const int max_EM_iter, //Maximum number of EM iterations
            const double EM_tol, //Tolerance for EM convergence
            const bool verbose, //Verbose output
            const double uniform_tail_factor = 1.1, //Uniform tail cutoff multiplier
            const bool flat = false //Use unit concentration parameters at every depth
){
    //Initialize the class object
    PT my_tree(X, prior, support, max_M, prior_params, allocated, max_EM_iter, EM_tol, verbose, uniform_tail_factor, flat);

    if (allocated == 1){
        my_tree.EM_step();
    }

    List out = my_tree.output();

    my_tree.clear();

    return out;
}

// [[Rcpp::export]]
List predict_survPT_cpp(const List model, const arma::vec x,
                       const int n_samples, const int burn_in, const int thin,
                       const bool density) {
    PT fitted(model);
    return fitted.predict(x, n_samples, burn_in, thin, density);
}
