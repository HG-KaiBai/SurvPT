// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "PT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
List SurvPT(const arma::mat X, //Data (n by d matrix)
            const arma::vec test_seq, //Test sequence (length m vector)
            const int prior,
            const arma::vec support, //Support domain (length 2 vector)
            const int max_M, //Maximum level of tree
            const arma::vec prior_params, //Prior parameters (only for Exponential and Weibull
            const int allocated, //0: not allocated (MW) method, 1: allocated (PT) method
            const int max_EM_iter, //Maximum number of EM iterations
            const double EM_tol, //Tolerance for EM convergence
            const bool verbose, //Verbose output
            Rcpp::Nullable<Rcpp::NumericVector> entry = R_NilValue //Left-truncation times
){
    arma::vec entry_times(X.n_rows, arma::fill::zeros);
    if (entry.isNotNull()) {
        const arma::vec supplied = Rcpp::as<arma::vec>(entry.get());
        if (supplied.n_elem == 1) entry_times.fill(supplied(0));
        else entry_times = supplied;
    }
    //Initialize the class object
    PT my_tree(X, test_seq, prior, support, max_M, prior_params, allocated, max_EM_iter, EM_tol, verbose, entry_times);

    if (allocated == 1){
        my_tree.EM_step();
    }

    List out = my_tree.output();

    my_tree.clear();

    return out;
}
