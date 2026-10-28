// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "OPT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
List OPT_right(const arma::vec X, //Data (n by d matrix)
            const arma::vec delta, //Censoring indicator (length n vector)
            const arma::vec test_seq, //Test sequence (length m vector)
            const double Rho, //Initial parameters (length p vector)
            const arma::vec support, //Support domain (length 2 vector)
            const int max_M, //Maximum level of tree
            const int max_iter, //Maximum number of EM iterations
            const double iter_tol, //Tolerance for EM convergence
            const bool verbose //Verbose output
){
    
    //Initialize the class object
    OPT my_tree(X, delta, test_seq, Rho, support, max_M, max_iter, iter_tol, verbose);

    my_tree.iter_step();
    
    List out = my_tree.output();

    my_tree.clear();

    return out;
}