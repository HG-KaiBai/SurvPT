// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include "PT.h"

// [[Rcpp::export(rng = false)]]
Rcpp::List SurvPT(const arma::mat X, const int prior, const int max_M, const arma::vec prior_params,
                 const int allocated, const int max_EM_iter, const double EM_tol,
                 const bool verbose, const double gamma = 0.5, const double kappa = 0.1,
                 const double s0 = -1, const double min_count = 1,
                 const double split_epsilon = 1e-6, const bool flat = false) {
    PT tree(X,prior,max_M,prior_params,allocated,max_EM_iter,EM_tol,
            verbose,gamma,kappa,s0,min_count,split_epsilon,flat);
    tree.fit();
    Rcpp::List result=tree.output();
    result["model_state"]=tree.state(result);
    return result;
}

// [[Rcpp::export]]
Rcpp::List SurvPTSample(const Rcpp::List state,const int n_samples,
                       const int burn_in,const int thin,const bool verbose) {
    if (n_samples<=0 || burn_in<0 || thin<1)
        Rcpp::stop("Require n_samples > 0, burn_in >= 0 and thin >= 1.");
    PT tree(state,verbose);
    return tree.sample(n_samples,burn_in,thin);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List SurvPTPredict(const Rcpp::List state,const arma::vec times,
                        const arma::mat leaf_draws) {
    PT tree(state,false);
    return tree.predict_at(times,leaf_draws);
}
