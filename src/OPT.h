#ifndef OPT_H
#define OPT_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

class OPT {

public:
    vec X; // data time
    vec delta; // data censoring indicator
    vec test_seq; // test sequence
    vec prior_rho; // probability of optional split
    vec support; // support domain
    int max_M; // maximum level of tree
    int NODES_NUM; // total number of nodes in the tree
    int n; // total sample size
    double largest_obs; // largest observation in the data
    int max_EM_iter; // maximum number of EM iterations
    double EM_tol; // tolerance for EM convergence
    bool verbose; // print log info if true

    // constructor
    OPT(const vec X, const vec delta, const vec test_seq, const double Rho, const vec support, const int max_M, const int max_EM_iter, const double EM_tol, const bool verbose);

    // functions
    void init(); // initialize the first tree structure
    void clear(); // clear all member variables
    double choose_split(int level, double left, double right); // choose split points for each node
    void init_allocations(int parent, int left_child, int right_child); // update number of observations allocated to each node
    void conditional_prob_WKM(const arma::vec& Y,
                             const arma::vec& Delta,
                             const arma::vec& w,
                             const double a,
                             const double s,
                             const double b); // compute conditional probabilities based on WKM
    double WKM_allocation(double t); // compute WKM-based allocation for a single observation
    void pos_rho(); // update rho values to ensure they are in (0,1)
    void pos_probs(); // initialize posterior probabilities
    void iter_step(); // perform iterations until convergence
    void re_allocation(); // re-allocate observations and update posterior probabilities
    void update_allocations(int parent, int left_child, int right_child); // allocate observations to children nodes
    double update_individual_allocation(double left_obs, double left_node, double right_node); // update individual allocation based on terminal nodes and previous posterior probabilities
    void update_tree_size(); // double the size of the tree when necessary
    List output(); // output results

private:
    int nodes_used; // current idx of the last used node
    vec left_bound;
    vec right_bound; // left and right bounds of each node
    vec node_levels; // levels of each node in the tree
    vec rho; // probability of optional split at each node
    vec alpha; // concentration parameters for beta distributions
    mat in_node_obs; // observations idx allocated to each node
    vec node_km; // Kaplan-Meier estimates at each node
    vec node_times; // unique event times at each node
    vec n_obs_allocated; // number of observations allocated to each node
    vec posterior_probs; // posterior probabilities of each node
    vec terminal_nodes; // indices of terminal nodes
    vec new_posterior_probs; // new posterior probabilities after update
    vec parent_nodes; // indices of parent nodes
    vec survival_probs; // survival probabilities at test_seq
};

#endif
