#ifndef PT_H
#define PT_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

class PT {

public:
    mat X; // data matrix
    int prior; // prior type (1: Uniform, 2: Exponential, 3: Weibull)
    vec support; // support domain
    int max_M; // maximum level of tree
    int NODES_NUM; // total number of nodes in the tree
    vec prior_params; // prior parameters (only for Exponential and Weibull)
    int n; // total sample size
    int allocated; // 0: not allocated (MW) method, 1: allocated (PT) method
    double largest_obs; // largest observation in the data

    // constructor
    PT(const mat X, const int prior, const vec support, const int max_M, const vec prior_params, const int allocated);

    // functions
    void init(); // initialize the first tree structure
    void clear(); // clear all member variables
    double choose_split(int level, double left, double right); // choose split points for each node
    double init_alpha(int level, double left, double right); // update concentration parameters alpha
    void init_allocations(int parent, int left_child, int right_child); // update number of observations allocated to each node
    double init_individual_allocation(double left_obs, double right_obs, double left_node, double right_node); // individual allocation based on prior
    void re_allocation(); // update allocations by the previous PT distribution
    void update_allocations(int parent, int left_child, int right_child); // update allocations for children nodes
    double update_individual_allocation(double left_obs, double right_obs, double left_node, double right_node); // individual allocation based on posterior
    List output(); // output results

private:
    vec left_bound;
    vec right_bound; // left and right bounds of each node
    vec node_levels; // levels of each node in the tree
    vec alpha; // concentration parameters for beta distributions
    mat in_node_obs; // observations idx allocated to each node
    vec n_obs_allocated; // number of observations allocated to each node
    vec posterior_probs; // posterior probabilities of each node
    vec terminal_nodes; // indices of terminal nodes
    vec new_posterior_probs; // new posterior probabilities after update
};

#endif
