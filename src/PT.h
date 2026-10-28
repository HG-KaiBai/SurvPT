#ifndef PT_H
#define PT_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <utility>
#include <vector>

using namespace Rcpp;
using namespace arma;
using namespace std;

class PT {

public:
    mat X; // data matrix
    vec test_seq; // test sequence
    int prior; // prior type (1: Uniform, 2: Exponential, 3: Weibull)
    vec support; // support domain
    int max_M; // maximum level of tree
    int NODES_NUM; // total number of nodes in the tree
    vec prior_params; // prior parameters (only for Exponential and Weibull)
    int n; // total sample size
    int allocated; // 0: not allocated (MW) method, 1: allocated (PT) method
    double largest_obs; // largest observation in the data
    int max_EM_iter; // maximum number of EM iterations
    double EM_tol; // tolerance for EM convergence
    bool verbose; // print log info if true

    // constructor
    PT(const mat X, const vec test_seq, const int prior, const vec support, const int max_M, const vec prior_params, const int allocated, const int max_EM_iter, const double EM_tol, const bool verbose);

    // functions
    void init(); // initialize the first tree structure
    void clear(); // clear all member variables
    double choose_split(int level, double left, double right); // choose split points for each node
    double init_alpha(int level, double left, double right); // update concentration parameters alpha
    void init_allocations(int parent, int left_child, int right_child); // update number of observations allocated to each node
    double init_individual_allocation(double left_obs, double right_obs, double left_node, double right_node); // individual allocation based on prior
    void re_allocation(); // update allocations by the previous PT distribution
    void update_tree_size(); // double the size of the tree when necessary
    void EM_step(); // perform EM algorithm to update allocations and posterior probabilities
    List output(); // output results

private:
    bool has_tail; // an unresolved tail is the right child (index 2) of the root
    void prepare_observation_cache();
    void update_expected_allocations();
    // Nonzero (leaf index, conditional observation probability/density) entries.
    std::vector<std::vector<std::pair<int, double>>> observation_kernels;
    std::vector<std::pair<int, int>> child_nodes;
    double log_likelihood() const;
    vec loglik_trace;
    int em_iterations = 0;
    bool converged = false;
    double base_survival(double x) const;
    double base_density(double x) const;
    double terminal_probability(double left, double right, double a, double b) const;
    double terminal_density(double t, double a, double b) const;
    int nodes_used; // current idx of the last used node
    vec left_bound;
    vec right_bound; // left and right bounds of each node
    vec node_levels; // levels of each node in the tree
    vec alpha; // concentration parameters for beta distributions
    mat in_node_obs; // observations idx allocated to each node
    vec n_obs_allocated; // number of observations allocated to each node
    vec posterior_probs; // posterior probabilities of each node
    vec terminal_nodes; // indices of terminal nodes
    vec new_posterior_probs; // new posterior probabilities after update
    vec parent_nodes; // indices of parent nodes
    vec survival_probs; // survival probabilities at test_seq
};

#endif
