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
    int prior; // 1: Uniform, 2: Exponential, 3: Weibull, 4: Uniform/exponential splice
    vec support; // support domain
    int max_M; // maximum level of tree
    int NODES_NUM; // total number of nodes in the tree
    vec prior_params; // base parameters (rate for priors 2/4; rate and shape for 3)
    int n; // total sample size
    int allocated; // 0: not allocated (MW) method, 1: allocated (PT) method
    double largest_obs; // largest observation in the data
    double uniform_tail_factor; // multiplier for the uniform finite tail cutoff
    int max_EM_iter; // maximum number of EM iterations
    double EM_tol; // tolerance for EM convergence
    bool verbose; // print log info if true

    // constructor
    PT(const mat X, const vec test_seq, const int prior, const vec support, const int max_M, const vec prior_params, const int allocated, const int max_EM_iter, const double EM_tol, const bool verbose, const double uniform_tail_factor);

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
    void update_tree_size(); // double the size of the tree when necessary
    void EM_step(); // perform EM algorithm to update allocations and posterior probabilities
    List output(); // output results

private:
    void prepare_em_cache();
    void update_observation_probabilities();
    double cached_individual_allocation(arma::uword obs, int node) const;
    bool allocation_cache_ready = false;
    arma::uvec terminal_indices;
    std::vector<std::pair<int, int>> child_nodes;
    // Leaf indices remain in the original terminal_nodes order.
    std::vector<std::vector<arma::uword>> node_terminal_indices;
    mat observation_leaf_probabilities;
    vec observation_probabilities;
    double terminal_right_bound(double right) const;
    double base_survival(double x) const;
    double base_density(double x) const;
    double hybrid_cdf(double x) const;
    double hybrid_quantile(double probability, bool lower_tail) const;
    double hybrid_probability(double left, double right) const;
    double hybrid_body_mass = 0.0;
    double hybrid_tail_mass = 0.0;
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
