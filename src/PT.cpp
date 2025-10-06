// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "PT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

#define SMALL_VALUE 1e-300
#define LARGE_VALUE 1e+300


// constructor
PT::PT(const mat X, const int prior, const vec support, const int max_M, const vec prior_params, const int allocated) :
    X(X), prior(prior), support(support), max_M(max_M), prior_params(prior_params), allocated(allocated) {

    n = X.n_rows;
    NODES_NUM = (1 << (max_M + 1)) - 1;

    largest_obs = -LARGE_VALUE;
    for (uword i = 0; i < X.n_rows; ++i) {
        for (uword j = 0; j < 2; ++j) {
            const double v = X(i, j);
            if (arma::is_finite(v) && v > largest_obs) {
                largest_obs = v;
            }
        }
    }

    // initialize other member variables
    left_bound = vec(NODES_NUM, fill::zeros);
    left_bound(0) = support(0);

    right_bound = vec(NODES_NUM, fill::zeros);
    right_bound(0) = support(1);

    node_levels = vec(NODES_NUM, fill::zeros);
    node_levels(0) = 0;

    alpha = vec(NODES_NUM, fill::zeros);
    alpha(0) = 1.0;

    in_node_obs = mat(NODES_NUM, n, fill::zeros);
    for (int i = 0; i < n; i++) {
        in_node_obs(0, i) = 1; // all observations start at the root
    }

    n_obs_allocated = vec(NODES_NUM, fill::zeros);
    n_obs_allocated(0) = n; // all observations start at the root

    posterior_probs = vec(NODES_NUM, fill::zeros);
    posterior_probs(0) = 1.0; // root node starts with probability 1

    new_posterior_probs = vec(NODES_NUM, fill::zeros);
    new_posterior_probs(0) = 1.0; // root node starts with probability 1

    init();
}

void PT::init() {
    for (int i = 0; i < NODES_NUM; i++) {
        int level = node_levels(i);
        if (level == max_M) {
            terminal_nodes.insert_rows(terminal_nodes.n_rows, vec({(double)i}));
        }

        if (level < max_M) {
            double split_point = choose_split(level, left_bound(i), right_bound(i));
            int left_child = 2 * i + 1;
            int right_child = 2 * i + 2;
            left_bound(left_child) = left_bound(i);
            right_bound(left_child) = split_point;

            left_bound(right_child) = split_point;
            right_bound(right_child) = right_bound(i);

            node_levels(left_child) = level + 1;
            node_levels(right_child) = level + 1;

            // prior alpha values
            alpha(left_child) = init_alpha(level + 1, left_bound(left_child), right_bound(left_child));
            alpha(right_child) = init_alpha(level + 1, left_bound(right_child), right_bound(right_child));

            // allocate observations to children
            init_allocations(i, left_child, right_child);

            // posterior probabilities (initially set to prior probabilities)
            double total = alpha(left_child) + alpha(right_child) + n_obs_allocated(left_child) + n_obs_allocated(right_child);
            if (total > SMALL_VALUE) {
                posterior_probs(left_child) = ((alpha(left_child) + n_obs_allocated(left_child)) / total) * posterior_probs(i);
                posterior_probs(right_child) = ((alpha(right_child) + n_obs_allocated(right_child)) / total) * posterior_probs(i);
            } else {
                posterior_probs(left_child) = 0.5 * posterior_probs(i);
                posterior_probs(right_child) = 0.5 * posterior_probs(i);
            }
        }
    }
}

double PT::choose_split(int level, double left, double right) {
    if (arma::is_finite(right)){
        return (left + right) / 2.0; // midpoint split
    } else {
        return left + (largest_obs - support(0)) / (1 << level); // exponential split
    }
}

double PT::init_alpha(int level, double left, double right) {
    if (level <= 0) {
        Rcpp::stop("Level must be positive.");
    }
    if (prior == 1) { // Uniform
        return pow(level, 2.0) * 0.5; // example: increasing with level
    } else if (prior == 2) { // Exponential
        const double lambda = prior_params(0);
        if (arma::is_finite(right)) {
            return pow(level, 2.0) * (exp(-lambda * left) - exp(-lambda * right));
        } else {
            return pow(level, 2.0) * exp(-lambda * left);
        }
    } else if (prior == 3) { // Weibull
        const double lambda = prior_params(0);
        const double k = prior_params(1);
        if (arma::is_finite(right)) {
            return pow(level, 2.0) * (exp(-pow(lambda * left, k)) - exp(-pow(lambda * right, k)));
        } else {
            return pow(level, 2.0) * exp(-pow(lambda * left, k));
        }
    } else {
        Rcpp::stop("Invalid prior type.");
    }
    return 0.0; // unreachable, quiets some compilers
}

void PT::init_allocations(int parent, int left_child, int right_child) {
    // indices of observations allocated to the parent node
    arma::uvec parent_idx = arma::find(in_node_obs.row(parent) == 1.0);

    for (arma::uword k = 0; k < parent_idx.n_elem; ++k) {
        const arma::uword obs = parent_idx(k);
        const double left_obs = X(obs, 0);
        const double right_obs = X(obs, 1);

        if (allocated == 1) {
            // update for left child
            if (max(left_obs, left_bound(left_child)) <= min(right_obs, right_bound(left_child))) {
                in_node_obs(left_child, obs) = 1.0;
                n_obs_allocated(left_child) += init_individual_allocation(left_obs, right_obs, left_bound(left_child), right_bound(left_child));
            }

            // update for right child
            if (max(left_obs, left_bound(right_child)) <= min(right_obs, right_bound(right_child))) {
                in_node_obs(right_child, obs) = 1.0;
                n_obs_allocated(right_child) += init_individual_allocation(left_obs, right_obs, left_bound(right_child), right_bound(right_child));
            }
        } else if (allocated == 0) {
            // hard allocation for MW method
            if (left_obs < right_bound(left_child) && right_obs > left_bound(left_child)) {
                in_node_obs(left_child, obs) = 1.0;
                n_obs_allocated(left_child) += 1.0;
            }
            if (left_obs < right_bound(right_child) && right_obs > left_bound(right_child)) {
                in_node_obs(right_child, obs) = 1.0;
                n_obs_allocated(right_child) += 1.0;
            }
        } else {
            Rcpp::stop("Invalid allocation type.");
        }
    }
    if (allocated == 1){
    // check for sum of allocations
        if (abs(n_obs_allocated(left_child) + n_obs_allocated(right_child) - n_obs_allocated(parent)) > 1e-5) {
            std::cout << "Parent node: " << n_obs_allocated(parent) << ", Left child: " << n_obs_allocated(left_child) << ", Right child: " << n_obs_allocated(right_child) << std::endl;
            Rcpp::stop("Sum of allocations exceeds parent sample size.");
        }
    }
}

double PT::init_individual_allocation(double left_obs, double right_obs, double left_node, double right_node) {
    double joint_left = max(left_obs, left_node);
    double joint_right = min(right_obs, right_node);
    if (prior == 1) { // Uniform
        double width_obs;
        if (arma::is_finite(right_obs)){
            width_obs = right_obs - left_obs;
        } else {
            width_obs = 2 * largest_obs - left_obs; // arbitrary large width for right-censored
        }
        double joint_width;
        if (arma::is_finite(joint_right)){
            joint_width = joint_right - joint_left;
        } else {
            joint_width = 2 * largest_obs - joint_left; // arbitrary large width for right-censored
        }
        if (width_obs > SMALL_VALUE) {
            return joint_width / width_obs;
        } else {
            return (right_node != left_obs) ? 1.0 : 0.0;
        }
    } else if (prior == 2) { // Exponential
        const double lambda = prior_params(0);
        double numerator, denominator;
        if (arma::is_finite(joint_right)) {
            numerator = exp(-lambda * joint_left) - exp(-lambda * joint_right);
        } else {
            numerator = exp(-lambda * joint_left);
        }
        if (arma::is_finite(right_obs)) {
            denominator = exp(-lambda * left_obs) - exp(-lambda * right_obs);
        } else {
            denominator = exp(-lambda * left_obs);
        }
        if (denominator > SMALL_VALUE) {
            return numerator / denominator;
        } else {
            return (right_node != left_obs) ? 1.0 : 0.0;
        }
    } else if (prior == 3) { // Weibull
        const double lambda = prior_params(0);
        const double k = prior_params(1);
        double numerator, denominator;
        if (arma::is_finite(joint_right)) {
            numerator = exp(-pow(lambda * joint_left, k)) - exp(-pow(lambda * joint_right, k));
        } else {
            numerator = exp(-pow(lambda * joint_left, k));
        }
        if (arma::is_finite(right_obs)) {
            denominator = exp(-pow(lambda * left_obs, k)) - exp(-pow(lambda * right_obs, k));
        } else {
            denominator = exp(-pow(lambda * left_obs, k));
        }
        if (denominator > SMALL_VALUE) {
            return numerator / denominator;
        } else {
            return (right_node != left_obs) ? 1.0 : 0.0;
        }
    } else {
        Rcpp::stop("Invalid prior type.");
    }
    return 0.0; // unreachable, quiets some compilers
}

void PT::re_allocation() {
    for (int i = 0; i < NODES_NUM; i++) {
        int level = node_levels(i);
        if (level < max_M) {
            int left_child = 2 * i + 1;
            int right_child = 2 * i + 2;

            // allocate observations to children
            update_allocations(i, left_child, right_child);

            // posterior probabilities (initially set to prior probabilities)
            double total = alpha(left_child) + alpha(right_child) + n_obs_allocated(left_child) + n_obs_allocated(right_child);
            if (total > SMALL_VALUE) {
                new_posterior_probs(left_child) = ((alpha(left_child) + n_obs_allocated(left_child)) / total) * new_posterior_probs(i);
                new_posterior_probs(right_child) = ((alpha(right_child) + n_obs_allocated(right_child)) / total) * new_posterior_probs(i);
            } else {
                new_posterior_probs(left_child) = 0.5 * new_posterior_probs(i);
                new_posterior_probs(right_child) = 0.5 * new_posterior_probs(i);
            }
        }
    }
    posterior_probs = new_posterior_probs;
}

void PT::update_allocations(int parent, int left_child, int right_child) {
    // reset allocations for children
    in_node_obs.row(left_child).zeros();
    in_node_obs.row(right_child).zeros();
    n_obs_allocated(left_child) = 0.0;
    n_obs_allocated(right_child) = 0.0;

    // indices of observations allocated to the parent node
    arma::uvec parent_idx = arma::find(in_node_obs.row(parent) == 1.0);

    for (arma::uword k = 0; k < parent_idx.n_elem; ++k) {
        const arma::uword obs = parent_idx(k);
        const double left_obs = X(obs, 0);
        const double right_obs = X(obs, 1);

        if (allocated == 1) {
            // update for left child
            if (max(left_obs, left_bound(left_child)) <= min(right_obs, right_bound(left_child))) {
                in_node_obs(left_child, obs) = 1.0;
                n_obs_allocated(left_child) += update_individual_allocation(left_obs, right_obs, left_bound(left_child), right_bound(left_child));
            }

            // update for right child
            if (max(left_obs, left_bound(right_child)) <= min(right_obs, right_bound(right_child))) {
                in_node_obs(right_child, obs) = 1.0;
                n_obs_allocated(right_child) += update_individual_allocation(left_obs, right_obs, left_bound(right_child), right_bound(right_child));
            }
        } else if (allocated == 0) {
            // hard allocation for MW method
            if (left_obs < right_bound(left_child) && right_obs > left_bound(left_child)) {
                in_node_obs(left_child, obs) = 1.0;
                n_obs_allocated(left_child) += 1.0;
            }
            if (left_obs < right_bound(right_child) && right_obs > left_bound(right_child)) {
                in_node_obs(right_child, obs) = 1.0;
                n_obs_allocated(right_child) += 1.0;
            }
        } else {
            Rcpp::stop("Invalid allocation type.");
        }
    }
    if (allocated == 1){
    // check for sum of allocations
        if (abs(n_obs_allocated(left_child) + n_obs_allocated(right_child) - n_obs_allocated(parent)) > 1e-5) {
            std::cout << "Parent node: " << n_obs_allocated(parent) << ", Left child: " << n_obs_allocated(left_child) << ", Right child: " << n_obs_allocated(right_child) << std::endl;
            Rcpp::stop("Sum of allocations exceeds parent sample size.");
        }
    }
}

double PT::update_individual_allocation(double left_obs, double right_obs, double left_node, double right_node) {
    // Update individual allocation based on the terminal nodes and previous posterior probabilities
    double joint_left = max(left_obs, left_node);
    double joint_right = min(right_obs, right_node);
    double p_joint = 0.0;
    double p_data = 0.0;

    // Loop over all terminal nodes
    for (arma::uword t = 0; t < terminal_nodes.n_rows; ++t) {
        int node = static_cast<int>(terminal_nodes(t));
        double node_left = left_bound(node);
        double node_right = (arma::is_finite(right_bound(node))) ? right_bound(node) : 2 * largest_obs; // handle infinite right bound

        // Compute overlap between [joint_left, joint_right] and [node_left, node_right]
        double overlap_left = std::max(joint_left, node_left);
        double overlap_right = std::min(joint_right, node_right);

        if (overlap_left < overlap_right) {
            p_joint += (overlap_right - overlap_left) / (node_right - node_left) * posterior_probs(node);
        }

        // Compute overlap between [left_obs, right_obs] and [node_left, node_right]
        double obs_overlap_left = std::max(left_obs, node_left);
        double obs_overlap_right = std::min(right_obs, node_right);

        if (obs_overlap_left < obs_overlap_right) {
            p_data += (obs_overlap_right - obs_overlap_left) / (node_right - node_left) * posterior_probs(node);
        }
    }

    // Avoid division by zero
    if (p_data > SMALL_VALUE) {
        return p_joint / p_data;
    } else {
        return (right_node != left_obs) ? 1.0 : 0.0;
    }

    return 0.0;
}

List PT::output() {
    return List::create(
        Named("left_bound") = left_bound,
        Named("right_bound") = right_bound,
        Named("node_levels") = node_levels,
        Named("alpha") = alpha,
        Named("n_obs_allocated") = n_obs_allocated,
        Named("posterior_probs") = posterior_probs,
        Named("terminal_nodes") = terminal_nodes
    );
}

void PT::clear(){
    // clear all member variables
    in_node_obs.clear();
    new_posterior_probs.clear();
}
