// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include "PT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

#define SMALL_VALUE 1e-300
#define LARGE_VALUE 1e+300

namespace {
constexpr double SPLIT_EPSILON = 1e-6;
}


// constructor
PT::PT(const mat X, const int prior, const vec support, const int max_M, const vec prior_params, const int allocated, const int max_EM_iter, const double EM_tol, const bool verbose, const double uniform_tail_factor, const bool flat) :
    X(X), prior(prior), support(support), max_M(max_M), prior_params(prior_params), allocated(allocated), uniform_tail_factor(uniform_tail_factor), max_EM_iter(max_EM_iter), EM_tol(EM_tol), verbose(verbose), flat(flat) {

    if (!arma::is_finite(uniform_tail_factor) || uniform_tail_factor <= 1.0) {
        Rcpp::stop("uniform_tail_factor must be a finite number greater than 1.");
    }
    if (prior == 4 && (X.n_cols < 2 || X.n_rows == 0)) {
        Rcpp::stop("prior = 4 requires nonempty data with two observation endpoint columns.");
    }

    n = X.n_rows;
    if (max_M <= 10){
        NODES_NUM = (1 << (max_M + 1)) - 1; // number of nodes for complete tree
    } else {
        NODES_NUM = (1 << (8 + 1)) - 1; // number of nodes first 8 complete levels
    }
    nodes_used = 0;

    largest_obs = -LARGE_VALUE;
    for (uword i = 0; i < X.n_rows; ++i) {
        for (uword j = 0; j < 2; ++j) {
            const double v = X(i, j);
            if (arma::is_finite(v) && v > largest_obs) {
                largest_obs = v;
            }
        }
    }

    if (prior == 4) {
        if (!arma::is_finite(largest_obs) || largest_obs <= 0.0) {
            Rcpp::stop("prior = 4 requires a positive finite largest observation endpoint.");
        }
        if (prior_params.n_elem != 1 || !arma::is_finite(prior_params(0)) || prior_params(0) <= 0.0) {
            Rcpp::stop("prior = 4 requires one finite positive exponential rate in prior_params.");
        }
        if (support.n_elem != 2 || !arma::is_finite(support(0)) || support(0) < 0.0 ||
            std::isnan(support(1)) || support(1) <= support(0)) {
            Rcpp::stop("prior = 4 requires nonnegative, ordered support bounds.");
        }
        for (uword obs = 0; obs < X.n_rows; ++obs) {
            if (!arma::is_finite(X(obs, 0)) || X(obs, 0) < 0.0 ||
                std::isnan(X(obs, 1)) || X(obs, 1) < X(obs, 0)) {
                Rcpp::stop("prior = 4 requires finite nonnegative left endpoints and ordered, nonmissing right endpoints.");
            }
        }
        hybrid_body_mass = -std::expm1(-prior_params(0) * largest_obs);
        hybrid_tail_mass = std::exp(-prior_params(0) * largest_obs);
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

    parent_nodes = vec(NODES_NUM, fill::zeros);
    for (int i = 0; i < NODES_NUM; i++) {
        parent_nodes(i) = -1; // root has no parent
    }

    in_node_obs = mat(NODES_NUM, n, fill::zeros);
    for (int i = 0; i < n; i++) {
        in_node_obs(0, i) = 1.0; // all observations start at the root
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

        double split_point = choose_split(level, left_bound(i), right_bound(i));
        if (split_point < 0 || n_obs_allocated(i) < 1) {
        // if (split_point < 0) {
            terminal_nodes.insert_rows(terminal_nodes.n_rows, vec({(double)i}));
            continue; // no split
        }
        int left_child = ++nodes_used;
        int right_child = ++nodes_used;
        if (left_child >= NODES_NUM || right_child >= NODES_NUM) {
            update_tree_size();
        }
        left_bound(left_child) = left_bound(i);
        right_bound(left_child) = split_point;

        left_bound(right_child) = split_point;
        right_bound(right_child) = right_bound(i);

        node_levels(left_child) = level + 1;
        node_levels(right_child) = level + 1;

        parent_nodes(left_child) = i;
        parent_nodes(right_child) = i;

        // prior alpha values
        alpha(left_child) = init_alpha(level + 1, left_bound(left_child), right_bound(left_child));
        alpha(right_child) = init_alpha(level + 1, left_bound(right_child), right_bound(right_child));

        // allocate observations to children
        init_allocations(i, left_child, right_child);

        // posterior probabilities (initially set to prior probabilities)
        double total = alpha(left_child) + alpha(right_child) + n_obs_allocated(left_child) + n_obs_allocated(right_child) - 2;
        double left_count = alpha(left_child) + n_obs_allocated(left_child) - 1;
        if (total > SMALL_VALUE) {
            double left_prob = std::clamp(left_count / total, SPLIT_EPSILON, 1.0 - SPLIT_EPSILON);
            posterior_probs(left_child) = left_prob * posterior_probs(i);
            posterior_probs(right_child) = (1.0 - left_prob) * posterior_probs(i);
        } else {
            posterior_probs(left_child) = 0.5 * posterior_probs(i);
            posterior_probs(right_child) = 0.5 * posterior_probs(i);
        }
    }

    // remove terminal node entries that reference nodes beyond the current nodes_used
    if (terminal_nodes.n_rows > 0) {
        arma::uvec keep = arma::find(terminal_nodes <= static_cast<double>(nodes_used));
        if (keep.n_elem < terminal_nodes.n_rows) {
            terminal_nodes = terminal_nodes.elem(keep);
        }
    }
}

double PT::choose_split(int level, double left, double right) {
    if (level == max_M) {
        return -1.0; // no split at the last level
    }
    if (prior == 1) { // Uniform
        if (arma::is_finite(right)){
            return (left + right) / 2.0; // midpoint split
        } else {
            return left + (uniform_tail_factor * largest_obs - support(0)) / (1 << (level + 1)); // exponential split
        }
    } else if (prior == 2) { // Exponential
        const double lambda = prior_params(0);
        if (arma::is_finite(right)){
            double cdf_left = 1 - exp(-lambda * left);
            double cdf_right = 1 - exp(-lambda * right);
            double split_cdf = 0.5 * (cdf_right + cdf_left);
            return -log(1 - split_cdf) / lambda;
        } else {
            double cdf_left = 1 - exp(-lambda * left);
            double split_cdf = 0.5 * (1 + cdf_left);
            return -log(1 - split_cdf) / lambda;
        }
    } else if (prior == 3) { // Weibull
        const double lambda = prior_params(0);
        const double k = prior_params(1);
        if (arma::is_finite(right)){
            double cdf_left = 1 - exp(-pow(lambda * left, k));
            double cdf_right = 1 - exp(-pow(lambda * right, k));
            double split_cdf =  0.5 * (cdf_right + cdf_left);
            return pow(-log(1 - split_cdf), 1.0 / k) / lambda;
        } else {
            double cdf_left = 1 - exp(-pow(lambda * left, k));
            double split_cdf = 0.5 * (1 + cdf_left);
            return pow(-log(1 - split_cdf), 1.0 / k) / lambda;
        }
    } else if (prior == 4) {
        const double midpoint = 0.5 * hybrid_cdf(left) + 0.5 * hybrid_cdf(right);
        // Use the survival scale in the upper half to avoid rounding the CDF
        // to one before inverting it in an exponential tail.
        const double split = midpoint <= 0.5 ? hybrid_quantile(midpoint, true) :
            hybrid_quantile(0.5 * base_survival(left) + 0.5 * base_survival(right), false);
        return arma::is_finite(split) && split > left && split < right ? split : -1.0;
    } else {
        Rcpp::stop("Invalid prior type.");
    }
}

double PT::init_alpha(int level, double left, double right) {
    if (level <= 0) {
        Rcpp::stop("Level must be positive.");
    }
    return flat ? 1.0 : std::max(1.0, pow(level, 2.0) * 0.5);
    // Uncomment below for prior-dependent alpha initialization
    // if (prior == 1) { // Uniform
    //     return std::max(1.0, pow(level, 2.0) * 0.5);
    // } else if (prior == 2) { // Exponential
    //     const double lambda = prior_params(0);
    //     if (arma::is_finite(right)) {
    //         return pow(level, 2.0) * (exp(-lambda * left) - exp(-lambda * right));
    //     } else {
    //         return pow(level, 2.0) * exp(-lambda * left);
    //     }
    // } else if (prior == 3) { // Weibull
    //     const double lambda = prior_params(0);
    //     const double k = prior_params(1);
    //     if (arma::is_finite(right)) {
    //         return pow(level, 2.0) * (exp(-pow(lambda * left, k)) - exp(-pow(lambda * right, k)));
    //     } else {
    //         return pow(level, 2.0) * exp(-pow(lambda * left, k));
    //     }
    // } else {
    //     Rcpp::stop("Invalid prior type.");
    // }
    // return 0.0; // unreachable, quiets some compilers
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
                if ((left_obs > left_bound(left_child) ||
                     (left_obs == 0.0 && left_bound(left_child) == 0.0)) &&
                    right_obs <= right_bound(left_child)) {
                    // point mass
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += 1.0;
                } else if (right_obs != left_obs) {
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += init_individual_allocation(left_obs, right_obs, left_bound(left_child), right_bound(left_child));
                }
            }

            // update for right child
            if (max(left_obs, left_bound(right_child)) <= min(right_obs, right_bound(right_child))) {
                if ((left_obs > left_bound(right_child) ||
                     (left_obs == 0.0 && left_bound(right_child) == 0.0)) &&
                    right_obs <= right_bound(right_child)) {
                    // point mass
                    in_node_obs(right_child, obs) = 1.0;
                    n_obs_allocated(right_child) += 1.0;
                } else if (right_obs != left_obs) {
                    in_node_obs(right_child, obs) = 1.0;
                    n_obs_allocated(right_child) += init_individual_allocation(left_obs, right_obs, left_bound(right_child), right_bound(right_child));
                }
            }
        } else if (allocated == 0) {
            // hard allocation for MW method
            if (left_obs == right_obs) {
                if ((left_obs > left_bound(left_child) ||
                     (left_obs == 0.0 && left_bound(left_child) == 0.0)) &&
                    right_obs <= right_bound(left_child)) {
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += 1.0;
                }
                if ((left_obs > left_bound(right_child) ||
                     (left_obs == 0.0 && left_bound(right_child) == 0.0)) &&
                    right_obs <= right_bound(right_child)) {
                    in_node_obs(right_child, obs) = 1.0;
                    n_obs_allocated(right_child) += 1.0;
                }
            } else {
                if (left_obs < right_bound(left_child) && right_obs > left_bound(left_child)) {
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += 1.0;
                }
                if (left_obs < right_bound(right_child) && right_obs > left_bound(right_child)) {
                    in_node_obs(right_child, obs) = 1.0;
                    n_obs_allocated(right_child) += 1.0;
                }
            }
        } else {
            Rcpp::stop("Invalid allocation type.");
        }
    }
    if (allocated == 1){
    // check for sum of allocations
        if (abs(n_obs_allocated(left_child) + n_obs_allocated(right_child) - n_obs_allocated(parent)) > 1e-5) {
            std::cout << "Parent node: " << n_obs_allocated(parent) << ", Left child: " << n_obs_allocated(left_child) << ", Right child: " << n_obs_allocated(right_child) << std::endl;
            std::cout << "Sum of allocations exceeds parent sample size." << std::endl;
        }
    }
}

double PT::init_individual_allocation(double left_obs, double right_obs, double left_node, double right_node) {
    double joint_left = max(left_obs, left_node);
    double joint_right = min(right_obs, right_node);
    if (prior == 4) {
        const double denominator = hybrid_probability(left_obs, right_obs);
        return denominator > SMALL_VALUE ?
            std::clamp(hybrid_probability(joint_left, joint_right) / denominator, 0.0, 1.0) : 0.0;
    }
    if (prior == 1) { // Uniform
        double width_obs;
        if (arma::is_finite(right_obs)){
            width_obs = right_obs - left_obs;
        } else {
            width_obs = uniform_tail_factor * largest_obs - left_obs; // finite width for right-censored
        }
        double joint_width;
        if (arma::is_finite(joint_right)){
            joint_width = joint_right - joint_left;
        } else {
            joint_width = uniform_tail_factor * largest_obs - joint_left; // finite width for right-censored
        }
        if (width_obs > SMALL_VALUE) {
            return joint_width / width_obs;
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
        }
    } else {
        Rcpp::stop("Invalid prior type.");
    }
    return 0.0; // unreachable, quiets some compilers
}


void PT::EM_step() {
    if (max_EM_iter > 0) prepare_em_cache();
    bool condition = true;
    int iter = 1;
    if (verbose) {
        std::cout << "Starting EM iterations..." << std::endl;
    }
    auto start = std::chrono::high_resolution_clock::now();
    while (condition && iter <= max_EM_iter)
    {
        vec old_terminal_probs = posterior_probs.elem(terminal_indices);
        re_allocation(); // EM step
        double diff = arma::norm(posterior_probs.elem(terminal_indices) - old_terminal_probs, 2);
        condition = (diff > EM_tol);
        if (verbose) {
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "EM iteration " << iter << ", diff = " << diff << ", time = " << elapsed.count() << " seconds" << std::endl;
        }
        iter++;
    }
}



void PT::re_allocation() {
    if (allocation_cache_ready) update_observation_probabilities();
    for (int i = 0; i < nodes_used; i++) {

        // skip terminal nodes
        if (child_nodes[i].first < 0) {
            continue;
        }

        const int left_child = child_nodes[i].first;
        const int right_child = child_nodes[i].second;

        // allocate observations to children
        update_allocations(i, left_child, right_child);

        // posterior probabilities (initially set to prior probabilities)
        double total = alpha(left_child) + alpha(right_child) + n_obs_allocated(left_child) + n_obs_allocated(right_child) - 2;
        double left_count = alpha(left_child) + n_obs_allocated(left_child) - 1;
        if (total > SMALL_VALUE) {
            double left_prob = std::clamp(left_count / total, SPLIT_EPSILON, 1.0 - SPLIT_EPSILON);
            new_posterior_probs(left_child) = left_prob * new_posterior_probs(i);
            new_posterior_probs(right_child) = (1 - left_prob) * new_posterior_probs(i);
        } else {
            new_posterior_probs(left_child) = 0.5 * new_posterior_probs(i);
            new_posterior_probs(right_child) = 0.5 * new_posterior_probs(i);
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

                if ((left_obs > left_bound(left_child) ||
                     (left_obs == 0.0 && left_bound(left_child) == 0.0)) &&
                    right_obs <= right_bound(left_child)) {
                    // point mass
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += 1.0;
                } else if (right_obs != left_obs) {
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += allocation_cache_ready ?
                        cached_individual_allocation(obs, left_child) :
                        update_individual_allocation(left_obs, right_obs, left_bound(left_child), right_bound(left_child));
                }
            }
            // update for right child
            if (max(left_obs, left_bound(right_child)) <= min(right_obs, right_bound(right_child))) {

                if ((left_obs > left_bound(right_child) ||
                     (left_obs == 0.0 && left_bound(right_child) == 0.0)) &&
                    right_obs <= right_bound(right_child)) {
                    // point mass
                    in_node_obs(right_child, obs) = 1.0;
                    n_obs_allocated(right_child) += 1.0;
                } else if (right_obs != left_obs) {
                    in_node_obs(right_child, obs) = 1.0;
                    n_obs_allocated(right_child) += allocation_cache_ready ?
                        cached_individual_allocation(obs, right_child) :
                        update_individual_allocation(left_obs, right_obs, left_bound(right_child), right_bound(right_child));
                }
            }
        } else if (allocated == 0) {
            // hard allocation for MW method
            if ((left_obs > left_bound(left_child) ||
                 (left_obs == 0.0 && left_bound(left_child) == 0.0)) &&
                right_obs <= right_bound(left_child)) {
                in_node_obs(left_child, obs) = 1.0;
                n_obs_allocated(left_child) += 1.0;
            }
            if ((left_obs > left_bound(right_child) ||
                 (left_obs == 0.0 && left_bound(right_child) == 0.0)) &&
                right_obs <= right_bound(right_child)) {
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
            std::cout << "Sum of allocations exceeds parent sample size." << std::endl;
        }
    }
}

void PT::prepare_em_cache() {
    terminal_indices = arma::conv_to<arma::uvec>::from(terminal_nodes);
    child_nodes.assign(nodes_used + 1, std::make_pair(-1, -1));
    for (int node = 1; node <= nodes_used; ++node) {
        const int parent = static_cast<int>(parent_nodes(node));
        if (child_nodes[parent].first < 0) child_nodes[parent].first = node;
        else child_nodes[parent].second = node;
    }

    allocation_cache_ready = false;
    if (nodes_used == 0 || n == 0 || !posterior_probs.is_finite()) return;

    bool has_censored = false;
    for (int obs = 0; obs < n; ++obs) {
        const double left = X(obs, 0), right = X(obs, 1);
        if (!arma::is_finite(left) || std::isnan(right) || right < left) return;
        has_censored = has_censored || left != right;
    }
    if (!has_censored) return;

    // Descendant leaves can replace interval intersections only for a proper
    // partition. Preserve the original path for degenerate numerical splits.
    for (int node = 0; node <= nodes_used; ++node) {
        if (!arma::is_finite(left_bound(node)) || std::isnan(right_bound(node)) ||
            right_bound(node) <= left_bound(node)) return;
        const int left = child_nodes[node].first;
        const int right = child_nodes[node].second;
        if (left >= 0 && (right < 0 || left_bound(left) != left_bound(node) ||
            right_bound(left) != left_bound(right) || right_bound(right) != right_bound(node) ||
            right_bound(left) <= left_bound(node) || right_bound(left) >= right_bound(node))) return;
    }
    for (uword k = 0; k < terminal_indices.n_elem; ++k) {
        const uword leaf = terminal_indices(k);
        if (terminal_right_bound(right_bound(leaf)) <= left_bound(leaf)) return;
    }

    // Bound the additional dense probability cache to 64 MiB per fit.
    const uword max_cache_elements = (64u * 1024u * 1024u) / sizeof(double);
    if (terminal_indices.n_elem > max_cache_elements / static_cast<uword>(n)) return;
    observation_leaf_probabilities.zeros(terminal_indices.n_elem, n);
    observation_probabilities.zeros(n);
    for (int obs = 0; obs < n; ++obs) {
        if (X(obs, 0) == X(obs, 1)) continue;
        for (uword k = 0; k < terminal_indices.n_elem; ++k) {
            const uword leaf = terminal_indices(k);
            const double probability = terminal_probability(X(obs, 0), X(obs, 1),
                                                            left_bound(leaf), right_bound(leaf));
            if (!arma::is_finite(probability)) {
                observation_leaf_probabilities.reset();
                observation_probabilities.reset();
                return;
            }
            observation_leaf_probabilities(k, obs) = probability;
        }
    }

    node_terminal_indices.assign(nodes_used + 1, std::vector<uword>());
    for (uword k = 0; k < terminal_indices.n_elem; ++k) {
        int node = static_cast<int>(terminal_indices(k));
        while (node > 0) {
            node_terminal_indices[node].push_back(k);
            node = static_cast<int>(parent_nodes(node));
        }
    }
    allocation_cache_ready = true;
}

void PT::update_observation_probabilities() {
    for (int obs = 0; obs < n; ++obs) {
        if (X(obs, 0) == X(obs, 1)) continue;
        double probability = 0.0;
        for (uword k = 0; k < terminal_indices.n_elem; ++k) {
            probability += posterior_probs(terminal_indices(k)) * observation_leaf_probabilities(k, obs);
        }
        observation_probabilities(obs) = probability;
    }
}

double PT::cached_individual_allocation(uword obs, int node) const {
    double joint_probability = 0.0;
    // Keep the original leaf summation order. Summing child totals or using
    // prefix-sum subtraction would change floating-point rounding.
    for (uword k : node_terminal_indices[node]) {
        joint_probability += posterior_probs(terminal_indices(k)) * observation_leaf_probabilities(k, obs);
    }
    const double probability = observation_probabilities(obs);
    return probability > SMALL_VALUE ? joint_probability / probability : 0.0;
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
        double node_right = right_bound(node);

        // Integrate the same conditional base distribution used by output().
        p_joint += posterior_probs(node) *
            terminal_probability(joint_left, joint_right, node_left, node_right);
        p_data += posterior_probs(node) *
            terminal_probability(left_obs, right_obs, node_left, node_right);
    }

    // Avoid division by zero
    if (p_data > SMALL_VALUE) {
        return p_joint / p_data;
    }

    return 0.0;
}


void PT::update_tree_size() {
    int old_size = NODES_NUM;
    NODES_NUM = 2 * NODES_NUM + 1; // double the size of the tree
    left_bound.resize(NODES_NUM);
    right_bound.resize(NODES_NUM);
    node_levels.resize(NODES_NUM);
    alpha.resize(NODES_NUM);
    in_node_obs.resize(NODES_NUM, n);
    n_obs_allocated.resize(NODES_NUM);
    posterior_probs.resize(NODES_NUM);
    new_posterior_probs.resize(NODES_NUM);
    parent_nodes.resize(NODES_NUM);

    // Initialize new nodes
    for (int i = old_size; i < NODES_NUM; i++) {
        left_bound(i) = 0.0;
        right_bound(i) = 0.0;
        node_levels(i) = 0;
        alpha(i) = 0.0;
        in_node_obs.row(i).zeros();
        n_obs_allocated(i) = 0.0;
        posterior_probs(i) = 0.0;
        new_posterior_probs(i) = 0.0;
        parent_nodes(i) = -1.0;
    }
}


double PT::terminal_right_bound(double right) const {
    if (prior == 1 && !arma::is_finite(right)) {
        return arma::is_finite(support(1)) ? support(1) : uniform_tail_factor * largest_obs;
    }
    return right;
}

double PT::base_survival(double x) const {
    if (prior == 2) {
        return exp(-prior_params(0) * x);
    } else if (prior == 3) {
        return exp(-pow(prior_params(0) * x, prior_params(1)));
    } else if (prior == 4) {
        if (x <= 0.0) return 1.0;
        if (x <= largest_obs) {
            return hybrid_tail_mass + hybrid_body_mass * ((largest_obs - x) / largest_obs);
        }
        return exp(-prior_params(0) * x);
    }
    Rcpp::stop("Invalid prior type.");
    return 0.0;
}

double PT::base_density(double x) const {
    const double lambda = prior_params(0);
    if (prior == 2) {
        return lambda * exp(-lambda * x);
    } else if (prior == 3) {
        const double k = prior_params(1);
        return k * lambda * pow(lambda * x, k - 1.0) * exp(-pow(lambda * x, k));
    } else if (prior == 4) {
        if (x < 0.0 || std::isinf(x)) return 0.0;
        return x <= largest_obs ? hybrid_body_mass / largest_obs : lambda * exp(-lambda * x);
    }
    Rcpp::stop("Invalid prior type.");
    return 0.0;
}

double PT::hybrid_cdf(double x) const {
    if (x <= 0.0) return 0.0;
    if (x <= largest_obs) return hybrid_body_mass * (x / largest_obs);
    return -std::expm1(-prior_params(0) * x);
}

double PT::hybrid_quantile(double probability, bool lower_tail) const {
    if (lower_tail) {
        if (probability <= 0.0) return 0.0;
        if (probability >= 1.0) return arma::datum::inf;
        if (probability <= hybrid_body_mass) return largest_obs * (probability / hybrid_body_mass);
        return -std::log1p(-probability) / prior_params(0);
    }
    if (probability >= 1.0) return 0.0;
    if (probability <= 0.0) return arma::datum::inf;
    if (probability >= hybrid_tail_mass) return largest_obs * ((1.0 - probability) / hybrid_body_mass);
    return -std::log(probability) / prior_params(0);
}

double PT::hybrid_probability(double left, double right) const {
    const double lo = std::max(0.0, left);
    if (right <= lo) return 0.0;
    double probability = 0.0;
    if (lo < largest_obs) {
        probability = hybrid_body_mass * ((std::min(right, largest_obs) - lo) / largest_obs);
    }
    if (right > largest_obs) {
        const double tail_left = std::max(lo, largest_obs);
        probability += std::exp(-prior_params(0) * tail_left) *
            (-std::expm1(-prior_params(0) * (right - tail_left)));
    }
    return probability;
}

double PT::terminal_probability(double left, double right, double a, double b) const {
    b = terminal_right_bound(b);
    const double lo = max(left, a);
    const double hi = min(right, b);
    if (lo >= hi) {
        return 0.0;
    }

    // Conditional continuous probability; endpoint inclusion does not change it.
    const double denom = prior == 4 ? hybrid_probability(a, b) :
        (prior == 1 ? b - a : base_survival(a) - base_survival(b));
    const double numerator = prior == 4 ? hybrid_probability(lo, hi) :
        (prior == 1 ? hi - lo : base_survival(lo) - base_survival(hi));
    if (denom > SMALL_VALUE) {
        return min(max(numerator / denom, 0.0), 1.0);
    }
    return 0.0;
}

double PT::terminal_density(double t, double a, double b) const {
    b = terminal_right_bound(b);
    // Match exact-event allocation: (a, b], or [0, b] when a is zero.
    // Infinity is not an event time, even for an unbounded terminal interval.
    if (std::isinf(t) || t < a || t > b || (t == a && a != 0.0)) {
        return 0.0;
    }
    const double denom = prior == 4 ? hybrid_probability(a, b) :
        (prior == 1 ? b - a : base_survival(a) - base_survival(b));
    if (denom > SMALL_VALUE) {
        return prior == 1 ? 1.0 / denom : base_density(t) / denom;
    }
    return 0.0;
}

List PT::output() {
    return List::create(
        Named("left_bound") = left_bound.head(nodes_used+1),
        Named("right_bound") = right_bound.head(nodes_used+1),
        Named("node_levels") = node_levels.head(nodes_used+1),
        Named("alpha") = alpha.head(nodes_used+1),
        Named("n_obs_allocated") = n_obs_allocated.head(nodes_used+1),
        Named("posterior_probs") = posterior_probs.head(nodes_used+1),
        Named("terminal_nodes") = terminal_nodes,
        Named("parent_nodes") = parent_nodes.head(nodes_used+1),
        Named("largest_obs") = largest_obs
    );
}





void PT::clear(){
    // clear all member variables
    in_node_obs.clear();
    new_posterior_probs.clear();
    parent_nodes.clear();
}
