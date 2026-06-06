// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>
#include "PT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

#define SMALL_VALUE 1e-300
#define LARGE_VALUE 1e+300


// constructor
PT::PT(const mat X, const vec test_seq, const int prior, const vec support, const int max_M, const vec prior_params, const int allocated, const int max_EM_iter, const double EM_tol, const bool verbose) :
    X(X), test_seq(test_seq), prior(prior), support(support), max_M(max_M), prior_params(prior_params), allocated(allocated), max_EM_iter(max_EM_iter), EM_tol(EM_tol), verbose(verbose) {

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
            double left_prob = min(max(left_count / total, 0.0), 1.0);
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
            return left + (1.1 * largest_obs - support(0)) / (1 << (level + 1)); // exponential split
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
    } else {
        Rcpp::stop("Invalid prior type.");
    }
}

double PT::init_alpha(int level, double left, double right) {
    if (level <= 0) {
        Rcpp::stop("Level must be positive.");
    }
    return pow(level, 2.0) * 0.5; // example: increasing with level
    // Uncomment below for prior-dependent alpha initialization
    // if (prior == 1) { // Uniform
    //     return pow(level, 2.0) * 0.5;
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
                if (left_obs > left_bound(left_child) && right_obs <= right_bound(left_child)) {
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
                if (left_obs > left_bound(right_child) && right_obs <= right_bound(right_child)) {
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
                if (left_obs > left_bound(left_child) && right_obs <= right_bound(left_child)) {
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += 1.0;
                }
                if (left_obs > left_bound(right_child) && right_obs <= right_bound(right_child)) {
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
    if (prior == 1) { // Uniform
        double width_obs;
        if (arma::is_finite(right_obs)){
            width_obs = right_obs - left_obs;
        } else {
            width_obs = 1.1 * largest_obs - left_obs; // arbitrary large width for right-censored
        }
        double joint_width;
        if (arma::is_finite(joint_right)){
            joint_width = joint_right - joint_left;
        } else {
            joint_width = 1.1 * largest_obs - joint_left; // arbitrary large width for right-censored
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
    bool condition = true;
    int iter = 1;
    if (verbose) {
        std::cout << "Starting EM iterations..." << std::endl;
    }
    auto start = std::chrono::high_resolution_clock::now();
    while (condition && iter <= max_EM_iter)
    {
        arma::uvec idx = arma::conv_to<arma::uvec>::from(terminal_nodes);
        vec old_terminal_probs = posterior_probs.elem(idx);
        re_allocation(); // EM step
        double diff = arma::norm(posterior_probs.elem(idx) - old_terminal_probs, 2);
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
    for (int i = 0; i < nodes_used; i++) {

        // skip terminal nodes
        if (arma::any(terminal_nodes == static_cast<double>(i))) {
            continue;
        }

        // find children by parent_nodes
        arma::uvec children = arma::find(parent_nodes == static_cast<double>(i));

        int left_child = min(children);
        int right_child = max(children);

        // allocate observations to children
        update_allocations(i, left_child, right_child);

        // posterior probabilities (initially set to prior probabilities)
        double total = alpha(left_child) + alpha(right_child) + n_obs_allocated(left_child) + n_obs_allocated(right_child) - 2;
        double left_count = alpha(left_child) + n_obs_allocated(left_child) - 1;
        if (total > SMALL_VALUE) {
            double left_prob = min(max(left_count / total, 0.0), 1.0);
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

                if (left_obs > left_bound(left_child) && right_obs <= right_bound(left_child)) {
                    // point mass
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += 1.0;
                } else if (right_obs != left_obs) {
                    in_node_obs(left_child, obs) = 1.0;
                    n_obs_allocated(left_child) += update_individual_allocation(left_obs, right_obs, left_bound(left_child), right_bound(left_child));
                }
            }
            // update for right child
            if (max(left_obs, left_bound(right_child)) <= min(right_obs, right_bound(right_child))) {

                if (left_obs > left_bound(right_child) && right_obs <= right_bound(right_child)) {
                    // point mass
                    in_node_obs(right_child, obs) = 1.0;
                    n_obs_allocated(right_child) += 1.0;
                } else if (right_obs != left_obs) {
                    in_node_obs(right_child, obs) = 1.0;
                    n_obs_allocated(right_child) += update_individual_allocation(left_obs, right_obs, left_bound(right_child), right_bound(right_child));
                }
            }
        } else if (allocated == 0) {
            // hard allocation for MW method
            if (left_obs > left_bound(left_child) && right_obs <= right_bound(left_child)) {
                in_node_obs(left_child, obs) = 1.0;
                n_obs_allocated(left_child) += 1.0;
            }
            if (left_obs > left_bound(right_child) && right_obs <= right_bound(right_child)) {
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

        // Compute overlap between [joint_left, joint_right] and [node_left, node_right]
        double overlap_left = std::max(joint_left, node_left);
        double overlap_right = std::min(joint_right, node_right);

        if (overlap_left < overlap_right) {
            if (arma::is_finite(overlap_right)) {
              p_joint += (overlap_right - overlap_left) / (node_right - node_left) * posterior_probs(node);
            } else {
              p_joint += posterior_probs(node);
            }
        }

        // Compute overlap between [left_obs, right_obs] and [node_left, node_right]
        double obs_overlap_left = std::max(left_obs, node_left);
        double obs_overlap_right = std::min(right_obs, node_right);

        if (obs_overlap_left < obs_overlap_right) {
            if (arma::is_finite(obs_overlap_right)) {
              p_data += (obs_overlap_right - obs_overlap_left) / (node_right - node_left) * posterior_probs(node);
            } else {
              p_data += posterior_probs(node);
            }
        }
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


List PT::output() {
    // compute survival probabilities and densities at test_seq by mixing over terminal nodes
    survival_probs = vec(test_seq.n_elem, fill::zeros);
    vec density_probs(test_seq.n_elem, fill::zeros);

    auto base_survival = [&](double x) -> double {
        if (prior == 2) {
            const double lambda = prior_params(0);
            return exp(-lambda * x);
        } else if (prior == 3) {
            const double lambda = prior_params(0);
            const double k = prior_params(1);
            return exp(-pow(lambda * x, k));
        }
        Rcpp::stop("Invalid prior type.");
        return 0.0;
    };

    auto base_density = [&](double x) -> double {
        if (prior == 2) {
            const double lambda = prior_params(0);
            return lambda * exp(-lambda * x);
        } else if (prior == 3) {
            const double lambda = prior_params(0);
            const double k = prior_params(1);
            return k * lambda * pow(lambda * x, k - 1.0) * exp(-pow(lambda * x, k));
        }
        Rcpp::stop("Invalid prior type.");
        return 0.0;
    };

    auto terminal_survival = [&](double t, double a, double b) -> double {
        if (prior == 1) {
            double right = b;
            if (!arma::is_finite(right)) {
                right = arma::is_finite(support(1)) ? support(1) : 1.1 * largest_obs;
            }
            const double width = right - a;
            if (width > SMALL_VALUE && t < right) {
                return min(max((right - t) / width, 0.0), 1.0);
            }
            return 0.0;
        }

        const double survival_a = base_survival(a);
        const double survival_t = base_survival(t);
        const double survival_b = arma::is_finite(b) ? base_survival(b) : 0.0;
        const double denom = survival_a - survival_b;
        if (denom > SMALL_VALUE) {
            return min(max((survival_t - survival_b) / denom, 0.0), 1.0);
        }
        return 0.0;
    };

    auto terminal_density = [&](double t, double a, double b) -> double {
        if (prior == 1) {
            double right = b;
            if (!arma::is_finite(right)) {
                right = arma::is_finite(support(1)) ? support(1) : 1.1 * largest_obs;
            }
            const double width = right - a;
            if (width > SMALL_VALUE && t >= a && t < right) {
                return 1.0 / width;
            }
            return 0.0;
        }

        if (t < a || (arma::is_finite(b) && t >= b)) {
            return 0.0;
        }

        const double survival_a = base_survival(a);
        const double survival_b = arma::is_finite(b) ? base_survival(b) : 0.0;
        const double denom = survival_a - survival_b;
        if (denom > SMALL_VALUE) {
            return base_density(t) / denom;
        }
        return 0.0;
    };

    for (arma::uword i = 0; i < test_seq.n_elem; ++i) {
        const double t = test_seq(i);
        double s = 0.0;
        double d = 0.0;

        for (arma::uword r = 0; r < terminal_nodes.n_rows; ++r) {
            const int node = static_cast<int>(terminal_nodes(r));
            const double a = left_bound(node);
            double b = right_bound(node);
            const double w = posterior_probs(node);

            d += w * terminal_density(t, a, b);

            if (t < a) {
                s += w;
            } else if (t >= b) {
                // contributes 0
            } else {
                s += w * terminal_survival(t, a, b);
            }
        }

        // clamp numeric noise
        survival_probs(i) = s;
        density_probs(i) = d;
    }

    return List::create(
        Named("left_bound") = left_bound.head(nodes_used+1),
        Named("right_bound") = right_bound.head(nodes_used+1),
        Named("node_levels") = node_levels.head(nodes_used+1),
        Named("alpha") = alpha.head(nodes_used+1),
        Named("n_obs_allocated") = n_obs_allocated.head(nodes_used+1),
        Named("posterior_probs") = posterior_probs.head(nodes_used+1),
        Named("terminal_nodes") = terminal_nodes,
        Named("survival_probs") = survival_probs,
        Named("density_probs") = density_probs
    );
}





void PT::clear(){
    // clear all member variables
    in_node_obs.clear();
    new_posterior_probs.clear();
    parent_nodes.clear();
}
