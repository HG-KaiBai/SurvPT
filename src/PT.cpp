// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>
#include "PT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

#define SMALL_VALUE 1e-300
#define LARGE_VALUE 1e+300

constexpr double SPLIT_EPSILON = 1e-6;

// constructor
PT::PT(const mat X, const vec test_seq, const int prior, const vec support, const int max_M, const vec prior_params, const int allocated, const int max_EM_iter, const double EM_tol, const bool verbose) :
    X(X), test_seq(test_seq), prior(prior), support(support), max_M(max_M), prior_params(prior_params), allocated(allocated), max_EM_iter(max_EM_iter), EM_tol(EM_tol), verbose(verbose) {

    n = X.n_rows;
    if (n == 0 || X.n_cols < 2 || support.n_elem != 2 ||
        !arma::is_finite(support(0)) || support(0) < 0 ||
        std::isnan(support(1)) || support(1) <= support(0) ||
        max_M < 0 || max_EM_iter < 0 || !arma::is_finite(EM_tol) || EM_tol < 0 ||
        (allocated != 0 && allocated != 1))
        Rcpp::stop("Invalid data, support, or iteration controls.");
    if (prior < 1 || prior > 3 ||
        (prior > 1 && (prior_params.n_elem != static_cast<uword>(prior - 1) ||
                       !prior_params.is_finite() || arma::any(prior_params <= 0))))
        Rcpp::stop("Invalid base distribution parameters.");
    has_tail = !arma::is_finite(support(1));
    for (int i = 0; i < n; ++i) {
        if (!arma::is_finite(X(i, 0)) || std::isnan(X(i, 1)) ||
            X(i, 0) < support(0) || X(i, 1) < X(i, 0) || X(i, 1) > support(1) ||
            (X(i, 0) == X(i, 1) && X(i, 0) == support(0)))
            Rcpp::stop("Invalid observation endpoints within support.");
    }
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

    if (has_tail && largest_obs <= support(0))
        Rcpp::stop("The largest finite observed time must exceed the lower support bound.");

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
    new_posterior_probs = posterior_probs;
    prepare_observation_cache();
}

void PT::init() {
    for (int i = 0; i <= nodes_used; i++) {
        int level = node_levels(i);

        // The first split is mandatory, and the tail is never subdivided.
        double split_point = has_tail && i == 0 ? largest_obs :
            (has_tail && i == 2 ? -1.0 :
             choose_split(level - (has_tail ? 1 : 0), left_bound(i), right_bound(i)));
        if (split_point < 0 || (!(has_tail && i == 0) &&
            (n_obs_allocated(i) < 1 || split_point <= left_bound(i) ||
             split_point >= right_bound(i)))) {
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
        alpha(left_child) = has_tail && i == 0 ? 1.0 :
            init_alpha(level + 1 - (has_tail ? 1 : 0), left_bound(left_child), right_bound(left_child));
        alpha(right_child) = has_tail && i == 0 ? 1.0 :
            init_alpha(level + 1 - (has_tail ? 1 : 0), left_bound(right_child), right_bound(right_child));

        // allocate observations to children
        init_allocations(i, left_child, right_child);

        // posterior probabilities (initially set to prior probabilities)
        double total = alpha(left_child) + alpha(right_child) + n_obs_allocated(left_child) + n_obs_allocated(right_child) - 2;
        double left_count = alpha(left_child) + n_obs_allocated(left_child) - 1;
        if (total > SMALL_VALUE) {
            double left_prob = min(max(left_count / total, SPLIT_EPSILON), 1.0 - SPLIT_EPSILON);
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
        return left + (right - left) / 2.0;
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
    return pow(level, 2.0); // alpha_l = gamma * l^2, with gamma = 1
    // Uncomment below for prior-dependent alpha initialization
    // if (prior == 1) { // Uniform
    //     return pow(level, 2.0);
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
    const double body_right = has_tail ? largest_obs : support(1);
    // The first layer has the symmetric level-one Beta prior (mean 1/2).
    const double body_weight = has_tail ? 0.5 : 1.0;
    const double denominator = body_weight *
        terminal_probability(left_obs, right_obs, support(0), body_right) +
        (has_tail && right_obs == arma::datum::inf ? 1.0 - body_weight : 0.0);
    const double numerator = has_tail && left_node == largest_obs &&
        right_node == arma::datum::inf ?
        (right_obs == arma::datum::inf ? 1.0 - body_weight : 0.0) :
        body_weight * terminal_probability(max(left_obs, left_node),
                                           min(right_obs, right_node),
                                           support(0), body_right);
    return denominator > SMALL_VALUE ? numerator / denominator : 0.0;
}


void PT::EM_step() {
    bool condition = true;
    int iter = 1;
    loglik_trace = vec({log_likelihood()});
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
        em_iterations = iter;
        converged = !condition;
        loglik_trace.insert_rows(loglik_trace.n_rows, vec({log_likelihood()}));
        if (verbose) {
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "EM iteration " << iter << ", diff = " << diff << ", time = " << elapsed.count() << " seconds" << std::endl;
        }
        iter++;
    }
}



void PT::re_allocation() {
    // All E-step counts use the same previous tree masses.
    update_expected_allocations();
    for (int i = 0; i < nodes_used; i++) {

        // skip terminal nodes
        if (child_nodes[i].first < 0) {
            continue;
        }

        const int left_child = child_nodes[i].first;
        const int right_child = child_nodes[i].second;

        // posterior probabilities (initially set to prior probabilities)
        double total = alpha(left_child) + alpha(right_child) + n_obs_allocated(left_child) + n_obs_allocated(right_child) - 2;
        double left_count = alpha(left_child) + n_obs_allocated(left_child) - 1;
        if (total > SMALL_VALUE) {
            double left_prob = min(max(left_count / total, SPLIT_EPSILON), 1.0 - SPLIT_EPSILON);
            new_posterior_probs(left_child) = left_prob * new_posterior_probs(i);
            new_posterior_probs(right_child) = (1 - left_prob) * new_posterior_probs(i);
        } else {
            new_posterior_probs(left_child) = 0.5 * new_posterior_probs(i);
            new_posterior_probs(right_child) = 0.5 * new_posterior_probs(i);
        }
    }
    posterior_probs = new_posterior_probs;
}

void PT::prepare_observation_cache() {
    child_nodes.assign(nodes_used + 1, std::make_pair(-1, -1));
    for (int node = 1; node <= nodes_used; ++node) {
        const int parent = static_cast<int>(parent_nodes(node));
        if (child_nodes[parent].first < 0) child_nodes[parent].first = node;
        else child_nodes[parent].second = node;
    }

    // The partition, observation intervals and base parameters stay fixed in EM.
    // Store only nonzero kernels, rather than a dense observations-by-leaves matrix.
    observation_kernels.resize(n);
    for (int obs = 0; obs < n; ++obs) {
        const double left = X(obs, 0), right = X(obs, 1);
        for (uword j = 0; j < terminal_nodes.n_elem; ++j) {
            const int node = static_cast<int>(terminal_nodes(j));
            if (left == right && has_tail && node == 2) continue;
            const double kernel = left == right ?
                terminal_density(left, left_bound(node), right_bound(node)) :
                terminal_probability(left, right, left_bound(node), right_bound(node));
            if (kernel != 0.0) observation_kernels[obs].emplace_back(node, kernel);
        }
    }
}

void PT::update_expected_allocations() {
    n_obs_allocated.zeros();
    n_obs_allocated(0) = n;
    vec node_likelihood(nodes_used + 1, fill::zeros);
    for (int obs = 0; obs < n; ++obs) {
        const double left = X(obs, 0), right = X(obs, 1);
        if (left == right) {
            // Exact-event membership, including split-point conventions, is fixed.
            for (int node = 1; node <= nodes_used; ++node)
                n_obs_allocated(node) += in_node_obs(node, obs);
            continue;
        }

        node_likelihood.zeros();
        double probability = 0.0;
        for (const auto& entry : observation_kernels[obs]) {
            const double joint = posterior_probs(entry.first) * entry.second;
            node_likelihood(entry.first) = joint;
            probability += joint;
        }
        // A node is a union of its descendant leaves. Aggregate once per
        // observation instead of reintegrating all leaves for every child.
        for (int node = nodes_used; node > 0; --node)
            node_likelihood(static_cast<int>(parent_nodes(node))) += node_likelihood(node);

        for (int node = 1; node <= nodes_used; ++node) {
            if (in_node_obs(node, obs) != 1.0) continue;
            if (left > left_bound(node) && right <= right_bound(node)) {
                // Preserve the original containment shortcut even at zero mass.
                n_obs_allocated(node) += 1.0;
            } else if (probability > SMALL_VALUE) {
                n_obs_allocated(node) += node_likelihood(node) / probability;
            }
        }
    }
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


double PT::base_survival(double x) const {
    if (prior == 2) {
        return exp(-prior_params(0) * x);
    } else if (prior == 3) {
        return exp(-pow(prior_params(0) * x, prior_params(1)));
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
    }
    Rcpp::stop("Invalid prior type.");
    return 0.0;
}

double PT::terminal_probability(double left, double right, double a, double b) const {
    if (has_tail && a == largest_obs && !arma::is_finite(b)) {
        if (prior != 2)
            return right == arma::datum::inf && left <= largest_obs ? 1.0 : 0.0;
        const double lo = max(left, a);
        if (lo >= right) return 0.0;
        const double lambda = prior_params(0);
        // Conditional exponential tail, avoiding division by exp(-lambda * tau).
        return exp(-lambda * (lo - a)) * -expm1(-lambda * (right - lo));
    }
    const double lo = max(left, a);
    const double hi = min(right, b);
    if (lo >= hi) {
        return 0.0;
    }

    // P(left < T < right | a < T < b) under the selected base distribution.
    const double denom = prior == 1 ? b - a : base_survival(a) - base_survival(b);
    const double numerator = prior == 1 ? hi - lo : base_survival(lo) - base_survival(hi);
    if (denom > SMALL_VALUE) {
        return min(max(numerator / denom, 0.0), 1.0);
    }
    return 0.0;
}

double PT::terminal_density(double t, double a, double b) const {
    // Exact events at a split point, including tau, belong to the left interval.
    if (t < a || t > b || (t == a && a > support(0))) {
        return 0.0;
    }
    if (has_tail && prior == 2 && a == largest_obs && !arma::is_finite(b))
        return prior_params(0) * exp(-prior_params(0) * (t - a));
    const double denom = prior == 1 ? b - a : base_survival(a) - base_survival(b);
    if (denom > SMALL_VALUE) {
        return prior == 1 ? 1.0 / denom : base_density(t) / denom;
    }
    return 0.0;
}

double PT::log_likelihood() const {
    double value = 0.0;
    for (int obs = 0; obs < n; ++obs) {
        double probability = 0.0;
        for (const auto& entry : observation_kernels[obs])
            probability += posterior_probs(entry.first) * entry.second;
        if (probability <= 0.0) return -arma::datum::inf;
        value += log(probability);
    }
    return value;
}
List PT::output() {
    // compute survival probabilities and densities at test_seq by mixing over terminal nodes
    survival_probs = vec(test_seq.n_elem, fill::zeros);
    vec density_probs(test_seq.n_elem, fill::zeros);

    for (arma::uword i = 0; i < test_seq.n_elem; ++i) {
        const double t = test_seq(i);
        if (std::isnan(t)) {
            survival_probs(i) = NA_REAL;
            density_probs(i) = NA_REAL;
            continue;
        }
        if (has_tail && prior != 2 && t > largest_obs && arma::is_finite(t)) {
            // Carry the uniform survival curve forward from tau.
            survival_probs(i) = prior == 1 ? posterior_probs(2) : NA_REAL;
            density_probs(i) = NA_REAL;
            continue;
        }
        double s = has_tail && prior != 2 && t <= largest_obs ? posterior_probs(2) : 0.0;
        double d = 0.0;

        for (arma::uword r = 0; r < terminal_nodes.n_rows; ++r) {
            const int node = static_cast<int>(terminal_nodes(r));
            if (has_tail && prior != 2 && node == 2) continue;
            const double a = left_bound(node);
            double b = right_bound(node);
            const double w = posterior_probs(node);

            d += w * terminal_density(t, a, b);

            s += w * terminal_probability(t, arma::datum::inf, a, b);
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
        Named("density_probs") = density_probs,
        Named("tau") = has_tail ? largest_obs : NA_REAL,
        Named("log_likelihood") = log_likelihood(),
        Named("loglik_trace") = loglik_trace,
        Named("em_iterations") = em_iterations,
        Named("converged") = converged
    );
}





void PT::clear(){
    // clear all member variables
    in_node_obs.clear();
    new_posterior_probs.clear();
    parent_nodes.clear();
}
