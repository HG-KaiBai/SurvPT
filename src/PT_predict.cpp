#include "PT.h"
#include <algorithm>
#include <cstdint>
#include <limits>

// A fitted model is ordinary R data, so prediction also works after saveRDS().
// This constructor deliberately does not call init() or EM_step().
PT::PT(const List& model) :
    X(as<mat>(model["X"])), prior(as<int>(model["Prior"])),
    support(as<vec>(model["Support"])), max_M(as<int>(model["Max_M"])),
    prior_params(as<vec>(model["Prior_Params"])),
    allocated(as<int>(model["Allocated"])),
    largest_obs(as<double>(model["Largest_Obs"])),
    uniform_tail_factor(as<double>(model["Uniform_Tail_Factor"])),
    max_EM_iter(0), EM_tol(0.0), verbose(false), flat(as<bool>(model["Flat"])) {
    left_bound = as<vec>(model["Left"]);
    right_bound = as<vec>(model["Right"]);
    alpha = as<vec>(model["Alpha"]);
    posterior_probs = as<vec>(model["Posterior_Probs"]);
    parent_nodes = as<vec>(model["Parent"]) - 1.0;
    terminal_nodes = as<vec>(model["Terminal_Nodes"]) - 1.0;
    n = X.n_rows;
    NODES_NUM = left_bound.n_elem;
    nodes_used = NODES_NUM - 1;
    if (n == 0 || X.n_cols != 2 || NODES_NUM == 0 ||
        right_bound.n_elem != left_bound.n_elem || alpha.n_elem != left_bound.n_elem ||
        posterior_probs.n_elem != left_bound.n_elem || parent_nodes.n_elem != left_bound.n_elem ||
        !left_bound.is_finite() || !alpha.is_finite() || arma::any(alpha < 1.0) ||
        !posterior_probs.is_finite() || arma::any(posterior_probs < 0.0) ||
        parent_nodes(0) != -1.0 || support.n_elem != 2 || prior < 1 || prior > 4 ||
        (prior != 1 && prior_params.n_elem < static_cast<uword>(prior == 3 ? 2 : 1))) {
        Rcpp::stop("Invalid fitted survPT model.");
    }
    child_nodes.assign(NODES_NUM, std::make_pair(-1, -1));
    for (int node = 1; node < NODES_NUM; ++node) {
        const double value = parent_nodes(node);
        if (!std::isfinite(value) || value != std::floor(value) || value < 0 || value >= node) {
            Rcpp::stop("Invalid parent indices in fitted survPT model.");
        }
        auto& children = child_nodes[static_cast<int>(value)];
        if (children.first < 0) children.first = node;
        else if (children.second < 0) children.second = node;
        else Rcpp::stop("Fitted survPT tree must be binary.");
    }
    std::vector<bool> is_leaf(NODES_NUM, false);
    for (double value : terminal_nodes) {
        if (!std::isfinite(value) || value != std::floor(value) || value < 0 || value >= NODES_NUM ||
            is_leaf[static_cast<int>(value)]) {
            Rcpp::stop("Invalid terminal indices in fitted survPT model.");
        }
        is_leaf[static_cast<int>(value)] = true;
    }
    for (int node = 0; node < NODES_NUM; ++node) {
        if (is_leaf[node] != (child_nodes[node].first < 0) ||
            (child_nodes[node].first >= 0 && child_nodes[node].second < 0)) {
            Rcpp::stop("Inconsistent fitted survPT tree.");
        }
    }
    terminal_indices = arma::conv_to<arma::uvec>::from(terminal_nodes);
    if (prior == 4) {
        hybrid_body_mass = -std::expm1(-prior_params(0) * largest_obs);
        hybrid_tail_mass = std::exp(-prior_params(0) * largest_obs);
    }
}

List PT::predict(const vec& x, int n_samples, int burn_in, int thin, bool density) {
    if (n_samples < 0 || burn_in < 0 || thin < 1) {
        Rcpp::stop("Invalid posterior sampling controls.");
    }
    const uword leaves = terminal_indices.n_elem;
    mat basis(leaves, x.n_elem);
    vec estimate(x.n_elem, fill::zeros);
    for (uword k = 0; k < leaves; ++k) {
        const uword node = terminal_indices(k);
        for (uword j = 0; j < x.n_elem; ++j) {
            basis(k, j) = density ? terminal_density(x(j), left_bound(node), right_bound(node)) :
                terminal_probability(x(j), arma::datum::inf, left_bound(node), right_bound(node));
            // Preserve the original prediction formula and leaf summation order.
            estimate(j) += posterior_probs(node) * basis(k, j);
        }
    }
    mat samples(n_samples, x.n_elem);
    if (n_samples == 0 || x.n_elem == 0) {
        return List::create(Named("estimate") = estimate, Named("samples") = samples);
    }

    // Coverage weights depend only on the retained tree and the observations.
    // Store nonzero entries per observation, avoiding a dense n-by-leaf matrix.
    struct Candidate { uword node; double log_coverage; };
    std::vector<std::vector<Candidate>> ambiguous;
    vec fixed_counts(NODES_NUM, fill::zeros);
    for (int obs = 0; obs < n; ++obs) {
        if (obs % 256 == 0) Rcpp::checkUserInterrupt();
        const double left = X(obs, 0), right = X(obs, 1);
        std::vector<Candidate> candidates;
        for (uword k = 0; k < leaves; ++k) {
            const uword node = terminal_indices(k);
            if (left == right) {
                // Exact observations have a known label, including the zero
                // endpoint convention used by the fitting algorithm.
                if ((left > left_bound(node) || (left == 0.0 && left_bound(node) == 0.0)) &&
                    left <= terminal_right_bound(right_bound(node))) {
                    candidates.push_back({node, 0.0});
                }
            } else {
                const double weight = terminal_probability(left, right, left_bound(node), right_bound(node));
                if (!std::isfinite(weight) || weight < 0.0) {
                    Rcpp::stop("Nonfinite coverage weight for observation %d.", obs + 1);
                }
                if (weight > 0.0) candidates.push_back({node, std::log(weight)});
            }
        }
        if (candidates.empty()) {
            Rcpp::stop("Observation %d has zero probability on the fitted tree.", obs + 1);
        }
        if (candidates.size() == 1) fixed_counts(candidates[0].node) += 1.0;
        else ambiguous.push_back(std::move(candidates));
    }

    // Start at the fitted MAP masses. Subsequent Beta draws are unrestricted:
    // the clipping used for MAP optimization is not part of this posterior.
    vec log_mass = arma::log(posterior_probs);
    vec counts(NODES_NUM), leaf_mass(leaves);
    std::vector<double> weights(leaves);
    const std::int64_t sweeps = static_cast<std::int64_t>(burn_in) +
        static_cast<std::int64_t>(n_samples) * thin;
    int saved = 0;
    for (std::int64_t sweep = 0; sweep < sweeps; ++sweep) {
        if (sweep % 32 == 0) Rcpp::checkUserInterrupt();
        counts = fixed_counts;
        for (const auto& candidates : ambiguous) {
            double largest = -std::numeric_limits<double>::infinity();
            for (uword k = 0; k < candidates.size(); ++k) {
                weights[k] = log_mass(candidates[k].node) + candidates[k].log_coverage;
                largest = std::max(largest, weights[k]);
            }
            if (!std::isfinite(largest)) Rcpp::stop("Posterior allocation probabilities vanished.");
            double total = 0.0;
            for (uword k = 0; k < candidates.size(); ++k) {
                weights[k] = std::exp(weights[k] - largest);
                total += weights[k];
            }
            const double draw = R::runif(0.0, total);
            double cumulative = 0.0;
            uword selected = candidates.size() - 1;
            for (uword k = 0; k < candidates.size(); ++k) {
                cumulative += weights[k];
                if (draw < cumulative) { selected = k; break; }
            }
            counts(candidates[selected].node) += 1.0;
        }
        // Each observation contributes one integer count along its entire path.
        for (int node = nodes_used; node > 0; --node) {
            counts(static_cast<uword>(parent_nodes(node))) += counts(node);
        }
        log_mass(0) = 0.0;
        for (int parent = 0; parent < NODES_NUM; ++parent) {
            const int left = child_nodes[parent].first, right = child_nodes[parent].second;
            if (left < 0) continue;
            const double theta = R::rbeta(alpha(left) + counts(left), alpha(right) + counts(right));
            log_mass(left) = log_mass(parent) + std::log(theta);
            log_mass(right) = log_mass(parent) + std::log1p(-theta);
        }
        if (sweep >= burn_in && (sweep - burn_in + 1) % thin == 0) {
            leaf_mass = arma::exp(log_mass.elem(terminal_indices));
            samples.row(saved++) = leaf_mass.t() * basis;
        }
    }
    if (!density) samples = arma::clamp(samples, 0.0, 1.0);
    return List::create(Named("estimate") = estimate, Named("samples") = samples);
}
