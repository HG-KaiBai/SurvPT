// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>
#include <vector>
#include "OPT.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

#define SMALL_VALUE 1e-300
#define LARGE_VALUE 1e+300


// constructor
OPT::OPT(const vec X, const vec delta, const vec test_seq, const double Rho, const vec support, const int max_M, const int max_EM_iter, const double EM_tol, const bool verbose) :
    X(X), delta(delta), test_seq(test_seq), support(support), max_M(max_M), max_EM_iter(max_EM_iter), EM_tol(EM_tol), verbose(verbose) {

    n = X.n_elem;
    if (max_M <= 10){
        NODES_NUM = (1 << (max_M + 1)) - 1; // number of nodes for complete tree
    } else {
        NODES_NUM = (1 << (8 + 1)) - 1; // number of nodes first 8 complete levels
    }
    nodes_used = 0;

    prior_rho = vec(NODES_NUM, fill::zeros);
    for (int i = 0; i < NODES_NUM; i++) {
        prior_rho(i) = Rho; // root has no parent
    }

    largest_obs = max(X);


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

void OPT::init() {
    for (int i = 0; i < NODES_NUM; i++) {
        int level = node_levels(i);

        double split_point = choose_split(level, left_bound(i), right_bound(i));
        if (split_point < 0) {
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
        alpha(left_child) = 0.5;
        alpha(right_child) = 0.5;

        // allocate observations to children
        init_allocations(i, left_child, right_child);
    }
    pos_rho();
    pos_probs();
}

double OPT::choose_split(int level, double left, double right) {
    if (level == max_M) {
        return -1.0; // no split at the last level
    }
    if (arma::is_finite(right)){
        return (left + right) / 2.0; // midpoint split
    } else {
        return left + (largest_obs - support(0)) / (1 << (level + 1)); // exponential split
    }
}

void OPT::init_allocations(int parent, int left_child, int right_child) {
    // indices of observations allocated to the parent node
    arma::uvec parent_idx = arma::find(in_node_obs.row(parent) == 1.0);

    arma::vec weight;
    weight.zeros(parent_idx.n_elem);
    for (arma::uword k = 0; k < parent_idx.n_elem; ++k) {
        const arma::uword obs = parent_idx(k);
        const double obs_time = X(obs);
        if (delta(obs) == 1.0) {
            weight(k) = 1.0;
        } else {
            double b = arma::is_finite(right_bound(parent)) ? right_bound(parent) : (1.1*largest_obs);
            weight(k) = (b - max(left_bound(parent), obs_time)) / (1.1*largest_obs - obs_time);
        }
    }

    conditional_prob_WKM(X.elem(parent_idx), delta.elem(parent_idx), weight, left_bound(parent), right_bound(left_child), right_bound(parent));

    for (arma::uword k = 0; k < parent_idx.n_elem; ++k) {
        const arma::uword obs = parent_idx(k);
        const double obs_time = X(obs);
        const double obs_delta = delta(obs);
        const double s = right_bound(left_child);
        // update for left child
        double temp = 0.0;
        if ((obs_delta == 1.0 && obs_time >= left_bound(parent) && obs_time < right_bound(parent)) ||
            (obs_delta == 0.0 && obs_time < right_bound(parent))) {
            if (obs_time >= s && obs_delta == 1.0) {
                in_node_obs(right_child, obs) = 1.0;
                temp = 0.0;
            } else if (obs_time < s && obs_delta == 1.0) {
                in_node_obs(left_child, obs) = 1.0;
                temp = 1.0;
            } else if (obs_time >= s && obs_delta == 0.0) {
                in_node_obs(right_child, obs) = 1.0;
                temp = 0.0;
            } else if (obs_time < s && obs_delta == 0.0) {
                in_node_obs(left_child, obs) = 1.0;
                in_node_obs(right_child, obs) = 1.0;
                double de = WKM_allocation(max(obs_time, left_bound(parent))) - WKM_allocation(right_bound(parent));
                if (de < SMALL_VALUE) {
                    continue;
                } else {
                    temp = (WKM_allocation(max(obs_time, left_bound(parent))) - WKM_allocation(s)) / de;
                }
            }
        }
        n_obs_allocated(left_child) += temp;
        n_obs_allocated(right_child) += (1.0 - temp);
    }
}

void OPT::conditional_prob_WKM(const arma::vec& Y,
                                const arma::vec& Delta,
                                const arma::vec& w,
                                const double a,
                                const double s,
                                const double b) { // 1-based index (R-style)

  const arma::uword n_inside = Y.n_elem;
  if (Delta.n_elem != n_inside || w.n_elem != n_inside) stop("Lengths mismatch.");

  // sort by time
  arma::uvec ord = arma::sort_index(Y, "ascend");
  arma::vec  ys  = Y.elem(ord);
  arma::vec  ds  = Delta.elem(ord);
  arma::vec  ws  = w.elem(ord);

  // group equal times, collect event weights dk and group-end index
  std::vector<double> times;   times.reserve(n_inside);
  std::vector<double> dks;     dks.reserve(n_inside);
  std::vector<double> rks;    rks.reserve(n_inside);

  arma::uword u = 0;
  while (u < n_inside) {
    if (ds[u] == 1.0) {
        if (!times.empty() && ys[u] == times.back()) {
            dks.back() += ws[u];
        } else {
            times.push_back(ys[u]);
            dks.push_back(ws[u]);
            rks.push_back(0.0); // initialize, will include ws[u] below
        }
    }
    for (auto& rk : rks) {
        rk += ws[u];
    }
    ++u;
  }

   const arma::uword G = times.size();

    // compute survival after each event
    std::vector<double> Se; Se.reserve(G);
    double Srun = 1.0;
    for (arma::uword g = 0; g < G; ++g) {
        double dk = dks[g];
        double rk = rks[g];
        double frac = min(max(dk / rk, 0.0), 1.0);
        Srun *= (1.0 - frac);
        Se.push_back(Srun);
    }
    node_km = arma::vec(Se);
    node_times = arma::vec(times);
}

// Return Kaplan–Meier survival S(t) at a given time.
// node_km holds survival values after each unique event time in ascending order.
// event_times must be the corresponding sorted unique event times.
double OPT::WKM_allocation(double t) {
    if (node_km.is_empty()) {
        return 1.0;
    }
    if (node_times.n_elem != node_km.n_elem) {
        Rcpp::stop("WKM_allocation: event_times length must equal node_km length.");
    }

    // Find the first event time strictly greater than 'time'
    arma::uvec idx = arma::find(node_times > t, 1, "first");

    if (idx.is_empty()) {
        // time is at or beyond the last event -> return last survival value
        return node_km(node_km.n_elem - 1);
    }

    arma::uword g = idx(0);
    if (g == 0) {
        // time is before the first event
        return 1.0;
    }

    // Otherwise, return survival right after the last event <= time
    return node_km(g - 1);
}


void OPT::pos_rho() {
    // Compute posterior stopping probabilities rho(A | n^) for all nodes
    // using the binary-split specialization of:
    // Phi(A) = rho(A)*Phi0(A) + (1-rho(A)) * [ D(n+alpha)/D(alpha) * prod_l Phi(A_l) ]
    //
    // where:
    // - Phi0(A) = mu(A)^(-n_A), with mu(A) the Lebesgue measure (interval width)
    // - D(t) = prod Gamma(t_i) / Gamma(sum t_i)
    // - For our tree, K=2 (left/right child), M(A)=1, lambda=1.
    //
    // We compute Phi(A) recursively (post-order), then set:
    // rho_post(A) = rho_prior(A) * Phi0(A) / Phi(A).

    // Snapshot prior rho; write back posterior into member rho.
    arma::vec new_rho(NODES_NUM, arma::fill::zeros);

    auto width_of = [&](int node) -> double {
        double a = left_bound(node);
        double b = right_bound(node);
        if (!arma::is_finite(b)) {
            b = 1.1 * largest_obs;  // consistent with rest of code
        }
        double w = b - a;
        if (!(w > SMALL_VALUE)) w = SMALL_VALUE; // guard
        return w;
    };

    auto children_of = [&](int node) -> std::pair<int,int> {
        arma::uvec ch = arma::find(parent_nodes.head(nodes_used + 1) == static_cast<double>(node));
        if (ch.n_elem < 2) return {-1, -1};
        int l = static_cast<int>(ch.min());
        int r = static_cast<int>(ch.max());
        return {l, r};
    };

    // Recursive lambda using log-domain to avoid overflow
    auto compute_phi = [&](auto&& self, int node) -> double {
        const double INF  = std::numeric_limits<double>::infinity();
        const double NINF = -INF;

        auto logsumexp2 = [&](double x, double y) -> double {
            if (std::isinf(x) && x < 0 && std::isinf(y) && y < 0) return NINF;
            if (std::isinf(x) && x > 0 && std::isinf(y) && y > 0) return INF;
            double d = std::fabs(x - y);
            if (d >= 100.0) return std::max(x, y);
            double m = std::min(x, y);
            return m + std::log1p(std::exp(d));
        };

        // log Phi0(A) = -n_A * log(width(A))
        const double nA = std::max(n_obs_allocated(node), 0.0);
        const double log_phi0 = -nA * std::log(width_of(node));

        // children (binary if present)
        auto [left_child, right_child] = children_of(node);
        const bool has_children = (left_child >= 0 && right_child >= 0);

        // prior rho at this node
        double p = prior_rho(node);
        p = std::min(std::max(p, 0.0), 1.0);
        const double lp = (p > 0.0) ? std::log(p) : NINF;
        const double lq = (p < 1.0) ? std::log1p(-p) : NINF;

        if (!has_children) {
            // cannot split further: Phi(A) = p * Phi0(A) and rho_post = 1
            const double log_phiA = lp + log_phi0; // may be -Inf if p==0
            new_rho(node) = 1.0;
            return log_phiA;
        }

        // Dirichlet integral ratio in log domain
        const double nL = std::max(n_obs_allocated(left_child), 0.0);
        const double nR = std::max(n_obs_allocated(right_child), 0.0);
        const double aL = std::max(alpha(left_child), SMALL_VALUE);
        const double aR = std::max(alpha(right_child), SMALL_VALUE);
        const double logD_alpha  = std::lgamma(aL) + std::lgamma(aR) - std::lgamma(aL + aR);
        const double logD_nalpha = std::lgamma(nL + aL) + std::lgamma(nR + aR) - std::lgamma(nL + aL + nR + aR);
        const double log_dir_ratio = logD_nalpha - logD_alpha;

        // recurse on children in log-space
        const double log_phi_left  = self(self, left_child);
        const double log_phi_right = self(self, right_child);
        const double log_split = log_dir_ratio + log_phi_left + log_phi_right;

        // log Phi(A) = logsumexp( log(p) + log_phi0, log(1-p) + log_split )
        const double log_num = lp + log_phi0;
        const double log_phiA = logsumexp2(log_num, lq + log_split);

        // posterior rho = (p * Phi0) / Phi(A) in log-space
        double rho_post;
        if (std::isinf(log_phiA) && log_phiA < 0) {
            // degenerate case: both terms -Inf, fall back to prior
            rho_post = p;
        } else {
            const double d = log_num - log_phiA;
            if      (d <= -50.0) rho_post = 0.0;
            else if (d >=  50.0) rho_post = 1.0;
            else                 rho_post = std::exp(d);
            rho_post = std::min(std::max(rho_post, 0.0), 1.0);
        }
        new_rho(node) = rho_post;
        return log_phiA;
    };

    // Start at root
    if (nodes_used >= 0) {
        compute_phi(compute_phi, 0);
    }
    rho = new_rho;
}

void OPT::pos_probs() {
    // Compute posterior mass for terminal nodes using the iterative recursions:
    // a_i = (mu(A_i)/mu(A_{i-1})) * a_{i-1} + [alpha_i/(alpha_i+alpha_i')] * rho_i * b_{i-1}
    // b_i = [alpha_i/(alpha_i+alpha_i')] * (1 - rho_i) * b_{i-1}
    //
    // where for each step we move from parent A_{i-1} to child A_i along the path from root to a terminal node.
    // We set a_0 = rho(root), b_0 = 1 - rho(root). The final posterior mass for a terminal node is a_last.

    // keep the working vector used by EM initialization unchanged
    new_posterior_probs.zeros();
    new_posterior_probs(0) = 1.0; // root node starts with probability 1 for EM

    // helper: interval width (Lebesgue measure), with infinite right bound approximated
    auto mu = [&](int node) -> double {
        double a = left_bound(node);
        double b = right_bound(node);
        if (!arma::is_finite(b)) {
            b = 1.1 * largest_obs;
        }
        double w = b - a;
        if (!(w > SMALL_VALUE)) w = SMALL_VALUE;
        return w;
    };

    // helper: sibling index of 'child' under 'parent' (binary tree). Returns -1 if unavailable.
    auto sibling_of = [&](int parent, int child) -> int {
        arma::uvec ch = arma::find(parent_nodes.head(nodes_used + 1) == static_cast<double>(parent));
        if (ch.n_elem < 2) return -1;
        int c0 = static_cast<int>(ch(0));
        int c1 = static_cast<int>(ch(1));
        return (c0 == child) ? c1 : c0;
    };

    posterior_probs.zeros();

    // iterate over terminal nodes
    arma::uvec terms = arma::conv_to<arma::uvec>::from(terminal_nodes);
    for (arma::uword tIdx = 0; tIdx < terms.n_elem; ++tIdx) {
        int leaf = static_cast<int>(terms(tIdx));

        // collect path from leaf up to root (reverse order): [leaf, parent, ..., root]
        std::vector<int> revPath;
        int cur = leaf;
        while (cur != -1) {
            revPath.push_back(cur);
            cur = static_cast<int>(parent_nodes(cur));
        }

        // initialize at root: a_0 = rho(root), b_0 = 1 - rho(root)
        int root = revPath.back();
        double a = std::min(std::max(rho(root), 0.0), 1.0);
        double b = 1.0 - a;

        // walk from root to leaf using revPath indices: parent = revPath[j], child = revPath[j-1]
        for (int j = static_cast<int>(revPath.size()) - 1; j >= 1; --j) {
            int parent = revPath[j];
            int child  = revPath[j - 1];

            // width ratio
            double w_ratio = mu(child) / mu(parent);

            // Dirichlet splitting weight to this child
            double a_child = std::max(n_obs_allocated(child) + alpha(child), SMALL_VALUE);
            int sib = sibling_of(parent, child);
            double a_sib = (sib >= 0) ? std::max(n_obs_allocated(sib) + alpha(sib), SMALL_VALUE) : a_child;
            double denom = a_child + a_sib;
            double pi = (denom > SMALL_VALUE) ? (a_child / denom) : 0.5;

            // stopping probability at the child
            double r = std::min(std::max(rho(child), 0.0), 1.0);

            // update
            double a_new = w_ratio * a + pi * r * b;
            double b_new = pi * (1.0 - r) * b;
            a = a_new;
            b = b_new;
        }

        // terminal node mass is a_last
        posterior_probs(leaf) = a + b;
    }
}


void OPT::iter_step() {
    bool condition = true;
    int iter = 0;
    if (verbose) {
        std::cout << "Starting OPT iterations..." << std::endl;
    }
    auto start = std::chrono::high_resolution_clock::now();
    while (condition && iter < max_EM_iter) {
        arma::uvec idx = arma::conv_to<arma::uvec>::from(terminal_nodes);
        vec old_terminal_probs = posterior_probs.elem(idx);
        re_allocation();
        double diff = arma::norm(posterior_probs.elem(idx) - old_terminal_probs, 2);
        // check convergence
        condition = (diff > EM_tol);
        if (verbose) {
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Iteration " << iter + 1 << ": diff = " << diff << ", time elapsed = " << elapsed.count() << " seconds." << std::endl;
        }
        iter++;
    }
}


void OPT::re_allocation() {
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
    }
    pos_rho();
    pos_probs();
}


void OPT::update_allocations(int parent, int left_child, int right_child) {
    // reset allocations for children
    in_node_obs.row(left_child).zeros();
    in_node_obs.row(right_child).zeros();
    n_obs_allocated(left_child) = 0.0;
    n_obs_allocated(right_child) = 0.0;

    // indices of observations allocated to the parent node
    arma::uvec parent_idx = arma::find(in_node_obs.row(parent) == 1.0);

    for (arma::uword k = 0; k < parent_idx.n_elem; ++k) {
        const arma::uword obs = parent_idx(k);
        const double obs_time = X(obs);
        const double obs_delta = delta(obs);
        const double s = right_bound(left_child);
        // update for left child
        if ((obs_delta == 1.0 && obs_time >= left_bound(parent) && obs_time < right_bound(parent)) ||
            (obs_delta == 0.0 && obs_time < right_bound(parent))) {

            if (obs_time >= s && obs_delta == 1.0) {
                in_node_obs(right_child, obs) = 1.0;
                n_obs_allocated(right_child) += 1.0;
            } else if (obs_time < s && obs_delta == 1.0) {
                in_node_obs(left_child, obs) = 1.0;
                n_obs_allocated(left_child) += 1.0;
            } else if (obs_time >= s && obs_delta == 0.0) {
                in_node_obs(right_child, obs) = 1.0;
                n_obs_allocated(right_child) += update_individual_allocation(obs_time, left_bound(right_child), right_bound(right_child));
            } else if (obs_time < s && obs_delta == 0.0) {
                in_node_obs(left_child, obs) = 1.0;
                in_node_obs(right_child, obs) = 1.0;
                n_obs_allocated(left_child) += update_individual_allocation(obs_time, left_bound(left_child), right_bound(left_child));
                n_obs_allocated(right_child) += update_individual_allocation(obs_time, left_bound(right_child), right_bound(right_child));
            }
        }
    }
    if (abs(n_obs_allocated(left_child) + n_obs_allocated(right_child) - n_obs_allocated(parent)) > 1e-5) {
        std::cout << "Parent node: " << n_obs_allocated(parent) << ", Left child: " << n_obs_allocated(left_child) << ", Right child: " << n_obs_allocated(right_child) << std::endl;
        std::cout << "Sum of allocations exceeds parent sample size." << std::endl;
    }
}

double OPT::update_individual_allocation(double left_obs, double left_node, double right_node) {
    // Update individual allocation based on the terminal nodes and previous posterior probabilities
    double joint_left = max(left_obs, left_node);
    double joint_right = right_node;
    double p_joint = 0.0;
    double p_data = 0.0;

    // Loop over all terminal nodes
    for (arma::uword t = 0; t < terminal_nodes.n_rows; ++t) {
        int node = static_cast<int>(terminal_nodes(t));
        double node_left = left_bound(node);
        double node_right = (arma::is_finite(right_bound(node))) ? right_bound(node) : 1.1 * largest_obs; // handle infinite right bound

        // Compute overlap between [joint_left, joint_right] and [node_left, node_right]
        double overlap_left = std::max(joint_left, node_left);
        double overlap_right = std::min(joint_right, node_right);

        if (overlap_left < overlap_right) {
            p_joint += (overlap_right - overlap_left) / (node_right - node_left) * posterior_probs(node);
        }

        // Compute overlap between [left_obs, right_obs] and [node_left, node_right]
        double obs_overlap_left = std::max(left_obs, node_left);
        double obs_overlap_right = node_right;

        if (obs_overlap_left < obs_overlap_right) {
            p_data += (obs_overlap_right - obs_overlap_left) / (node_right - node_left) * posterior_probs(node);
        }
    }

    // Avoid division by zero
    if (p_data > SMALL_VALUE) {
        return p_joint / p_data;
    }

    return 0.0;
}


void OPT::update_tree_size() {
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


List OPT::output() {
    // compute survival probabilities at test_seq by mixing over terminal nodes
    survival_probs = vec(test_seq.n_elem, fill::zeros);

    for (arma::uword i = 0; i < test_seq.n_elem; ++i) {
        const double t = test_seq(i);
        double s = 0.0;

        for (arma::uword r = 0; r < terminal_nodes.n_rows; ++r) {
            const int node = static_cast<int>(terminal_nodes(r));
            const double a = left_bound(node);
            double b = right_bound(node);
            const double w = posterior_probs(node);

            // approximate infinite right bound using 1.1 * largest_obs
            if (!arma::is_finite(b)) {
                b = 1.1 * largest_obs;
            }

            if (t < a) {
                s += w;
            } else if (t >= b) {
                // contributes 0
            } else {
                const double width = b - a;
                if (width > SMALL_VALUE) {
                    s += w * (b - t) / width;
                }
            }
        }

        // clamp numeric noise
        survival_probs(i) = s;
    }

    return List::create(
        Named("left_bound") = left_bound.head(nodes_used+1),
        Named("right_bound") = right_bound.head(nodes_used+1),
        Named("node_levels") = node_levels.head(nodes_used+1),
        Named("rho") = rho.head(nodes_used+1),
        Named("n_obs_allocated") = n_obs_allocated.head(nodes_used+1),
        Named("posterior_probs") = posterior_probs.head(nodes_used+1),
        Named("terminal_nodes") = terminal_nodes,
        Named("survival_probs") = survival_probs
    );
}





void OPT::clear(){
    // clear all member variables
    in_node_obs.clear();
    new_posterior_probs.clear();
    parent_nodes.clear();
    node_km.clear();
    node_times.clear();
}
