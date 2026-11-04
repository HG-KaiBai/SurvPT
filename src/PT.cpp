// [[Rcpp::depends(RcppArmadillo)]]
#include "PT.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace {
const double INF = std::numeric_limits<double>::infinity();
double log_add(double a, double b) {
    if (a == -INF) return b;
    if (b == -INF) return a;
    const double hi = std::max(a, b);
    return hi + std::log1p(std::exp(std::min(a, b) - hi));
}
bool contains(double t, double a, double b) {
    return (t > a || (t == 0 && a == 0)) && t <= b;
}
}

PT::PT(const arma::mat& X, int prior, int max_depth, const arma::vec& parameters,
       int allocated, int max_iter, double tolerance, bool verbose,
       double gamma, double kappa, double s0, double min_count,
       double epsilon, bool flat) :
    data(X), parameters(parameters), prior(prior),
    max_depth(max_depth), allocated(allocated), max_iter(max_iter),
    tolerance(tolerance), gamma(gamma), kappa(kappa), s0(s0),
    min_count(min_count), epsilon(epsilon), tau(0), h_tau(0),
    verbose(verbose), flat(flat) {
    if (data.n_rows == 0 || data.n_cols < 2)
        Rcpp::stop("X must have observations and at least two endpoint columns.");
    if (prior < 1 || prior > 3 || allocated < 0 || allocated > 1)
        Rcpp::stop("Invalid prior or allocation mode.");
    if (max_depth < 0 || max_depth > 30 || max_iter < 0 ||
        !std::isfinite(tolerance) || tolerance < 0)
        Rcpp::stop("Invalid depth, iteration cap, or EM tolerance.");
    if (!std::isfinite(gamma) || gamma < 0 || !std::isfinite(kappa) || kappa <= 0 ||
        !std::isfinite(min_count) || min_count < 0 ||
        !std::isfinite(epsilon) || epsilon <= 0 || epsilon >= 0.5)
        Rcpp::stop("Invalid gamma, kappa, min_count, or split_epsilon.");
    if ((prior == 2 && parameters.n_elem != 1) || (prior == 3 && parameters.n_elem != 2))
        Rcpp::stop("Exponential needs a rate; Weibull needs c(rate, shape).");
    if (prior != 1 && (!parameters.is_finite() || arma::any(parameters <= 0)))
        Rcpp::stop("Base parameters must be finite and positive.");
    for (arma::uword i = 0; i < data.n_rows; ++i) {
        double l = data(i, 0), r = data(i, 1);
        if (!std::isfinite(l) || l < 0 || std::isnan(r) || r < l)
            Rcpp::stop("Endpoints must satisfy 0 <= L <= R, with finite L and no missing values.");
        tau = std::max(tau, l);
        if (std::isfinite(r)) tau = std::max(tau, r);
    }
    if (tau <= 0) Rcpp::stop("At least one positive finite endpoint is required.");
    if (prior != 1) {
        h_tau = hazard(tau);
        if (!std::isfinite(h_tau) || h_tau <= 0)
            Rcpp::stop("The base cumulative hazard at the largest endpoint is not numerically representable.");
    }
    if (s0 == -1) {
        // Equation (4); Weibull extends the manuscript's uniform/exponential choices.
        double log_extra = std::log(kappa);
        if (prior != 1) {
            const double shape = prior == 3 ? parameters(1) : 1.0;
            log_extra += std::log(shape) + std::log(h_tau) - h_tau - std::log(-std::expm1(-h_tau));
        }
        this->s0 = log_extra >= 0 ? 1.0 / (1.0 + std::exp(-log_extra)) :
            std::exp(log_extra) / (1.0 + std::exp(log_extra));
    }
    if (!std::isfinite(this->s0) || this->s0 <= 0 || this->s0 >= 1)
        Rcpp::stop("s0 must be strictly between zero and one; if automatic calibration underflows, supply s0 explicitly.");
    build_tree();
    prepare_coverage();
    mass.set_size(nodes.size()); log_mass.set_size(nodes.size());
    counts.zeros(nodes.size()); theta.zeros(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        log_mass(i) = nodes[i].log_center; mass(i) = std::exp(log_mass(i));
    }
    for (int parent : splits)
        theta(parent) = std::exp(nodes[nodes[parent].left_child].log_center - nodes[parent].log_center);
    allocate(counts, false);
}

double PT::hazard(double t) const {
    return prior == 2 ? parameters(0) * t : std::pow(parameters(0) * t, parameters(1));
}
double PT::log_body_probability(double a, double b) const {
    a = std::max(0.0, a); b = std::min(tau, b);
    if (b <= a) return -INF;
    if (prior == 1) return std::log(b-a) - std::log(tau);
    const double ha = hazard(a), width = hazard(b) - ha;
    if (width <= 0) return -INF;
    return -ha + std::log(-std::expm1(-width)) - std::log(-std::expm1(-h_tau));
}
double PT::log_reference_probability(double a, double b) const {
    const double body = std::log1p(-s0) + log_body_probability(a, b);
    return b == INF ? log_add(body, std::log(s0)) : body;
}
double PT::reference_count(double a, double b) const {
    double total = 0;
    for (arma::uword i = 0; i < data.n_rows; ++i) {
        const double l = data(i,0), r = data(i,1);
        if (l == r) { total += contains(l,a,b); continue; }
        const double lo = std::max(l,a), hi = std::min(r,b);
        if (hi <= lo) continue;
        total += std::exp(log_reference_probability(lo,hi) - log_reference_probability(l,r));
    }
    return total;
}
double PT::split_point(const Node& node) const {
    if (node.depth >= max_depth || node.initial_count < min_count) return -1;
    const double a = node.left, b = std::min(tau,node.right);
    if (b <= a) return -1;
    double cut;
    if (prior == 1) cut = a + (b-a)/2;
    else {
        const double ha = hazard(a), width = hazard(b)-ha;
        const double increment = -std::log1p(0.5*std::expm1(-width));
        cut = prior == 2 ? a + increment/parameters(0) :
            std::pow(ha+increment, 1/parameters(1))/parameters(0);
    }
    return std::isfinite(cut) && cut > a && cut < b ? cut : -1;
}
void PT::add_children(int parent, double cut) {
    if (nodes.size() >= 2000000) Rcpp::stop("Tree too large; reduce max_M or increase min_count.");
    const Node p = nodes[parent];
    const int left = nodes.size(), right = left+1, depth = p.depth+1;
    const double log_left = log_reference_probability(p.left, cut);
    const double log_right = log_reference_probability(cut, p.right);
    const double scale = flat ? 0.0 : 2*gamma*depth*depth;
    nodes.push_back({p.left, cut, log_left, reference_count(p.left,cut),
        std::max(1.0,scale*std::exp(log_left-p.log_center)), depth,parent,-1,-1});
    nodes.push_back({cut,p.right,log_right,reference_count(cut,p.right),
        std::max(1.0,scale*std::exp(log_right-p.log_center)),depth,parent,-1,-1});
    nodes[parent].left_child = left; nodes[parent].right_child = right;
}
void PT::build_tree() {
    nodes.push_back({0,INF,0,static_cast<double>(data.n_rows),1,0,-1,-1,-1});
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (i % 1024 == 0) Rcpp::checkUserInterrupt();
        const double cut = split_point(nodes[i]);
        if (cut > 0) add_children(i,cut);
    }
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].left_child < 0 && nodes[i].right == INF) {
            add_children(i,tau); break;
        }
    }
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].left_child < 0) leaves.push_back(i); else splits.push_back(i);
    }
}
double PT::log_conditional_probability(double l, double r, double a, double b) const {
    if (b == INF) return r == INF ? 0.0 : -INF;
    const double lo = std::max(l,a), hi = std::min(r,b);
    if (hi <= lo) return -INF;
    if (prior == 1) return std::log(hi-lo)-std::log(b-a);
    const double ha = hazard(a), hlo = hazard(lo);
    const double width = hazard(b)-ha, overlap = hazard(hi)-hlo;
    if (width <= 0 || overlap <= 0) return -INF;
    return std::min(0.0, ha-hlo + std::log(-std::expm1(-overlap)) - std::log(-std::expm1(-width)));
}
double PT::log_weight(int obs, int leaf) const {
    if (exact_leaf[obs] >= 0) return exact_leaf[obs] == leaf ? 0.0 : -INF;
    const Node& node = nodes[leaves[leaf]];
    return log_conditional_probability(data(obs,0),data(obs,1),node.left,node.right);
}
double PT::weight(int obs, int leaf) const {
    return cached ? log_coverage(leaf,obs) : log_weight(obs,leaf);
}
void PT::prepare_coverage() {
    exact_leaf.assign(data.n_rows,-1);
    for (arma::uword i=0; i<data.n_rows; ++i) {
        if (data(i,0) != data(i,1)) continue;
        for (size_t j=0; j<leaves.size(); ++j) {
            const Node& node = nodes[leaves[j]];
            if (contains(data(i,0),node.left,node.right)) { exact_leaf[i]=j; break; }
        }
        if (exact_leaf[i] < 0) Rcpp::stop("An exact observation is outside the retained tree.");
    }
    cached = leaves.size() <= (64u*1024u*1024u/sizeof(double))/data.n_rows;
    if (cached) {
        log_coverage.set_size(leaves.size(),data.n_rows);
        for (arma::uword i=0; i<data.n_rows; ++i)
            for (size_t j=0; j<leaves.size(); ++j) log_coverage(j,i)=log_weight(i,j);
    }
    coverage_prepared = true;
}
void PT::aggregate(arma::vec& values) const {
    for (int i=static_cast<int>(nodes.size())-1; i>0; --i) values(nodes[i].parent) += values(i);
}
double PT::allocate(arma::vec& values, bool draw_labels) const {
    values.zeros(nodes.size());
    std::vector<double> probabilities(leaves.size());
    double likelihood = 0;
    for (arma::uword i=0; i<data.n_rows; ++i) {
        if (i % 1024 == 0) Rcpp::checkUserInterrupt();
        if (exact_leaf[i] >= 0) {
            const int node=leaves[exact_leaf[i]];
            values(node) += 1; likelihood += log_mass(node); continue;
        }
        double largest = -INF;
        for (size_t j=0; j<leaves.size(); ++j) {
            probabilities[j] = log_mass(leaves[j]) + weight(i,j);
            largest = std::max(largest,probabilities[j]);
        }
        if (!std::isfinite(largest)) Rcpp::stop("An observation has zero probability under the current tree.");
        double total=0;
        for (double& p : probabilities) { p=std::exp(p-largest); total+=p; }
        likelihood += largest+std::log(total);
        if (draw_labels) {
            const double u=R::runif(0,total);
            double cumulative=0; int chosen=-1;
            for (size_t j=0; j<leaves.size(); ++j) {
                if (probabilities[j] <= 0) continue;
                chosen=j; cumulative+=probabilities[j];
                if (u < cumulative) break;
            }
            if (chosen < 0) Rcpp::stop("No compatible leaf for a sampled label.");
            values(leaves[chosen]) += 1;
        } else {
            for (size_t j=0; j<leaves.size(); ++j) values(leaves[j])+=probabilities[j]/total;
        }
    }
    aggregate(values);
    return likelihood; // Omits exact-density factors constant in the splits.
}
void PT::propagate() {
    mass(0)=1; log_mass(0)=0;
    for (int parent : splits) {
        const double q=theta(parent);
        const int left=nodes[parent].left_child, right=nodes[parent].right_child;
        log_mass(left)=log_mass(parent)+std::log(q);
        log_mass(right)=log_mass(parent)+std::log1p(-q);
        mass(left)=std::exp(log_mass(left)); mass(right)=std::exp(log_mass(right));
    }
}
void PT::update_splits(const arma::vec& values) {
    for (int parent : splits) {
        const int left=nodes[parent].left_child, right=nodes[parent].right_child;
        const double a=nodes[left].alpha+values(left);
        const double b=nodes[right].alpha+values(right);
        const double q=a+b>2 ? (a-1)/(a+b-2) : 0.5;
        // Both children share level m; bound their conditional probabilities
        // by epsilon * 2^(-m), including the final body/tail split.
        const double bound=std::ldexp(epsilon,-nodes[left].depth);
        theta(parent)=std::min(1-bound,std::max(bound,q));
    }
    propagate();
}
double PT::objective(double log_likelihood) const {
    double result=log_likelihood;
    for (int parent : splits) {
        const double a=nodes[nodes[parent].left_child].alpha;
        const double b=nodes[nodes[parent].right_child].alpha;
        if (a>1) result+=(a-1)*std::log(theta(parent));
        if (b>1) result+=(b-1)*std::log1p(-theta(parent));
    }
    return result;
}
PT::Stage PT::run_em() {
    Stage result;
    double likelihood=allocate(counts,false);
    result.objective.push_back(objective(likelihood));
    for (int iter=1; iter<=max_iter; ++iter) {
        Rcpp::checkUserInterrupt();
        const arma::vec previous=mass;
        update_splits(counts);
        double delta=0;
        for (int leaf : leaves) delta+=std::pow(mass(leaf)-previous(leaf),2);
        result.difference=std::sqrt(delta);
        likelihood=allocate(counts,false); // Refresh final counts at final masses.
        result.objective.push_back(objective(likelihood));
        result.iterations=iter;
        if (verbose && (iter==1 || iter%100==0 || result.difference<=tolerance))
            Rcpp::Rcout << (flat ? "Flat-prior" : "MAP") << " iteration " << iter
                        << ", terminal-mass change = " << result.difference << "\n";
        if (result.difference<=tolerance) { result.converged=true; break; }
    }
    return result;
}
void PT::fit() {
    if (allocated == 0) {
        // Historical overlap-count option, separate from the paper's MAP-EM.
        counts.zeros(nodes.size()); counts(0)=data.n_rows;
        for (size_t j=1; j<nodes.size(); ++j) {
            const Node& node=nodes[j];
            for (arma::uword i=0; i<data.n_rows; ++i) {
                const double l=data(i,0),r=data(i,1);
                if (prior==1 && r==INF) counts(j)+=(node.right==INF);
                else counts(j)+=l==r ? contains(l,node.left,node.right) :
                    std::max(l,node.left)<std::min(r,node.right);
            }
        }
        update_splits(counts); return;
    }
    if (max_iter == 0) return;
    // The first E-step uses the centering reference masses. The selected
    // prior is fixed for the entire fit and the subsequent Gibbs sampler.
    map_stage=run_em();
}
double PT::conditional_density(double t,double a,double b) const {
    if (b==INF || !contains(t,a,b)) return 0;
    if (prior==1) return 1/(b-a);
    const double ha=hazard(a),width=hazard(b)-ha;
    if (width<=0) return 0;
    double log_hazard=std::log(parameters(0));
    if (prior==3 && parameters(1)!=1) {
        const double shape=parameters(1);
        log_hazard+=std::log(shape)+(shape-1)*std::log(parameters(0)*t);
    }
    return std::exp(log_hazard+ha-hazard(t)-std::log(-std::expm1(-width)));
}
void PT::predict(arma::vec& survival,arma::vec& density) const {
    if (predictions_cached) {
        arma::vec leaf_mass(leaves.size());
        for (size_t j=0; j<leaves.size(); ++j) leaf_mass(j)=mass(leaves[j]);
        survival=survival_kernel*leaf_mass;
        density=density_kernel*leaf_mass;
        return;
    }
    survival.zeros(grid.n_elem); density.zeros(grid.n_elem);
    for (arma::uword i=0; i<grid.n_elem; ++i) {
        const double t=grid(i);
        if (t>tau) { survival(i)=NA_REAL; density(i)=NA_REAL; continue; }
        for (int leaf : leaves) {
            const Node& node=nodes[leaf];
            survival(i)+=mass(leaf)*std::exp(log_conditional_probability(t,INF,node.left,node.right));
            density(i)+=mass(leaf)*conditional_density(t,node.left,node.right);
        }
    }
}
void PT::prepare_predictions() {
    predictions_cached = grid.n_elem > 0 &&
        leaves.size() <= (64u*1024u*1024u/(2*sizeof(double)))/grid.n_elem;
    if (!predictions_cached) return;
    survival_kernel.set_size(grid.n_elem,leaves.size());
    density_kernel.set_size(grid.n_elem,leaves.size());
    for (arma::uword i=0; i<grid.n_elem; ++i) {
        for (size_t j=0; j<leaves.size(); ++j) {
            const Node& node=nodes[leaves[j]];
            survival_kernel(i,j)=grid(i)>tau ? NA_REAL :
                std::exp(log_conditional_probability(grid(i),INF,node.left,node.right));
            density_kernel(i,j)=grid(i)>tau ? NA_REAL :
                conditional_density(grid(i),node.left,node.right);
        }
    }
}
Rcpp::List PT::sample(int n_samples,int burn_in,int thin) {
    if (allocated!=1) Rcpp::stop("DA-Gibbs requires allocated = 1.");
    if (n_samples<=0 || burn_in<0 || thin<1)
        Rcpp::stop("Require n_samples > 0, burn_in >= 0 and thin >= 1.");
    if (!coverage_prepared) prepare_coverage();
    const int kept=n_samples;
    const std::int64_t sweeps=static_cast<std::int64_t>(burn_in)+
        static_cast<std::int64_t>(n_samples)*thin;
    const double cells=static_cast<double>(kept)*(leaves.size()+splits.size());
    if (cells>512.0*1024*1024/sizeof(double))
        Rcpp::stop("Requested draws exceed 512 MiB; reduce n_samples or tree depth.");
    arma::mat leaf_draws(kept,leaves.size()),branch_draws(kept,splits.size());
    Rcpp::NumericVector iterations(kept);
    const arma::vec saved_mass=mass,saved_log_mass=log_mass,saved_theta=theta;
    arma::vec integer_counts;
    int saved=0;
    for (std::int64_t g=1; g<=sweeps; ++g) {
        Rcpp::checkUserInterrupt();
        allocate(integer_counts,true);
        for (int parent : splits) {
            const int left=nodes[parent].left_child,right=nodes[parent].right_child;
            // Unrestricted conjugate draws: never apply the MAP clipping bound.
            theta(parent)=R::rbeta(nodes[left].alpha+integer_counts(left),
                                  nodes[right].alpha+integer_counts(right));
        }
        propagate();
        if (g>burn_in && (g-burn_in)%thin==0) {
            for (size_t j=0; j<leaves.size(); ++j) leaf_draws(saved,j)=mass(leaves[j]);
            for (size_t j=0; j<splits.size(); ++j) branch_draws(saved,j)=theta(splits[j]);
            iterations[saved]=g; ++saved;
        }
        if (verbose && g%500==0) Rcpp::Rcout << "DA-Gibbs sweep " << g << "/" << sweeps << "\n";
    }
    mass=saved_mass; log_mass=saved_log_mass; theta=saved_theta;
    return Rcpp::List::create(Rcpp::Named("terminal_probs")=leaf_draws,
        Rcpp::Named("branch_probs")=branch_draws,Rcpp::Named("iterations")=iterations);
}
Rcpp::List PT::stage_output(const Stage& stage) const {
    return Rcpp::List::create(Rcpp::Named("iterations")=stage.iterations,
        Rcpp::Named("converged")=stage.converged,Rcpp::Named("difference")=stage.difference,
        Rcpp::Named("log_posterior_kernel")=stage.objective);
}
Rcpp::List PT::output() const {
    const size_t n=nodes.size();
    arma::vec left(n),right(n),depth(n),alpha(n),center(n),center_branch(n),branch(n),initial(n);
    Rcpp::IntegerVector parent(n),left_child(n),right_child(n);
    double prior_tail_mean=1; int tail=-1;
    for (size_t i=0; i<n; ++i) {
        const Node& node=nodes[i];
        left(i)=node.left; right(i)=node.right; depth(i)=node.depth; alpha(i)=node.alpha;
        center(i)=std::exp(node.log_center); initial(i)=node.initial_count;
        parent[i]=node.parent; left_child[i]=node.left_child; right_child[i]=node.right_child;
        center_branch(i)=i==0 ? 1 : std::exp(node.log_center-nodes[node.parent].log_center);
        branch(i)=i==0 ? 1 : (nodes[node.parent].left_child==static_cast<int>(i) ? theta(node.parent) : 1-theta(node.parent));
        if (node.right==INF && node.left_child>=0) {
            const double a=nodes[node.left_child].alpha,b=nodes[node.right_child].alpha;
            prior_tail_mean*=b/(a+b);
        }
        if (node.left_child<0 && node.right==INF) tail=i;
    }
    return Rcpp::List::create(Rcpp::Named("left_bound")=left,Rcpp::Named("right_bound")=right,
        Rcpp::Named("node_levels")=depth,Rcpp::Named("alpha")=alpha,
        Rcpp::Named("n_obs_allocated")=counts,Rcpp::Named("posterior_probs")=mass,
        Rcpp::Named("terminal_nodes")=leaves,Rcpp::Named("parent_nodes")=parent,
        Rcpp::Named("left_children")=left_child,Rcpp::Named("right_children")=right_child,
        Rcpp::Named("split_nodes")=splits,Rcpp::Named("branch_probs")=branch,
        Rcpp::Named("reference_probs")=center,Rcpp::Named("reference_branch_probs")=center_branch,
        Rcpp::Named("initial_allocations")=initial,Rcpp::Named("largest_obs")=tau,
        Rcpp::Named("tail_node")=tail,Rcpp::Named("tail_prob")=mass(tail),
        Rcpp::Named("tail_prior")=s0,Rcpp::Named("prior_tail_mean")=prior_tail_mean,
        Rcpp::Named("map_fit")=stage_output(map_stage));
}

// Plain R data are used instead of an external pointer, so saved model objects
// remain usable after saveRDS/readRDS and in a new R process.
Rcpp::List PT::state(const Rcpp::List& tree_output) const {
    arma::vec log_center(nodes.size());
    for (size_t i=0; i<nodes.size(); ++i) log_center(i)=nodes[i].log_center;
    return Rcpp::List::create(Rcpp::Named("version")=1,
        Rcpp::Named("X")=data,Rcpp::Named("prior")=prior,
        Rcpp::Named("parameters")=parameters,Rcpp::Named("tau")=tau,
        Rcpp::Named("gamma")=gamma,Rcpp::Named("kappa")=kappa,Rcpp::Named("s0")=s0,
        Rcpp::Named("min_count")=min_count,Rcpp::Named("epsilon")=epsilon,
        Rcpp::Named("allocated")=allocated,Rcpp::Named("max_depth")=max_depth,
        Rcpp::Named("flat")=flat,Rcpp::Named("tree")=tree_output,
        Rcpp::Named("log_mass")=log_mass,Rcpp::Named("log_center")=log_center,
        Rcpp::Named("theta")=theta);
}

PT::PT(const Rcpp::List& saved, bool verbose) : verbose(verbose) {
    if (!saved.containsElementNamed("version") || Rcpp::as<int>(saved["version"])!=1)
        Rcpp::stop("Unsupported saved PT model; refit with the current package.");
    data=Rcpp::as<arma::mat>(saved["X"]);
    prior=Rcpp::as<int>(saved["prior"]);
    parameters=Rcpp::as<arma::vec>(saved["parameters"]);
    tau=Rcpp::as<double>(saved["tau"]);
    gamma=Rcpp::as<double>(saved["gamma"]);
    kappa=Rcpp::as<double>(saved["kappa"]);
    s0=Rcpp::as<double>(saved["s0"]);
    min_count=Rcpp::as<double>(saved["min_count"]);
    epsilon=Rcpp::as<double>(saved["epsilon"]);
    allocated=Rcpp::as<int>(saved["allocated"]);
    max_depth=Rcpp::as<int>(saved["max_depth"]);
    // Older 0.4.0 states already contain the actual Alpha values; those remain
    // sufficient for prediction even though the old warm-start flag is gone.
    flat=saved.containsElementNamed("flat") ? Rcpp::as<bool>(saved["flat"]) : false;
    max_iter=0; tolerance=0; h_tau=prior==1 ? 0 : hazard(tau);
    const Rcpp::List tree=saved["tree"];
    const arma::vec left=Rcpp::as<arma::vec>(tree["left_bound"]);
    const arma::vec right=Rcpp::as<arma::vec>(tree["right_bound"]);
    const arma::vec depth=Rcpp::as<arma::vec>(tree["node_levels"]);
    const arma::vec alpha=Rcpp::as<arma::vec>(tree["alpha"]);
    const arma::vec initial=Rcpp::as<arma::vec>(tree["initial_allocations"]);
    const arma::vec center=Rcpp::as<arma::vec>(saved["log_center"]);
    const Rcpp::IntegerVector parent=tree["parent_nodes"], lc=tree["left_children"], rc=tree["right_children"];
    mass=Rcpp::as<arma::vec>(tree["posterior_probs"]);
    log_mass=Rcpp::as<arma::vec>(saved["log_mass"]);
    theta=Rcpp::as<arma::vec>(saved["theta"]);
    counts=Rcpp::as<arma::vec>(tree["n_obs_allocated"]);
    const size_t n=left.n_elem;
    if (n<3 || right.n_elem!=n || depth.n_elem!=n || alpha.n_elem!=n ||
        initial.n_elem!=n || center.n_elem!=n || mass.n_elem!=n || log_mass.n_elem!=n ||
        theta.n_elem!=n || counts.n_elem!=n || static_cast<size_t>(parent.size())!=n ||
        static_cast<size_t>(lc.size())!=n || static_cast<size_t>(rc.size())!=n)
        Rcpp::stop("Invalid saved PT model dimensions.");
    for (size_t i=0; i<n; ++i) {
        if ((i==0 && parent[i]!=-1) || (i>0 && (parent[i]<0 || parent[i]>=static_cast<int>(i))) ||
            ((lc[i]<0)!=(rc[i]<0)) || (lc[i]>=0 &&
             (lc[i]<=static_cast<int>(i) || rc[i]<=static_cast<int>(i) ||
              lc[i]>=static_cast<int>(n) || rc[i]>=static_cast<int>(n))))
            Rcpp::stop("Invalid saved PT tree links.");
        nodes.push_back({left(i),right(i),center(i),initial(i),alpha(i),
                         static_cast<int>(depth(i)),parent[i],lc[i],rc[i]});
        if (lc[i]<0) leaves.push_back(i); else splits.push_back(i);
    }
    // No tree building, rate estimation, EM, or grid evaluation occurs here.
}

Rcpp::List PT::predict_at(const arma::vec& times,const arma::mat& leaf_draws) {
    for (double t : times)
        if (std::isnan(t) || t<0) Rcpp::stop("Prediction times must be nonnegative and nonmissing.");
    if (leaf_draws.n_rows>0 && (leaf_draws.n_cols!=leaves.size() ||
        !leaf_draws.is_finite() || arma::any(arma::vectorise(leaf_draws)<0)))
        Rcpp::stop("Invalid saved posterior leaf draws.");
    const double cells=2.0*times.n_elem*leaf_draws.n_rows;
    if (cells>512.0*1024*1024/sizeof(double))
        Rcpp::stop("Prediction draws exceed 512 MiB; predict in smaller time batches.");
    grid=times;
    prepare_predictions();
    arma::vec survival,density;
    predict(survival,density);
    Rcpp::List result=Rcpp::List::create(Rcpp::Named("survival_probs")=survival,
                                        Rcpp::Named("density_probs")=density);
    if (leaf_draws.n_rows>0) {
        arma::mat survival_draws(leaf_draws.n_rows,grid.n_elem),density_draws(leaf_draws.n_rows,grid.n_elem);
        if (predictions_cached) {
            survival_draws=leaf_draws*survival_kernel.t();
            density_draws=leaf_draws*density_kernel.t();
        } else {
            const arma::vec saved_mass=mass;
            for (arma::uword g=0; g<leaf_draws.n_rows; ++g) {
                Rcpp::checkUserInterrupt();
                for (size_t j=0; j<leaves.size(); ++j) mass(leaves[j])=leaf_draws(g,j);
                predict(survival,density);
                survival_draws.row(g)=survival.t(); density_draws.row(g)=density.t();
            }
            mass=saved_mass;
        }
        result["survival_draws"]=survival_draws;
        result["density_draws"]=density_draws;
    }
    return result;
}
