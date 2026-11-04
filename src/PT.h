#ifndef SURVPT_PT_H
#define SURVPT_PT_H
#include <RcppArmadillo.h>
#include <vector>

// The final leaf has an unspecified tail shape, for every centering family.
class PT {
public:
    PT(const arma::mat& X, int prior, int max_depth, const arma::vec& parameters,
       int allocated, int max_iter, double tolerance, bool verbose,
       double gamma, double kappa, double s0, double min_count,
       double epsilon, bool flat);
    PT(const Rcpp::List& state, bool verbose);
    void fit();
    Rcpp::List sample(int n_samples, int burn_in, int thin);
    Rcpp::List output() const;
    Rcpp::List state(const Rcpp::List& tree_output) const;
    Rcpp::List predict_at(const arma::vec& times, const arma::mat& leaf_draws);
private:
    struct Node {
        double left, right, log_center, initial_count, alpha;
        int depth, parent, left_child, right_child;
    };
    struct Stage {
        int iterations = 0;
        bool converged = false;
        double difference = NA_REAL;
        std::vector<double> objective;
    };
    arma::mat data;
    arma::vec grid, parameters;
    int prior, max_depth, allocated, max_iter;
    double tolerance, gamma, kappa, s0, min_count, epsilon, tau, h_tau;
    bool verbose, flat, cached = false, predictions_cached = false, coverage_prepared = false;
    std::vector<Node> nodes;
    std::vector<int> leaves, splits, exact_leaf;
    arma::mat log_coverage, survival_kernel, density_kernel;
    arma::vec mass, log_mass, counts, theta;
    Stage map_stage;
    double hazard(double t) const;
    double log_body_probability(double a, double b) const;
    double log_reference_probability(double a, double b) const;
    double reference_count(double a, double b) const;
    double split_point(const Node& node) const;
    void add_children(int parent, double cut);
    void build_tree();
    void prepare_coverage();
    void prepare_predictions();
    double log_conditional_probability(double l, double r, double a, double b) const;
    double log_weight(int obs, int leaf) const;
    double weight(int obs, int leaf) const;
    double conditional_density(double t, double a, double b) const;
    void propagate();
    void aggregate(arma::vec& values) const;
    double allocate(arma::vec& values, bool draw_labels) const;
    void update_splits(const arma::vec& values);
    double objective(double log_likelihood) const;
    Stage run_em();
    void predict(arma::vec& survival, arma::vec& density) const;
    Rcpp::List stage_output(const Stage& stage) const;
};
#endif
