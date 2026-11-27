// Multivariate Hawkes process with exponential kernels: exact simulation, likelihood, EM calibration.
#pragma once
#include <span>
#include <vector>

#include <qc/linalg.hpp>
#include <qc/rng.hpp>

namespace ofl {

// lambda_i(t) = mu_i + sum_j sum_{t^j_k<t} alpha_ij exp(-beta_ij (t - t^j_k)); rates in 1/s, beta in 1/s.
struct HawkesParams {
    int dim = 0;
    std::vector<double> mu;
    qc::Mat alpha, beta;
    // Throws std::invalid_argument on size mismatch, mu<=0, alpha<0, beta<=0, or rho(G)>=1.
    void validate() const;
};
// Branching matrix G_ij = alpha_ij / beta_ij (dimensionless).
[[nodiscard]] qc::Mat branching_matrix(const HawkesParams& p);
// Stationary mean rates (I-G)^{-1} mu (events/s).
[[nodiscard]] std::vector<double> stationary_rates(const HawkesParams& p);
// mu = (I-G) lambda_bar, alpha = G .* beta. Throws if rho(G)>=1 or any mu_i <= 0.
[[nodiscard]] HawkesParams params_from_rates(const qc::Mat& G, const qc::Mat& beta, std::span<const double> lambda_bar);
// G = G_shape * (target / rho(G_shape)), then params_from_rates. Requires 0 <= target < 1.
[[nodiscard]] HawkesParams with_branching_ratio(const qc::Mat& G_shape, const qc::Mat& beta,
                                                std::span<const double> lambda_bar, double target);

// Recursive intensity state (time in seconds). S_ij = sum over past type-j events of alpha_ij exp(-beta_ij (t - t_k)).
class HawkesState {
  public:
    explicit HawkesState(const HawkesParams& p);
    // Decays state to time t (precondition t >= time()).
    void decay_to(double t) noexcept;
    // Registers a type-j event at the current time().
    void add_event(int j) noexcept;
    [[nodiscard]] double intensity(int i) const noexcept;
    [[nodiscard]] double total_intensity() const noexcept;
    [[nodiscard]] double time() const noexcept { return t_; }
    // Integral of lambda_i over (t, t+dt] assuming no events in between.
    [[nodiscard]] double compensator(int i, double dt) const noexcept;
    [[nodiscard]] const HawkesParams& params() const noexcept { return p_; }

  private:
    HawkesParams p_;
    std::vector<double> S_;
    double t_ = 0.0;
};

struct HawkesEvent {
    double t;
    int type;
};

// Exact Ogata thinning on [0, T]; draw order per proposal: exponential, uniform (accept), uniform (type).
[[nodiscard]] std::vector<HawkesEvent> simulate_hawkes(const HawkesParams& p, double T, qc::RngStream& rng);
// Log-likelihood of events (sorted, in (0,T]) on [0,T].
[[nodiscard]] double hawkes_loglik(const HawkesParams& p, std::span<const HawkesEvent> events, double T);
// Per-dimension compensator increments between consecutive same-type events (first from 0); Exp(1) under the model.
[[nodiscard]] std::vector<std::vector<double>> rescaled_interarrivals(const HawkesParams& p,
                                                                      std::span<const HawkesEvent> events, double T);

struct EmOptions {
    int max_iter = 500;
    double rel_tol = 1e-10;
};
struct EmResult {
    HawkesParams params;
    std::vector<double> loglik_trace;
    int iterations = 0;
    bool converged = false;
};
// One EM step with beta fixed (no precondition checks): returns updated params and the log-likelihood at `p`.
struct EmStep {
    HawkesParams params;
    double loglik;
};
[[nodiscard]] EmStep em_step(std::span<const HawkesEvent> events, double T, const HawkesParams& p);
// EM with beta known (Veen & Schoenberg 2008, recursive form). Throws if events unsorted/out of (0,T] or any
// dimension has < 10 events.
[[nodiscard]] EmResult fit_em(std::span<const HawkesEvent> events, double T, const qc::Mat& beta,
                              const HawkesParams& init, const EmOptions& opt = {});
// Default EM initial guess: mu_i = 0.5 N_i/T, alpha_ij = 0.1 beta_ij.
[[nodiscard]] HawkesParams default_em_init(std::span<const HawkesEvent> events, double T, const qc::Mat& beta);

}  // namespace ofl
