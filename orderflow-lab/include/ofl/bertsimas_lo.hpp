// Bertsimas-Lo (1998) dynamic execution with serially correlated order-flow information.
#pragma once
#include <cstdint>
#include <vector>

#include <qc/stats.hpp>

namespace ofl {

// Buy S_bar shares over T periods. sigma_eps, sigma_eta are standard deviations of the price and information shocks.
struct BlParams {
    double S_bar;
    int T;
    double P0, theta, gamma, rho, sigma_eps, sigma_eta, X0;
};
// Value function V_t(P,X,W) = P W + a_t W^2 + b_t X W + c_t X^2 + d_t; vectors indexed 1..T (index 0 unused).
struct BlSolution {
    BlParams p;
    std::vector<double> a, b, c, d;
};
// Throws std::invalid_argument if theta <= 0 or T < 1.
[[nodiscard]] BlSolution bertsimas_lo_solve(const BlParams& p);
// Optimal shares S*_t to buy in period t given remaining W_t and previous information state X_{t-1}.
[[nodiscard]] double bl_policy(const BlSolution& s, int t, double W, double X_prev);
// Expected remaining cost V_t(P_{t-1}, X_{t-1}, W_t).
[[nodiscard]] double bl_value(const BlSolution& s, int t, double P_prev, double X_prev, double W);
// Analytic expected cost of the TWAP policy in the same model.
[[nodiscard]] double bl_twap_expected_cost(const BlParams& p);
struct BlMcResult {
    qc::RunningStats opt, twap, diff;  // cost of optimal, TWAP, and opt - twap per path (common random numbers)
};
// Monte Carlo (1024 paths per block, stream_id = block, purpose 1: per period eta then eps normal draws).
[[nodiscard]] BlMcResult bl_monte_carlo(const BlParams& p, std::int64_t paths, std::uint64_t seed, unsigned threads);

}  // namespace ofl
