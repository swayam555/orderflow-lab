// Avellaneda-Stoikov (2008) market-making formulas and the paper's own-model simulation.
#pragma once
#include <cstdint>

#include <qc/stats.hpp>

namespace ofl {

struct AsQuote {
    double bid, ask;  // price units
};
// Inventory strategy: reservation price r = s - q*gamma*sigma^2*(T-t), total spread
// gamma*sigma^2*(T-t) + (2/gamma) ln(1+gamma/k); quotes r -/+ spread/2. t, T in time units, sigma per sqrt(time).
[[nodiscard]] AsQuote as_quotes(double s, double q, double t, double T, double gamma, double sigma, double k);
// Time average over [0,T] of the AS spread: (2/gamma) ln(1+gamma/k) + gamma sigma^2 T / 2.
[[nodiscard]] double as_symmetric_spread(double gamma, double sigma, double k, double T);

// Own-model parameters (paper defaults).
struct AsModelParams {
    double s0 = 100.0, T = 1.0, sigma = 2.0, dt = 0.005, k = 1.5, A = 140.0, gamma = 0.1;
};
struct AsStrategyStats {
    qc::RunningStats pnl, q_final, spread;  // spread: per-path average quoted spread
};
struct AsOwnResult {
    AsStrategyStats inventory, symmetric;
    qc::RunningStats pnl_diff;  // symmetric - inventory, same path
};
// Simulates `paths` paths in blocks of 1024 (stream_id = block, purpose 1, draws u_a, u_b, u_s per step), both
// strategies on identical draws; block results are merged in block order (bitwise independent of `threads`).
[[nodiscard]] AsOwnResult as_own_model(const AsModelParams& p, std::int64_t paths, std::uint64_t seed, unsigned threads);

}  // namespace ofl
