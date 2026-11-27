#include "ofl/avellaneda_stoikov.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <qc/parallel.hpp>
#include <qc/rng.hpp>

namespace ofl {

AsQuote as_quotes(double s, double q, double t, double T, double gamma, double sigma, double k) {
    if (!(gamma > 0.0) || !(k > 0.0))
        throw std::invalid_argument("as_quotes: gamma = " + std::to_string(gamma) + " and k = " + std::to_string(k) + " must be > 0");
    const double r = s - q * gamma * sigma * sigma * (T - t);
    const double spread = gamma * sigma * sigma * (T - t) + (2.0 / gamma) * std::log1p(gamma / k);
    return {r - 0.5 * spread, r + 0.5 * spread};
}

double as_symmetric_spread(double gamma, double sigma, double k, double T) {
    return (2.0 / gamma) * std::log1p(gamma / k) + 0.5 * gamma * sigma * sigma * T;
}

AsOwnResult as_own_model(const AsModelParams& p, std::int64_t paths, std::uint64_t seed, unsigned threads) {
    if (paths < 1) throw std::invalid_argument("as_own_model: paths must be >= 1");
    if (!(p.dt > 0.0) || !(p.T > 0.0) || !(p.k > 0.0) || !(p.A >= 0.0))
        throw std::invalid_argument("as_own_model: need dt > 0, T > 0, k > 0, A >= 0");
    constexpr std::int64_t kBlock = 1024;
    const std::int64_t nblocks = (paths + kBlock - 1) / kBlock;
    const int steps = static_cast<int>(std::llround(p.T / p.dt));
    if (steps < 1) throw std::invalid_argument("as_own_model: T/dt rounds to zero steps");
    const double sym_spread = as_symmetric_spread(p.gamma, p.sigma, p.k, p.T);
    std::vector<AsOwnResult> parts(static_cast<std::size_t>(nblocks));
    qc::parallel_blocks(nblocks, threads, [&](std::int64_t blk) {
        qc::RngStream rng(seed, static_cast<std::uint64_t>(blk), 1);
        AsOwnResult& out = parts[static_cast<std::size_t>(blk)];
        const std::int64_t lo = blk * kBlock, hi = std::min(paths, lo + kBlock);
        const double sq = p.sigma * std::sqrt(p.dt);
        for (std::int64_t i = lo; i < hi; ++i) {
            double s = p.s0;
            double q_i = 0, cash_i = 0, q_s = 0, cash_s = 0, spread_i = 0;
            for (int n = 0; n < steps; ++n) {
                const double t = n * p.dt;
                const AsQuote qi = as_quotes(s, q_i, t, p.T, p.gamma, p.sigma, p.k);
                const double bs = s - 0.5 * sym_spread, as = s + 0.5 * sym_spread;
                const double ua = rng.uniform(), ub = rng.uniform(), us = rng.uniform();
                spread_i += qi.ask - qi.bid;
                if (ua < std::min(1.0, p.A * std::exp(-p.k * (qi.ask - s)) * p.dt)) {
                    q_i -= 1;
                    cash_i += qi.ask;
                }
                if (ub < std::min(1.0, p.A * std::exp(-p.k * (s - qi.bid)) * p.dt)) {
                    q_i += 1;
                    cash_i -= qi.bid;
                }
                if (ua < std::min(1.0, p.A * std::exp(-p.k * (as - s)) * p.dt)) {
                    q_s -= 1;
                    cash_s += as;
                }
                if (ub < std::min(1.0, p.A * std::exp(-p.k * (s - bs)) * p.dt)) {
                    q_s += 1;
                    cash_s -= bs;
                }
                s += sq * (us < 0.5 ? 1.0 : -1.0);
            }
            const double pnl_i = cash_i + q_i * s, pnl_s = cash_s + q_s * s;
            out.inventory.pnl.add(pnl_i);
            out.inventory.q_final.add(q_i);
            out.inventory.spread.add(spread_i / steps);
            out.symmetric.pnl.add(pnl_s);
            out.symmetric.q_final.add(q_s);
            out.symmetric.spread.add(sym_spread);
            out.pnl_diff.add(pnl_s - pnl_i);
        }
    });
    AsOwnResult tot;
    for (const AsOwnResult& r : parts) {
        tot.inventory.pnl.merge(r.inventory.pnl);
        tot.inventory.q_final.merge(r.inventory.q_final);
        tot.inventory.spread.merge(r.inventory.spread);
        tot.symmetric.pnl.merge(r.symmetric.pnl);
        tot.symmetric.q_final.merge(r.symmetric.q_final);
        tot.symmetric.spread.merge(r.symmetric.spread);
        tot.pnl_diff.merge(r.pnl_diff);
    }
    return tot;
}

}  // namespace ofl
