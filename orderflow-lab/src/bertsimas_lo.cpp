#include "ofl/bertsimas_lo.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include <qc/parallel.hpp>
#include <qc/rng.hpp>

namespace ofl {

BlSolution bertsimas_lo_solve(const BlParams& p) {
    if (p.T < 1) throw std::invalid_argument("BlParams: T = " + std::to_string(p.T) + " must be >= 1");
    if (!(p.theta > 0.0)) throw std::invalid_argument("BlParams: theta = " + std::to_string(p.theta) + " must be > 0");
    BlSolution s;
    s.p = p;
    const auto n = static_cast<std::size_t>(p.T) + 1;
    s.a.assign(n, 0.0);
    s.b.assign(n, 0.0);
    s.c.assign(n, 0.0);
    s.d.assign(n, 0.0);
    const auto T = static_cast<std::size_t>(p.T);
    s.a[T] = p.theta;
    s.b[T] = p.gamma * p.rho;
    const double var_eta = p.sigma_eta * p.sigma_eta;
    for (std::size_t t = T; t-- > 1;) {
        const double an = s.a[t + 1], bn = s.b[t + 1], cn = s.c[t + 1], dn = s.d[t + 1];
        s.a[t] = p.theta - p.theta * p.theta / (4.0 * an);
        s.b[t] = p.gamma * p.rho + p.theta * p.rho * bn / (2.0 * an);
        s.c[t] = p.rho * p.rho * cn - p.rho * p.rho * bn * bn / (4.0 * an);
        s.d[t] = dn + cn * var_eta;
    }
    return s;
}

double bl_policy(const BlSolution& s, int t, double W, double X_prev) {
    if (t < 1 || t > s.p.T) throw std::invalid_argument("bl_policy: t = " + std::to_string(t) + " out of range");
    if (t == s.p.T) return W;
    const double an = s.a[static_cast<std::size_t>(t) + 1], bn = s.b[static_cast<std::size_t>(t) + 1];
    return (1.0 - s.p.theta / (2.0 * an)) * W + (s.p.rho * bn / (2.0 * an)) * X_prev;
}

double bl_value(const BlSolution& s, int t, double P_prev, double X_prev, double W) {
    if (t < 1 || t > s.p.T) throw std::invalid_argument("bl_value: t = " + std::to_string(t) + " out of range");
    const auto u = static_cast<std::size_t>(t);
    return P_prev * W + s.a[u] * W * W + s.b[u] * X_prev * W + s.c[u] * X_prev * X_prev + s.d[u];
}

double bl_twap_expected_cost(const BlParams& p) {
    double sum = 0.0, geo = 0.0, rp = 1.0;
    for (int t = 1; t <= p.T; ++t) {
        rp *= p.rho;
        geo += rp;
        sum += p.P0 + p.theta * p.S_bar * t / p.T + p.gamma * p.X0 * geo;
    }
    return (p.S_bar / p.T) * sum;
}

BlMcResult bl_monte_carlo(const BlParams& p, std::int64_t paths, std::uint64_t seed, unsigned threads) {
    if (paths < 1) throw std::invalid_argument("bl_monte_carlo: paths must be >= 1");
    const BlSolution sol = bertsimas_lo_solve(p);
    constexpr std::int64_t kBlock = 1024;
    const std::int64_t nblocks = (paths + kBlock - 1) / kBlock;
    std::vector<BlMcResult> parts(static_cast<std::size_t>(nblocks));
    qc::parallel_blocks(nblocks, threads, [&](std::int64_t blk) {
        qc::RngStream rng(seed, static_cast<std::uint64_t>(blk), 1);
        BlMcResult& out = parts[static_cast<std::size_t>(blk)];
        const std::int64_t lo = blk * kBlock, hi = std::min(paths, lo + kBlock);
        for (std::int64_t i = lo; i < hi; ++i) {
            double X = p.X0, Po = p.P0, Pt = p.P0, Wo = p.S_bar, co = 0.0, ct = 0.0;
            const double s_twap = p.S_bar / p.T;
            for (int t = 1; t <= p.T; ++t) {
                const double eta = p.sigma_eta * rng.normal();
                const double eps = p.sigma_eps * rng.normal();
                const double So = bl_policy(sol, t, Wo, X);
                X = p.rho * X + eta;
                Po += p.theta * So + p.gamma * X + eps;
                Pt += p.theta * s_twap + p.gamma * X + eps;
                co += Po * So;
                ct += Pt * s_twap;
                Wo -= So;
            }
            out.opt.add(co);
            out.twap.add(ct);
            out.diff.add(co - ct);
        }
    });
    BlMcResult tot;
    for (const BlMcResult& r : parts) {
        tot.opt.merge(r.opt);
        tot.twap.merge(r.twap);
        tot.diff.merge(r.diff);
    }
    return tot;
}

}  // namespace ofl
