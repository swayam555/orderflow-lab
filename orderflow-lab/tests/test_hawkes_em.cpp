#include <qc/test.hpp>

#include <cmath>
#include <vector>

#include <qc/linalg.hpp>
#include <qc/rng.hpp>

#include "ofl/hawkes.hpp"

using namespace ofl;
using qc::test::close;

namespace {
[[maybe_unused]] HawkesParams model2d() {
    HawkesParams p;
    p.dim = 2;
    p.mu = {0.5, 0.3};
    p.alpha = qc::Mat(2, 2);
    p.alpha(0, 0) = 0.8;
    p.alpha(0, 1) = 0.3;
    p.alpha(1, 0) = 0.4;
    p.alpha(1, 1) = 0.6;
    p.beta = qc::Mat(2, 2);
    for (double& b : p.beta.a) b = 2.0;
    return p;
}
// Default 4-D market-simulator flow (DESIGN 5.1): G rows = affected, cols = source.
[[maybe_unused]] HawkesParams model4d() {
    qc::Mat G(4, 4), beta(4, 4);
    G(0, 0) = 0.3;
    G(0, 3) = 0.3;
    G(1, 1) = 0.3;
    G(1, 2) = 0.3;
    G(2, 2) = 0.5;
    G(3, 3) = 0.5;
    for (double& b : beta.a) b = 2.0;
    const double lam[4] = {1.0, 1.0, 0.25, 0.25};
    return params_from_rates(G, beta, lam);
}
}  // namespace

#ifndef OFL_SLOW_EM
TEST_CASE("em 1: one iteration equals O(K^2) hand computation (2-D, 8 events)") {
    const HawkesParams p = model2d();
    const double T = 5.0;
    const std::vector<HawkesEvent> ev = {{0.3, 0}, {0.9, 1}, {1.1, 1}, {1.8, 0}, {2.5, 0}, {3.0, 1}, {3.2, 0}, {4.4, 1}};
    double num_mu[2] = {0, 0}, num_a[2][2] = {{0, 0}, {0, 0}}, C[2][2] = {{0, 0}, {0, 0}};
    double ll = 0.0;
    for (std::size_t k = 0; k < ev.size(); ++k) {
        const int i = ev[k].type;
        double lam = p.mu[static_cast<std::size_t>(i)];
        double contrib[2] = {0, 0};
        for (std::size_t l = 0; l < k; ++l) {
            const int j = ev[l].type;
            const double v = p.alpha(i, j) * std::exp(-p.beta(i, j) * (ev[k].t - ev[l].t));
            lam += v;
            contrib[j] += v;
        }
        ll += std::log(lam);
        num_mu[i] += p.mu[static_cast<std::size_t>(i)] / lam;
        for (int j = 0; j < 2; ++j) num_a[i][j] += contrib[j] / lam;
        for (int r = 0; r < 2; ++r) C[r][i] += (1.0 - std::exp(-p.beta(r, i) * (T - ev[k].t))) / p.beta(r, i);
    }
    double comp = 0.0;
    for (int i = 0; i < 2; ++i) {
        comp += p.mu[static_cast<std::size_t>(i)] * T;
        for (int j = 0; j < 2; ++j) comp += p.alpha(i, j) * C[i][j];
    }
    const EmStep s = em_step(ev, T, p);
    CHECK(close(s.loglik, ll - comp, 1e-12, 1e-12));
    for (int i = 0; i < 2; ++i) {
        CHECK(close(s.params.mu[static_cast<std::size_t>(i)], num_mu[i] / T, 1e-12, 1e-12));
        for (int j = 0; j < 2; ++j) CHECK(close(s.params.alpha(i, j), num_a[i][j] / C[i][j], 1e-12, 1e-12));
    }
    // 1-D 3-event example from the likelihood hand check
    HawkesParams q;
    q.dim = 1;
    q.mu = {1.0};
    q.alpha = qc::Mat(1, 1);
    q.alpha(0, 0) = 0.5;
    q.beta = qc::Mat(1, 1);
    q.beta(0, 0) = 1.0;
    const HawkesEvent e3[3] = {{1.0, 0}, {2.0, 0}, {2.5, 0}};
    const double l2 = 1.0 + 0.5 * std::exp(-1.0);
    const double A3 = std::exp(-1.5) + std::exp(-0.5), l3 = 1.0 + 0.5 * A3;
    const double nm = 1.0 + 1.0 / l2 + 1.0 / l3;
    const double na = 0.5 * std::exp(-1.0) / l2 + 0.5 * A3 / l3;
    const double c = (1 - std::exp(-3.0)) + (1 - std::exp(-2.0)) + (1 - std::exp(-1.5));
    const EmStep s1 = em_step(e3, 4.0, q);
    CHECK(close(s1.params.mu[0], nm / 4.0, 1e-12));
    CHECK(close(s1.params.alpha(0, 0), na / c, 1e-12));
}

TEST_CASE("em 2: preconditions throw") {
    const HawkesParams p = model2d();
    std::vector<HawkesEvent> few = {{0.3, 0}, {0.9, 1}};
    bool t1 = false, t2 = false;
    try {
        (void)fit_em(few, 5.0, p.beta, p);
    } catch (const std::invalid_argument&) {
        t1 = true;
    }
    std::vector<HawkesEvent> uns;
    for (int i = 0; i < 40; ++i) uns.push_back({1.0 + 0.1 * (i % 2 ? 40 - i : i), i % 2});
    try {
        (void)fit_em(uns, 100.0, p.beta, p);
    } catch (const std::invalid_argument&) {
        t2 = true;
    }
    CHECK(t1);
    CHECK(t2);
}

TEST_CASE("em 3-6: recovery on 2-D model (T=1e5)") {
    const HawkesParams truth = model2d();
    const double T = 1e5;
    qc::RngStream rng(1002, 0, 1);
    const auto ev = simulate_hawkes(truth, T, rng);
    const EmResult fit = fit_em(ev, T, truth.beta, default_em_init(ev, T, truth.beta));
    MESSAGE("iterations " << fit.iterations << " converged " << fit.converged);
    bool mono = true;
    for (std::size_t k = 1; k < fit.loglik_trace.size(); ++k)
        mono = mono && fit.loglik_trace[k] >= fit.loglik_trace[k - 1] - 1e-9 * std::fabs(fit.loglik_trace[k]);
    CHECK(mono);
    CHECK(fit.converged);
    CHECK(fit.iterations <= 500);
    for (int i = 0; i < 2; ++i) {
        const double m = fit.params.mu[static_cast<std::size_t>(i)], mt = truth.mu[static_cast<std::size_t>(i)];
        MESSAGE("mu" << i << " fit " << m << " true " << mt << " rel " << (m - mt) / mt);
        CHECK(std::fabs(m - mt) / mt < 0.10);
        for (int j = 0; j < 2; ++j) {
            const double a = fit.params.alpha(i, j), at = truth.alpha(i, j);
            MESSAGE("alpha" << i << j << " fit " << a << " true " << at << " rel " << (a - at) / at);
            CHECK(std::fabs(a - at) / at < 0.10);
        }
    }
    const double rho = qc::spectral_radius_nonneg(branching_matrix(fit.params));
    MESSAGE("rho fit " << rho << " vs 0.53028");
    CHECK(std::fabs(rho - 0.53027756377) < 0.03);
    CHECK(hawkes_loglik(fit.params, ev, T) >= hawkes_loglik(truth, ev, T) - 1e-6);
}

#else
TEST_CASE("em 7 (slow): recovery on 4-D market-simulator flow (T=2e5)") {
    const HawkesParams truth = model4d();
    const double T = 2e5;
    qc::RngStream rng(1002, 1, 1);
    const auto ev = simulate_hawkes(truth, T, rng);
    const EmResult fit = fit_em(ev, T, truth.beta, default_em_init(ev, T, truth.beta));
    MESSAGE("events " << ev.size() << " iterations " << fit.iterations << " converged " << fit.converged);
    for (int i = 0; i < 4; ++i) {
        const double m = fit.params.mu[static_cast<std::size_t>(i)], mt = truth.mu[static_cast<std::size_t>(i)];
        MESSAGE("mu" << i << " fit " << m << " true " << mt);
        CHECK(std::fabs(m - mt) / mt < 0.10);
        for (int j = 0; j < 4; ++j) {
            const double a = fit.params.alpha(i, j), at = truth.alpha(i, j);
            if (at > 0.0) {
                MESSAGE("alpha" << i << j << " fit " << a << " true " << at);
                CHECK(std::fabs(a - at) / at < 0.15);
            } else {
                MESSAGE("alpha" << i << j << " fit " << a << " true 0");
                CHECK(a < 0.05);
            }
        }
    }
}
#endif
