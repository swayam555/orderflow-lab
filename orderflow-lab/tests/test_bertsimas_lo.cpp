#include <qc/test.hpp>

#include <cmath>
#include <vector>

#include <qc/parallel.hpp>
#include <qc/rng.hpp>

#include "ofl/bertsimas_lo.hpp"

using namespace ofl;
using qc::test::close;

namespace {
BlParams test_params(double X0 = 0.0) {
    return BlParams{1e5, 20, 50.0, 5e-5, 5.0, 0.5, 0.125, std::sqrt(0.001), X0};
}

// Gauss-Hermite nodes/weights (physicists', weight exp(-x^2)) by Newton iteration on H_n (Numerical Recipes gauher).
void gauss_hermite(int n, std::vector<double>& x, std::vector<double>& w) {
    x.assign(static_cast<std::size_t>(n), 0.0);
    w.assign(static_cast<std::size_t>(n), 0.0);
    const double pim4 = 0.7511255444649425;  // pi^(-1/4)
    double z = 0.0, pp = 0.0;
    for (int i = 0; i < (n + 1) / 2; ++i) {
        if (i == 0)
            z = std::sqrt(2.0 * n + 1.0) - 1.85575 * std::pow(2.0 * n + 1.0, -0.16667);
        else if (i == 1)
            z -= 1.14 * std::pow(static_cast<double>(n), 0.426) / z;
        else if (i == 2)
            z = 1.86 * z - 0.86 * x[0];
        else if (i == 3)
            z = 1.91 * z - 0.91 * x[1];
        else
            z = 2.0 * z - x[static_cast<std::size_t>(i - 2)];
        for (int its = 0; its < 100; ++its) {
            double p1 = pim4, p2 = 0.0, p3;
            for (int j = 0; j < n; ++j) {
                p3 = p2;
                p2 = p1;
                p1 = z * std::sqrt(2.0 / (j + 1)) * p2 - std::sqrt(static_cast<double>(j) / (j + 1)) * p3;
            }
            pp = std::sqrt(2.0 * n) * p2;
            const double z1 = z;
            z = z1 - p1 / pp;
            if (std::fabs(z - z1) < 1e-14) break;
        }
        x[static_cast<std::size_t>(i)] = z;
        x[static_cast<std::size_t>(n - 1 - i)] = -z;
        w[static_cast<std::size_t>(i)] = 2.0 / (pp * pp);
        w[static_cast<std::size_t>(n - 1 - i)] = w[static_cast<std::size_t>(i)];
    }
}
}  // namespace

TEST_CASE("bl 1: gamma = 0 closed forms") {
    for (int T : {1, 2, 20}) {
        BlParams p = test_params();
        p.gamma = 0.0;
        p.T = T;
        const BlSolution s = bertsimas_lo_solve(p);
        for (int t = 1; t <= T; ++t) {
            const double W = 7e4 - 100.0 * t;
            CHECK(close(bl_policy(s, t, W, 0.37), W / (T - t + 1), 0.0, 1e-12));
        }
        const double v1 = p.P0 * p.S_bar + p.theta * p.S_bar * p.S_bar * (T + 1.0) / (2.0 * T);
        CHECK(close(bl_value(s, 1, p.P0, p.X0, p.S_bar), v1, 0.0, 1e-12));
    }
}

TEST_CASE("bl 2: T = 1") {
    BlParams p = test_params(0.3);
    p.T = 1;
    const BlSolution s = bertsimas_lo_solve(p);
    CHECK(close(bl_policy(s, 1, p.S_bar, p.X0), p.S_bar, 0.0, 1e-15));
    const double v = p.P0 * p.S_bar + p.theta * p.S_bar * p.S_bar + p.gamma * p.rho * p.X0 * p.S_bar;
    CHECK(close(bl_value(s, 1, p.P0, p.X0, p.S_bar), v, 0.0, 1e-12));
}

TEST_CASE("bl 3: Bellman consistency by Gauss-Hermite") {
    const BlParams p = test_params();
    const BlSolution s = bertsimas_lo_solve(p);
    std::vector<double> gx, gw;
    gauss_hermite(20, gx, gw);
    double wsum = 0;
    for (double w : gw) wsum += w;
    CHECK(close(wsum, std::sqrt(M_PI), 0.0, 1e-12));
    qc::RngStream rng(3000, 0, 1);
    for (int trial = 0; trial < 5; ++trial) {
        const int t = static_cast<int>(rng.uniform_int(1, p.T - 1));
        const double P = 40.0 + 20.0 * rng.uniform(), X = (rng.uniform() - 0.5) * 2.0, W = 1e4 + 8e4 * rng.uniform();
        auto Q = [&](double S) {
            double acc = 0.0;
            for (std::size_t i = 0; i < gx.size(); ++i) {
                const double eta = std::sqrt(2.0) * p.sigma_eta * gx[i];
                const double Xn = p.rho * X + eta;
                const double Pn = P + p.theta * S + p.gamma * Xn;
                acc += gw[i] * (Pn * S + bl_value(s, t + 1, Pn, Xn, W - S));
            }
            return acc / std::sqrt(M_PI);
        };
        const double Sstar = bl_policy(s, t, W, X);
        const double v = bl_value(s, t, P, X, W);
        CHECK_MESSAGE(close(Q(Sstar), v, 0.0, 1e-9), "t=" << t);
        CHECK(Q(Sstar + 1.0) > Q(Sstar));
        CHECK(Q(Sstar - 1.0) > Q(Sstar));
    }
}

TEST_CASE("bl 4: Monte Carlo vs analytic (X0 = 0)") {
    const BlParams p = test_params(0.0);
    const BlSolution s = bertsimas_lo_solve(p);
    const BlMcResult r = bl_monte_carlo(p, 200000, 3001, qc::default_threads());
    const double v1 = bl_value(s, 1, p.P0, p.X0, p.S_bar);
    const double tw = bl_twap_expected_cost(p);
    MESSAGE("optimal: MC " << r.opt.mean() << " +- " << r.opt.se() << " vs V1 " << v1);
    MESSAGE("TWAP: MC " << r.twap.mean() << " +- " << r.twap.se() << " vs analytic " << tw);
    MESSAGE("diff: MC " << r.diff.mean() << " +- " << r.diff.se() << " vs " << v1 - tw);
    CHECK(std::fabs(r.opt.mean() - v1) < 4.0 * r.opt.se());
    CHECK(std::fabs(r.twap.mean() - tw) < 4.0 * r.twap.se());
    CHECK(r.diff.mean() <= 0.0);
    CHECK(std::fabs(r.diff.mean() - (v1 - tw)) < 4.0 * r.diff.se());
}

TEST_CASE("bl 5: information has value (X0 = 0.1)") {
    const BlParams p = test_params(0.1);
    const BlMcResult r = bl_monte_carlo(p, 200000, 3001, qc::default_threads());
    MESSAGE("paired diff (opt - TWAP): " << r.diff.mean() << " +- " << r.diff.se());
    CHECK(r.diff.mean() < -4.0 * r.diff.se());
}
