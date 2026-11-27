#include <qc/test.hpp>

#include <cmath>
#include <vector>

#include <qc/parallel.hpp>
#include <qc/rng.hpp>
#include <qc/stats.hpp>

#include "ofl/almgren_chriss.hpp"

using namespace ofl;
using qc::test::close;

namespace {
AcParams base(double lambda, int N = 5) { return AcParams{1e6, 5.0, N, 0.95, 2.5e-7, 2.5e-6, 0.0625, lambda}; }
const double kLambdas[4] = {0.0, 1e-6, 2e-6, 1e-5};
}  // namespace

TEST_CASE("ac 1: lambda = 0 is TWAP") {
    const AcSchedule s = almgren_chriss(base(0.0));
    for (int k = 0; k <= 5; ++k) CHECK(std::fabs(s.x[static_cast<std::size_t>(k)] - 1e6 * (1.0 - k / 5.0)) < 1e-9 * 1e6);
}

TEST_CASE("ac 2: closed form vs tridiagonal reference") {
    for (int N : {1, 2, 10, 250})
        for (double lam : kLambdas) {
            const AcSchedule a = almgren_chriss(base(lam, N)), b = almgren_chriss_tridiagonal(base(lam, N));
            double worst = 0;
            for (std::size_t k = 0; k < a.x.size(); ++k) worst = std::max(worst, std::fabs(a.x[k] - b.x[k]));
            CHECK_MESSAGE(worst <= 1e-10 * 1e6, "N=" << N << " lambda=" << lam << " worst=" << worst);
        }
}

TEST_CASE("ac 3: trades sum to X, monotone holdings") {
    for (double lam : kLambdas) {
        const AcSchedule s = almgren_chriss(base(lam));
        double sum = 0;
        for (double n : s.n) sum += n;
        CHECK(std::fabs(sum - 1e6) < 1e-9 * 1e6);
        if (lam > 0)
            for (std::size_t k = 1; k < s.x.size(); ++k) CHECK(s.x[k] < s.x[k - 1]);
    }
}

TEST_CASE("ac 4: E/Var equal direct evaluation and Monte Carlo (stat)") {
    for (double lam : kLambdas) {
        const AcParams p = base(lam);
        const AcSchedule s = almgren_chriss(p);
        const double tau = p.T / p.N;
        // direct: trade k executes at S_{k-1} - h(n_k/tau); S_k = S_{k-1} - gamma n_k + sigma sqrt(tau) xi
        double e = 0.0;
        for (int k = 1; k <= p.N; ++k) {
            const double nk = s.n[static_cast<std::size_t>(k - 1)];
            e += nk * (p.gamma * (p.X - s.x[static_cast<std::size_t>(k - 1)]) + p.epsilon + p.eta * nk / tau);
        }
        double v = 0.0;
        for (int j = 1; j <= p.N; ++j) v += p.sigma * p.sigma * tau * s.x[static_cast<std::size_t>(j)] * s.x[static_cast<std::size_t>(j)];
        CHECK(close(s.expected_cost, e, 0.0, 1e-10));
        CHECK(close(s.variance, v, 0.0, 1e-10));
        // Monte Carlo on price paths
        const std::int64_t n_paths = 200000;
        qc::RngStream rng(2001, static_cast<std::uint64_t>(lam * 1e7), 1);
        qc::RunningStats cost;
        const double S0 = 50.0;
        for (std::int64_t i = 0; i < n_paths; ++i) {
            double S = S0, proceeds = 0.0;
            for (int k = 1; k <= p.N; ++k) {
                const double nk = s.n[static_cast<std::size_t>(k - 1)];
                proceeds += nk * (S - p.epsilon - p.eta * nk / tau);
                S += -p.gamma * nk + p.sigma * std::sqrt(tau) * rng.normal();
            }
            cost.add(p.X * S0 - proceeds);
        }
        const double se_mean = cost.se();
        const double se_var = s.variance * std::sqrt(2.0 / static_cast<double>(n_paths - 1));
        MESSAGE("lambda " << lam << ": E model " << s.expected_cost << " MC " << cost.mean() << " se " << se_mean
                          << "; Var model " << s.variance << " MC " << cost.variance() << " se " << se_var);
        CHECK(std::fabs(cost.mean() - s.expected_cost) < 4.0 * se_mean);
        CHECK(std::fabs(cost.variance() - s.variance) < 4.0 * se_var);
    }
}

TEST_CASE("ac 5: efficient frontier ordering") {
    double pe = -1, pv = 1e300;
    for (double lam : kLambdas) {
        const AcSchedule s = almgren_chriss(base(lam));
        if (lam > 0) {
            CHECK(s.variance < pv);
            CHECK(s.expected_cost > pe);
        }
        pe = s.expected_cost;
        pv = s.variance;
    }
}

TEST_CASE("ac 6: optimality under perturbation") {
    const AcParams p = base(2e-6);
    const AcSchedule s = almgren_chriss(p);
    auto obj = [&](const std::vector<double>& x) {
        double e, v;
        ac_evaluate(p, x, e, v);
        return e + p.lambda * v;
    };
    const double f0 = obj(s.x);
    for (int j = 1; j < p.N; ++j)
        for (double sgn : {-1.0, 1.0}) {
            std::vector<double> x = s.x;
            x[static_cast<std::size_t>(j)] += sgn * 1e-3 * p.X;
            CHECK(obj(x) > f0);
        }
}

TEST_CASE("ac 7-8: large kappa*T and lambda inversion") {
    const AcParams b = base(0.0);
    for (double kT : {50.0, 700.0}) {
        AcParams p = b;
        p.lambda = ac_lambda_for_kappaT(b, kT);
        const AcSchedule s = almgren_chriss(p);
        for (std::size_t k = 0; k < s.x.size(); ++k) {
            CHECK(std::isfinite(s.x[k]));
            CHECK(s.x[k] >= 0.0);
            CHECK(s.x[k] <= p.X);
            if (k > 0) CHECK(s.x[k] < s.x[k - 1]);
        }
        CHECK(close(s.kappa * p.T, kT, 0.0, 1e-12));
        if (kT == 50.0) {
            const long double kappa = 10.0L, T = 5.0L, tau = 1.0L;
            const long double ref = sinhl(kappa * (T - tau)) / sinhl(kappa * T);
            CHECK(std::fabs(static_cast<double>(s.x[1] / p.X) / static_cast<double>(ref) - 1.0) < 1e-12);
        }
    }
    for (double kT : {0.5, 1.0, 3.0, 50.0, 700.0}) {
        AcParams p = b;
        p.lambda = ac_lambda_for_kappaT(b, kT);
        CHECK(close(almgren_chriss(p).kappa * p.T, kT, 0.0, 1e-12));
    }
}

TEST_CASE("ac 9: eta_tilde <= 0 throws") {
    AcParams p = base(1e-6);
    p.eta = 0.5 * p.gamma * (p.T / p.N);  // eta_tilde == 0
    bool t1 = false, t2 = false;
    try {
        (void)almgren_chriss(p);
    } catch (const std::invalid_argument&) {
        t1 = true;
    }
    p.eta *= 0.5;
    try {
        (void)almgren_chriss(p);
    } catch (const std::invalid_argument&) {
        t2 = true;
    }
    CHECK(t1);
    CHECK(t2);
}

TEST_CASE("ac 10: round_schedule") {
    std::vector<double> n(10, 6.4);
    const auto r = round_schedule(n, 64);
    const std::vector<Qty> expect = {6, 7, 6, 7, 6, 6, 7, 6, 7, 6};
    CHECK(r == expect);
    Qty sum = 0;
    for (Qty q : r) sum += q;
    CHECK(sum == 64);
    const AcSchedule s = almgren_chriss(AcParams{64, 300, 10, 0.5, 0.0, 1e-3, 0.0, 1e-3});
    const auto rr = round_schedule(s.n, 64);
    Qty sr = 0;
    for (std::size_t k = 0; k < rr.size(); ++k) {
        sr += rr[k];
        CHECK(std::fabs(static_cast<double>(rr[k]) - s.n[k]) <= 1.0);
    }
    CHECK(sr == 64);
}
