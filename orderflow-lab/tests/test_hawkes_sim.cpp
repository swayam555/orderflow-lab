#include <qc/test.hpp>

#include <cmath>
#include <string>

#include <qc/linalg.hpp>
#include <qc/rng.hpp>
#include <qc/stats.hpp>

#include "ofl/hawkes.hpp"

using namespace ofl;
using qc::test::close;

namespace {
HawkesParams model1d() {
    HawkesParams p;
    p.dim = 1;
    p.mu = {1.0};
    p.alpha = qc::Mat(1, 1);
    p.alpha(0, 0) = 1.5;
    p.beta = qc::Mat(1, 1);
    p.beta(0, 0) = 2.0;
    return p;
}
HawkesParams model2d() {
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
}  // namespace

TEST_CASE("hawkes 1: spectral radius") {
    qc::Mat d(2, 2);
    d(0, 0) = 0.3;
    d(1, 1) = 0.7;
    CHECK(close(qc::spectral_radius_nonneg(d), 0.7, 1e-9));
    qc::Mat nil(2, 2);
    nil(0, 1) = 1.0;
    CHECK(std::fabs(qc::spectral_radius_nonneg(nil)) < 1e-6);
    qc::Mat jor(2, 2);
    jor(0, 0) = 0.5;
    jor(0, 1) = 1.0;
    jor(1, 1) = 0.5;
    CHECK(close(qc::spectral_radius_nonneg(jor), 0.5, 1e-8));
    CHECK(close(qc::spectral_radius_nonneg(branching_matrix(model2d())), 0.53027756377, 1e-9));
    qc::Mat neg(2, 2);
    neg(0, 1) = -0.1;
    bool threw = false;
    try {
        (void)qc::spectral_radius_nonneg(neg);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK(threw);
    CHECK(qc::spectral_radius_nonneg(qc::Mat(3, 3)) == 0.0);
}

TEST_CASE("hawkes 2: rates, params_from_rates, branching ratio") {
    const HawkesParams p = model2d();
    const auto lam = stationary_rates(p);
    CHECK(close(lam[0], 1.012820513, 1e-8));
    CHECK(close(lam[1], 0.717948718, 1e-8));
    const HawkesParams q = params_from_rates(branching_matrix(p), p.beta, lam);
    CHECK(close(q.mu[0], 0.5, 1e-12));
    CHECK(close(q.mu[1], 0.3, 1e-12));
    CHECK(close(q.alpha(0, 1), 0.3, 1e-12));
    const auto lam2 = stationary_rates(q);
    CHECK(close(lam2[0], lam[0], 1e-12));
    CHECK(close(lam2[1], lam[1], 1e-12));
    const double target_rates[2] = {2.0, 1.5};
    const HawkesParams r = with_branching_ratio(branching_matrix(p), p.beta, target_rates, 0.8);
    CHECK(close(qc::spectral_radius_nonneg(branching_matrix(r)), 0.8, 1e-10));
    const auto lam3 = stationary_rates(r);
    CHECK(close(lam3[0], 2.0, 1e-12));
    CHECK(close(lam3[1], 1.5, 1e-12));
    HawkesParams bad = p;
    bad.alpha(0, 0) = 2.5;  // rho >= 1
    bool t1 = false, t2 = false;
    try {
        bad.validate();
    } catch (const std::invalid_argument&) {
        t1 = true;
    }
    bad = p;
    bad.mu[1] = 0.0;
    try {
        bad.validate();
    } catch (const std::invalid_argument&) {
        t2 = true;
    }
    CHECK(t1);
    CHECK(t2);
}

TEST_CASE("hawkes 3: state intensities and compensator") {
    const HawkesParams p = model2d();
    HawkesState st(p);
    st.decay_to(0.5);
    st.add_event(0);
    st.decay_to(1.0);
    st.add_event(1);
    st.decay_to(1.7);
    const double l0 = 0.5 + 0.8 * std::exp(-2.4) + 0.3 * std::exp(-1.4);
    const double l1 = 0.3 + 0.4 * std::exp(-2.4) + 0.6 * std::exp(-1.4);
    CHECK(close(st.intensity(0), l0, 0.0, 1e-14));
    CHECK(close(st.intensity(1), l1, 0.0, 1e-14));
    CHECK(close(st.total_intensity(), l0 + l1, 0.0, 1e-14));
    // compensator over (1.0, 1.7] vs Simpson
    HawkesState base(p);
    base.decay_to(0.5);
    base.add_event(0);
    base.decay_to(1.0);
    base.add_event(1);
    for (int i = 0; i < 2; ++i) {
        const int n = 10000;
        const double h = 0.7 / n;
        double s = 0.0;
        for (int k = 0; k <= n; ++k) {
            HawkesState c = base;
            c.decay_to(1.0 + k * h);
            const double w = (k == 0 || k == n) ? 1.0 : (k % 2 ? 4.0 : 2.0);
            s += w * c.intensity(i);
        }
        s *= h / 3.0;
        CHECK(close(base.compensator(i, 0.7), s, 0.0, 1e-10));
    }
}

TEST_CASE("hawkes 4: hand log-likelihood") {
    HawkesParams p;
    p.dim = 1;
    p.mu = {1.0};
    p.alpha = qc::Mat(1, 1);
    p.alpha(0, 0) = 0.5;
    p.beta = qc::Mat(1, 1);
    p.beta(0, 0) = 1.0;
    const HawkesEvent ev[3] = {{1.0, 0}, {2.0, 0}, {2.5, 0}};
    const double l1 = 1.0, l2 = 1.0 + 0.5 * std::exp(-1.0), l3 = 1.0 + 0.5 * (std::exp(-1.5) + std::exp(-0.5));
    const double comp = 4.0 + 0.5 * ((1 - std::exp(-3.0)) + (1 - std::exp(-2.0)) + (1 - std::exp(-1.5)));
    const double expect = std::log(l1) + std::log(l2) + std::log(l3) - comp;
    CHECK(close(hawkes_loglik(p, ev, 4.0), expect, 1e-12));
}

TEST_CASE("hawkes 5: determinism") {
    const HawkesParams p = model2d();
    qc::RngStream a(1001, 5, 1), b(1001, 5, 1), c(1001, 6, 1);
    const auto ea = simulate_hawkes(p, 200.0, a);
    const auto eb = simulate_hawkes(p, 200.0, b);
    const auto ec = simulate_hawkes(p, 200.0, c);
    REQUIRE(ea.size() == eb.size());
    bool same = true;
    for (std::size_t i = 0; i < ea.size(); ++i) same = same && ea[i].t == eb[i].t && ea[i].type == eb[i].type;
    CHECK(same);
    bool differs = ea.size() != ec.size();
    for (std::size_t i = 0; !differs && i < ea.size(); ++i) differs = ea[i].t != ec[i].t;
    CHECK(differs);
}

namespace {
void check_ks(const HawkesParams& p, const std::vector<HawkesEvent>& ev, double T, const std::string& tag) {
    const auto res = rescaled_interarrivals(p, ev, T);
    for (std::size_t i = 0; i < res.size(); ++i) {
        const auto n = static_cast<std::int64_t>(res[i].size());
        const double ks = qc::ks_statistic(res[i], [](double x) { return 1.0 - std::exp(-x); });
        qc::RunningStats rs;
        for (double r : res[i]) rs.add(r);
        MESSAGE(tag << " dim " << i << ": n=" << n << " KS=" << ks << " crit=" << qc::ks_critical_001(n)
                    << " mean resid=" << rs.mean() << " se=" << rs.se());
        CHECK(ks < qc::ks_critical_001(n));
        CHECK(std::fabs(rs.mean() - 1.0) < 4.0 * rs.se());
    }
}
}  // namespace

TEST_CASE("hawkes 6-9: stationary rates, time rescaling, power") {
    const double T = 1e5;
    const HawkesParams p1 = model1d();
    qc::RngStream r1(1001, 0, 1);
    const auto e1 = simulate_hawkes(p1, T, r1);
    const double rate1 = static_cast<double>(e1.size()) / T;
    MESSAGE("1-D rate " << rate1 << " vs 4");
    CHECK(std::fabs(rate1 - 4.0) < 0.03 * 4.0);

    const HawkesParams p2 = model2d();
    qc::RngStream r2(1001, 1, 1);
    const auto e2 = simulate_hawkes(p2, T, r2);
    double n0 = 0, n1 = 0;
    for (const auto& e : e2) (e.type == 0 ? n0 : n1) += 1.0;
    MESSAGE("2-D rates " << n0 / T << ", " << n1 / T << " vs 1.012820513, 0.717948718");
    CHECK(std::fabs(n0 / T - 1.012820513) < 0.03 * 1.012820513);
    CHECK(std::fabs(n1 / T - 0.717948718) < 0.03 * 0.717948718);

    check_ks(p1, e1, T, "1-D");
    check_ks(p2, e2, T, "2-D");

    HawkesParams wrong = p1;
    wrong.alpha(0, 0) = 1.2;
    const auto res = rescaled_interarrivals(wrong, e1, T);
    const auto n = static_cast<std::int64_t>(res[0].size());
    const double ks = qc::ks_statistic(res[0], [](double x) { return 1.0 - std::exp(-x); });
    MESSAGE("wrong-model KS=" << ks << " vs 3*crit=" << 3.0 * qc::ks_critical_001(n));
    CHECK(ks > 3.0 * qc::ks_critical_001(n));
}
