#include <qc/test.hpp>

#include <cmath>
#include <vector>

#include <qc/stats.hpp>

#include "ofl/bertsimas_lo.hpp"
#include "ofl/market_sim.hpp"

using namespace ofl;

namespace {
bool same(const qc::RunningStats& a, const qc::RunningStats& b) {
    return a.n() == b.n() && a.mean() == b.mean() && a.variance() == b.variance();
}
std::vector<Fill> fills_of(std::uint64_t path) {
    MarketSim sim(SimConfig{}, 1003, path);
    std::vector<Fill> v;
    sim.set_fill_observer([&](const Fill& f, double) { v.push_back(f); });
    sim.run_until(900.0);
    return v;
}
}  // namespace

TEST_CASE("det 3: MarketSim fill streams") {
    const auto a = fills_of(5), b = fills_of(5), c = fills_of(6);
    CHECK(a.size() > 500);
    CHECK(a == b);
    CHECK(!(a == c));
}

TEST_CASE("det 2a: Bertsimas-Lo Monte Carlo is thread-count independent") {
    const BlParams p{1e5, 20, 50.0, 5e-5, 5.0, 0.5, 0.125, std::sqrt(0.001), 0.1};
    const auto r1 = bl_monte_carlo(p, 20000, 3001, 1);
    const auto r3 = bl_monte_carlo(p, 20000, 3001, 3);
    const auto r8 = bl_monte_carlo(p, 20000, 3001, 8);
    CHECK(same(r1.opt, r8.opt));
    CHECK(same(r1.twap, r8.twap));
    CHECK(same(r1.diff, r8.diff));
    CHECK(same(r1.opt, r3.opt));
}
