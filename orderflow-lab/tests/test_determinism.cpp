#include <qc/test.hpp>

#include <cmath>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include <qc/stats.hpp>

#include "ofl/avellaneda_stoikov.hpp"
#include "ofl/bertsimas_lo.hpp"
#include "ofl/experiments.hpp"
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

TEST_CASE("det 2b: Avellaneda-Stoikov own-model MC is thread-count independent") {
    const AsModelParams p;
    const auto r1 = as_own_model(p, 5000, 8001, 1);
    const auto r8 = as_own_model(p, 5000, 8001, 8);
    CHECK(same(r1.inventory.pnl, r8.inventory.pnl));
    CHECK(same(r1.inventory.q_final, r8.inventory.q_final));
    CHECK(same(r1.symmetric.pnl, r8.symmetric.pnl));
    CHECK(same(r1.pnl_diff, r8.pnl_diff));
}

namespace {
std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}
}  // namespace

TEST_CASE("det 1: impact experiment CSV bytes are thread-count independent") {
    ImpactConfig ic;
    ic.R = 16;
    ic.Qs = {8};
    std::string ref[3];
    int idx = 0;
    for (unsigned th : {1u, 3u, 8u}) {
        const ImpactResult r = run_impact_experiment(SimConfig{}, ic, 9001, th);
        const std::string dir = "ofl_test_det_" + std::to_string(th);
        write_impact_csvs(r, dir);
        ref[idx++] = slurp(dir + "/impact.csv") + slurp(dir + "/impact_profile.csv") + slurp(dir + "/impact_fit.csv");
    }
    CHECK(ref[0].size() > 200);
    CHECK(ref[0] == ref[1]);
    CHECK(ref[0] == ref[2]);
}
