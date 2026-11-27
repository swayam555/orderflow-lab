#include <qc/test.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include <qc/parallel.hpp>
#include <qc/stats.hpp>

#include "ofl/agents.hpp"
#include "ofl/market_sim.hpp"

using namespace ofl;

namespace {
// Scripted agent: market buy every `period` s from t0, plus a resting bid, records fills.
class ScriptAgent final : public Agent {
  public:
    ScriptAgent(double t0, double period) : next_(t0), period_(period) {}
    [[nodiscard]] double next_wakeup() const override { return next_; }
    void on_wakeup(MarketSim& sim) override {
        const OwnerId me = owner;
        (void)sim.agent_market(me, Side::Buy, 2);
        if (const auto b = sim.best_bid()) {
            const OrderResult r = sim.agent_limit(me, Side::Buy, *b, 1);
            if (r.resting > 0) last_limit = r.id;
        }
        next_ += period_;
    }
    void on_fill(const Fill& f, double) override { fills.push_back(f); }
    OwnerId owner = 0;
    OrderId last_limit = kNoOrder;
    std::vector<Fill> fills;

  private:
    double next_, period_;
};

struct Snapshot {
    std::vector<LevelView> bids, asks;
    SimConfig dummy;
};
}  // namespace

TEST_CASE("sim 1: initial book") {
    MarketSim sim(SimConfig{}, 1003, 0);
    const auto bids = sim.book().levels(Side::Buy, 100);
    const auto asks = sim.book().levels(Side::Sell, 100);
    REQUIRE(bids.size() == 20);
    REQUIRE(asks.size() == 20);
    for (const auto& l : bids) {
        CHECK(l.qty == 9);
        CHECK(l.orders == 3);
    }
    CHECK(bids.front().price == 9999);
    CHECK(bids.back().price == 9980);
    CHECK(asks.front().price == 10001);
    CHECK(asks.back().price == 10020);
    CHECK(sim.bg_live_count() == 120);
    CHECK(sim.mid() == 10000.0);
    sim.book().check_invariants();
}

TEST_CASE("sim 2: determinism of identical runs") {
    auto run = [](std::vector<Fill>& fills, std::vector<LevelView>& b, std::vector<LevelView>& a, MarketSim::Counters& c,
                  double& tm) {
        MarketSim sim(SimConfig{}, 1003, 7);
        ScriptAgent ag(620.0, 13.0);
        ag.owner = sim.add_agent(&ag);
        sim.set_fill_observer([&](const Fill& f, double) { fills.push_back(f); });
        for (double t : {100.0, 400.0, 700.0}) {
            sim.run_until(t);
            tm = sim.time();
            CHECK(sim.time() == t);
        }
        b = sim.book().levels(Side::Buy, 1000);
        a = sim.book().levels(Side::Sell, 1000);
        c = sim.counters();
    };
    std::vector<Fill> f1, f2;
    std::vector<LevelView> b1, a1, b2, a2;
    MarketSim::Counters c1, c2;
    double t1 = 0, t2 = 0;
    run(f1, b1, a1, c1, t1);
    run(f2, b2, a2, c2, t2);
    CHECK(f1.size() > 100);
    CHECK(f1 == f2);
    CHECK(b1 == b2);
    CHECK(a1 == a2);
    CHECK(c1.hawkes_events == c2.hawkes_events);
    CHECK(c1.cancels == c2.cancels);
    CHECK(c1.bg_live_integral == c2.bg_live_integral);
    CHECK(t1 == 700.0);
}

TEST_CASE("sim 3: background orders never cross") {
    MarketSim sim(SimConfig{}, 1003, 1);
    sim.run_until(1000.0);
    CHECK(sim.counters().bg_crossing_attempts == 0);
    CHECK(sim.counters().hawkes_events[0] > 500);
    sim.book().check_invariants();
}

TEST_CASE("sim 4: agent ownership and on_fill") {
    MarketSim sim(SimConfig{}, 1003, 2);
    ScriptAgent ag(50.0, 5.0);
    ag.owner = sim.add_agent(&ag);
    CHECK(ag.owner == 1);
    std::vector<Fill> all;
    sim.set_fill_observer([&](const Fill& f, double) { all.push_back(f); });
    sim.run_until(49.0);
    // cancel of a live background order must be rejected
    const auto bids = sim.book().levels(Side::Buy, 1);
    REQUIRE(!bids.empty());
    OrderId bg_id = kNoOrder;
    for (OrderId id = 1; id < sim.book().next_id(); ++id)
        if (auto v = sim.book().find(id)) {
            bg_id = id;
            break;
        }
    REQUIRE(bg_id != kNoOrder);
    const OrderResult r = sim.agent_cancel(ag.owner, bg_id);
    CHECK(r.status == Status::Rejected);
    CHECK(r.reject == Reject::UnknownId);
    CHECK(sim.book().find(bg_id).has_value());
    sim.run_until(400.0);
    std::vector<Fill> mine;
    for (const Fill& f : all)
        if (f.taker_owner == ag.owner || f.maker_owner == ag.owner) mine.push_back(f);
    CHECK(mine.size() > 20);
    CHECK(mine == ag.fills);
    // own resting order can be cancelled by the owner
    const OrderResult l = sim.agent_limit(ag.owner, Side::Buy, 9000, 1);
    CHECK(l.status == Status::Resting);
    CHECK(sim.agent_cancel(ag.owner, l.id).status == Status::Cancelled);
    CHECK(sim.agent_cancel(ag.owner, l.id).reject == Reject::UnknownId);
}

namespace {
class OneShot final : public Agent {
  public:
    explicit OneShot(double t) : t_(t) {}
    [[nodiscard]] double next_wakeup() const override { return t_; }
    void on_wakeup(MarketSim& sim) override {
        before = sim.hawkes_intensity(2);
        (void)sim.agent_market(owner, Side::Buy, 1);
        after = sim.hawkes_intensity(2);
        t_ = std::numeric_limits<double>::infinity();
    }
    OwnerId owner = 0;
    double before = 0, after = 0;

  private:
    double t_;
};
}  // namespace

TEST_CASE("sim 5: agent market orders excite the Hawkes state") {
    for (bool excite : {true, false}) {
        SimConfig cfg;
        cfg.agent_market_orders_excite = excite;
        MarketSim sim(cfg, 1003, 3);
        OneShot a(30.0);
        a.owner = sim.add_agent(&a);
        sim.run_until(60.0);
        if (excite)
            CHECK(qc::test::close(a.after - a.before, 1.0, 1e-12));
        else
            CHECK(a.after == a.before);
    }
}

namespace {
struct PathStats {
    std::array<std::int64_t, 4> events{};
    std::int64_t cancels = 0;
    double integral = 0.0;
    double sum_bid_orders = 0, sum_ask_orders = 0, sum_spread = 0;
    std::int64_t samples = 0, spread_samples = 0, both_sides = 0;
    std::vector<double> mids;  // every 10 s on [600,3000]
};

class Sampler final : public Agent {
  public:
    explicit Sampler(PathStats& ps) : ps_(ps) {}
    [[nodiscard]] double next_wakeup() const override { return t_; }
    void on_wakeup(MarketSim& sim) override {
        if (!started_) {
            c0_ = sim.counters();
            started_ = true;
        }
        const auto& bk = sim.book();
        std::int64_t nb = 0, na = 0;
        for (const auto& l : bk.levels(Side::Buy, 100000)) nb += l.orders;
        for (const auto& l : bk.levels(Side::Sell, 100000)) na += l.orders;
        ps_.sum_bid_orders += static_cast<double>(nb);
        ps_.sum_ask_orders += static_cast<double>(na);
        ++ps_.samples;
        if (bk.best_bid() && bk.best_ask()) {
            ++ps_.both_sides;
            ps_.sum_spread += *bk.best_ask() - *bk.best_bid();
        }
        const auto k = static_cast<std::int64_t>(std::llround(t_ - 600.0));
        if (k % 10 == 0) ps_.mids.push_back(sim.mid());
        t_ += 1.0;
        last_ = sim.counters();
    }
    PathStats& ps_;
    double t_ = 600.0;
    bool started_ = false;
    MarketSim::Counters c0_, last_;
};

std::vector<PathStats> run_paths(int first, int n, std::uint64_t seed) {
    std::vector<PathStats> out(static_cast<std::size_t>(n));
    qc::parallel_blocks(n, qc::default_threads(), [&](std::int64_t p) {
        MarketSim sim(SimConfig{}, seed, static_cast<std::uint64_t>(first + p));
        PathStats& ps = out[static_cast<std::size_t>(p)];
        Sampler s(ps);
        sim.add_agent(&s);
        sim.run_until(3000.0);
        // final sampler call at t=3000 happened (wakeups are <= t_end); counters at the last sample vs first
        for (int i = 0; i < 4; ++i) ps.events[static_cast<std::size_t>(i)] = s.last_.hawkes_events[static_cast<std::size_t>(i)] - s.c0_.hawkes_events[static_cast<std::size_t>(i)];
        ps.cancels = s.last_.cancels - s.c0_.cancels;
        ps.integral = s.last_.bg_live_integral - s.c0_.bg_live_integral;
    });
    return out;
}
}  // namespace

TEST_CASE("sim 6-9: stationary statistics (40 paths x 3000 s)") {
    const auto ps = run_paths(0, 40, 1003);
    std::array<double, 4> ev{};
    double cancels = 0, integral = 0, bid = 0, ask = 0, spread = 0, samples = 0, both = 0;
    for (const auto& p : ps) {
        for (int i = 0; i < 4; ++i) ev[static_cast<std::size_t>(i)] += static_cast<double>(p.events[static_cast<std::size_t>(i)]);
        cancels += static_cast<double>(p.cancels);
        integral += p.integral;
        bid += p.sum_bid_orders;
        ask += p.sum_ask_orders;
        spread += p.sum_spread;
        samples += static_cast<double>(p.samples);
        both += static_cast<double>(p.both_sides);
    }
    const double dur = 2400.0 * 40;  // wakeups at 600..3000 inclusive: 2400 s between first and last sample
    const double lam[4] = {1.0, 1.0, 0.25, 0.25};
    for (int i = 0; i < 4; ++i) {
        const double rate = ev[static_cast<std::size_t>(i)] / dur;
        MESSAGE("type " << i << " rate " << rate << " target " << lam[i]);
        CHECK(std::fabs(rate - lam[i]) < 0.05 * lam[i]);
    }
    const double expect = 0.01 * integral;
    MESSAGE("cancels " << cancels << " expected " << expect << " z=" << (cancels - expect) / std::sqrt(expect));
    CHECK(std::fabs(cancels - expect) < 4.0 * std::sqrt(expect));
    MESSAGE("mean live per side: bid " << bid / samples << " ask " << ask / samples << "; mean spread " << spread / both
            << "; fraction both sides non-empty " << both / samples);
    CHECK(bid / samples >= 30.0);
    CHECK(bid / samples <= 150.0);
    CHECK(ask / samples >= 30.0);
    CHECK(ask / samples <= 150.0);
    CHECK(spread / both >= 1.0);
    CHECK(spread / both <= 15.0);
    CHECK(both / samples >= 0.99);

    // 9: variance of 10-s mid changes, two disjoint path sets
    auto var10 = [&](std::size_t lo, std::size_t hi) {
        qc::RunningStats rs;
        for (std::size_t p = lo; p < hi; ++p)
            for (std::size_t k = 1; k < ps[p].mids.size(); ++k) rs.add(ps[p].mids[k] - ps[p].mids[k - 1]);
        return rs.variance();
    };
    const double v1 = var10(0, 20), v2 = var10(20, 40);
    MESSAGE("var(mid change 10 s): set A " << v1 << " set B " << v2);
    CHECK(std::isfinite(v1));
    CHECK(v1 > 0.0);
    CHECK(v2 > 0.0);
    CHECK(std::fabs(v1 - v2) / std::max(v1, v2) < 0.20);
}
