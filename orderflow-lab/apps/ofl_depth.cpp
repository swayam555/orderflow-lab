// Mean order-book depth profile of the background market (distance in ticks from each side's best quote).
#include <cstdio>
#include <vector>

#include <qc/cli.hpp>
#include <qc/csv.hpp>
#include <qc/parallel.hpp>
#include <qc/stats.hpp>

#include "ofl/agents.hpp"
#include "ofl/market_sim.hpp"

#ifndef OFL_RESULTS_DIR
#define OFL_RESULTS_DIR "results"
#endif

using namespace ofl;

namespace {
constexpr int kDist = 40;
struct Acc {
    std::vector<double> bid, ask;
    double spread_sum = 0;
    std::int64_t spread_n = 0;
    std::int64_t samples = 0;
    double orders_sum = 0;
    Acc() : bid(kDist, 0.0), ask(kDist, 0.0) {}
};

class Sampler final : public Agent {
  public:
    explicit Sampler(Acc& a) : a_(a) {}
    [[nodiscard]] double next_wakeup() const override { return t_; }
    void on_wakeup(MarketSim& sim) override {
        const auto& b = sim.book();
        const auto bb = b.best_bid(), ba = b.best_ask();
        for (int d = 0; d < kDist; ++d) {
            if (bb) a_.bid[static_cast<std::size_t>(d)] += static_cast<double>(b.depth_at(Side::Buy, *bb - d));
            if (ba) a_.ask[static_cast<std::size_t>(d)] += static_cast<double>(b.depth_at(Side::Sell, *ba + d));
        }
        if (bb && ba) {
            a_.spread_sum += *ba - *bb;
            ++a_.spread_n;
        }
        a_.orders_sum += static_cast<double>(sim.bg_live_count());
        ++a_.samples;
        t_ += 10.0;
    }
    Acc& a_;
    double t_ = 600.0;
};
}  // namespace

int main(int argc, char** argv) {
    try {
        qc::Args args(argc, argv, {"paths", "seed", "out"});
        const auto paths = args.get_int("paths", 200);
        const auto seed = args.get_u64("seed", 9100);
        const std::string out = args.get_str("out", OFL_RESULTS_DIR);
        std::vector<Acc> acc(static_cast<std::size_t>(paths));
        qc::parallel_blocks(paths, qc::default_threads(), [&](std::int64_t p) {
            MarketSim sim(SimConfig{}, seed, static_cast<std::uint64_t>(p));
            Sampler s(acc[static_cast<std::size_t>(p)]);
            sim.add_agent(&s);
            sim.run_until(1200.0);
        });
        Acc tot;
        for (const Acc& a : acc) {
            for (int d = 0; d < kDist; ++d) {
                tot.bid[static_cast<std::size_t>(d)] += a.bid[static_cast<std::size_t>(d)];
                tot.ask[static_cast<std::size_t>(d)] += a.ask[static_cast<std::size_t>(d)];
            }
            tot.spread_sum += a.spread_sum;
            tot.spread_n += a.spread_n;
            tot.samples += a.samples;
            tot.orders_sum += a.orders_sum;
        }
        qc::CsvWriter csv(out + "/depth_profile.csv", {"distance_ticks", "mean_bid_depth", "mean_ask_depth"});
        double cum = 0;
        for (int d = 0; d < kDist; ++d) {
            const double bd = tot.bid[static_cast<std::size_t>(d)] / static_cast<double>(tot.samples);
            const double ad = tot.ask[static_cast<std::size_t>(d)] / static_cast<double>(tot.samples);
            csv.row({std::int64_t{d}, bd, ad});
            cum += 0.5 * (bd + ad);
        }
        std::printf("paths %lld, samples %lld\n", static_cast<long long>(paths), static_cast<long long>(tot.samples));
        std::printf("mean spread %.3f ticks, mean live background orders (both sides) %.1f\n",
                    tot.spread_sum / static_cast<double>(tot.spread_n), tot.orders_sum / static_cast<double>(tot.samples));
        std::printf("depth at touch: bid %.2f ask %.2f lots; sum over first %d levels (avg side) %.1f lots\n",
                    tot.bid[0] / static_cast<double>(tot.samples), tot.ask[0] / static_cast<double>(tot.samples), kDist, cum);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
