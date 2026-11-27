// Engine throughput / latency benchmark (DESIGN 9). Single thread, in-memory, pre-generated op stream.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <qc/cli.hpp>
#include <qc/csv.hpp>
#include <qc/rng.hpp>
#include <qc/timer.hpp>

#include "ofl/order_book.hpp"
#include "ofl/reference_book.hpp"

#ifndef OFL_RESULTS_DIR
#define OFL_RESULTS_DIR "results"
#endif

using namespace ofl;

namespace {
enum class OpType : std::uint8_t { Limit = 0, Cancel = 1, Market = 2, Modify = 3 };
const char* const kOpNames[] = {"limit", "cancel", "market", "modify"};

struct Op {
    OrderId id;  // cancel / modify target
    Tick price;
    std::int32_t qty;
    OpType type;
    Side side;
};

constexpr Tick kMinTick = 0, kMaxTick = 65535, kStartMid = 32768;

struct Stream {
    std::vector<Op> ops;          // prefill + warmup + measured
    std::size_t n_prefill = 0;
    std::size_t max_live = 0;
    double mean_levels = 0.0;     // mean (bid levels + ask levels) sampled during generation
    double fills_per_op = 0.0;
};

OrderResult apply(OrderBook& b, const Op& op) {
    switch (op.type) {
        case OpType::Limit: return b.submit_limit(op.side, op.price, op.qty, Tif::GTC, 0);
        case OpType::Cancel: return b.cancel(op.id);
        case OpType::Market: return b.submit_market(op.side, op.qty, 0);
        case OpType::Modify: return b.modify(op.id, op.price, op.qty);
    }
    return {};
}

Stream generate(std::uint64_t seed, std::size_t total_ops) {
    qc::RngStream rng(seed, 0, 1);
    OrderBook book(kMinTick, kMaxTick, 1u << 20, total_ops + 4000);
    Stream s;
    s.ops.reserve(total_ops + 2000);
    std::vector<OrderId> live;
    std::uint64_t fills = 0;
    double level_sum = 0;
    std::size_t level_samples = 0;
    auto clamp = [](std::int64_t p) { return static_cast<Tick>(std::clamp<std::int64_t>(p, kMinTick, kMaxTick)); };
    auto run = [&](const Op& op) {
        const OrderId nid = book.next_id();
        OrderResult r = apply(book, op);
        s.ops.push_back(op);
        fills += book.last_fills().size();
        if ((op.type == OpType::Limit) && (r.status == Status::Resting || r.status == Status::PartialResting))
            live.push_back(nid);
        s.max_live = std::max(s.max_live, book.live_orders());
    };
    for (int side = 0; side < 2; ++side)
        for (int i = 0; i < 1000; ++i) {
            const auto off = static_cast<Tick>(rng.uniform_int(1, 20));
            const auto q = static_cast<std::int32_t>(rng.uniform_int(1, 10));
            run(Op{0, side == 0 ? kStartMid - off : kStartMid + off, q, OpType::Limit,
                   side == 0 ? Side::Buy : Side::Sell});
        }
    s.n_prefill = s.ops.size();
    fills = 0;
    auto pick_live = [&]() -> OrderId {
        while (!live.empty()) {
            const auto i = static_cast<std::size_t>(rng.uniform_int(0, static_cast<std::int64_t>(live.size()) - 1));
            const OrderId id = live[i];
            if (book.find(id)) return id;
            live[i] = live.back();
            live.pop_back();
        }
        return kNoOrder;
    };
    auto passive_limit = [&]() {
        const Side side = rng.uniform() < 0.5 ? Side::Buy : Side::Sell;
        const auto u = rng.uniform_int(1, 20);
        const auto q = static_cast<std::int32_t>(rng.uniform_int(1, 10));
        Tick px;
        if (side == Side::Buy) px = clamp((book.best_ask() ? *book.best_ask() : kStartMid + 1) - u);
        else px = clamp((book.best_bid() ? *book.best_bid() : kStartMid - 1) + u);
        return Op{0, px, q, OpType::Limit, side};
    };
    for (std::size_t i = 0; i < total_ops; ++i) {
        const double u = rng.uniform();
        Op op{};
        if (u < 0.55) {
            if (rng.uniform() < 0.9) {
                op = passive_limit();
            } else {
                const Side side = rng.uniform() < 0.5 ? Side::Buy : Side::Sell;
                const auto c = rng.uniform_int(0, 2);
                const auto q = static_cast<std::int32_t>(rng.uniform_int(1, 10));
                const Tick px = side == Side::Buy ? clamp((book.best_ask() ? *book.best_ask() : kStartMid + 1) + c)
                                                  : clamp((book.best_bid() ? *book.best_bid() : kStartMid - 1) - c);
                op = Op{0, px, q, OpType::Limit, side};
            }
        } else if (u < 0.85) {
            const OrderId id = pick_live();
            op = id == kNoOrder ? passive_limit() : Op{id, 0, 0, OpType::Cancel, Side::Buy};
        } else if (u < 0.95) {
            const Side side = rng.uniform() < 0.5 ? Side::Buy : Side::Sell;
            op = Op{0, 0, static_cast<std::int32_t>(rng.uniform_int(1, 4)), OpType::Market, side};
        } else {
            const OrderId id = pick_live();
            if (id == kNoOrder) {
                op = passive_limit();
            } else {
                const auto ov = *book.find(id);
                if (rng.uniform() < 0.5) {
                    const Qty nq = ov.remaining > 1 ? rng.uniform_int(1, ov.remaining) : 1;
                    op = Op{id, ov.price, static_cast<std::int32_t>(nq), OpType::Modify, ov.side};
                } else {
                    const auto d = rng.uniform_int(1, 3) * (rng.uniform() < 0.5 ? -1 : 1);
                    op = Op{id, clamp(ov.price + d), static_cast<std::int32_t>(ov.remaining), OpType::Modify, ov.side};
                }
            }
        }
        run(op);
        if (i % 10000 == 0) {
            level_sum += static_cast<double>(book.levels(Side::Buy, 100000).size() + book.levels(Side::Sell, 100000).size());
            ++level_samples;
        }
    }
    s.fills_per_op = static_cast<double>(fills) / static_cast<double>(total_ops);
    s.mean_levels = level_sum / static_cast<double>(std::max<std::size_t>(level_samples, 1));
    return s;
}

std::string cpu_model() {
    std::ifstream f("/proc/cpuinfo");
    std::string line;
    while (std::getline(f, line))
        if (line.rfind("model name", 0) == 0) return line.substr(line.find(':') + 2);
    return "unknown";
}

double pct(std::vector<std::uint32_t>& v, double p) {
    if (v.empty()) return 0;
    const auto k = static_cast<std::size_t>(p * static_cast<double>(v.size() - 1));
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k), v.end());
    return v[k];
}

bool verify(const Stream& s, std::size_t n_ops) {
    OrderBook a(kMinTick, kMaxTick, 1u << 20, s.ops.size() + 10);
    ReferenceBook r(kMinTick, kMaxTick);
    const std::size_t n = std::min(n_ops, s.ops.size());
    for (std::size_t i = 0; i < n; ++i) {
        const Op& op = s.ops[i];
        OrderResult ra = apply(a, op), rr;
        switch (op.type) {
            case OpType::Limit: rr = r.submit_limit(op.side, op.price, op.qty, Tif::GTC, 0); break;
            case OpType::Cancel: rr = r.cancel(op.id); break;
            case OpType::Market: rr = r.submit_market(op.side, op.qty, 0); break;
            case OpType::Modify: rr = r.modify(op.id, op.price, op.qty); break;
        }
        const auto fa = a.last_fills();
        const auto& fr = r.last_fills();
        if (!(ra == rr) || fa.size() != fr.size() || !std::equal(fa.begin(), fa.end(), fr.begin())) {
            std::fprintf(stderr, "verify: mismatch at op %zu\n", i);
            return false;
        }
        if (i % 20000 == 0) {
            a.check_invariants();
            if (a.levels(Side::Buy, 100000) != r.levels(Side::Buy, 100000) ||
                a.levels(Side::Sell, 100000) != r.levels(Side::Sell, 100000)) {
                std::fprintf(stderr, "verify: level mismatch at op %zu\n", i);
                return false;
            }
        }
    }
    a.check_invariants();
    return a.levels(Side::Buy, 100000) == r.levels(Side::Buy, 100000) &&
           a.levels(Side::Sell, 100000) == r.levels(Side::Sell, 100000);
}
}  // namespace

int main(int argc, char** argv) {
    try {
        qc::Args args(argc, argv, {"ops", "warmup", "reps", "verify", "seed", "out"});
        const auto ops = static_cast<std::size_t>(args.get_int("ops", 5000000));
        const auto warmup = static_cast<std::size_t>(args.get_int("warmup", 1000000));
        const auto reps = static_cast<int>(args.get_int("reps", 5));
        const std::string out = args.get_str("out", OFL_RESULTS_DIR);
        std::printf("generating op stream (%zu + %zu ops) ...\n", warmup, ops);
        Stream s = generate(args.get_u64("seed", 42), warmup + ops);
        std::printf("max live orders %zu, mean book levels %.1f, fills/op %.3f\n", s.max_live, s.mean_levels, s.fills_per_op);
        if (args.has("verify")) {
            const bool ok = verify(s, 200000);
            std::printf("verify (first 200000 ops, OrderBook vs ReferenceBook): %s\n", ok ? "IDENTICAL" : "MISMATCH");
            if (!ok) return 1;
        }
        const std::size_t cap = s.max_live + s.max_live / 4 + 1024;
        qc::CsvWriter csv(out + "/bench_engine.csv", {"rep", "ops", "seconds", "ops_per_s", "ns_per_op"});
        volatile std::int64_t sink = 0;
        std::vector<double> rates;
        bool cap_ok = true;
        for (int rep = 0; rep < reps; ++rep) {
            OrderBook b(kMinTick, kMaxTick, cap, s.ops.size() + 4000);
            const std::size_t pc = b.pool_capacity();
            std::size_t i = 0;
            for (; i < s.n_prefill + warmup; ++i) {
                const OrderResult r = apply(b, s.ops[i]);
                sink = sink + r.filled + r.resting;
            }
            qc::Stopwatch sw;
            sw.start();
            for (; i < s.ops.size(); ++i) {
                const OrderResult r = apply(b, s.ops[i]);
                sink = sink + r.filled + r.resting;
            }
            const double sec = sw.elapsed_s();
            cap_ok = cap_ok && (b.pool_capacity() == pc);
            const double rate = static_cast<double>(ops) / sec;
            rates.push_back(rate);
            csv.row({std::int64_t{rep}, static_cast<std::int64_t>(ops), sec, rate, 1e9 / rate});
            std::printf("rep %d: %.3f s, %.3f Mops/s, %.1f ns/op\n", rep, sec, rate / 1e6, 1e9 / rate);
        }
        std::sort(rates.begin(), rates.end());
        const double med = rates[rates.size() / 2];
        std::printf("throughput median %.3f Mops/s (%.1f ns/op), min %.3f, max %.3f; pool_capacity unchanged: %s\n",
                    med / 1e6, 1e9 / med, rates.front() / 1e6, rates.back() / 1e6, cap_ok ? "yes" : "NO");
        // latency rep
        std::vector<std::uint32_t> lat(ops);
        {
            OrderBook b(kMinTick, kMaxTick, cap, s.ops.size() + 4000);
            std::size_t i = 0;
            for (; i < s.n_prefill + warmup; ++i) {
                const OrderResult r = apply(b, s.ops[i]);
                sink = sink + r.filled + r.resting;
            }
            for (std::size_t k = 0; i < s.ops.size(); ++i, ++k) {
                const std::uint64_t t0 = qc::now_ns();
                const OrderResult r = apply(b, s.ops[i]);
                const std::uint64_t t1 = qc::now_ns();
                lat[k] = static_cast<std::uint32_t>(std::min<std::uint64_t>(t1 - t0, 0xFFFFFFFFu));
                sink = sink + r.filled + r.resting;
            }
        }
        // Batched latency: time blocks of 64 consecutive ops (timer overhead amortised over the block), report ns per op.
        std::vector<std::uint32_t> batch;
        {
            OrderBook b(kMinTick, kMaxTick, cap, s.ops.size() + 4000);
            std::size_t i = 0;
            for (; i < s.n_prefill + warmup; ++i) {
                const OrderResult r = apply(b, s.ops[i]);
                sink = sink + r.filled + r.resting;
            }
            constexpr std::size_t kB = 64;
            for (; i + kB <= s.ops.size(); i += kB) {
                const std::uint64_t t0 = qc::now_ns();
                for (std::size_t k = 0; k < kB; ++k) {
                    const OrderResult r = apply(b, s.ops[i + k]);
                    sink = sink + r.filled + r.resting;
                }
                const std::uint64_t t1 = qc::now_ns();
                batch.push_back(static_cast<std::uint32_t>((t1 - t0) / kB));
            }
        }
        std::vector<std::uint64_t> ov(1000000);
        for (auto& o : ov) {
            const std::uint64_t t0 = qc::now_ns();
            const std::uint64_t t1 = qc::now_ns();
            o = t1 - t0;
        }
        std::nth_element(ov.begin(), ov.begin() + 500000, ov.end());
        const double timer_ns = static_cast<double>(ov[500000]);
        qc::CsvWriter lc(out + "/bench_latency.csv", {"op_type", "p50", "p90", "p99", "p999", "max", "count"});
        auto emit = [&](const std::string& name, std::vector<std::uint32_t> v) {
            if (v.empty()) return;
            const double mx = *std::max_element(v.begin(), v.end());
            const double p50 = pct(v, 0.5), p90 = pct(v, 0.9), p99 = pct(v, 0.99), p999 = pct(v, 0.999);
            lc.row({name, p50, p90, p99, p999, mx, static_cast<std::int64_t>(v.size())});
            std::printf("latency %-7s p50 %6.0f p90 %6.0f p99 %7.0f p99.9 %8.0f max %9.0f ns (n=%zu)\n", name.c_str(),
                        p50, p90, p99, p999, mx, v.size());
        };
        emit("all", lat);
        emit("batch64_per_op", batch);
        for (int t = 0; t < 4; ++t) {
            std::vector<std::uint32_t> v;
            for (std::size_t k = 0; k < ops; ++k)
                if (static_cast<int>(s.ops[s.n_prefill + warmup + k].type) == t) v.push_back(lat[k]);
            emit(kOpNames[t], std::move(v));
        }
        std::printf("timer overhead (median of 1e6 back-to-back now_ns pairs): %.0f ns (NOT subtracted)\n", timer_ns);
        std::printf("CPU: %s\ncompiler: %s\nflags: -O3 -DNDEBUG (no -march=native), single thread, not pinned\n",
                    cpu_model().c_str(), __clang_version__);
        (void)sink;
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
