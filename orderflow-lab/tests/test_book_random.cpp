#include <qc/test.hpp>

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include <qc/rng.hpp>

#include "ofl/order_book.hpp"
#include "ofl/reference_book.hpp"

using namespace ofl;

namespace {
constexpr Tick kLo = 0, kHi = 300;

struct SeqResult {
    bool ok = true;
    std::string msg;
    std::size_t pool_before = 0, pool_after = 0;
};

SeqResult run_sequence(std::uint64_t seq, int n_ops, int invariant_every) {
    qc::RngStream rng(11, seq, 1);
    OrderBook book(kLo, kHi, 10000, static_cast<std::size_t>(n_ops) + 16);
    ReferenceBook ref(kLo, kHi);
    SeqResult res;
    res.pool_before = book.pool_capacity();
    std::vector<OrderId> issued;
    std::int64_t entered = 0, cancelled = 0, taker_filled = 0, maker_filled = 0;
    auto fail = [&](const std::string& m, int i) {
        if (res.ok) res.msg = "seq " + std::to_string(seq) + " op " + std::to_string(i) + ": " + m;
        res.ok = false;
    };
    auto pick_id = [&](bool prefer_live) -> OrderId {
        if (issued.empty()) return static_cast<OrderId>(rng.uniform_int(1, 50));
        if (prefer_live) {
            for (int tries = 0; tries < 8; ++tries) {
                const OrderId id = issued[static_cast<std::size_t>(rng.uniform_int(0, static_cast<std::int64_t>(issued.size()) - 1))];
                if (ref.find(id)) return id;
            }
        }
        return static_cast<OrderId>(rng.uniform_int(0, static_cast<std::int64_t>(issued.size()) + 20));
    };
    auto rand_qty = [&]() -> Qty { return rng.uniform() < 0.02 ? rng.uniform_int(-2, 0) : rng.uniform_int(1, 10); };
    for (int i = 0; i < n_ops && res.ok; ++i) {
        const double u = rng.uniform();
        OrderResult a, b;
        OrderId touched = kNoOrder;
        Qty old_remaining = 0;
        bool in_place = false;
        enum { Sub, Can, Mod } kind = Sub;
        Qty sub_qty = 0;
        if (u < 0.40 || (u >= 0.40 && u < 0.50)) {
            const Side s = rng.uniform() < 0.5 ? Side::Buy : Side::Sell;
            const Tick p = static_cast<Tick>(rng.uniform_int(kLo, kHi));
            const Qty q = rand_qty();
            const Tif t = u < 0.40 ? Tif::GTC : Tif::IOC;
            a = book.submit_limit(s, p, q, t, 1);
            b = ref.submit_limit(s, p, q, t, 1);
            sub_qty = q;
        } else if (u < 0.60) {
            const Side s = rng.uniform() < 0.5 ? Side::Buy : Side::Sell;
            const Qty q = rand_qty();
            a = book.submit_market(s, q, 2);
            b = ref.submit_market(s, q, 2);
            sub_qty = q;
        } else if (u < 0.85) {
            kind = Can;
            const OrderId id = pick_id(rng.uniform() < 0.7);
            a = book.cancel(id);
            b = ref.cancel(id);
        } else {
            kind = Mod;
            const OrderId id = pick_id(rng.uniform() < 0.85);
            touched = id;
            const auto before = ref.find(id);
            if (before) old_remaining = before->remaining;
            const Tick cur = before ? before->price : 150;
            const Tick np = rng.uniform() < 0.5 ? cur : static_cast<Tick>(rng.uniform_int(kLo - 5, kHi + 5));
            const Qty nq = rng.uniform() < 0.5 ? rng.uniform_int(-1, 10) : (before ? rng.uniform_int(1, before->remaining + 3) : 3);
            in_place = before && np == before->price && nq <= before->remaining;
            a = book.modify(id, np, nq);
            b = ref.modify(id, np, nq);
            sub_qty = nq;
        }
        if (!(a == b)) fail("OrderResult mismatch", i);
        const auto fa = book.last_fills();
        const auto& fr = ref.last_fills();
        if (fa.size() != fr.size() || !std::equal(fa.begin(), fa.end(), fr.begin())) fail("fills mismatch", i);
        // ledger
        if (a.reject == Reject::None) {
            if (kind == Sub) {
                entered += sub_qty;
                issued.push_back(a.id);
                cancelled += a.cancelled;
            } else if (kind == Can) {
                cancelled += a.cancelled;
            } else {
                if (in_place) {
                    cancelled += a.cancelled;  // in-place reduction (or no-op)
                } else {
                    cancelled += old_remaining;  // removed and re-submitted
                    entered += sub_qty;
                }
            }
        }
        for (const Fill& f : fa) {
            taker_filled += f.qty;
            maker_filled += f.qty;
        }
        if (kind == Mod && a.reject == Reject::None && !ref.find(touched) != !book.find(touched)) fail("find mismatch", i);
        if (invariant_every > 0 && i % invariant_every == 0) {
            try {
                book.check_invariants();
            } catch (const std::exception& e) {
                fail(e.what(), i);
            }
        }
        if (i % 50 == 0) {
            if (book.levels(Side::Buy, 1000) != ref.levels(Side::Buy, 1000) ||
                book.levels(Side::Sell, 1000) != ref.levels(Side::Sell, 1000))
                fail("levels mismatch", i);
            if (book.best_bid() != ref.best_bid() || book.best_ask() != ref.best_ask()) fail("best mismatch", i);
            if (book.live_orders() != ref.live_orders()) fail("live count mismatch", i);
        }
    }
    book.check_invariants();
    if (book.levels(Side::Buy, 1000) != ref.levels(Side::Buy, 1000) ||
        book.levels(Side::Sell, 1000) != ref.levels(Side::Sell, 1000))
        fail("final levels mismatch", n_ops);
    std::int64_t resting = 0;
    for (Side s : {Side::Buy, Side::Sell})
        for (const auto& l : book.levels(s, 1000)) resting += l.qty;
    // Each submitted lot ends as: taker-filled, maker-filled (later), cancelled, or still resting.
    if (entered != taker_filled + maker_filled + cancelled + resting) fail("conservation violated", n_ops);
    res.pool_after = book.pool_capacity();
    return res;
}
}  // namespace

#ifndef OFL_SLOW_RANDOM
TEST_CASE("randomized differential vs reference book (unit: 200 x 2000 ops)") {
    for (std::uint64_t s = 0; s < 200; ++s) {
        const SeqResult r = run_sequence(s, 2000, 1);
        CHECK_MESSAGE(r.ok, r.msg);
        CHECK(r.pool_before == r.pool_after);
        if (!r.ok) break;
    }
}
#endif
