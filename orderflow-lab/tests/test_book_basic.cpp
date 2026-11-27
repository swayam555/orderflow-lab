#include <qc/test.hpp>

#include <vector>

#include "ofl/order_book.hpp"

using namespace ofl;

namespace {
// Wrapper that runs check_invariants() after every call.
struct B {
    OrderBook b{0, 1000, 64, 256};
    OrderResult lim(Side s, Tick p, Qty q, Tif t = Tif::GTC, OwnerId o = 0) {
        auto r = b.submit_limit(s, p, q, t, o);
        b.check_invariants();
        return r;
    }
    OrderResult mkt(Side s, Qty q, OwnerId o = 0) {
        auto r = b.submit_market(s, q, o);
        b.check_invariants();
        return r;
    }
    OrderResult can(OrderId id) {
        auto r = b.cancel(id);
        b.check_invariants();
        return r;
    }
    OrderResult mod(OrderId id, Tick p, Qty q) {
        auto r = b.modify(id, p, q);
        b.check_invariants();
        return r;
    }
};
}  // namespace

TEST_CASE("book 1: empty book market order") {
    B x;
    CHECK_FALSE(x.b.best_bid().has_value());
    CHECK_FALSE(x.b.best_ask().has_value());
    auto r = x.mkt(Side::Buy, 5);
    CHECK(r.status == Status::Cancelled);
    CHECK(r.filled == 0);
    CHECK(r.cancelled == 5);
    CHECK(x.b.next_id() == 2);
}

TEST_CASE("book 2: rejections") {
    B x;
    auto r = x.lim(Side::Buy, 100, 0);
    CHECK(r.status == Status::Rejected);
    CHECK(r.reject == Reject::BadQty);
    CHECK(r.id == kNoOrder);
    CHECK(x.b.next_id() == 1);
    CHECK(x.lim(Side::Buy, -1, 5).reject == Reject::BadPrice);
    CHECK(x.lim(Side::Buy, 1001, 5).reject == Reject::BadPrice);
    CHECK(x.b.next_id() == 1);
    CHECK(x.can(999).reject == Reject::UnknownId);
    CHECK(x.mod(999, 100, 5).reject == Reject::UnknownId);
    auto live = x.lim(Side::Buy, 100, 5);
    const auto before = x.b.find(live.id);
    CHECK(x.mod(live.id, 100, 0).reject == Reject::BadQty);
    CHECK(x.mod(live.id, 2000, 5).reject == Reject::BadPrice);
    CHECK(x.b.find(live.id) == before);
}

TEST_CASE("book 3: resting") {
    B x;
    auto r = x.lim(Side::Buy, 100, 10);
    CHECK(r.status == Status::Resting);
    CHECK(r.resting == 10);
    CHECK(x.b.best_bid() == 100);
    CHECK(x.b.depth_at(Side::Buy, 100) == 10);
}

TEST_CASE("book 4: price priority") {
    B x;
    auto a = x.lim(Side::Sell, 102, 5);
    auto b = x.lim(Side::Sell, 101, 5);
    auto r = x.mkt(Side::Buy, 7);
    CHECK(r.status == Status::Filled);
    auto f = x.b.last_fills();
    REQUIRE(f.size() == 2);
    CHECK(f[0].maker == b.id);
    CHECK(f[0].qty == 5);
    CHECK(f[0].price == 101);
    CHECK(f[0].maker_done);
    CHECK(f[1].maker == a.id);
    CHECK(f[1].qty == 2);
    CHECK(f[1].price == 102);
    CHECK_FALSE(f[1].maker_done);
}

TEST_CASE("book 5: time priority") {
    B x;
    auto A = x.lim(Side::Sell, 101, 3);
    auto Bo = x.lim(Side::Sell, 101, 4);
    auto C = x.lim(Side::Sell, 101, 5);
    auto r = x.lim(Side::Buy, 101, 6);
    CHECK(r.status == Status::Filled);
    auto f = x.b.last_fills();
    REQUIRE(f.size() == 2);
    CHECK(f[0].maker == A.id);
    CHECK(f[0].qty == 3);
    CHECK(f[1].maker == Bo.id);
    CHECK(f[1].qty == 3);
    CHECK(x.b.find(Bo.id)->remaining == 1);
    CHECK(x.b.find(C.id)->remaining == 5);
}

TEST_CASE("book 6/7: limit price protection and IOC") {
    {
        B x;
        (void)x.lim(Side::Sell, 101, 5);
        (void)x.lim(Side::Sell, 103, 5);
        auto r = x.lim(Side::Buy, 102, 8);
        CHECK(r.status == Status::PartialResting);
        CHECK(r.filled == 5);
        CHECK(r.resting == 3);
        CHECK(x.b.best_bid() == 102);
        CHECK(x.b.best_ask() == 103);
    }
    {
        B x;
        (void)x.lim(Side::Sell, 101, 5);
        (void)x.lim(Side::Sell, 103, 5);
        auto r = x.lim(Side::Buy, 102, 8, Tif::IOC);
        CHECK(r.status == Status::PartialCancelled);
        CHECK(r.filled == 5);
        CHECK(r.cancelled == 3);
        CHECK_FALSE(x.b.best_bid().has_value());
        auto r2 = x.lim(Side::Buy, 100, 4, Tif::IOC);
        CHECK(r2.status == Status::Cancelled);
        CHECK(r2.cancelled == 4);
        CHECK_FALSE(x.b.best_bid().has_value());
    }
}

TEST_CASE("book 8: sell-side mirror") {
    B x;
    auto a = x.lim(Side::Buy, 98, 5);
    auto b = x.lim(Side::Buy, 99, 5);
    auto r = x.mkt(Side::Sell, 7);
    CHECK(r.status == Status::Filled);
    auto f = x.b.last_fills();
    REQUIRE(f.size() == 2);
    CHECK(f[0].maker == b.id);
    CHECK(f[0].price == 99);
    CHECK(f[1].maker == a.id);
    CHECK(f[1].qty == 2);
    CHECK(f[1].taker_side == Side::Sell);
    // time priority
    auto A = x.lim(Side::Buy, 98, 3);
    auto Bo = x.lim(Side::Buy, 98, 4);
    auto r2 = x.lim(Side::Sell, 98, 6);  // takes 3 remaining of a first, then A3
    CHECK(r2.status == Status::Filled);
    auto f2 = x.b.last_fills();
    REQUIRE(f2.size() == 2);
    CHECK(f2[0].maker == a.id);
    CHECK(f2[0].qty == 3);
    CHECK(f2[1].maker == A.id);
    CHECK(f2[1].qty == 3);
    CHECK(x.b.find(Bo.id)->remaining == 4);
    // protection
    B y;
    (void)y.lim(Side::Buy, 99, 5);
    (void)y.lim(Side::Buy, 97, 5);
    auto r3 = y.lim(Side::Sell, 98, 8);
    CHECK(r3.status == Status::PartialResting);
    CHECK(r3.filled == 5);
    CHECK(y.b.best_ask() == 98);
    CHECK(y.b.best_bid() == 97);
    auto r4 = y.lim(Side::Sell, 90, 10, Tif::IOC);
    CHECK(r4.status == Status::PartialCancelled);
    CHECK(r4.filled == 5);
    CHECK(r4.cancelled == 5);
}

TEST_CASE("book 9: cancel") {
    B x;
    auto A = x.lim(Side::Sell, 101, 1);
    auto Bo = x.lim(Side::Sell, 101, 2);
    auto C = x.lim(Side::Sell, 101, 3);
    auto r = x.can(Bo.id);
    CHECK(r.status == Status::Cancelled);
    CHECK(r.cancelled == 2);
    CHECK(x.b.depth_at(Side::Sell, 101) == 4);
    CHECK(x.b.orders_at(Side::Sell, 101) == 2);
    CHECK(x.can(Bo.id).reject == Reject::UnknownId);
    auto m = x.mkt(Side::Buy, 4);
    auto f = x.b.last_fills();
    REQUIRE(f.size() == 2);
    CHECK(f[0].maker == A.id);
    CHECK(f[1].maker == C.id);
    (void)m;
    // best moves across a gap larger than 64 ticks (and word boundary)
    B y;
    auto near = y.lim(Side::Sell, 200, 1);
    (void)y.lim(Side::Sell, 700, 1);
    CHECK(y.b.best_ask() == 200);
    (void)y.can(near.id);
    CHECK(y.b.best_ask() == 700);
    auto nb = y.lim(Side::Buy, 600, 1);
    (void)y.lim(Side::Buy, 100, 1);
    (void)y.lim(Side::Buy, 3, 1);
    (void)y.can(nb.id);
    CHECK(y.b.best_bid() == 100);
    auto lv = y.b.levels(Side::Buy, 10);
    REQUIRE(lv.size() == 2);
    CHECK(lv[0].price == 100);
    CHECK(lv[1].price == 3);
}

TEST_CASE("book 10-12: modify") {
    {
        B x;
        auto A = x.lim(Side::Sell, 101, 5);
        auto Bo = x.lim(Side::Sell, 101, 5);
        auto r = x.mod(A.id, 101, 2);
        CHECK(r.status == Status::Resting);
        CHECK(r.cancelled == 3);
        (void)x.mkt(Side::Buy, 3);
        auto f = x.b.last_fills();
        REQUIRE(f.size() == 2);
        CHECK(f[0].maker == A.id);
        CHECK(f[0].qty == 2);
        CHECK(f[1].maker == Bo.id);
        CHECK(f[1].qty == 1);
    }
    {
        B x;
        auto A = x.lim(Side::Sell, 101, 5);
        auto Bo = x.lim(Side::Sell, 101, 5);
        auto r = x.mod(A.id, 101, 8);
        CHECK(r.id == A.id);
        (void)x.mkt(Side::Buy, 3);
        auto f = x.b.last_fills();
        REQUIRE(f.size() == 1);
        CHECK(f[0].maker == Bo.id);
    }
    {
        B x;
        auto X = x.lim(Side::Buy, 99, 5);
        (void)x.lim(Side::Sell, 101, 3);
        auto r = x.mod(X.id, 101, 5);
        CHECK(r.status == Status::PartialResting);
        CHECK(r.id == X.id);
        CHECK(r.filled == 3);
        CHECK(r.resting == 2);
        auto f = x.b.last_fills();
        REQUIRE(f.size() == 1);
        CHECK(f[0].taker == X.id);
        CHECK(x.b.best_bid() == 101);
        CHECK(x.b.depth_at(Side::Buy, 101) == 2);
        CHECK_FALSE(x.b.best_ask().has_value());
    }
}

TEST_CASE("book 13: owners and seq") {
    B x;
    (void)x.lim(Side::Sell, 101, 2, Tif::GTC, 7);
    (void)x.lim(Side::Sell, 102, 2, Tif::GTC, 8);
    (void)x.mkt(Side::Buy, 3, 9);
    auto f = x.b.last_fills();
    REQUIRE(f.size() == 2);
    CHECK(f[0].taker_owner == 9);
    CHECK(f[0].maker_owner == 7);
    CHECK(f[1].maker_owner == 8);
    CHECK(f[1].seq > f[0].seq);
    const auto last = f[1].seq;
    (void)x.mkt(Side::Buy, 1, 9);
    CHECK(x.b.last_fills()[0].seq > last);
}

TEST_CASE("book 14: large sweep") {
    B x;
    for (Tick p = 500; p < 600; ++p) (void)x.lim(Side::Sell, p, 1);
    auto r = x.mkt(Side::Buy, 100);
    CHECK(r.status == Status::Filled);
    CHECK(x.b.last_fills().size() == 100);
    CHECK_FALSE(x.b.best_ask().has_value());
    CHECK(x.b.live_orders() == 0);
}

TEST_CASE("book 15: band edges") {
    B x;
    (void)x.lim(Side::Buy, 0, 2);
    (void)x.lim(Side::Sell, 1000, 2);
    CHECK(x.b.best_bid() == 0);
    CHECK(x.b.best_ask() == 1000);
    auto r = x.mkt(Side::Sell, 1);
    CHECK(x.b.last_fills()[0].price == 0);
    CHECK(r.status == Status::Filled);
    auto r2 = x.mkt(Side::Buy, 2);
    CHECK(x.b.last_fills()[0].price == 1000);
    (void)r2;
    auto r3 = x.lim(Side::Buy, 1000, 1);
    CHECK(r3.status == Status::Resting);
    auto r4 = x.lim(Side::Sell, 0, 1);
    CHECK(r4.status == Status::Filled);  // sell @0 takes the bid at 1000 first (price priority)
    CHECK(x.b.last_fills()[0].price == 1000);
    CHECK(x.b.best_bid() == 0);
}
