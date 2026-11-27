#include "ofl/order_book.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <stdexcept>
#include <string>

namespace ofl {

namespace {
[[noreturn]] void fail(const std::string& what) { throw std::logic_error("OrderBook invariant: " + what); }

OrderResult rejected(Reject r, OrderId id = kNoOrder) {
    return OrderResult{id, Status::Rejected, r, 0, 0, 0};
}

// Highest set bit index <= idx, or -1.
std::int64_t find_at_or_below(const std::vector<std::uint64_t>& bits, std::int64_t idx) noexcept {
    if (idx < 0) return -1;
    std::int64_t w = idx >> 6;
    const unsigned off = static_cast<unsigned>(idx & 63);
    std::uint64_t word = bits[static_cast<std::size_t>(w)];
    if (off != 63) word &= (std::uint64_t{1} << (off + 1)) - 1;
    while (true) {
        if (word != 0) return (w << 6) + 63 - std::countl_zero(word);
        if (--w < 0) return -1;
        word = bits[static_cast<std::size_t>(w)];
    }
}
// Lowest set bit index >= idx, or -1.
std::int64_t find_at_or_above(const std::vector<std::uint64_t>& bits, std::int64_t idx) noexcept {
    const auto nw = static_cast<std::int64_t>(bits.size());
    std::int64_t w = idx >> 6;
    if (w >= nw) return -1;
    std::uint64_t word = bits[static_cast<std::size_t>(w)] & (~std::uint64_t{0} << (idx & 63));
    while (true) {
        if (word != 0) return (w << 6) + std::countr_zero(word);
        if (++w >= nw) return -1;
        word = bits[static_cast<std::size_t>(w)];
    }
}
}  // namespace

OrderBook::OrderBook(Tick min_tick, Tick max_tick, std::size_t expected_live_orders, std::size_t expected_total_ids)
    : min_tick_(min_tick), max_tick_(max_tick) {
    if (max_tick < min_tick)
        throw std::invalid_argument("OrderBook: max_tick " + std::to_string(max_tick) + " < min_tick " +
                                    std::to_string(min_tick));
    if (expected_live_orders == 0) throw std::invalid_argument("OrderBook: expected_live_orders must be > 0");
    const auto n = static_cast<std::size_t>(static_cast<std::int64_t>(max_tick) - min_tick + 1);
    for (SideBook* b : {&bid_, &ask_}) {
        b->levels.assign(n, Level{});
        b->bits.assign((n + 63) / 64, 0);
        b->best = -1;
    }
    pool_.resize(expected_live_orders);
    for (std::size_t i = 0; i < pool_.size(); ++i) {
        pool_[i] = Node{0, 0, 0, 0, Side::Buy, kNil, i + 1 < pool_.size() ? static_cast<std::uint32_t>(i + 1) : kNil};
    }
    free_head_ = 0;
    slot_of_id_.assign(std::max<std::size_t>(expected_total_ids, 2) + 1, kNil);
    fills_.reserve(1024);
}

std::uint32_t OrderBook::alloc_node() {
    if (free_head_ == kNil) {  // grow by doubling (allowed, not expected in benchmarks)
        const std::size_t old = pool_.size();
        pool_.resize(old * 2);
        for (std::size_t i = old; i < pool_.size(); ++i) {
            pool_[i] = Node{0, 0, 0, 0, Side::Buy, kNil,
                            i + 1 < pool_.size() ? static_cast<std::uint32_t>(i + 1) : kNil};
        }
        free_head_ = static_cast<std::uint32_t>(old);
    }
    const std::uint32_t s = free_head_;
    free_head_ = pool_[s].next;
    return s;
}

void OrderBook::free_node(std::uint32_t slot) noexcept {
    Node& n = pool_[slot];
    slot_of_id_[n.id] = kNil;
    n.next = free_head_;
    n.prev = kNil;
    free_head_ = slot;
}

void OrderBook::ensure_id_capacity(OrderId id) {
    if (id >= slot_of_id_.size()) slot_of_id_.resize(std::max<std::size_t>(slot_of_id_.size() * 2, id + 1), kNil);
}

void OrderBook::link_tail(SideBook& b, std::uint32_t slot) {
    Node& n = pool_[slot];
    const auto idx = static_cast<std::size_t>(static_cast<std::int64_t>(n.price) - min_tick_);
    Level& L = b.levels[idx];
    n.prev = L.tail;
    n.next = kNil;
    if (L.tail != kNil)
        pool_[L.tail].next = slot;
    else
        L.head = slot;
    L.tail = slot;
    L.qty += n.remaining;
    ++L.count;
    b.bits[idx >> 6] |= std::uint64_t{1} << (idx & 63);
    const auto sidx = static_cast<std::int64_t>(idx);
    if (b.best < 0 || (n.side == Side::Buy ? sidx > b.best : sidx < b.best)) b.best = sidx;
}

// Removes the node from its level (qty/count/bitmap/best updated); does not free the slot.
void OrderBook::unlink(SideBook& b, std::uint32_t slot) {
    Node& n = pool_[slot];
    const auto idx = static_cast<std::size_t>(static_cast<std::int64_t>(n.price) - min_tick_);
    Level& L = b.levels[idx];
    if (n.prev != kNil)
        pool_[n.prev].next = n.next;
    else
        L.head = n.next;
    if (n.next != kNil)
        pool_[n.next].prev = n.prev;
    else
        L.tail = n.prev;
    L.qty -= n.remaining;
    --L.count;
    if (L.count == 0) {
        b.bits[idx >> 6] &= ~(std::uint64_t{1} << (idx & 63));
        if (b.best == static_cast<std::int64_t>(idx)) recompute_best(b, n.side);
    }
}

void OrderBook::recompute_best(SideBook& b, Side s) noexcept {
    if (b.best < 0) return;
    b.best = (s == Side::Buy) ? find_at_or_below(b.bits, b.best) : find_at_or_above(b.bits, b.best);
}

Qty OrderBook::match(Side side, bool has_limit, Tick limit, Qty qty, OrderId taker, OwnerId owner) {
    SideBook& opp = sb(side == Side::Buy ? Side::Sell : Side::Buy);
    const Side opp_side = side == Side::Buy ? Side::Sell : Side::Buy;
    Qty remaining = qty;
    while (remaining > 0 && opp.best >= 0) {
        const auto px = static_cast<Tick>(opp.best + min_tick_);
        if (has_limit && (side == Side::Buy ? px > limit : px < limit)) break;
        Level& L = opp.levels[static_cast<std::size_t>(opp.best)];
        while (remaining > 0 && L.head != kNil) {
            const std::uint32_t slot = L.head;
            Node& m = pool_[slot];
            const Qty q = std::min(remaining, m.remaining);
            remaining -= q;
            const bool done = (m.remaining == q);
            fills_.push_back(Fill{seq_++, taker, m.id, owner, m.owner, px, q, side, done});
            if (done) {
                // unlink the head (it is the level head, so no prev)
                L.head = m.next;
                if (m.next != kNil)
                    pool_[m.next].prev = kNil;
                else
                    L.tail = kNil;
                L.qty -= q;
                --L.count;
                --live_;
                free_node(slot);
            } else {
                m.remaining -= q;
                L.qty -= q;
            }
        }
        if (L.count == 0) {
            const auto idx = static_cast<std::size_t>(opp.best);
            opp.bits[idx >> 6] &= ~(std::uint64_t{1} << (idx & 63));
            recompute_best(opp, opp_side);
        }
    }
    return qty - remaining;
}

OrderResult OrderBook::finish_limit(OrderId id, Side side, Tick price, Qty qty, Tif tif, OwnerId owner) {
    const Qty filled = match(side, true, price, qty, id, owner);
    const Qty rem = qty - filled;
    if (rem == 0) return {id, Status::Filled, Reject::None, filled, 0, 0};
    if (tif == Tif::IOC)
        return {id, filled > 0 ? Status::PartialCancelled : Status::Cancelled, Reject::None, filled, 0, rem};
    const std::uint32_t slot = alloc_node();
    Node& n = pool_[slot];
    n.id = id;
    n.remaining = rem;
    n.price = price;
    n.owner = owner;
    n.side = side;
    ensure_id_capacity(id);
    slot_of_id_[id] = slot;
    link_tail(sb(side), slot);
    ++live_;
    return {id, filled > 0 ? Status::PartialResting : Status::Resting, Reject::None, filled, rem, 0};
}

OrderResult OrderBook::submit_limit(Side side, Tick price, Qty qty, Tif tif, OwnerId owner) {
    fills_.clear();
    if (qty <= 0) return rejected(Reject::BadQty);
    if (price < min_tick_ || price > max_tick_) return rejected(Reject::BadPrice);
    const OrderId id = next_id_++;
    return finish_limit(id, side, price, qty, tif, owner);
}

OrderResult OrderBook::submit_market(Side side, Qty qty, OwnerId owner) {
    fills_.clear();
    if (qty <= 0) return rejected(Reject::BadQty);
    const OrderId id = next_id_++;
    const Qty filled = match(side, false, 0, qty, id, owner);
    const Qty rem = qty - filled;
    if (rem == 0) return {id, Status::Filled, Reject::None, filled, 0, 0};
    return {id, filled > 0 ? Status::PartialCancelled : Status::Cancelled, Reject::None, filled, 0, rem};
}

OrderResult OrderBook::cancel(OrderId id) {
    fills_.clear();
    if (id == kNoOrder || id >= slot_of_id_.size() || slot_of_id_[id] == kNil) return rejected(Reject::UnknownId, id);
    const std::uint32_t slot = slot_of_id_[id];
    const Qty rem = pool_[slot].remaining;
    unlink(sb(pool_[slot].side), slot);
    free_node(slot);
    --live_;
    return {id, Status::Cancelled, Reject::None, 0, 0, rem};
}

OrderResult OrderBook::modify(OrderId id, Tick new_price, Qty new_qty) {
    fills_.clear();
    if (id == kNoOrder || id >= slot_of_id_.size() || slot_of_id_[id] == kNil) return rejected(Reject::UnknownId, id);
    if (new_qty <= 0) return rejected(Reject::BadQty, id);
    if (new_price < min_tick_ || new_price > max_tick_) return rejected(Reject::BadPrice, id);
    const std::uint32_t slot = slot_of_id_[id];
    Node& n = pool_[slot];
    if (new_price == n.price && new_qty <= n.remaining) {
        const Qty diff = n.remaining - new_qty;
        const auto idx = static_cast<std::size_t>(static_cast<std::int64_t>(n.price) - min_tick_);
        sb(n.side).levels[idx].qty -= diff;
        n.remaining = new_qty;
        return {id, Status::Resting, Reject::None, 0, new_qty, diff};
    }
    const Side side = n.side;
    const OwnerId owner = n.owner;
    unlink(sb(side), slot);
    free_node(slot);
    --live_;
    return finish_limit(id, side, new_price, new_qty, Tif::GTC, owner);
}

std::optional<Tick> OrderBook::best_bid() const noexcept {
    if (bid_.best < 0) return std::nullopt;
    return static_cast<Tick>(bid_.best + min_tick_);
}
std::optional<Tick> OrderBook::best_ask() const noexcept {
    if (ask_.best < 0) return std::nullopt;
    return static_cast<Tick>(ask_.best + min_tick_);
}

Qty OrderBook::depth_at(Side side, Tick price) const {
    if (price < min_tick_ || price > max_tick_) return 0;
    return sb(side).levels[static_cast<std::size_t>(static_cast<std::int64_t>(price) - min_tick_)].qty;
}
std::int32_t OrderBook::orders_at(Side side, Tick price) const {
    if (price < min_tick_ || price > max_tick_) return 0;
    return sb(side).levels[static_cast<std::size_t>(static_cast<std::int64_t>(price) - min_tick_)].count;
}

std::vector<LevelView> OrderBook::levels(Side side, std::size_t max_levels) const {
    std::vector<LevelView> out;
    const SideBook& b = sb(side);
    std::int64_t idx = b.best;
    while (idx >= 0 && out.size() < max_levels) {
        const Level& L = b.levels[static_cast<std::size_t>(idx)];
        out.push_back({static_cast<Tick>(idx + min_tick_), L.qty, L.count});
        idx = side == Side::Buy ? find_at_or_below(b.bits, idx - 1) : find_at_or_above(b.bits, idx + 1);
    }
    return out;
}

std::optional<OrderView> OrderBook::find(OrderId id) const {
    if (id == kNoOrder || id >= slot_of_id_.size() || slot_of_id_[id] == kNil) return std::nullopt;
    const Node& n = pool_[slot_of_id_[id]];
    return OrderView{n.id, n.owner, n.side, n.price, n.remaining};
}

void OrderBook::check_invariants() const {
    std::size_t total_live = 0;
    for (Side s : {Side::Buy, Side::Sell}) {
        const SideBook& b = sb(s);
        std::int64_t actual_best = -1;
        for (std::size_t idx = 0; idx < b.levels.size(); ++idx) {
            const Level& L = b.levels[idx];
            const bool bit = (b.bits[idx >> 6] >> (idx & 63)) & 1u;
            if (bit != (L.count > 0)) fail("bitmap bit mismatch at index " + std::to_string(idx));
            if (L.count == 0) {
                if (L.head != kNil || L.tail != kNil || L.qty != 0) fail("empty level not clean at " + std::to_string(idx));
                continue;
            }
            if (s == Side::Buy || actual_best < 0) actual_best = static_cast<std::int64_t>(idx);
            Qty sum = 0;
            std::int32_t cnt = 0;
            std::uint32_t prev = kNil;
            for (std::uint32_t cur = L.head; cur != kNil; cur = pool_[cur].next) {
                const Node& n = pool_[cur];
                if (n.prev != prev) fail("prev link mismatch at level " + std::to_string(idx));
                if (n.remaining <= 0) fail("non-positive remaining");
                if (static_cast<std::int64_t>(n.price) - min_tick_ != static_cast<std::int64_t>(idx))
                    fail("node price != level price");
                if (n.side != s) fail("node side mismatch");
                if (n.id == kNoOrder || n.id >= slot_of_id_.size() || slot_of_id_[n.id] != cur)
                    fail("id table does not point back to node " + std::to_string(n.id));
                sum += n.remaining;
                ++cnt;
                prev = cur;
                if (cnt > static_cast<std::int32_t>(pool_.size())) fail("cycle in level list");
            }
            if (prev != L.tail) fail("tail mismatch at level " + std::to_string(idx));
            if (sum != L.qty) fail("level qty != sum at " + std::to_string(idx));
            if (cnt != L.count) fail("level count != nodes at " + std::to_string(idx));
            total_live += static_cast<std::size_t>(cnt);
        }
        if (actual_best != b.best) fail("cached best != actual best");
    }
    if (total_live != live_) fail("live count mismatch");
    std::size_t in_table = 0;
    for (std::size_t id = 0; id < slot_of_id_.size(); ++id) {
        const std::uint32_t s = slot_of_id_[id];
        if (s == kNil) continue;
        ++in_table;
        if (s >= pool_.size() || pool_[s].id != id) fail("id table entry points to wrong node");
    }
    if (in_table != live_) fail("id table size != live orders");
    if (bid_.best >= 0 && ask_.best >= 0 && bid_.best + min_tick_ >= ask_.best + min_tick_) fail("book is crossed");
}

}  // namespace ofl
