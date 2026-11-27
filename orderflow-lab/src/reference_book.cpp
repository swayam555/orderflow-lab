#include "ofl/reference_book.hpp"

#include <algorithm>

namespace ofl {

namespace {
OrderResult rej(Reject r, OrderId id = kNoOrder) { return OrderResult{id, Status::Rejected, r, 0, 0, 0}; }
}  // namespace

OrderResult ReferenceBook::run_limit(OrderId id, Side side, Tick price, Qty qty, Tif tif, OwnerId owner) {
    Qty rem = qty;
    const bool is_market = (price == kMarketSentinel_());
    Book& opp = book(side == Side::Buy ? Side::Sell : Side::Buy);
    while (rem > 0 && !opp.empty()) {
        auto it = side == Side::Buy ? opp.begin() : std::prev(opp.end());
        const Tick px = it->first;
        if (!is_market && (side == Side::Buy ? px > price : px < price)) break;
        auto& dq = it->second;
        while (rem > 0 && !dq.empty()) {
            Ord& m = dq.front();
            const Qty q = std::min(rem, m.remaining);
            rem -= q;
            m.remaining -= q;
            const bool done = m.remaining == 0;
            fills_.push_back(Fill{seq_++, id, m.id, owner, m.owner, px, q, side, done});
            if (done) {
                loc_.erase(m.id);
                dq.pop_front();
            }
        }
        if (dq.empty()) opp.erase(it);
    }
    const Qty filled = qty - rem;
    if (rem == 0) return {id, Status::Filled, Reject::None, filled, 0, 0};
    if (is_market || tif == Tif::IOC)
        return {id, filled > 0 ? Status::PartialCancelled : Status::Cancelled, Reject::None, filled, 0, rem};
    book(side)[price].push_back(Ord{id, owner, rem});
    loc_[id] = Loc{side, price};
    return {id, filled > 0 ? Status::PartialResting : Status::Resting, Reject::None, filled, rem, 0};
}

OrderResult ReferenceBook::submit_limit(Side side, Tick price, Qty qty, Tif tif, OwnerId owner) {
    fills_.clear();
    if (qty <= 0) return rej(Reject::BadQty);
    if (price < lo_ || price > hi_) return rej(Reject::BadPrice);
    return run_limit(next_id_++, side, price, qty, tif, owner);
}

OrderResult ReferenceBook::submit_market(Side side, Qty qty, OwnerId owner) {
    fills_.clear();
    if (qty <= 0) return rej(Reject::BadQty);
    return run_limit(next_id_++, side, kMarketSentinel_(), qty, Tif::IOC, owner);
}

OrderResult ReferenceBook::cancel(OrderId id) {
    fills_.clear();
    auto it = loc_.find(id);
    if (it == loc_.end()) return rej(Reject::UnknownId, id);
    auto& bk = book(it->second.side);
    auto lv = bk.find(it->second.price);
    auto& dq = lv->second;
    auto o = std::find_if(dq.begin(), dq.end(), [&](const Ord& x) { return x.id == id; });
    const Qty rem = o->remaining;
    dq.erase(o);
    if (dq.empty()) bk.erase(lv);
    loc_.erase(it);
    return {id, Status::Cancelled, Reject::None, 0, 0, rem};
}

OrderResult ReferenceBook::modify(OrderId id, Tick new_price, Qty new_qty) {
    fills_.clear();
    auto it = loc_.find(id);
    if (it == loc_.end()) return rej(Reject::UnknownId, id);
    if (new_qty <= 0) return rej(Reject::BadQty, id);
    if (new_price < lo_ || new_price > hi_) return rej(Reject::BadPrice, id);
    const Loc l = it->second;
    auto& bk = book(l.side);
    auto lv = bk.find(l.price);
    auto& dq = lv->second;
    auto o = std::find_if(dq.begin(), dq.end(), [&](const Ord& x) { return x.id == id; });
    if (new_price == l.price && new_qty <= o->remaining) {
        const Qty diff = o->remaining - new_qty;
        o->remaining = new_qty;
        return {id, Status::Resting, Reject::None, 0, new_qty, diff};
    }
    const OwnerId owner = o->owner;
    dq.erase(o);
    if (dq.empty()) bk.erase(lv);
    loc_.erase(id);
    return run_limit(id, l.side, new_price, new_qty, Tif::GTC, owner);
}

std::optional<Tick> ReferenceBook::best_bid() const {
    if (bids_.empty()) return std::nullopt;
    return bids_.rbegin()->first;
}
std::optional<Tick> ReferenceBook::best_ask() const {
    if (asks_.empty()) return std::nullopt;
    return asks_.begin()->first;
}

std::vector<LevelView> ReferenceBook::levels(Side side, std::size_t max_levels) const {
    std::vector<LevelView> out;
    auto add = [&](Tick p, const std::deque<Ord>& dq) {
        Qty q = 0;
        for (const Ord& o : dq) q += o.remaining;
        out.push_back({p, q, static_cast<std::int32_t>(dq.size())});
    };
    if (side == Side::Buy) {
        for (auto it = bids_.rbegin(); it != bids_.rend() && out.size() < max_levels; ++it) add(it->first, it->second);
    } else {
        for (auto it = asks_.begin(); it != asks_.end() && out.size() < max_levels; ++it) add(it->first, it->second);
    }
    return out;
}

std::optional<OrderView> ReferenceBook::find(OrderId id) const {
    auto it = loc_.find(id);
    if (it == loc_.end()) return std::nullopt;
    const Book& bk = it->second.side == Side::Buy ? bids_ : asks_;
    for (const Ord& o : bk.at(it->second.price))
        if (o.id == id) return OrderView{o.id, o.owner, it->second.side, it->second.price, o.remaining};
    return std::nullopt;
}

std::size_t ReferenceBook::live_orders() const { return loc_.size(); }

}  // namespace ofl
