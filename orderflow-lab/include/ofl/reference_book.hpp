// Naive, obviously-correct reference order book (std::map of deques). Same semantics as OrderBook; tests only.
#pragma once
#include <cstddef>
#include <deque>
#include <limits>
#include <map>
#include <optional>
#include <vector>

#include "ofl/types.hpp"

namespace ofl {

class ReferenceBook {
  public:
    // Valid prices are [min_tick, max_tick] (ticks).
    ReferenceBook(Tick min_tick, Tick max_tick) : lo_(min_tick), hi_(max_tick) {}
    [[nodiscard]] OrderResult submit_limit(Side side, Tick price, Qty qty, Tif tif, OwnerId owner);
    [[nodiscard]] OrderResult submit_market(Side side, Qty qty, OwnerId owner);
    [[nodiscard]] OrderResult cancel(OrderId id);
    [[nodiscard]] OrderResult modify(OrderId id, Tick new_price, Qty new_qty);
    [[nodiscard]] const std::vector<Fill>& last_fills() const noexcept { return fills_; }
    [[nodiscard]] std::optional<Tick> best_bid() const;
    [[nodiscard]] std::optional<Tick> best_ask() const;
    // Up to max_levels non-empty levels, best first.
    [[nodiscard]] std::vector<LevelView> levels(Side side, std::size_t max_levels) const;
    [[nodiscard]] std::optional<OrderView> find(OrderId id) const;
    [[nodiscard]] std::size_t live_orders() const;
    [[nodiscard]] OrderId next_id() const noexcept { return next_id_; }

  private:
    struct Ord {
        OrderId id;
        OwnerId owner;
        Qty remaining;
    };
    using Book = std::map<Tick, std::deque<Ord>>;
    struct Loc {
        Side side;
        Tick price;
    };
    Book& book(Side s) { return s == Side::Buy ? bids_ : asks_; }
    OrderResult run_limit(OrderId id, Side side, Tick price, Qty qty, Tif tif, OwnerId owner);
    static constexpr Tick kMarketSentinel_() { return std::numeric_limits<Tick>::min(); }
    Tick lo_, hi_;
    Book bids_, asks_;
    std::map<OrderId, Loc> loc_;
    std::vector<Fill> fills_;
    OrderId next_id_ = 1;
    std::uint64_t seq_ = 0;
};

}  // namespace ofl
