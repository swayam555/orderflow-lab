// Single-instrument limit order book: price-time priority, GTC/IOC limit, market, cancel, modify.
#pragma once
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "ofl/types.hpp"

namespace ofl {

class OrderBook {
  public:
    // Valid prices are [min_tick, max_tick] (ticks); reserves pool for expected_live_orders (lots-orders)
    // and the id table for expected_total_ids.
    OrderBook(Tick min_tick, Tick max_tick, std::size_t expected_live_orders, std::size_t expected_total_ids);
    // Submits a limit order (price in ticks, qty in lots).
    [[nodiscard]] OrderResult submit_limit(Side side, Tick price, Qty qty, Tif tif, OwnerId owner);
    // Submits a market order (qty in lots); any unfilled remainder is cancelled.
    [[nodiscard]] OrderResult submit_market(Side side, Qty qty, OwnerId owner);
    // Cancels a live order.
    [[nodiscard]] OrderResult cancel(OrderId id);
    // Modifies price (ticks) and remaining quantity (lots) of a live order; see DESIGN 3.3.
    [[nodiscard]] OrderResult modify(OrderId id, Tick new_price, Qty new_qty);
    // Fills of the most recent call only; invalidated by the next call.
    [[nodiscard]] std::span<const Fill> last_fills() const noexcept { return fills_; }
    [[nodiscard]] std::optional<Tick> best_bid() const noexcept;
    [[nodiscard]] std::optional<Tick> best_ask() const noexcept;
    // Total resting lots at a price.
    [[nodiscard]] Qty depth_at(Side side, Tick price) const;
    // Number of resting orders at a price.
    [[nodiscard]] std::int32_t orders_at(Side side, Tick price) const;
    // Up to max_levels non-empty levels, best first.
    [[nodiscard]] std::vector<LevelView> levels(Side side, std::size_t max_levels) const;
    [[nodiscard]] std::optional<OrderView> find(OrderId id) const;
    [[nodiscard]] std::size_t live_orders() const noexcept { return live_; }
    // Id the next accepted submit will receive.
    [[nodiscard]] OrderId next_id() const noexcept { return next_id_; }
    // O(book) consistency check; throws std::logic_error describing the first violation.
    void check_invariants() const;
    [[nodiscard]] std::size_t pool_capacity() const noexcept { return pool_.size(); }
    [[nodiscard]] Tick min_tick() const noexcept { return min_tick_; }
    [[nodiscard]] Tick max_tick() const noexcept { return max_tick_; }

  private:
    static constexpr std::uint32_t kNil = 0xFFFFFFFFu;
    struct Node {
        OrderId id;
        Qty remaining;
        Tick price;
        OwnerId owner;
        Side side;
        std::uint32_t prev, next;
    };
    struct Level {
        std::uint32_t head = kNil, tail = kNil;
        Qty qty = 0;
        std::int32_t count = 0;
    };
    struct SideBook {
        std::vector<Level> levels;
        std::vector<std::uint64_t> bits;
        std::int64_t best = -1;  // index into levels, -1 if empty
    };

    SideBook& sb(Side s) noexcept { return s == Side::Buy ? bid_ : ask_; }
    const SideBook& sb(Side s) const noexcept { return s == Side::Buy ? bid_ : ask_; }
    std::uint32_t alloc_node();
    void free_node(std::uint32_t slot) noexcept;
    void ensure_id_capacity(OrderId id);
    void link_tail(SideBook& b, std::uint32_t slot);
    void unlink(SideBook& b, std::uint32_t slot);
    void recompute_best(SideBook& b, Side s) noexcept;
    // Matches an incoming order against the opposite side; returns the filled quantity.
    Qty match(Side side, bool has_limit, Tick limit, Qty qty, OrderId taker, OwnerId owner);
    OrderResult finish_limit(OrderId id, Side side, Tick price, Qty qty, Tif tif, OwnerId owner);

    Tick min_tick_, max_tick_;
    SideBook bid_, ask_;
    std::vector<Node> pool_;
    std::uint32_t free_head_ = kNil;
    std::vector<std::uint32_t> slot_of_id_;
    std::vector<Fill> fills_;
    OrderId next_id_ = 1;
    std::uint64_t seq_ = 0;
    std::size_t live_ = 0;
};

}  // namespace ofl
