// Core types of the order-flow lab: time in seconds, price in integer ticks, quantity in integer lots.
#pragma once
#include <cstdint>

namespace ofl {
using Tick = std::int32_t;
using Qty = std::int64_t;
using OrderId = std::uint64_t;
using OwnerId = std::uint32_t;
inline constexpr OrderId kNoOrder = 0;     // ids start at 1
inline constexpr OwnerId kBackground = 0;  // agents get 1,2,...
enum class Side : std::uint8_t { Buy = 0, Sell = 1 };
enum class Tif : std::uint8_t { GTC, IOC };
enum class Status : std::uint8_t { Resting, Filled, PartialResting, PartialCancelled, Cancelled, Rejected };
enum class Reject : std::uint8_t { None, BadQty, BadPrice, UnknownId };

// One trade; price = maker's resting price (ticks), qty in lots.
struct Fill {
    std::uint64_t seq;
    OrderId taker, maker;
    OwnerId taker_owner, maker_owner;
    Tick price;
    Qty qty;
    Side taker_side;
    bool maker_done;
    bool operator==(const Fill&) const = default;
};
struct OrderResult {
    OrderId id;
    Status status;
    Reject reject;
    Qty filled, resting, cancelled;
    bool operator==(const OrderResult&) const = default;
};
struct LevelView {
    Tick price;
    Qty qty;
    std::int32_t orders;
    bool operator==(const LevelView&) const = default;
};
struct OrderView {
    OrderId id;
    OwnerId owner;
    Side side;
    Tick price;
    Qty remaining;
    bool operator==(const OrderView&) const = default;
};
}  // namespace ofl
