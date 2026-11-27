// Discrete-event market simulator: Hawkes background order flow + per-order cancellation hazard + agents.
#pragma once
#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <vector>

#include <qc/linalg.hpp>
#include <qc/rng.hpp>

#include "ofl/hawkes.hpp"
#include "ofl/order_book.hpp"
#include "ofl/types.hpp"

namespace ofl {

// Default branching matrix (rows = affected type LB,LS,MB,MS; cols = source type), dimensionless.
[[nodiscard]] qc::Mat default_G();

struct SimConfig {
    qc::Mat G = default_G();
    double beta = 2.0;                                       // kernel decay, 1/s (same for every pair)
    std::array<double, 4> lambda_bar{1.0, 1.0, 0.25, 0.25};  // stationary Hawkes rates, 1/s
    double cancel_hazard = 0.01;                             // nu, 1/s per resting background order
    int placement_K = 20;                                    // limit orders placed 1..K ticks behind the opposite best
    Qty size_min = 1, size_max = 4;                          // lots
    Tick min_tick = 0, max_tick = 20000;
    Tick init_best_bid = 9999, init_best_ask = 10001;
    int init_levels = 20, init_orders_per_level = 3;
    Qty init_order_size = 3;
    bool agent_market_orders_excite = true;  // agent MB/MS are injected into the Hawkes state as type 2/3 events
    // Throws std::invalid_argument naming the offending parameter.
    void validate() const;
};

class MarketSim;

class Agent {
  public:
    virtual ~Agent() = default;
    // Absolute time (s) of the next wakeup; +inf = never.
    [[nodiscard]] virtual double next_wakeup() const = 0;
    // Called at the wakeup time; may place/cancel orders via sim.agent_*; must advance next_wakeup().
    virtual void on_wakeup(MarketSim& sim) = 0;
    // Fills where the agent is taker or maker (time in s).
    virtual void on_fill(const Fill& f, double t) {
        (void)f;
        (void)t;
    }
};

class MarketSim {
  public:
    MarketSim(const SimConfig& cfg, std::uint64_t seed, std::uint64_t path_id);
    // Registers a non-owning agent; returns its OwnerId (1 + index).
    OwnerId add_agent(Agent* a);
    // Processes all events with time <= t_end (seconds), then sets time = t_end.
    void run_until(double t_end);
    [[nodiscard]] double time() const noexcept { return t_; }
    [[nodiscard]] const OrderBook& book() const noexcept { return book_; }
    // (bid+ask)/2 in ticks if both sides exist, else the last valid mid (initially 10000.0 for the default config).
    [[nodiscard]] double mid() const noexcept { return mid_; }
    [[nodiscard]] std::optional<Tick> best_bid() const noexcept { return book_.best_bid(); }
    [[nodiscard]] std::optional<Tick> best_ask() const noexcept { return book_.best_ask(); }
    // Current Hawkes intensity (1/s) of type 0..3.
    [[nodiscard]] double hawkes_intensity(int type) const noexcept { return hawkes_.intensity(type); }
    [[nodiscard]] OrderResult agent_market(OwnerId owner, Side side, Qty qty);
    // GTC limit order (price in ticks, qty in lots).
    [[nodiscard]] OrderResult agent_limit(OwnerId owner, Side side, Tick price, Qty qty);
    // Cancels an order owned by `owner`; Rejected/UnknownId if unknown or owned by someone else.
    [[nodiscard]] OrderResult agent_cancel(OwnerId owner, OrderId id);
    struct Counters {
        std::array<std::int64_t, 4> hawkes_events{};  // accepted background Hawkes events per type
        std::int64_t cancels = 0, agent_orders = 0, fills = 0;
        std::int64_t bg_crossing_attempts = 0;  // background limit orders that traded on arrival (must stay 0)
        double bg_live_integral = 0.0;          // integral of the number of live background orders over time (s)
    };
    [[nodiscard]] const Counters& counters() const noexcept { return counters_; }
    [[nodiscard]] std::size_t bg_live_count() const noexcept { return bg_live_.size(); }
    // Optional observer called for every fill (all owners) in engine order.
    void set_fill_observer(std::function<void(const Fill&, double)> f) { observer_ = std::move(f); }

  private:
    static constexpr std::uint32_t kNil = 0xFFFFFFFFu;
    void process_fills();
    void after_engine_call();
    void bg_add(OrderId id);
    void bg_remove(OrderId id);
    void background_event(int type);

    SimConfig cfg_;
    OrderBook book_;
    HawkesState hawkes_;
    qc::RngStream rng_h_exp_, rng_h_acc_, rng_attr_, rng_c_exp_, rng_c_sel_;
    std::vector<Agent*> agents_;
    std::vector<OrderId> bg_live_;
    std::vector<std::uint32_t> bg_pos_;
    double t_ = 0.0;
    double mid_;
    Tick last_bid_, last_ask_;
    Counters counters_;
    std::function<void(const Fill&, double)> observer_;
};

}  // namespace ofl
