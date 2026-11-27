// Agents for the market simulator.
#pragma once
#include <limits>
#include <vector>

#include "ofl/market_sim.hpp"

namespace ofl {

// Records the mid price (ticks) at t0 + j*dt, j = 0..n-1.
class ObserverAgent final : public Agent {
  public:
    ObserverAgent(double t0, double dt, int n) : t0_(t0), dt_(dt), n_(n) {}
    [[nodiscard]] double next_wakeup() const override {
        return j_ < n_ ? t0_ + static_cast<double>(j_) * dt_ : std::numeric_limits<double>::infinity();
    }
    void on_wakeup(MarketSim& sim) override {
        mids_.push_back(sim.mid());
        ++j_;
    }
    [[nodiscard]] const std::vector<double>& mids() const noexcept { return mids_; }

  private:
    double t0_, dt_;
    int n_, j_ = 0;
    std::vector<double> mids_;
};

// One scheduled market child order: time in seconds, qty in lots.
struct ScheduledOrder {
    double time;
    Qty qty;
};

// Executes a parent order with scheduled market orders. A child that is not fully filled carries over to the next
// child; at final_time (if finite) any carried remainder is sent as one last market order. Records its fills.
// RNG: uses none.
class ScheduleAgent final : public Agent {
  public:
    // ref_time: time (s) at which the arrival mid is recorded (inf = at the first child time).
    ScheduleAgent(Side side, std::vector<ScheduledOrder> children, double final_time, double ref_time);
    // Registers with the sim and remembers the OwnerId.
    OwnerId bind(MarketSim& sim);
    [[nodiscard]] double next_wakeup() const override;
    void on_wakeup(MarketSim& sim) override;
    void on_fill(const Fill& f, double t) override;
    [[nodiscard]] Qty filled() const noexcept { return filled_; }
    [[nodiscard]] Qty unfilled() const noexcept { return carry_ + unsent_; }
    [[nodiscard]] double arrival_mid() const noexcept { return arrival_mid_; }
    // Notional of fills (ticks*lots).
    [[nodiscard]] double notional() const noexcept { return notional_; }
    // Implementation shortfall (ticks*lots): signed so that a worse price is positive; unfilled lots excluded.
    [[nodiscard]] double implementation_shortfall() const noexcept;
    [[nodiscard]] const std::vector<Fill>& fills() const noexcept { return fills_; }

  private:
    Side side_;
    std::vector<ScheduledOrder> children_;
    double final_time_, ref_time_;
    std::size_t idx_ = 0;
    bool ref_done_ = false, final_done_ = false;
    Qty carry_ = 0, unsent_ = 0, filled_ = 0;
    double notional_ = 0.0, arrival_mid_ = 0.0;
    OwnerId owner_ = 0;
    std::vector<Fill> fills_;
};

// Sends alternating buy/sell market orders at t0 + period*m, sizes 1 + (m mod 8) lots, m = 0..n-1, and records the
// per-lot slippage side*(VWAP - mid_before) (ticks).
class ProbeAgent final : public Agent {
  public:
    ProbeAgent(double t0, double period, int n) : t0_(t0), period_(period), n_(n) {}
    OwnerId bind(MarketSim& sim) {
        owner_ = sim.add_agent(this);
        return owner_;
    }
    [[nodiscard]] double next_wakeup() const override {
        return m_ < n_ ? t0_ + period_ * m_ : std::numeric_limits<double>::infinity();
    }
    void on_wakeup(MarketSim& sim) override;
    void on_fill(const Fill& f, double) override {
        qty_ += f.qty;
        notional_ += static_cast<double>(f.qty) * f.price;
    }
    struct Obs {
        double size, slippage;
    };
    [[nodiscard]] const std::vector<Obs>& observations() const noexcept { return obs_; }

  private:
    double t0_, period_;
    int n_, m_ = 0;
    OwnerId owner_ = 0;
    Qty qty_ = 0;
    double notional_ = 0.0;
    std::vector<Obs> obs_;
};

}  // namespace ofl
