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

}  // namespace ofl
