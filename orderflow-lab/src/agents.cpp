#include "ofl/agents.hpp"

#include <algorithm>

namespace ofl {

ScheduleAgent::ScheduleAgent(Side side, std::vector<ScheduledOrder> children, double final_time, double ref_time)
    : side_(side), children_(std::move(children)), final_time_(final_time), ref_time_(ref_time) {
    for (const auto& c : children_) unsent_ += c.qty;
}

OwnerId ScheduleAgent::bind(MarketSim& sim) {
    owner_ = sim.add_agent(this);
    return owner_;
}

double ScheduleAgent::next_wakeup() const {
    constexpr double kInf = std::numeric_limits<double>::infinity();
    double t = kInf;
    if (!ref_done_ && ref_time_ < kInf) t = std::min(t, ref_time_);
    if (idx_ < children_.size()) t = std::min(t, children_[idx_].time);
    if (!final_done_ && final_time_ < kInf) t = std::min(t, final_time_);
    return t;
}

void ScheduleAgent::on_wakeup(MarketSim& sim) {
    const double now = sim.time();
    if (!ref_done_ && (ref_time_ <= now || (ref_time_ == std::numeric_limits<double>::infinity() && idx_ < children_.size() && children_[idx_].time <= now))) {
        arrival_mid_ = sim.mid();
        ref_done_ = true;
    }
    while (idx_ < children_.size() && children_[idx_].time <= now) {
        const Qty want = children_[idx_].qty + carry_;
        unsent_ -= children_[idx_].qty;
        ++idx_;
        if (want > 0) {
            const OrderResult r = sim.agent_market(owner_, side_, want);
            carry_ = r.status == Status::Rejected ? want : r.cancelled;
        }
    }
    if (!final_done_ && final_time_ <= now) {
        final_done_ = true;
        if (carry_ > 0) {
            const OrderResult r = sim.agent_market(owner_, side_, carry_);
            carry_ = r.status == Status::Rejected ? carry_ : r.cancelled;
        }
    }
}

void ScheduleAgent::on_fill(const Fill& f, double) {
    if (f.taker_owner != owner_) return;
    filled_ += f.qty;
    notional_ += static_cast<double>(f.qty) * f.price;
    fills_.push_back(f);
}

double ScheduleAgent::implementation_shortfall() const noexcept {
    const double sgn = side_ == Side::Buy ? 1.0 : -1.0;
    return sgn * (notional_ - static_cast<double>(filled_) * arrival_mid_);
}

void ProbeAgent::on_wakeup(MarketSim& sim) {
    const bool buy = (m_ % 2) == 0;
    const Qty v = 1 + (m_ % 8);
    const double mid_before = sim.mid();
    qty_ = 0;
    notional_ = 0.0;
    (void)sim.agent_market(owner_, buy ? Side::Buy : Side::Sell, v);
    if (qty_ > 0) {
        const double vwap = notional_ / static_cast<double>(qty_);
        obs_.push_back({static_cast<double>(v), (buy ? 1.0 : -1.0) * (vwap - mid_before)});
    }
    ++m_;
}

}  // namespace ofl
