#include "ofl/agents.hpp"

#include <algorithm>
#include <cmath>

#include "ofl/avellaneda_stoikov.hpp"

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

void AsAgent::on_wakeup(MarketSim& sim) {
    if (bid_id_ != kNoOrder) (void)sim.agent_cancel(owner_, bid_id_);
    if (ask_id_ != kNoOrder) (void)sim.agent_cancel(owner_, ask_id_);
    bid_id_ = ask_id_ = kNoOrder;
    abs_q_sum_ += std::fabs(static_cast<double>(q_));
    const double s = sim.mid();
    const double t = sim.time() - c_.t_start, T = c_.session;
    const double sigma = std::sqrt(c_.sigma2);
    double r, spread;
    if (c_.symmetric) {
        r = s;
        spread = as_symmetric_spread(c_.gamma, sigma, c_.k, T);
    } else {
        const AsQuote q = as_quotes(s, static_cast<double>(q_), t, T, c_.gamma, sigma, c_.k);
        r = 0.5 * (q.bid + q.ask);
        spread = q.ask - q.bid;
    }
    auto bid_px = static_cast<Tick>(std::floor(r - 0.5 * spread));
    auto ask_px = static_cast<Tick>(std::ceil(r + 0.5 * spread));
    if (const auto a = sim.best_ask()) bid_px = std::min<Tick>(bid_px, *a - 1);
    if (const auto b = sim.best_bid()) ask_px = std::max<Tick>(ask_px, *b + 1);
    if (q_ < c_.q_limit) {
        const OrderResult r1 = sim.agent_limit(owner_, Side::Buy, bid_px, 1);
        if (r1.resting > 0) bid_id_ = r1.id;
    }
    if (q_ > -c_.q_limit) {
        const OrderResult r2 = sim.agent_limit(owner_, Side::Sell, ask_px, 1);
        if (r2.resting > 0) ask_id_ = r2.id;
    }
    ++step_;
}

void AsAgent::on_fill(const Fill& f, double) {
    ++n_fills_;
    // The agent only rests limit orders, so it is the maker: taker_side Sell means the agent bought.
    if (f.taker_side == Side::Sell) {
        q_ += f.qty;
        cash_ -= static_cast<double>(f.qty) * f.price;
    } else {
        q_ -= f.qty;
        cash_ += static_cast<double>(f.qty) * f.price;
    }
}

void QuoteProbeAgent::on_wakeup(MarketSim& sim) {
    if (bid_id_ != kNoOrder) (void)sim.agent_cancel(owner_, bid_id_);
    if (ask_id_ != kNoOrder) (void)sim.agent_cancel(owner_, ask_id_);
    bid_id_ = ask_id_ = kNoOrder;
    const double m = sim.mid();
    const auto bid_px = static_cast<Tick>(std::floor(m - delta_));
    const auto ask_px = static_cast<Tick>(std::ceil(m + delta_));
    const OrderResult r1 = sim.agent_limit(owner_, Side::Buy, bid_px, 1);
    if (r1.resting > 0) bid_id_ = r1.id;
    const OrderResult r2 = sim.agent_limit(owner_, Side::Sell, ask_px, 1);
    if (r2.resting > 0) ask_id_ = r2.id;
    ++step_;
}

}  // namespace ofl
