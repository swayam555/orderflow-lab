#include "ofl/market_sim.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace ofl {

qc::Mat default_G() {
    qc::Mat G(4, 4);
    G(0, 0) = 0.3;
    G(0, 3) = 0.3;
    G(1, 1) = 0.3;
    G(1, 2) = 0.3;
    G(2, 2) = 0.5;
    G(3, 3) = 0.5;
    return G;
}

void SimConfig::validate() const {
    if (G.rows != 4 || G.cols != 4) throw std::invalid_argument("SimConfig: G must be 4x4");
    if (!(beta > 0.0)) throw std::invalid_argument("SimConfig: beta = " + std::to_string(beta) + " must be > 0");
    if (!(cancel_hazard >= 0.0))
        throw std::invalid_argument("SimConfig: cancel_hazard = " + std::to_string(cancel_hazard) + " must be >= 0");
    if (placement_K < 1) throw std::invalid_argument("SimConfig: placement_K = " + std::to_string(placement_K) + " must be >= 1");
    if (size_min < 1 || size_max < size_min)
        throw std::invalid_argument("SimConfig: bad size range [" + std::to_string(size_min) + "," + std::to_string(size_max) + "]");
    if (max_tick <= min_tick) throw std::invalid_argument("SimConfig: max_tick must exceed min_tick");
    if (init_best_bid < min_tick + init_levels || init_best_ask + init_levels > max_tick || init_best_bid >= init_best_ask)
        throw std::invalid_argument("SimConfig: initial book does not fit in the band or is crossed");
    if (init_levels < 1 || init_orders_per_level < 1 || init_order_size < 1)
        throw std::invalid_argument("SimConfig: initial book sizes must be >= 1");
    for (double l : lambda_bar)
        if (!(l > 0.0)) throw std::invalid_argument("SimConfig: lambda_bar entries must be > 0");
}

namespace {
HawkesParams make_params(const SimConfig& c) {
    c.validate();
    qc::Mat beta(4, 4);
    for (double& b : beta.a) b = c.beta;
    return params_from_rates(c.G, beta, std::span<const double>(c.lambda_bar.data(), 4));
}
}  // namespace

MarketSim::MarketSim(const SimConfig& cfg, std::uint64_t seed, std::uint64_t path_id)
    : cfg_((cfg.validate(), cfg)),  // validate first so errors name the offending SimConfig field
      book_(cfg.min_tick, cfg.max_tick, 1u << 14, 1u << 16),
      hawkes_(make_params(cfg)),
      rng_h_exp_(seed, path_id, 1),
      rng_h_acc_(seed, path_id, 2),
      rng_attr_(seed, path_id, 3),
      rng_c_exp_(seed, path_id, 4),
      rng_c_sel_(seed, path_id, 5),
      rng_f_(seed, path_id, 6),
      mid_(0.5 * (static_cast<double>(cfg.init_best_bid) + static_cast<double>(cfg.init_best_ask))),
      last_bid_(cfg.init_best_bid),
      last_ask_(cfg.init_best_ask) {
    fund_ = mid_;
    bg_pos_.assign(1u << 12, kNil);
    for (int side = 0; side < 2; ++side)
        for (int l = 0; l < cfg.init_levels; ++l) {
            const Tick px = side == 0 ? cfg.init_best_bid - l : cfg.init_best_ask + l;
            for (int o = 0; o < cfg.init_orders_per_level; ++o) {
                const OrderResult r = book_.submit_limit(side == 0 ? Side::Buy : Side::Sell, px, cfg.init_order_size,
                                                         Tif::GTC, kBackground);
                bg_add(r.id);
            }
        }
}

OwnerId MarketSim::add_agent(Agent* a) {
    agents_.push_back(a);
    return static_cast<OwnerId>(agents_.size());
}

void MarketSim::bg_add(OrderId id) {
    if (id >= bg_pos_.size()) bg_pos_.resize(std::max<std::size_t>(bg_pos_.size() * 2, id + 1), kNil);
    bg_pos_[id] = static_cast<std::uint32_t>(bg_live_.size());
    bg_live_.push_back(id);
}

void MarketSim::bg_remove(OrderId id) {
    const std::uint32_t pos = bg_pos_[id];
    const OrderId last = bg_live_.back();
    bg_live_[pos] = last;
    bg_pos_[last] = pos;
    bg_live_.pop_back();
    bg_pos_[id] = kNil;
}

void MarketSim::process_fills() {
    for (const Fill& f : book_.last_fills()) {
        ++counters_.fills;
        if (f.maker_owner == kBackground && f.maker_done) bg_remove(f.maker);
        if (observer_) observer_(f, t_);
        if (f.taker_owner != kBackground) agents_[f.taker_owner - 1]->on_fill(f, t_);
        if (f.maker_owner != kBackground && f.maker_owner != f.taker_owner) agents_[f.maker_owner - 1]->on_fill(f, t_);
    }
}

void MarketSim::after_engine_call() {
    process_fills();
    const auto b = book_.best_bid();
    const auto a = book_.best_ask();
    if (b) last_bid_ = *b;
    if (a) last_ask_ = *a;
    if (b && a) mid_ = 0.5 * (static_cast<double>(*b) + static_cast<double>(*a));
}

void MarketSim::advance_fund() {
    const double dt = t_ - t_fund_;
    if (dt > 0.0) fund_ += cfg_.fund_sigma * std::sqrt(dt) * rng_f_.normal();
    t_fund_ = t_;
}

void MarketSim::background_event(int type) {
    ++counters_.hawkes_events[static_cast<std::size_t>(type)];
    const auto size_of = [&]() { return static_cast<Qty>(rng_attr_.uniform_int(cfg_.size_min, cfg_.size_max)); };
    if (type == 0 || type == 1) {
        const auto k = rng_attr_.uniform_int(1, cfg_.placement_K);
        const Qty sz = size_of();
        const Side side = type == 0 ? Side::Buy : Side::Sell;
        std::int64_t px;
        if (cfg_.fundamental_anchor) {
            advance_fund();
            if (type == 0) {
                std::int64_t p = static_cast<std::int64_t>(std::floor(fund_)) - k;
                if (book_.best_ask()) p = std::min<std::int64_t>(p, static_cast<std::int64_t>(*book_.best_ask()) - 1);
                px = std::max<std::int64_t>(p, cfg_.min_tick);
            } else {
                std::int64_t p = static_cast<std::int64_t>(std::ceil(fund_)) + k;
                if (book_.best_bid()) p = std::max<std::int64_t>(p, static_cast<std::int64_t>(*book_.best_bid()) + 1);
                px = std::min<std::int64_t>(p, cfg_.max_tick);
            }
        } else if (type == 0) {
            const Tick ref = book_.best_ask() ? *book_.best_ask() : last_ask_;
            px = std::max<std::int64_t>(static_cast<std::int64_t>(ref) - k, cfg_.min_tick);
        } else {
            const Tick ref = book_.best_bid() ? *book_.best_bid() : last_bid_;
            px = std::min<std::int64_t>(static_cast<std::int64_t>(ref) + k, cfg_.max_tick);
        }
        const OrderResult r = book_.submit_limit(side, static_cast<Tick>(px), sz, Tif::GTC, kBackground);
        if (r.filled > 0) ++counters_.bg_crossing_attempts;
        if (r.resting > 0) bg_add(r.id);
    } else {
        const Qty sz = size_of();
        (void)book_.submit_market(type == 2 ? Side::Buy : Side::Sell, sz, kBackground);
    }
    after_engine_call();
    hawkes_.add_event(type);
}

void MarketSim::run_until(double t_end) {
    if (t_end < t_) throw std::invalid_argument("run_until: t_end = " + std::to_string(t_end) + " < time()");
    constexpr double kInf = std::numeric_limits<double>::infinity();
    while (true) {
        const double lbar = hawkes_.total_intensity();
        const double t_h = t_ + rng_h_exp_.exponential() / lbar;
        const std::size_t n_c = bg_live_.size();
        const double t_c = (n_c > 0 && cfg_.cancel_hazard > 0.0)
                               ? t_ + rng_c_exp_.exponential() / (cfg_.cancel_hazard * static_cast<double>(n_c))
                               : kInf;
        double t_a = kInf;
        for (const Agent* a : agents_) t_a = std::min(t_a, std::max(a->next_wakeup(), t_));
        const double t_n = std::min({t_h, t_c, t_a});
        if (t_n > t_end) {
            counters_.bg_live_integral += static_cast<double>(n_c) * (t_end - t_);
            hawkes_.decay_to(t_end);
            t_ = t_end;
            return;
        }
        counters_.bg_live_integral += static_cast<double>(n_c) * (t_n - t_);
        hawkes_.decay_to(t_n);
        t_ = t_n;
        if (t_n == t_a) {
            for (Agent* a : agents_)
                if (std::max(a->next_wakeup(), t_) <= t_n) a->on_wakeup(*this);
        } else if (t_n == t_c) {
            const auto idx = static_cast<std::size_t>(rng_c_sel_.uniform_int(0, static_cast<std::int64_t>(n_c) - 1));
            const OrderId id = bg_live_[idx];
            (void)book_.cancel(id);
            bg_remove(id);
            ++counters_.cancels;
            after_engine_call();
        } else {
            const double u = rng_h_acc_.uniform();
            const double lam = hawkes_.total_intensity();
            if (u * lbar <= lam) {
                const double v = rng_h_acc_.uniform() * lam;
                double c = 0.0;
                int type = 3;
                for (int i = 0; i < 4; ++i) {
                    c += hawkes_.intensity(i);
                    if (c >= v) {
                        type = i;
                        break;
                    }
                }
                background_event(type);
            }
        }
    }
}

void MarketSim::check_owner(OwnerId owner) const {
    if (owner < 1 || owner > agents_.size())
        throw std::invalid_argument("MarketSim: owner id " + std::to_string(owner) + " is not a registered agent (bind() the agent first)");
}

OrderResult MarketSim::agent_market(OwnerId owner, Side side, Qty qty) {
    check_owner(owner);
    ++counters_.agent_orders;
    const OrderResult r = book_.submit_market(side, qty, owner);
    after_engine_call();
    if (cfg_.fundamental_anchor && cfg_.fund_feedback != 0.0) {
        advance_fund();
        fund_ += (side == Side::Buy ? 1.0 : -1.0) * cfg_.fund_feedback * static_cast<double>(r.filled);
    }
    if (cfg_.agent_market_orders_excite && r.status != Status::Rejected) hawkes_.add_event(side == Side::Buy ? 2 : 3);
    return r;
}

OrderResult MarketSim::agent_limit(OwnerId owner, Side side, Tick price, Qty qty) {
    check_owner(owner);
    ++counters_.agent_orders;
    const OrderResult r = book_.submit_limit(side, price, qty, Tif::GTC, owner);
    after_engine_call();
    return r;
}

OrderResult MarketSim::agent_cancel(OwnerId owner, OrderId id) {
    check_owner(owner);
    ++counters_.agent_orders;
    const auto v = book_.find(id);
    if (!v || v->owner != owner) return OrderResult{id, Status::Rejected, Reject::UnknownId, 0, 0, 0};
    const OrderResult r = book_.cancel(id);
    after_engine_call();
    return r;
}

}  // namespace ofl
