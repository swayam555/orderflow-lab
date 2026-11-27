#include "ofl/hawkes.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace ofl {

void HawkesParams::validate() const {
    if (dim < 1) throw std::invalid_argument("HawkesParams: dim = " + std::to_string(dim) + " must be >= 1");
    const auto d = static_cast<std::size_t>(dim);
    if (mu.size() != d) throw std::invalid_argument("HawkesParams: mu has wrong size");
    if (alpha.rows != dim || alpha.cols != dim) throw std::invalid_argument("HawkesParams: alpha has wrong shape");
    if (beta.rows != dim || beta.cols != dim) throw std::invalid_argument("HawkesParams: beta has wrong shape");
    for (int i = 0; i < dim; ++i) {
        if (!(mu[static_cast<std::size_t>(i)] > 0.0))
            throw std::invalid_argument("HawkesParams: mu[" + std::to_string(i) + "] = " +
                                        std::to_string(mu[static_cast<std::size_t>(i)]) + " must be > 0");
        for (int j = 0; j < dim; ++j) {
            if (!(alpha(i, j) >= 0.0))
                throw std::invalid_argument("HawkesParams: alpha(" + std::to_string(i) + "," + std::to_string(j) +
                                            ") = " + std::to_string(alpha(i, j)) + " must be >= 0");
            if (!(beta(i, j) > 0.0))
                throw std::invalid_argument("HawkesParams: beta(" + std::to_string(i) + "," + std::to_string(j) +
                                            ") = " + std::to_string(beta(i, j)) + " must be > 0");
        }
    }
    const double rho = qc::spectral_radius_nonneg(branching_matrix(*this));
    if (!(rho < 1.0))
        throw std::invalid_argument("HawkesParams: spectral radius of G = " + std::to_string(rho) + " must be < 1");
}

qc::Mat branching_matrix(const HawkesParams& p) {
    qc::Mat G(p.dim, p.dim);
    for (int i = 0; i < p.dim; ++i)
        for (int j = 0; j < p.dim; ++j) G(i, j) = p.alpha(i, j) / p.beta(i, j);
    return G;
}

std::vector<double> stationary_rates(const HawkesParams& p) {
    qc::Mat A = qc::Mat::identity(p.dim);
    const qc::Mat G = branching_matrix(p);
    for (int i = 0; i < p.dim; ++i)
        for (int j = 0; j < p.dim; ++j) A(i, j) -= G(i, j);
    return qc::solve_lu(A, p.mu);
}

HawkesParams params_from_rates(const qc::Mat& G, const qc::Mat& beta, std::span<const double> lambda_bar) {
    const int d = G.rows;
    if (G.cols != d || beta.rows != d || beta.cols != d || static_cast<int>(lambda_bar.size()) != d)
        throw std::invalid_argument("params_from_rates: inconsistent dimensions");
    const double rho = qc::spectral_radius_nonneg(G);
    if (!(rho < 1.0))
        throw std::invalid_argument("params_from_rates: spectral radius " + std::to_string(rho) + " must be < 1");
    HawkesParams p;
    p.dim = d;
    p.mu.assign(static_cast<std::size_t>(d), 0.0);
    p.alpha = qc::Mat(d, d);
    p.beta = beta;
    for (int i = 0; i < d; ++i) {
        double m = lambda_bar[static_cast<std::size_t>(i)];
        for (int j = 0; j < d; ++j) {
            m -= G(i, j) * lambda_bar[static_cast<std::size_t>(j)];
            p.alpha(i, j) = G(i, j) * beta(i, j);
        }
        if (!(m > 0.0))
            throw std::invalid_argument("params_from_rates: mu[" + std::to_string(i) + "] = " + std::to_string(m) +
                                        " must be > 0");
        p.mu[static_cast<std::size_t>(i)] = m;
    }
    return p;
}

HawkesParams with_branching_ratio(const qc::Mat& G_shape, const qc::Mat& beta, std::span<const double> lambda_bar,
                                  double target) {
    if (!(target >= 0.0 && target < 1.0))
        throw std::invalid_argument("with_branching_ratio: target = " + std::to_string(target) + " must be in [0,1)");
    const double rho = qc::spectral_radius_nonneg(G_shape);
    qc::Mat G = G_shape;
    if (rho > 0.0)
        for (double& g : G.a) g *= target / rho;
    return params_from_rates(G, beta, lambda_bar);
}

HawkesState::HawkesState(const HawkesParams& p)
    : p_(p), S_(static_cast<std::size_t>(p.dim) * static_cast<std::size_t>(p.dim), 0.0) {}

void HawkesState::decay_to(double t) noexcept {
    const double dt = t - t_;
    if (dt > 0.0) {
        const int d = p_.dim;
        for (int i = 0; i < d; ++i)
            for (int j = 0; j < d; ++j) S_[static_cast<std::size_t>(i * d + j)] *= std::exp(-p_.beta(i, j) * dt);
    }
    t_ = t;
}

void HawkesState::add_event(int j) noexcept {
    const int d = p_.dim;
    for (int i = 0; i < d; ++i) S_[static_cast<std::size_t>(i * d + j)] += p_.alpha(i, j);
}

double HawkesState::intensity(int i) const noexcept {
    double s = p_.mu[static_cast<std::size_t>(i)];
    for (int j = 0; j < p_.dim; ++j) s += S_[static_cast<std::size_t>(i * p_.dim + j)];
    return s;
}

double HawkesState::total_intensity() const noexcept {
    double s = 0.0;
    for (int i = 0; i < p_.dim; ++i) s += intensity(i);
    return s;
}

double HawkesState::compensator(int i, double dt) const noexcept {
    double s = p_.mu[static_cast<std::size_t>(i)] * dt;
    for (int j = 0; j < p_.dim; ++j)
        s += S_[static_cast<std::size_t>(i * p_.dim + j)] * (-std::expm1(-p_.beta(i, j) * dt)) / p_.beta(i, j);
    return s;
}

std::vector<HawkesEvent> simulate_hawkes(const HawkesParams& p, double T, qc::RngStream& rng) {
    p.validate();
    HawkesState st(p);
    std::vector<HawkesEvent> ev;
    while (true) {
        const double lbar = st.total_intensity();
        const double t_new = st.time() + rng.exponential() / lbar;
        if (t_new > T) break;
        st.decay_to(t_new);
        const double lam = st.total_intensity();
        if (rng.uniform() * lbar <= lam) {
            const double u = rng.uniform() * lam;
            double c = 0.0;
            int idx = p.dim - 1;
            for (int i = 0; i < p.dim; ++i) {
                c += st.intensity(i);
                if (c >= u) {
                    idx = i;
                    break;
                }
            }
            ev.push_back({t_new, idx});
            st.add_event(idx);
        }
    }
    return ev;
}

double hawkes_loglik(const HawkesParams& p, std::span<const HawkesEvent> events, double T) {
    HawkesState st(p);
    double ll = 0.0;
    for (const HawkesEvent& e : events) {
        const double dt = e.t - st.time();
        for (int i = 0; i < p.dim; ++i) ll -= st.compensator(i, dt);
        st.decay_to(e.t);
        ll += std::log(st.intensity(e.type));
        st.add_event(e.type);
    }
    const double dt = T - st.time();
    for (int i = 0; i < p.dim; ++i) ll -= st.compensator(i, dt);
    return ll;
}

std::vector<std::vector<double>> rescaled_interarrivals(const HawkesParams& p, std::span<const HawkesEvent> events,
                                                        double /*T*/) {
    HawkesState st(p);
    std::vector<double> acc(static_cast<std::size_t>(p.dim), 0.0);
    std::vector<std::vector<double>> out(static_cast<std::size_t>(p.dim));
    for (const HawkesEvent& e : events) {
        const double dt = e.t - st.time();
        for (int i = 0; i < p.dim; ++i) acc[static_cast<std::size_t>(i)] += st.compensator(i, dt);
        st.decay_to(e.t);
        out[static_cast<std::size_t>(e.type)].push_back(acc[static_cast<std::size_t>(e.type)]);
        acc[static_cast<std::size_t>(e.type)] = 0.0;
        st.add_event(e.type);
    }
    return out;
}

EmStep em_step(std::span<const HawkesEvent> events, double T, const HawkesParams& p) {
    const int d = p.dim;
    const auto ud = static_cast<std::size_t>(d);
    bool uniform_beta = true;
    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j)
            if (p.beta(i, j) != p.beta(0, 0)) uniform_beta = false;
    std::vector<double> A(ud * ud, 0.0), num_mu(ud, 0.0), num_a(ud * ud, 0.0), C(ud * ud, 0.0);
    double ll = 0.0, t_prev = 0.0;
    for (const HawkesEvent& e : events) {
        const double dt = e.t - t_prev;
        if (uniform_beta) {
            const double f = std::exp(-p.beta(0, 0) * dt);
            for (double& a : A) a *= f;
        } else {
            for (int i = 0; i < d; ++i)
                for (int j = 0; j < d; ++j) A[static_cast<std::size_t>(i * d + j)] *= std::exp(-p.beta(i, j) * dt);
        }
        t_prev = e.t;
        const auto i = static_cast<std::size_t>(e.type);
        double lam = p.mu[i];
        for (std::size_t j = 0; j < ud; ++j) lam += p.alpha(e.type, static_cast<int>(j)) * A[i * ud + j];
        ll += std::log(lam);
        num_mu[i] += p.mu[i] / lam;
        for (std::size_t j = 0; j < ud; ++j) num_a[i * ud + j] += p.alpha(e.type, static_cast<int>(j)) * A[i * ud + j] / lam;
        for (std::size_t r = 0; r < ud; ++r) A[r * ud + i] += 1.0;
        // C_rj contribution of this event (type i = source j): (1-exp(-beta_ri (T-t)))/beta_ri
        for (int r = 0; r < d; ++r)
            C[static_cast<std::size_t>(r) * ud + i] +=
                -std::expm1(-p.beta(r, e.type) * (T - e.t)) / p.beta(r, e.type);
    }
    EmStep out;
    out.params = p;
    double comp = 0.0;
    for (std::size_t i = 0; i < ud; ++i) {
        comp += p.mu[i] * T;
        for (std::size_t j = 0; j < ud; ++j) comp += p.alpha(static_cast<int>(i), static_cast<int>(j)) * C[i * ud + j];
    }
    out.loglik = ll - comp;
    for (std::size_t i = 0; i < ud; ++i) {
        out.params.mu[i] = num_mu[i] / T;
        for (std::size_t j = 0; j < ud; ++j)
            out.params.alpha(static_cast<int>(i), static_cast<int>(j)) = C[i * ud + j] > 0.0 ? num_a[i * ud + j] / C[i * ud + j] : 0.0;
    }
    return out;
}

HawkesParams default_em_init(std::span<const HawkesEvent> events, double T, const qc::Mat& beta) {
    HawkesParams p;
    p.dim = beta.rows;
    p.beta = beta;
    p.alpha = qc::Mat(p.dim, p.dim);
    std::vector<double> n(static_cast<std::size_t>(p.dim), 0.0);
    for (const HawkesEvent& e : events) n[static_cast<std::size_t>(e.type)] += 1.0;
    p.mu.resize(n.size());
    for (int i = 0; i < p.dim; ++i) {
        p.mu[static_cast<std::size_t>(i)] = 0.5 * n[static_cast<std::size_t>(i)] / T;
        for (int j = 0; j < p.dim; ++j) p.alpha(i, j) = 0.1 * beta(i, j);
    }
    return p;
}

EmResult fit_em(std::span<const HawkesEvent> events, double T, const qc::Mat& beta, const HawkesParams& init,
                const EmOptions& opt) {
    const int d = beta.rows;
    if (init.dim != d) throw std::invalid_argument("fit_em: init dimension != beta dimension");
    if (!(T > 0.0)) throw std::invalid_argument("fit_em: T must be > 0");
    std::vector<int> cnt(static_cast<std::size_t>(d), 0);
    double prev = 0.0;
    for (const HawkesEvent& e : events) {
        if (e.type < 0 || e.type >= d) throw std::invalid_argument("fit_em: event type out of range");
        if (!(e.t > 0.0 && e.t <= T)) throw std::invalid_argument("fit_em: event time outside (0,T]");
        if (e.t < prev) throw std::invalid_argument("fit_em: events not sorted by time");
        prev = e.t;
        ++cnt[static_cast<std::size_t>(e.type)];
    }
    for (int i = 0; i < d; ++i)
        if (cnt[static_cast<std::size_t>(i)] < 10)
            throw std::invalid_argument("fit_em: dimension " + std::to_string(i) + " has only " +
                                        std::to_string(cnt[static_cast<std::size_t>(i)]) + " events (< 10)");
    EmResult res;
    HawkesParams cur = init;
    cur.beta = beta;
    for (int it = 0; it < opt.max_iter; ++it) {
        EmStep s = em_step(events, T, cur);
        res.loglik_trace.push_back(s.loglik);
        res.iterations = it + 1;
        if (it > 0) {
            const double prevl = res.loglik_trace[res.loglik_trace.size() - 2];
            if (std::fabs(s.loglik - prevl) <= opt.rel_tol * std::fabs(s.loglik)) {
                res.converged = true;
                break;
            }
        }
        cur = std::move(s.params);
    }
    res.params = std::move(cur);
    return res;
}

}  // namespace ofl
