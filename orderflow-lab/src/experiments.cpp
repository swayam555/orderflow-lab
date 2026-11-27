#include "ofl/experiments.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include <qc/parallel.hpp>

#include "ofl/agents.hpp"
#include "ofl/bertsimas_lo.hpp"

namespace ofl {

// ------------------------------------------------------------------ execution

AcParams ac_params_from(const ExecConfig& ec, const AcCalibration& cal, double kappaT) {
    AcParams p{static_cast<double>(ec.X), ec.H, ec.N, std::sqrt(cal.sigma2), cal.gamma, cal.eta, cal.epsilon, 0.0};
    if (kappaT > 0.0) p.lambda = ac_lambda_for_kappaT(p, kappaT);
    return p;
}

std::vector<Qty> exec_child_schedule(const ExecStrategy& s, const ExecConfig& ec, const AcCalibration& cal) {
    const auto n = static_cast<std::size_t>(ec.N);
    switch (s.kind) {
        case ExecKind::Twap: {
            std::vector<double> d(n, static_cast<double>(ec.X) / ec.N);
            return round_schedule(d, ec.X);
        }
        case ExecKind::Ac: {
            // The schedule depends on kappa and tau only; use a well-posed eta_tilde so it is defined even for
            // degenerate calibrations.
            AcParams p = ac_params_from(ec, cal, 0.0);
            if (!(p.sigma > 0.0)) p.sigma = 1.0;
            if (!(p.eta - 0.5 * p.gamma * (p.T / p.N) > 0.0)) {
                p.eta = 1.0;
                p.gamma = 0.0;
            }
            p.lambda = ac_lambda_for_kappaT(p, s.kappaT);
            return round_schedule(almgren_chriss(p).n, ec.X);
        }
        case ExecKind::BlGamma0: {
            BlParams bp{static_cast<double>(ec.X), ec.N, 0.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0};
            const BlSolution sol = bertsimas_lo_solve(bp);
            std::vector<double> d(n);
            double W = bp.S_bar;
            for (int t = 1; t <= ec.N; ++t) {
                d[static_cast<std::size_t>(t - 1)] = bl_policy(sol, t, W, 0.0);
                W -= d[static_cast<std::size_t>(t - 1)];
            }
            return round_schedule(d, ec.X);
        }
    }
    return {};
}

ExecPathResult run_exec_path(const SimConfig& cfg, const ExecConfig& ec, const std::vector<Qty>& children,
                             std::uint64_t seed, std::uint64_t path_id) {
    MarketSim sim(cfg, seed, path_id);
    std::vector<ScheduledOrder> sched;
    const double tau = ec.H / ec.N;
    for (std::size_t k = 0; k < children.size(); ++k)
        sched.push_back({ec.t0 + static_cast<double>(k) * tau, children[k]});
    ScheduleAgent ag(Side::Buy, std::move(sched), ec.t0 + ec.H, ec.t0);
    ag.bind(sim);
    sim.run_until(ec.t0 + ec.H);
    return {ag.implementation_shortfall() / static_cast<double>(ec.X), ag.filled(), ag.unfilled()};
}

AcCalibration calibrate_ac_params(const SimConfig& cfg, const ExecConfig& ec, std::uint64_t seed, int n_sigma_paths,
                                  int n_probe_paths, int n_gamma_paths, unsigned threads) {
    AcCalibration cal;
    const double tau = ec.H / ec.N;
    // sigma^2: agent-free paths, mid sampled every 10 s on [600,1200]
    {
        std::vector<std::vector<double>> mids(static_cast<std::size_t>(n_sigma_paths));
        qc::parallel_blocks(n_sigma_paths, threads, [&](std::int64_t p) {
            MarketSim sim(cfg, seed, static_cast<std::uint64_t>(p));
            ObserverAgent ob(600.0, 10.0, 61);
            sim.add_agent(&ob);
            sim.run_until(1200.0);
            mids[static_cast<std::size_t>(p)] = ob.mids();
        });
        qc::RunningStats d;
        for (const auto& m : mids)
            for (std::size_t k = 1; k < m.size(); ++k) d.add(m[k] - m[k - 1]);
        cal.sigma2 = d.variance() / 10.0;
    }
    // epsilon, eta: probe orders
    {
        std::vector<std::vector<ProbeAgent::Obs>> obs(static_cast<std::size_t>(n_probe_paths));
        qc::parallel_blocks(n_probe_paths, threads, [&](std::int64_t p) {
            MarketSim sim(cfg, seed + 1, static_cast<std::uint64_t>(p));
            ProbeAgent pr(600.0, 15.0, 40);
            pr.bind(sim);
            sim.run_until(600.0 + 15.0 * 40);
            obs[static_cast<std::size_t>(p)] = pr.observations();
        });
        std::vector<double> x, y;
        for (const auto& v : obs)
            for (const auto& o : v) {
                x.push_back(o.size);
                y.push_back(o.slippage);
            }
        const qc::OlsResult r = qc::ols(x, y);
        cal.epsilon = r.intercept;
        cal.se_epsilon = r.se_intercept;
        cal.se_c = r.se_slope;
        cal.eta = r.slope * tau;
        cal.n_probe_obs = r.n;
    }
    // gamma: TWAP program, mid(t0 + 2H) - mid(t0)
    {
        std::vector<double> dm(static_cast<std::size_t>(n_gamma_paths));
        const std::vector<Qty> twap = exec_child_schedule({"TWAP", ExecKind::Twap, 0.0}, ec, cal);
        qc::parallel_blocks(n_gamma_paths, threads, [&](std::int64_t p) {
            MarketSim sim(cfg, seed + 2, static_cast<std::uint64_t>(p));
            std::vector<ScheduledOrder> sched;
            for (std::size_t k = 0; k < twap.size(); ++k) sched.push_back({ec.t0 + static_cast<double>(k) * tau, twap[k]});
            ScheduleAgent ag(Side::Buy, std::move(sched), ec.t0 + ec.H, ec.t0);
            ag.bind(sim);
            ObserverAgent ob(ec.t0, 2.0 * ec.H, 2);
            sim.add_agent(&ob);
            sim.run_until(ec.t0 + 2.0 * ec.H);
            dm[static_cast<std::size_t>(p)] = ob.mids()[1] - ob.mids()[0];
        });
        qc::RunningStats rs;
        for (double v : dm) rs.add(v);
        cal.gamma = rs.mean() / static_cast<double>(ec.X);
        cal.se_gamma = rs.se() / static_cast<double>(ec.X);
    }
    const double eta_t = cal.eta - 0.5 * cal.gamma * tau;
    if (!(eta_t > 0.0) || cal.epsilon < 0.0 || !std::isfinite(cal.sigma2) || !(cal.sigma2 > 0.0)) {
        std::ostringstream os;
        os << "calibrate_ac_params: unusable estimates: sigma2=" << cal.sigma2 << " epsilon=" << cal.epsilon
           << " eta=" << cal.eta << " gamma=" << cal.gamma << " eta_tilde=" << eta_t;
        throw std::runtime_error(os.str());
    }
    return cal;
}

std::vector<ExecStrategyResult> run_execution_experiment(const SimConfig& cfg, const ExecConfig& ec,
                                                         const AcCalibration& cal,
                                                         const std::vector<ExecStrategy>& strategies, std::int64_t R,
                                                         std::uint64_t seed, unsigned threads) {
    if (strategies.empty() || strategies[0].kind != ExecKind::Twap)
        throw std::invalid_argument("run_execution_experiment: strategies[0] must be TWAP");
    const std::size_t S = strategies.size();
    std::vector<ExecStrategyResult> out(S);
    for (std::size_t s = 0; s < S; ++s) {
        out[s].name = strategies[s].name;
        out[s].children = exec_child_schedule(strategies[s], ec, cal);
        const AcParams p = ac_params_from(ec, cal, strategies[s].kind == ExecKind::Ac ? strategies[s].kappaT : 0.0);
        const AcSchedule ac = almgren_chriss(p);
        out[s].model_E_per_lot = ac.expected_cost / static_cast<double>(ec.X);
        out[s].model_sd_per_lot = std::sqrt(ac.variance) / static_cast<double>(ec.X);
    }
    std::vector<ExecPathResult> res(static_cast<std::size_t>(R) * S);
    qc::parallel_blocks(R, threads, [&](std::int64_t p) {
        for (std::size_t s = 0; s < S; ++s)
            res[static_cast<std::size_t>(p) * S + s] =
                run_exec_path(cfg, ec, out[s].children, seed, static_cast<std::uint64_t>(p));
    });
    for (std::int64_t p = 0; p < R; ++p)
        for (std::size_t s = 0; s < S; ++s) {
            const ExecPathResult& r = res[static_cast<std::size_t>(p) * S + s];
            out[s].is.add(r.is_per_lot);
            out[s].diff.add(r.is_per_lot - res[static_cast<std::size_t>(p) * S].is_per_lot);
            out[s].unfilled_total += r.unfilled;
        }
    return out;
}

}  // namespace ofl
