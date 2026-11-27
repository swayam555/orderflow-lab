#include "ofl/experiments.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include <qc/csv.hpp>
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

double estimate_sigma2(const SimConfig& cfg, std::uint64_t seed, int n_paths, unsigned threads) {
    std::vector<std::vector<double>> mids(static_cast<std::size_t>(n_paths));
    qc::parallel_blocks(n_paths, threads, [&](std::int64_t p) {
        MarketSim sim(cfg, seed, static_cast<std::uint64_t>(p));
        ObserverAgent ob(600.0, 10.0, 61);
        sim.add_agent(&ob);
        sim.run_until(1200.0);
        mids[static_cast<std::size_t>(p)] = ob.mids();
    });
    qc::RunningStats d;
    for (const auto& m : mids)
        for (std::size_t k = 1; k < m.size(); ++k) d.add(m[k] - m[k - 1]);
    return d.variance() / 10.0;
}

AcCalibration calibrate_ac_params(const SimConfig& cfg, const ExecConfig& ec, std::uint64_t seed, int n_sigma_paths,
                                  int n_probe_paths, int n_gamma_paths, unsigned threads) {
    AcCalibration cal;
    const double tau = ec.H / ec.N;
    cal.sigma2 = estimate_sigma2(cfg, seed, n_sigma_paths, threads);
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

// ------------------------------------------------------------------ Avellaneda-Stoikov in the simulator

AsCalibration calibrate_as_params(const SimConfig& cfg, std::uint64_t seed, unsigned threads) {
    AsCalibration cal;
    cal.sigma2 = estimate_sigma2(cfg, seed, 400, threads);
    constexpr int kPaths = 100;
    std::vector<double> x, y;
    for (int d = 1; d <= 6; ++d) {
        std::vector<std::int64_t> fills(kPaths, 0);
        qc::parallel_blocks(kPaths, threads, [&](std::int64_t p) {
            MarketSim sim(cfg, seed + 1, static_cast<std::uint64_t>(d * 1000 + p));
            QuoteProbeAgent pr(static_cast<double>(d), 600.0, 600.0, 1.0);
            pr.bind(sim);
            sim.run_until(1200.0);
            fills[static_cast<std::size_t>(p)] = pr.fills();
        });
        std::int64_t tot = 0;
        for (auto f : fills) tot += f;
        const double lam = static_cast<double>(tot) / (static_cast<double>(kPaths) * 2.0 * 600.0);
        cal.lambda_by_delta[d - 1] = lam;
        if (!(lam > 0.0)) throw std::runtime_error("calibrate_as_params: no fills at delta = " + std::to_string(d));
        x.push_back(d);
        y.push_back(std::log(lam));
    }
    const qc::OlsResult r = qc::ols(x, y);
    cal.A = std::exp(r.intercept);
    cal.k = -r.slope;
    cal.se_k = r.se_slope;
    if (!(cal.k > 0.0)) {
        std::ostringstream os;
        os << "calibrate_as_params: fitted k = " << cal.k << " <= 0 (A = " << cal.A << ")";
        throw std::runtime_error(os.str());
    }
    cal.gamma = 0.4 / (cal.sigma2 * 600.0);
    return cal;
}

AsSimResult run_as_sim_experiment(const SimConfig& cfg, const AsCalibration& cal, std::int64_t R, std::uint64_t seed,
                                  unsigned threads) {
    struct One {
        double pnl[2], q[2], absq[2], fills[2];
    };
    std::vector<One> res(static_cast<std::size_t>(R));
    qc::parallel_blocks(R, threads, [&](std::int64_t p) {
        for (int sym = 0; sym < 2; ++sym) {
            AsAgentConfig ac;
            ac.gamma = cal.gamma;
            ac.sigma2 = cal.sigma2;
            ac.k = cal.k;
            ac.symmetric = sym == 1;
            MarketSim sim(cfg, seed, static_cast<std::uint64_t>(p));
            AsAgent ag(ac);
            ag.bind(sim);
            sim.run_until(1200.0);
            One& o = res[static_cast<std::size_t>(p)];
            o.pnl[sym] = ag.pnl(sim.mid());
            o.q[sym] = static_cast<double>(ag.inventory());
            o.absq[sym] = ag.mean_abs_inventory();
            o.fills[sym] = static_cast<double>(ag.n_fills());
        }
    });
    AsSimResult out;
    for (const One& o : res) {
        out.as.pnl.add(o.pnl[0]);
        out.as.q_final.add(o.q[0]);
        out.as.abs_q.add(o.absq[0]);
        out.as.fills.add(o.fills[0]);
        out.symmetric.pnl.add(o.pnl[1]);
        out.symmetric.q_final.add(o.q[1]);
        out.symmetric.abs_q.add(o.absq[1]);
        out.symmetric.fills.add(o.fills[1]);
        out.pnl_diff.add(o.pnl[1] - o.pnl[0]);
    }
    return out;
}

// ------------------------------------------------------------------ impact validation

ImpactResult run_impact_experiment(const SimConfig& cfg, const ImpactConfig& ic, std::uint64_t seed, unsigned threads) {
    std::vector<int> qs{0};
    for (int q : ic.Qs) qs.push_back(q);
    const std::size_t nq = qs.size();
    const auto np = static_cast<std::size_t>(ic.profile_points);
    const std::int64_t R = ic.R;
    std::vector<double> mids(static_cast<std::size_t>(R) * nq * np);
    std::vector<double> volume(static_cast<std::size_t>(R), 0.0);
    qc::parallel_blocks(R, threads, [&](std::int64_t p) {
        for (std::size_t qi = 0; qi < nq; ++qi) {
            MarketSim sim(cfg, seed, static_cast<std::uint64_t>(p));
            ObserverAgent ob(ic.t0, ic.H / 10.0, ic.profile_points);
            sim.add_agent(&ob);
            std::vector<ScheduledOrder> sched;
            const int Q = qs[qi];
            for (int m = 0; m < Q; ++m) sched.push_back({ic.t0 + (m + 0.5) * ic.H / Q, 1});
            ScheduleAgent meta(Side::Buy, std::move(sched), std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity());
            if (Q > 0) meta.bind(sim);
            double vol = 0.0;
            if (qi == 0)
                sim.set_fill_observer([&](const Fill& f, double t) {
                    if (f.taker_owner == kBackground && t >= ic.t0 && t <= ic.t0 + ic.H) vol += static_cast<double>(f.qty);
                });
            sim.run_until(ic.t0 + ic.H / 10.0 * (ic.profile_points - 1));
            if (qi == 0) volume[static_cast<std::size_t>(p)] = vol;
            const auto& m = ob.mids();
            for (std::size_t j = 0; j < np; ++j) mids[(static_cast<std::size_t>(p) * nq + qi) * np + j] = m[j];
        }
    });
    auto mid_at = [&](std::int64_t p, std::size_t qi, std::size_t j) { return mids[(static_cast<std::size_t>(p) * nq + qi) * np + j]; };
    const std::size_t jH = 10;  // t0 + H
    ImpactResult out;
    {
        qc::RunningStats sh, vh;
        for (std::int64_t p = 0; p < R; ++p) {
            sh.add(mid_at(p, 0, jH) - mid_at(p, 0, 0));
            vh.add(volume[static_cast<std::size_t>(p)]);
        }
        out.sigma_H = sh.stddev();
        out.V_H = vh.mean();
    }
    for (std::size_t qi = 0; qi < nq; ++qi) {
        qc::RunningStats imp, step;
        for (std::int64_t p = 0; p < R; ++p) {
            const double d = (mid_at(p, qi, jH) - mid_at(p, qi, 0)) - (mid_at(p, 0, jH) - mid_at(p, 0, 0));
            imp.add(d);
            if (qi > 0) {
                const double dprev = (mid_at(p, qi - 1, jH) - mid_at(p, qi - 1, 0)) - (mid_at(p, 0, jH) - mid_at(p, 0, 0));
                step.add(d - dprev);
            } else {
                step.add(0.0);
            }
        }
        out.points.push_back({qs[qi], imp.mean(), imp.se(), imp.n(), step.mean(), step.se()});
        for (std::size_t j = 0; j < np; ++j) {
            qc::RunningStats pr;
            for (std::int64_t p = 0; p < R; ++p)
                pr.add((mid_at(p, qi, j) - mid_at(p, qi, 0)) - (mid_at(p, 0, j) - mid_at(p, 0, 0)));
            out.profile.push_back({qs[qi], static_cast<int>(j), ic.t0 + static_cast<double>(j) * ic.H / 10.0, pr.mean(), pr.se()});
        }
    }
    // fit: log I on log Q over points with I > 4 SE
    std::vector<double> lx, ly;
    double ysum = 0.0;
    for (const ImpactPoint& pt : out.points) {
        if (pt.Q <= 0 || !(pt.I > 4.0 * pt.se)) continue;
        lx.push_back(std::log(static_cast<double>(pt.Q)));
        ly.push_back(std::log(pt.I));
        ysum += pt.I / (out.sigma_H * std::sqrt(static_cast<double>(pt.Q) / out.V_H));
    }
    if (lx.size() >= 4) {
        const qc::OlsResult f = qc::ols(lx, ly);
        out.fit = {true, f.slope, f.se_slope, f.r2, ysum / static_cast<double>(lx.size()), static_cast<int>(lx.size())};
    } else {
        out.fit.n_points = static_cast<int>(lx.size());
    }
    return out;
}

void write_impact_csvs(const ImpactResult& r, const std::string& dir) {
    {
        qc::CsvWriter c(dir + "/impact.csv", {"Q", "I", "se", "n"});
        for (const auto& p : r.points) c.row({std::int64_t{p.Q}, p.I, p.se, p.n});
    }
    {
        qc::CsvWriter c(dir + "/impact_profile.csv", {"Q", "j", "time", "mean", "se"});
        for (const auto& p : r.profile) c.row({std::int64_t{p.Q}, std::int64_t{p.j}, p.time, p.mean, p.se});
    }
    {
        qc::CsvWriter c(dir + "/impact_fit.csv", {"delta", "se_delta", "r2", "Y", "n_points"});
        c.row({r.fit.delta, r.fit.se_delta, r.fit.r2, r.fit.Y, std::int64_t{r.fit.n_points}});
    }
}

}  // namespace ofl
