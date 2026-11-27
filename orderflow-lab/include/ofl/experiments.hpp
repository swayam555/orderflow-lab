// Experiments built on MarketSim: execution (AC/BL vs TWAP), Avellaneda-Stoikov in the simulator, impact validation.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include <qc/stats.hpp>

#include "ofl/almgren_chriss.hpp"
#include "ofl/market_sim.hpp"

namespace ofl {

// ---------------------------------------------------------------- execution (DESIGN 6.3)
struct ExecConfig {
    Qty X = 64;       // lots to buy
    double H = 300;   // horizon, s
    double t0 = 600;  // start time, s
    int N = 10;       // number of intervals
};
// Estimated Almgren-Chriss inputs (ticks, lots, seconds).
struct AcCalibration {
    double sigma2 = 0;    // ticks^2 / s
    double epsilon = 0;   // ticks
    double eta = 0;       // ticks per (lot/s)
    double gamma = 0;     // ticks per lot
    double se_epsilon = 0, se_c = 0, se_gamma = 0;
    std::int64_t n_probe_obs = 0;
};
// Estimates AC parameters from simulated paths (seed 7001 in the spec). Throws std::runtime_error (with the estimates
// in the message) if eta - gamma*tau/2 <= 0 or epsilon < 0.
[[nodiscard]] AcCalibration calibrate_ac_params(const SimConfig& cfg, const ExecConfig& ec, std::uint64_t seed,
                                                int n_sigma_paths, int n_probe_paths, int n_gamma_paths, unsigned threads);

enum class ExecKind { Twap, Ac, BlGamma0 };
struct ExecStrategy {
    std::string name;
    ExecKind kind;
    double kappaT = 0.0;  // Ac only
};
// AC parameter set (lambda chosen for the given kappa*T) from a calibration.
[[nodiscard]] AcParams ac_params_from(const ExecConfig& ec, const AcCalibration& cal, double kappaT);
// Integer child orders (lots) of a strategy.
[[nodiscard]] std::vector<Qty> exec_child_schedule(const ExecStrategy& s, const ExecConfig& ec, const AcCalibration& cal);

struct ExecPathResult {
    double is_per_lot;  // ticks per lot
    Qty filled, unfilled;
};
// Runs one buy program (children at t0 + (k-1) H/N, final sweep at t0+H) on path `path_id`.
[[nodiscard]] ExecPathResult run_exec_path(const SimConfig& cfg, const ExecConfig& ec, const std::vector<Qty>& children,
                                           std::uint64_t seed, std::uint64_t path_id);
struct ExecStrategyResult {
    std::string name;
    std::vector<Qty> children;
    qc::RunningStats is;    // IS per lot
    qc::RunningStats diff;  // IS per lot minus TWAP's, same path
    std::int64_t unfilled_total = 0;
    double model_E_per_lot = 0, model_sd_per_lot = 0;  // AC model prediction (ticks per lot)
};
// Evaluates all strategies on paths 0..R-1 with common random numbers; strategies[0] must be TWAP.
[[nodiscard]] std::vector<ExecStrategyResult> run_execution_experiment(const SimConfig& cfg, const ExecConfig& ec,
                                                                       const AcCalibration& cal,
                                                                       const std::vector<ExecStrategy>& strategies,
                                                                       std::int64_t R, std::uint64_t seed, unsigned threads);

// ---------------------------------------------------------------- Avellaneda-Stoikov in the simulator (DESIGN 7.2)
struct AsCalibration {
    double sigma2 = 0;  // ticks^2 / s
    double A = 0, k = 0;  // fill intensity lambda(delta) = A exp(-k delta), 1/s and 1/tick
    double gamma = 0;   // 0.4 / (sigma2 * T)
    double lambda_by_delta[6] = {};  // fills per second per side at delta = 1..6 ticks
    double se_k = 0;
};
// Estimates sigma^2 (seed), and A, k from quote probes (100 paths per delta in 1..6). Throws std::runtime_error if some
// delta has no fills or the fitted k <= 0.
[[nodiscard]] AsCalibration calibrate_as_params(const SimConfig& cfg, std::uint64_t seed, unsigned threads);
struct AsSimStats {
    qc::RunningStats pnl, q_final, abs_q, fills;
};
struct AsSimResult {
    AsSimStats as, symmetric;
    qc::RunningStats pnl_diff;  // symmetric - AS on the same path
};
// Runs the AS and symmetric agents on paths 0..R-1 (same stream per path for both: common random numbers).
[[nodiscard]] AsSimResult run_as_sim_experiment(const SimConfig& cfg, const AsCalibration& cal, std::int64_t R,
                                                std::uint64_t seed, unsigned threads);
// sigma^2 (ticks^2/s) from the pooled variance of 10-s mid changes on [600,1200] over agent-free paths.
[[nodiscard]] double estimate_sigma2(const SimConfig& cfg, std::uint64_t seed, int n_paths, unsigned threads);

// ---------------------------------------------------------------- impact validation (DESIGN 8)
struct ImpactConfig {
    std::vector<int> Qs{8, 16, 32, 64, 128};  // metaorder sizes (lots); a control Q = 0 is always run first
    double t0 = 600.0, H = 300.0;
    std::int64_t R = 2000;
    int profile_points = 21;  // mid sampled at t0 + j*H/10, j = 0..20
};
struct ImpactPoint {
    int Q;
    double I, se;  // paired impact over [t0, t0+H] (ticks) and its standard error
    std::int64_t n;
    double step_diff, step_se;  // I(Q) - I(previous Q) and the SE of that paired difference (0 for Q = 0)
};
struct ImpactProfilePoint {
    int Q, j;
    double time, mean, se;  // paired mid change since t0 (ticks)
};
struct ImpactFit {
    bool ok = false;  // false = "insufficient" (fewer than 4 points with I > 4 SE)
    double delta = 0, se_delta = 0, r2 = 0, Y = 0;
    int n_points = 0;
};
struct ImpactResult {
    std::vector<ImpactPoint> points;           // Q = 0 first
    std::vector<ImpactProfilePoint> profile;   // for every Q including the control (all zeros)
    ImpactFit fit;
    double sigma_H = 0;  // sd of the control mid change over H (ticks)
    double V_H = 0;      // mean background traded volume in [t0, t0+H] on control paths (lots)
};
// Runs control + metaorders on paths 0..R-1 (common random numbers across Q; path p uses stream p).
[[nodiscard]] ImpactResult run_impact_experiment(const SimConfig& cfg, const ImpactConfig& ic, std::uint64_t seed,
                                                 unsigned threads);
// Writes impact.csv (Q,I,se,n), impact_profile.csv (Q,j,time,mean,se) and impact_fit.csv (delta,se_delta,r2,Y,n_points).
void write_impact_csvs(const ImpactResult& r, const std::string& dir);

}  // namespace ofl
