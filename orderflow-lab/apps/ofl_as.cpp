// Avellaneda-Stoikov: own-model reproduction (paper parameters) and in-simulator market making.
#include <cmath>
#include <cstdio>
#include <string>

#include <qc/cli.hpp>
#include <qc/csv.hpp>
#include <qc/parallel.hpp>

#include "ofl/avellaneda_stoikov.hpp"
#include "ofl/experiments.hpp"

#ifndef OFL_RESULTS_DIR
#define OFL_RESULTS_DIR "results"
#endif

using namespace ofl;

int main(int argc, char** argv) {
    try {
        qc::Args args(argc, argv, {"paths", "model_paths", "seed", "model_seed", "calib_seed", "out"});
        const auto R = args.get_int("paths", 1000);
        const auto mpaths = args.get_int("model_paths", 10000);
        const auto seed = args.get_u64("seed", 8003);
        const auto mseed = args.get_u64("model_seed", 8001);
        const auto cseed = args.get_u64("calib_seed", 8002);
        const std::string out = args.get_str("out", OFL_RESULTS_DIR);
        const unsigned th = qc::default_threads();

        // --- own model
        qc::CsvWriter mc(out + "/as_model.csv", {"gamma", "strategy", "mean_pnl", "sd_pnl", "mean_q", "sd_q", "mean_spread"});
        std::printf("Own model (paper parameters s0=100 T=1 sigma=2 dt=0.005 k=1.5 A=140; %lld paths, seed %llu)\n",
                    static_cast<long long>(mpaths), static_cast<unsigned long long>(mseed));
        std::printf("%6s %-10s %9s %9s %8s %8s %8s\n", "gamma", "strategy", "mean PnL", "sd PnL", "mean q", "sd q", "spread");
        for (double g : {0.01, 0.1, 0.5}) {
            AsModelParams p;
            p.gamma = g;
            const AsOwnResult r = as_own_model(p, mpaths, mseed, th);
            const std::pair<const char*, const AsStrategyStats*> rows[2] = {{"inventory", &r.inventory}, {"symmetric", &r.symmetric}};
            for (const auto& [name, st] : rows) {
                mc.row({g, std::string(name), st->pnl.mean(), st->pnl.stddev(), st->q_final.mean(), st->q_final.stddev(),
                        st->spread.mean()});
                std::printf("%6.2f %-10s %9.3f %9.3f %8.3f %8.3f %8.4f\n", g, name, st->pnl.mean(), st->pnl.stddev(),
                            st->q_final.mean(), st->q_final.stddev(), st->spread.mean());
            }
            std::printf("       paired PnL diff (symmetric - inventory) %.3f +- %.3f\n", r.pnl_diff.mean(), r.pnl_diff.se());
        }

        // --- in simulator
        const SimConfig cfg;
        const AsCalibration cal = calibrate_as_params(cfg, cseed, th);
        std::printf("\nSimulator calibration (seed %llu): sigma^2 %.4f ticks^2/s, A %.4f /s, k %.4f /tick (se %.4f), gamma %.6f\n",
                    static_cast<unsigned long long>(cseed), cal.sigma2, cal.A, cal.k, cal.se_k, cal.gamma);
        std::printf("  lambda(delta=1..6) per side per second:");
        for (double l : cal.lambda_by_delta) std::printf(" %.4f", l);
        std::printf("\n");
        {
            qc::CsvWriter pc(out + "/as_params.csv", {"param", "value"});
            pc.row({std::string("sigma2"), cal.sigma2});
            pc.row({std::string("A"), cal.A});
            pc.row({std::string("k"), cal.k});
            pc.row({std::string("se_k"), cal.se_k});
            pc.row({std::string("gamma"), cal.gamma});
            for (int d = 0; d < 6; ++d) pc.row({"lambda_delta" + std::to_string(d + 1), cal.lambda_by_delta[d]});
        }
        const AsSimResult r = run_as_sim_experiment(cfg, cal, R, seed, th);
        qc::CsvWriter sc(out + "/as_sim.csv", {"strategy", "paths", "mean_pnl", "sd_pnl", "se_pnl", "mean_q_T", "sd_q_T",
                                                "mean_abs_q", "mean_fills", "paired_diff_pnl", "se_paired_diff"});
        std::printf("In simulator (%lld paths, seed %llu, session 600-1200 s):\n", static_cast<long long>(R),
                    static_cast<unsigned long long>(seed));
        const std::pair<const char*, const AsSimStats*> rows[2] = {{"AS", &r.as}, {"symmetric", &r.symmetric}};
        for (const auto& [name, st] : rows) {
            const bool is_sym = std::string(name) == "symmetric";
            sc.row({std::string(name), R, st->pnl.mean(), st->pnl.stddev(), st->pnl.se(), st->q_final.mean(),
                    st->q_final.stddev(), st->abs_q.mean(), st->fills.mean(), is_sym ? r.pnl_diff.mean() : 0.0,
                    is_sym ? r.pnl_diff.se() : 0.0});
            std::printf("  %-10s PnL %9.2f sd %8.2f (se %.2f)  q_T mean %6.2f sd %6.2f  mean|q| %6.2f  fills %6.1f\n", name,
                        st->pnl.mean(), st->pnl.stddev(), st->pnl.se(), st->q_final.mean(), st->q_final.stddev(),
                        st->abs_q.mean(), st->fills.mean());
        }
        std::printf("  paired PnL diff (symmetric - AS) %.2f +- %.2f\n", r.pnl_diff.mean(), r.pnl_diff.se());
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
