// Execution experiment: calibrates Almgren-Chriss inputs from the simulator, compares TWAP / AC / BL(gamma=0) by IS.
#include <cmath>
#include <cstdio>
#include <string>

#include <qc/cli.hpp>
#include <qc/csv.hpp>
#include <qc/parallel.hpp>

#include "ofl/bertsimas_lo.hpp"
#include "ofl/experiments.hpp"

#ifndef OFL_RESULTS_DIR
#define OFL_RESULTS_DIR "results"
#endif

using namespace ofl;

int main(int argc, char** argv) {
    try {
        qc::Args args(argc, argv, {"paths", "seed", "calib_seed", "out"});
        const auto R = args.get_int("paths", 2000);
        const auto seed = args.get_u64("seed", 7002);
        const auto cseed = args.get_u64("calib_seed", 7001);
        const std::string out = args.get_str("out", OFL_RESULTS_DIR);
        const unsigned th = qc::default_threads();
        const SimConfig cfg;
        const ExecConfig ec;
        const AcCalibration cal = calibrate_ac_params(cfg, ec, cseed, 400, 200, 500, th);
        {
            qc::CsvWriter c(out + "/ac_params.csv", {"param", "value", "se"});
            c.row({std::string("sigma2_ticks2_per_s"), cal.sigma2, 0.0});
            c.row({std::string("epsilon_ticks"), cal.epsilon, cal.se_epsilon});
            c.row({std::string("eta_ticks_per_lot_per_s"), cal.eta, cal.se_c * (ec.H / ec.N)});
            c.row({std::string("gamma_ticks_per_lot"), cal.gamma, cal.se_gamma});
        }
        std::printf("AC calibration: sigma^2 %.4f ticks^2/s, epsilon %.4f (se %.4f) ticks, eta %.4f ticks/(lot/s), gamma %.5f (se %.5f) ticks/lot\n",
                    cal.sigma2, cal.epsilon, cal.se_epsilon, cal.eta, cal.gamma, cal.se_gamma);
        const std::vector<ExecStrategy> strat = {{"TWAP", ExecKind::Twap, 0.0},
                                                 {"AC_kappaT1", ExecKind::Ac, 1.0},
                                                 {"AC_kappaT3", ExecKind::Ac, 3.0},
                                                 {"BL_gamma0", ExecKind::BlGamma0, 0.0}};
        const auto res = run_execution_experiment(cfg, ec, cal, strat, R, seed, th);
        qc::CsvWriter csv(out + "/execution.csv", {"strategy", "mean_is_per_lot", "sd", "se", "mean_diff_vs_twap", "se_diff",
                                                    "model_E", "model_sd", "unfilled_total"});
        std::printf("%-12s %12s %9s %8s %14s %9s %9s %9s %8s\n", "strategy", "IS/lot", "sd", "se", "diff vs TWAP", "se_diff",
                    "model E", "model sd", "unfilled");
        for (const auto& r : res) {
            csv.row({r.name, r.is.mean(), r.is.stddev(), r.is.se(), r.diff.mean(), r.diff.se(), r.model_E_per_lot,
                     r.model_sd_per_lot, r.unfilled_total});
            std::printf("%-12s %12.4f %9.4f %8.4f %14.4f %9.4f %9.4f %9.4f %8lld\n", r.name.c_str(), r.is.mean(), r.is.stddev(),
                        r.is.se(), r.diff.mean(), r.diff.se(), r.model_E_per_lot, r.model_sd_per_lot,
                        static_cast<long long>(r.unfilled_total));
        }
        std::printf("children (lots): ");
        for (const auto& r : res) {
            std::printf("%s=[", r.name.c_str());
            for (Qty q : r.children) std::printf("%lld ", static_cast<long long>(q));
            std::printf("] ");
        }
        std::printf("\n");
        // Bertsimas-Lo in its own linear model (not in the simulator: see README).
        qc::CsvWriter bl(out + "/bl_model.csv", {"x0", "paths", "opt_mc", "opt_se", "opt_analytic", "twap_mc", "twap_se",
                                                  "twap_analytic", "paired_diff_mc", "paired_diff_se", "paired_diff_analytic"});
        for (double x0 : {0.0, 0.1}) {
            const BlParams p{1e5, 20, 50.0, 5e-5, 5.0, 0.5, 0.125, std::sqrt(0.001), x0};
            const BlSolution sol = bertsimas_lo_solve(p);
            const auto mc = bl_monte_carlo(p, 200000, 3001, th);
            const double v1 = bl_value(sol, 1, p.P0, p.X0, p.S_bar), tw = bl_twap_expected_cost(p);
            bl.row({x0, std::int64_t{200000}, mc.opt.mean(), mc.opt.se(), v1, mc.twap.mean(), mc.twap.se(), tw, mc.diff.mean(),
                    mc.diff.se(), v1 - tw});
            std::printf("BL own model X0=%.1f: opt %.1f +- %.1f (analytic %.1f); TWAP %.1f (analytic %.1f); paired diff %.1f +- %.1f (analytic %.1f)\n",
                        x0, mc.opt.mean(), mc.opt.se(), v1, mc.twap.mean(), tw, mc.diff.mean(), mc.diff.se(), v1 - tw);
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
