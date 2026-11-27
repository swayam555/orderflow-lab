// Metaorder impact experiment: I(Q), impact profile, power-law fit, optional parameter sweep.
#include <cmath>
#include <cstdio>
#include <string>

#include <qc/cli.hpp>
#include <qc/csv.hpp>
#include <qc/parallel.hpp>

#include "ofl/experiments.hpp"

#ifndef OFL_RESULTS_DIR
#define OFL_RESULTS_DIR "results"
#endif

using namespace ofl;

namespace {
void print_result(const ImpactResult& r) {
    std::printf("  control: sigma_H (sd of mid change over H) %.3f ticks, V_H (background volume in H) %.1f lots\n", r.sigma_H, r.V_H);
    std::printf("  %5s %10s %8s %8s %10s\n", "Q", "I (ticks)", "se", "I/se", "step/se");
    for (const auto& p : r.points)
        std::printf("  %5d %10.4f %8.4f %8.1f %10.1f\n", p.Q, p.I, p.se, p.se > 0 ? p.I / p.se : 0.0,
                    p.step_se > 0 ? p.step_diff / p.step_se : 0.0);
    if (r.fit.ok) {
        std::printf("  fit log I = a + delta log Q over %d points: delta = %.4f +- %.4f (R^2 %.4f), Y = %.3f\n", r.fit.n_points,
                    r.fit.delta, r.fit.se_delta, r.fit.r2, r.fit.Y);
        const bool sr = r.fit.delta >= 0.4 && r.fit.delta <= 0.6;
        std::printf("  -> %s\n", sr ? "square-root-like (delta in [0.4, 0.6])"
                                     : "NOT square-root-like (delta outside [0.4, 0.6]); report as concave/convex with exponent delta");
    } else {
        std::printf("  fit insufficient (%d points with I > 4 SE, need 4)\n", r.fit.n_points);
    }
}
}  // namespace

int main(int argc, char** argv) {
    try {
        qc::Args args(argc, argv, {"paths", "seed", "sweep", "out"});
        const auto R = args.get_int("paths", 2000);
        const auto seed = args.get_u64("seed", 9001);
        const std::string out = args.get_str("out", OFL_RESULTS_DIR);
        const unsigned th = qc::default_threads();
        ImpactConfig ic;
        ic.R = R;
        std::printf("Impact experiment: %lld paths per Q, seed %llu, Q in {0,8,16,32,64,128}, agent_market_orders_excite = true\n",
                    static_cast<long long>(R), static_cast<unsigned long long>(seed));
        const SimConfig cfg;
        const ImpactResult r = run_impact_experiment(cfg, ic, seed, th);
        write_impact_csvs(r, out);
        print_result(r);
        if (args.get_int("sweep", 0) != 0) {
            ImpactConfig sc = ic;
            sc.R = 1000;
            qc::CsvWriter csv(out + "/impact_sweep.csv", {"nu", "K", "delta", "se_delta", "r2", "Y", "n_points", "fit_ok", "I_8",
                                                           "I_16", "I_32", "I_64", "I_128", "sigma_H", "V_H"});
            for (double nu : {0.005, 0.01, 0.02})
                for (int K : {10, 20, 40}) {
                    SimConfig c2 = cfg;
                    c2.cancel_hazard = nu;
                    c2.placement_K = K;
                    std::printf("Sweep nu = %.3f, K = %d (%lld paths):\n", nu, K, static_cast<long long>(sc.R));
                    const ImpactResult s = run_impact_experiment(c2, sc, seed, th);
                    print_result(s);
                    csv.row({nu, std::int64_t{K}, s.fit.delta, s.fit.se_delta, s.fit.r2, s.fit.Y, std::int64_t{s.fit.n_points},
                             std::int64_t{s.fit.ok ? 1 : 0}, s.points[1].I, s.points[2].I, s.points[3].I, s.points[4].I,
                             s.points[5].I, s.sigma_H, s.V_H});
                }
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
