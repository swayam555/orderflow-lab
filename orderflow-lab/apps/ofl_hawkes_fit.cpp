// Simulates the 2-D test Hawkes model, fits it by EM (beta known) and writes results/hawkes_fit.csv.
#include <cmath>
#include <cstdio>
#include <string>

#include <qc/cli.hpp>
#include <qc/csv.hpp>
#include <qc/rng.hpp>

#include "ofl/hawkes.hpp"

#ifndef OFL_RESULTS_DIR
#define OFL_RESULTS_DIR "results"
#endif

int main(int argc, char** argv) {
    try {
        qc::Args args(argc, argv, {"T", "seed", "out"});
        const double T = args.get_double("T", 100000.0);
        const std::uint64_t seed = args.get_u64("seed", 1001);
        const std::string out = args.get_str("out", OFL_RESULTS_DIR);
        ofl::HawkesParams truth;
        truth.dim = 2;
        truth.mu = {0.5, 0.3};
        truth.alpha = qc::Mat(2, 2);
        truth.alpha(0, 0) = 0.8;
        truth.alpha(0, 1) = 0.3;
        truth.alpha(1, 0) = 0.4;
        truth.alpha(1, 1) = 0.6;
        truth.beta = qc::Mat(2, 2);
        for (double& b : truth.beta.a) b = 2.0;
        qc::RngStream rng(seed, 0, 1);
        const auto ev = ofl::simulate_hawkes(truth, T, rng);
        const auto fit = ofl::fit_em(ev, T, truth.beta, ofl::default_em_init(ev, T, truth.beta));
        qc::CsvWriter csv(out + "/hawkes_fit.csv", {"param", "true", "fitted", "rel_err"});
        std::printf("events %zu, EM iterations %d, converged %d\n", ev.size(), fit.iterations, fit.converged ? 1 : 0);
        auto emit = [&](const std::string& name, double t, double f) {
            const double rel = std::fabs(f - t) / std::fabs(t);
            csv.row({name, t, f, rel});
            std::printf("  %-8s true %.4f fitted %.4f rel_err %.2f%%\n", name.c_str(), t, f, 100.0 * rel);
        };
        for (int i = 0; i < 2; ++i) emit("mu" + std::to_string(i), truth.mu[static_cast<std::size_t>(i)], fit.params.mu[static_cast<std::size_t>(i)]);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) emit("alpha" + std::to_string(i) + std::to_string(j), truth.alpha(i, j), fit.params.alpha(i, j));
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
