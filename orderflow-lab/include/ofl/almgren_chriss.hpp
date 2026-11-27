// Almgren-Chriss (2000) optimal execution, discrete time.
#pragma once
#include <span>
#include <vector>

#include "ofl/types.hpp"

namespace ofl {

// X shares over horizon T (s) in N intervals; sigma in price-units/sqrt(s); gamma permanent impact (price per share),
// eta temporary impact slope (price per share/s), epsilon fixed half-spread cost (price), lambda risk aversion.
struct AcParams {
    double X, T;
    int N;
    double sigma, gamma, eta, epsilon, lambda;
};
struct AcSchedule {
    std::vector<double> x;  // holdings x_0..x_N
    std::vector<double> n;  // trades n_1..n_N
    double kappa;           // 1/s
    double expected_cost;   // E[cost] (price*shares)
    double variance;        // Var[cost]
};
// Closed-form schedule. Throws std::invalid_argument if eta - gamma*tau/2 <= 0 or other parameters are invalid.
[[nodiscard]] AcSchedule almgren_chriss(const AcParams& p);
// Reference: solves the first-order conditions with the Thomas algorithm.
[[nodiscard]] AcSchedule almgren_chriss_tridiagonal(const AcParams& p);
// Risk aversion lambda for which kappa*T = kappaT (other fields of base are used).
[[nodiscard]] double ac_lambda_for_kappaT(const AcParams& base, double kappaT);
// Cumulative rounding of child orders (shares) so that the result sums to X exactly.
[[nodiscard]] std::vector<Qty> round_schedule(std::span<const double> n, Qty X);
// E[cost] and Var[cost] of a given trajectory x_0..x_N under the AC model (direct evaluation).
void ac_evaluate(const AcParams& p, std::span<const double> x, double& expected_cost, double& variance);

}  // namespace ofl
