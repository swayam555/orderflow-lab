#include "ofl/almgren_chriss.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace ofl {

namespace {
void check(const AcParams& p) {
    if (p.N < 1) throw std::invalid_argument("AcParams: N = " + std::to_string(p.N) + " must be >= 1");
    if (!(p.T > 0.0)) throw std::invalid_argument("AcParams: T = " + std::to_string(p.T) + " must be > 0");
    if (!(p.X >= 0.0)) throw std::invalid_argument("AcParams: X = " + std::to_string(p.X) + " must be >= 0");
    if (!(p.sigma >= 0.0)) throw std::invalid_argument("AcParams: sigma must be >= 0");
    if (!(p.lambda >= 0.0)) throw std::invalid_argument("AcParams: lambda = " + std::to_string(p.lambda) + " must be >= 0");
    const double tau = p.T / p.N;
    const double eta_t = p.eta - 0.5 * p.gamma * tau;
    if (!(eta_t > 0.0))
        throw std::invalid_argument("AcParams: eta_tilde = eta - gamma*tau/2 = " + std::to_string(eta_t) + " must be > 0");
}
double eta_tilde(const AcParams& p) { return p.eta - 0.5 * p.gamma * (p.T / p.N); }

AcSchedule finish(const AcParams& p, std::vector<double> x, double kappa) {
    AcSchedule s;
    s.kappa = kappa;
    s.x = std::move(x);
    s.n.resize(static_cast<std::size_t>(p.N));
    for (int k = 1; k <= p.N; ++k)
        s.n[static_cast<std::size_t>(k - 1)] = s.x[static_cast<std::size_t>(k - 1)] - s.x[static_cast<std::size_t>(k)];
    ac_evaluate(p, s.x, s.expected_cost, s.variance);
    return s;
}
}  // namespace

void ac_evaluate(const AcParams& p, std::span<const double> x, double& expected_cost, double& variance) {
    const double tau = p.T / p.N, et = eta_tilde(p);
    double abs_sum = 0.0, sq_sum = 0.0, xsq = 0.0;
    for (int k = 1; k <= p.N; ++k) {
        const double nk = x[static_cast<std::size_t>(k - 1)] - x[static_cast<std::size_t>(k)];
        abs_sum += std::fabs(nk);
        sq_sum += nk * nk;
        xsq += x[static_cast<std::size_t>(k)] * x[static_cast<std::size_t>(k)];
    }
    expected_cost = 0.5 * p.gamma * p.X * p.X + p.epsilon * abs_sum + (et / tau) * sq_sum;
    variance = p.sigma * p.sigma * tau * xsq;
}

AcSchedule almgren_chriss(const AcParams& p) {
    check(p);
    const double tau = p.T / p.N, et = eta_tilde(p);
    const double kt2 = p.lambda * p.sigma * p.sigma / et;  // kappa_tilde^2
    const double z = 0.5 * kt2 * tau * tau;
    const double kappa = std::log1p(z + std::sqrt(z * (2.0 + z))) / tau;  // acosh(1+z)/tau
    std::vector<double> x(static_cast<std::size_t>(p.N) + 1);
    const double kT = kappa * p.T;
    for (int k = 0; k <= p.N; ++k) {
        const double tk = k * tau;
        double v;
        if (kT < 1e-8)
            v = p.X * (1.0 - tk / p.T);
        else
            v = p.X * std::exp(-kappa * tk) * (-std::expm1(-2.0 * kappa * (p.T - tk))) / (-std::expm1(-2.0 * kT));
        x[static_cast<std::size_t>(k)] = v;
    }
    x[0] = p.X;
    x[static_cast<std::size_t>(p.N)] = 0.0;
    return finish(p, std::move(x), kappa);
}

AcSchedule almgren_chriss_tridiagonal(const AcParams& p) {
    check(p);
    const double tau = p.T / p.N, et = eta_tilde(p);
    const double kt2 = p.lambda * p.sigma * p.sigma / et;
    const double diag = 2.0 + kt2 * tau * tau;
    const int m = p.N - 1;  // unknowns x_1..x_{N-1}
    std::vector<double> x(static_cast<std::size_t>(p.N) + 1, 0.0);
    x[0] = p.X;
    if (m > 0) {
        std::vector<double> c(static_cast<std::size_t>(m)), d(static_cast<std::size_t>(m));
        // -x_{j-1} + diag x_j - x_{j+1} = 0; x_0 = X on the right-hand side of j=1.
        for (int j = 0; j < m; ++j) {
            const double rhs = (j == 0) ? p.X : 0.0;
            const double sub = -1.0;
            const double denom = (j == 0) ? diag : diag - sub * c[static_cast<std::size_t>(j - 1)];
            c[static_cast<std::size_t>(j)] = -1.0 / denom;
            d[static_cast<std::size_t>(j)] =
                (j == 0) ? rhs / denom : (rhs - sub * d[static_cast<std::size_t>(j - 1)]) / denom;
        }
        x[static_cast<std::size_t>(m)] = d[static_cast<std::size_t>(m - 1)];
        for (int j = m - 2; j >= 0; --j)
            x[static_cast<std::size_t>(j + 1)] = d[static_cast<std::size_t>(j)] - c[static_cast<std::size_t>(j)] * x[static_cast<std::size_t>(j + 2)];
    }
    const double z = 0.5 * kt2 * tau * tau;
    return finish(p, std::move(x), std::log1p(z + std::sqrt(z * (2.0 + z))) / tau);
}

double ac_lambda_for_kappaT(const AcParams& base, double kappaT) {
    check(base);
    if (!(kappaT >= 0.0)) throw std::invalid_argument("ac_lambda_for_kappaT: kappaT must be >= 0");
    if (!(base.sigma > 0.0)) throw std::invalid_argument("ac_lambda_for_kappaT: sigma must be > 0");
    const double tau = base.T / base.N, kappa = kappaT / base.T;
    const double s = std::sinh(0.5 * kappa * tau);
    const double kt2 = 4.0 * s * s / (tau * tau);  // 2(cosh(kappa tau)-1)/tau^2
    return kt2 * eta_tilde(base) / (base.sigma * base.sigma);
}

std::vector<Qty> round_schedule(std::span<const double> n, Qty X) {
    std::vector<Qty> out(n.size());
    double cum = 0.0;
    Qty prev = 0;
    for (std::size_t k = 0; k < n.size(); ++k) {
        cum += n[k];
        Qty c = (k + 1 == n.size()) ? X : static_cast<Qty>(std::llround(cum));
        out[k] = c - prev;
        prev = c;
    }
    return out;
}

}  // namespace ofl
