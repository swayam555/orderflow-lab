#include <qc/test.hpp>

#include <cmath>
#include <string>

#include <qc/parallel.hpp>

#include "ofl/avellaneda_stoikov.hpp"

using namespace ofl;
using qc::test::close;

TEST_CASE("as 1: quote formulas") {
    const double s = 100, T = 1, g = 0.1, sg = 2, k = 1.5;
    const AsQuote q0 = as_quotes(s, 0, 0.3, T, g, sg, k);
    CHECK(close(0.5 * (q0.bid + q0.ask), s, 1e-12));
    const AsQuote q3 = as_quotes(s, 3, 0.3, T, g, sg, k);
    CHECK(close(0.5 * (q3.bid + q3.ask), s - 3 * g * sg * sg * (T - 0.3), 1e-12));
    const AsQuote qT = as_quotes(s, 2, T, T, g, sg, k);
    CHECK(close(qT.ask - qT.bid, (2.0 / g) * std::log(1.0 + g / k), 1e-12));
}

TEST_CASE("as 2: gamma -> 0 spread is 2/k") {
    const AsQuote q = as_quotes(100, 0, 1.0, 1.0, 1e-8, 2.0, 1.5);
    CHECK(close(q.ask - q.bid, 2.0 / 1.5, 0.0, 1e-6));
}

TEST_CASE("as 3: symmetric spread is the time average of the AS spread") {
    const double T = 1, g = 0.1, sg = 2, k = 1.5;
    const int n = 2000;
    double acc = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double t = T * i / n;
        const AsQuote q = as_quotes(100, 0, t, T, g, sg, k);
        const double w = (i == 0 || i == n) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        acc += w * (q.ask - q.bid);
    }
    const double avg = acc * (T / n) / 3.0 / T;
    CHECK(close(as_symmetric_spread(g, sg, k, T), avg, 0.0, 1e-12));
}

TEST_CASE("as 4-6: own model, gamma = 0.1, 10000 paths") {
    const AsModelParams p;
    const auto r = as_own_model(p, 10000, 8001, qc::default_threads());
    const auto &I = r.inventory, &S = r.symmetric;
    MESSAGE("inventory: profit " << I.pnl.mean() << " sd " << I.pnl.stddev() << " q_T sd " << I.q_final.stddev() << " mean q "
                                 << I.q_final.mean() << " mean spread " << I.spread.mean());
    MESSAGE("symmetric: profit " << S.pnl.mean() << " sd " << S.pnl.stddev() << " q_T sd " << S.q_final.stddev() << " mean q "
                                 << S.q_final.mean() << " spread " << S.spread.mean());
    MESSAGE("paired diff (sym - inv) " << r.pnl_diff.mean() << " +- " << r.pnl_diff.se());
    CHECK(S.pnl.variance() / I.pnl.variance() > 1.5);
    CHECK(I.q_final.stddev() / S.q_final.stddev() < 0.6);
    CHECK(r.pnl_diff.mean() > 4.0 * r.pnl_diff.se());
    // soft comparison with recalled paper values (never fails)
    struct Ref {
        const char* what;
        double paper, mine;
    } refs[] = {{"inv profit", 62.94, I.pnl.mean()}, {"inv profit sd", 5.89, I.pnl.stddev()}, {"inv q sd", 2.80, I.q_final.stddev()},
                {"sym profit", 67.21, S.pnl.mean()}, {"sym profit sd", 13.43, S.pnl.stddev()}, {"sym q sd", 8.66, S.q_final.stddev()}};
    for (const Ref& x : refs)
        MESSAGE("soft " << std::string(x.what) << ": recalled " << x.paper << " measured " << x.mine << " dev "
                        << 100.0 * (x.mine - x.paper) / x.paper << "%" << std::string(std::fabs(x.mine - x.paper) / x.paper > 0.15 ? " FLAG>15%" : ""));
}

TEST_CASE("as 7: sd(q_T) decreases with gamma") {
    double prev = 1e300;
    for (double g : {0.01, 0.1, 0.5}) {
        AsModelParams p;
        p.gamma = g;
        const auto r = as_own_model(p, 10000, 8001, qc::default_threads());
        MESSAGE("gamma " << g << ": inventory sd(q_T) " << r.inventory.q_final.stddev() << " profit " << r.inventory.pnl.mean()
                         << " sd " << r.inventory.pnl.stddev());
        CHECK(r.inventory.q_final.stddev() < prev);
        prev = r.inventory.q_final.stddev();
    }
}
