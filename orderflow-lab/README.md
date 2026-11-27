# orderflow-lab

A C++20 limit-order-book (LOB) matching engine, an exact multivariate Hawkes order-flow simulator, a discrete-event market
simulator built on both, and four experiments on top of it: optimal execution (Almgren-Chriss, Bertsimas-Lo), Avellaneda-Stoikov
market making, Hawkes calibration and a metaorder-impact measurement. Only `qc` (the shared support library: Philox RNG,
statistics, parallel loop) and the vendored doctest header are used; no other third-party code.

Every number below comes from a CSV or log produced by the apps / tests in this repository (`results/`, git-ignored), with the
seeds stated. Where a hypothesis failed, it is reported as failed.

## What is modelled, and what is not

Modelled: one instrument; price-time priority; GTC limit, IOC limit, market, cancel, modify (reduce in place keeps priority, anything
else loses it); integer ticks and lots; background flow = 4-dimensional exponential-kernel Hawkes process (limit buy/sell, market
buy/sell) plus a per-order cancellation hazard (zero-intelligence style); pluggable agents.

Not modelled (non-goals): multiple instruments, networking/FIX/ITCH, persistence, self-trade prevention, iceberg/hidden/stop/post-only
orders, auctions, fees, multi-threaded matching, non-exponential Hawkes kernels, estimation of the Hawkes decay (beta is fixed and
known), real market data, RL agents. It does **not** reproduce or claim compatibility with the Oxford LOB-Impact benchmark.

References: Hawkes (1971); Ogata (1981) thinning; Veen and Schoenberg (2008) EM; Brown et al. (2002) time rescaling; Smith, Farmer,
Gillemot, Krishnamurthy (2003); Almgren and Chriss (2000); Bertsimas and Lo (1998); Avellaneda and Stoikov (2008); Toth et al. (2011).

## Reproduce

```
./build.sh release                 # builds libofl, apps, bench
./build.sh release --test          # unit + stat + slow tests (20 ctest tests incl. common)
./build.sh ubsan --test            # unit + stat (17 tests); ./build.sh tsan --test  (6 parallel-labelled tests)
build/release/orderflow-lab/ofl_bench_engine --verify=1
build/release/orderflow-lab/ofl_hawkes_fit      # -> orderflow-lab/results/hawkes_fit.csv
build/release/orderflow-lab/ofl_depth           # -> depth_profile.csv
build/release/orderflow-lab/ofl_impact [--sweep=1]
build/release/orderflow-lab/ofl_execution
build/release/orderflow-lab/ofl_as
```
Results are bitwise identical for any `QC_THREADS` (tested for the impact CSVs, Bertsimas-Lo MC, AS own-model MC and the simulator fill stream). Cross-machine bitwise identity is not claimed.

## Matching engine benchmark

`ofl_bench_engine` (seed 42; 5 000 000 measured ops after 1 000 000 warm-up; 5 reps; band [0, 65535]; op mix 55% limit / 30% cancel /
10% market / 5% modify; 0.219 fills per op; mean book depth 51.9 price levels; up to 618 525 live orders):

- Methodology: single thread, in-memory, pre-generated op stream, no I/O or networking, timer overhead 1396 ns not subtracted, CPU AMD Ryzen 7 5700U with Radeon Graphics, clang 22.1.8 (zig c++), -O3 (no -march=native), not pinned; not comparable to exchange tick-to-trade latency.
- Throughput: median **8.66 Mops/s = 115.5 ns/op** (5 reps: 8.48 to 8.75 Mops/s). Machine was not fully idle (load average about 3.5 from another build job).
- Per-op latency, every op timed individually: p50 1397 ns, p99 1956 ns, p99.9 10 966 ns, max 105 885 ns. **These are dominated by the timer**: a back-to-back pair of `steady_clock` reads costs 1396 ns on this machine, equal to the p50. They bound the latency from above and are not an engine latency. The less distorted estimate is the 64-op batch mean (timer cost amortised): p50 127 ns/op, p99 261 ns/op, p99.9 410 ns/op (`bench_latency.csv`, row `batch64_per_op`); batch means hide tails inside a batch.
- `--verify`: the first 200 000 ops replayed on `OrderBook` and on the naive `ReferenceBook` give identical results and fills; `pool_capacity()` does not change across timed reps.

## Hawkes simulation and calibration

Exact Ogata thinning; time-rescaled residuals pass a KS test at the 1% level on 1-D and 2-D test models (e.g. KS 6.9e-4 vs critical 2.6e-3 for n = 401 213), and the same test rejects a deliberately wrong model (KS 0.065 vs 3 x critical 0.0077). EM with known beta (O(K D^2) recursion): on the 2-D model (T = 1e5, 173 005 events, seed 1001) converged in 41 iterations with relative errors mu 0.86% / 1.33%, alpha 0.56% / 0.96% / 0.32% / 1.17%; on a 4-D market-simulator flow (T = 2e5) non-zero parameters within 2% and zero parameters below 0.006 (not converged to 1e-10 in 500 iterations).

## Market simulator

Default background flow: Lambda = (1, 1, 0.25, 0.25) per second, branching radius 0.5, beta = 2/s, cancellation hazard 0.01/s per order, limit orders placed 1 to 20 ticks behind the opposite best. Measured (40 paths x 2400 s): rates 1.001 / 1.000 / 0.249 / 0.252; cancels match nu x integral of live orders (z = -1.3); 75.4 / 75.8 live orders per side; mean spread 2.56 ticks; both sides always non-empty in the samples; Var(10 s mid change) about 4.2 ticks^2 (two disjoint path sets 4.14 and 4.37). `ofl_depth` (200 paths): depth at the touch 4.5 lots, 8.4 at 4 ticks, 11.2 at 9 ticks, about 0 beyond 25 ticks.

## Impact validation (hypothesis rejected)

Metaorder: Q one-lot market buys spread evenly over H = 300 s (t0 = 600 s), paired against a control on the same random streams (R = 2000 paths per Q, seed 9001). `agent_market_orders_excite = true` (agent market orders feed the Hawkes state), which is a modelling choice. This is inspired by, not an implementation of, the Oxford LOB-Impact benchmark idea.

| Q | I(Q) ticks | SE |
|---|---|---|
| 8 | 1.919 | 0.184 |
| 16 | 3.764 | 0.189 |
| 32 | 7.980 | 0.191 |
| 64 | 15.480 | 0.217 |
| 128 | 29.997 | 0.257 |

Fit of log I on log Q: **delta = 0.997 +- 0.013** (R^2 0.9995, 5 points; plot data `results/impact.csv`, `impact_profile.csv`, `impact_fit.csv`). delta is outside [0.4, 0.6] and even outside the concave range: **in this model impact is essentially linear in Q (about 0.23 ticks per lot) and permanent** (for Q = 64 the profile is 15.5 ticks at t0 + H and 14.9 +- 0.3 at t0 + 2H). It is therefore **not** square-root-like, and the DESIGN hypothesis that a zero-intelligence + Hawkes book would produce delta near 0.5 is rejected here. The documented sweep (nu in {0.005, 0.01, 0.02} x K in {10, 20, 40}, R = 1000, `results/impact_sweep.csv`) gave delta between 0.92 and 1.17 in all nine configurations. A plausible but untested reason is that background orders are priced relative to the current opposite best, so a displaced price is never restored. The test `exp 5b` records the failed concavity check as an expected failure (see NOTES.md).

## Execution (Almgren-Chriss, Bertsimas-Lo)

Buy 64 lots over 300 s in 10 market-order slices. AC inputs estimated from the simulator (seed 7001): sigma^2 = 0.416 ticks^2/s, epsilon = 1.135 +- 0.031 ticks, eta = 5.89 ticks per (lot/s), gamma = 0.089 +- 0.007 ticks/lot. Evaluation: 2000 paths per strategy, common random numbers, seed 7002; IS in ticks per lot (positive = worse than the arrival mid):

| strategy | mean IS | SE | sd | mean diff vs TWAP (SE) | AC model E / sd | unfilled |
|---|---|---|---|---|---|---|
| TWAP | 5.458 | 0.091 | 4.089 | 0 | 4.95 / 5.96 | 0 |
| AC kappa T = 1 | 5.368 | 0.087 | 3.904 | -0.090 (0.070) | 4.97 / 5.55 | 0 |
| AC kappa T = 3 | 5.112 | 0.071 | 3.166 | -0.346 (0.079) | 5.48 / 3.83 | 0 |
| BL, gamma = 0 | identical to TWAP | | | | | 0 |

Front-loading (kappa T = 3) lowered the IS standard deviation (variance ratio 0.60 vs TWAP) and, in this simulator, also the mean; the AC model's predicted E and sd are only roughly in line with the simulated values (it over-predicts the sd of TWAP and kappa T = 1). The ordering of means is a result for this simulator, not a general claim.
Bertsimas-Lo with information (gamma > 0) is evaluated **only in its own linear model** (`bl_model.csv`): in the simulator the information signal (the Hawkes intensity) has a half-life of about 0.35 s and has decorrelated by the next 30-s slice, so an in-simulator BL-with-information strategy would be meaningless. In its own model (S = 1e5, T = 20, gamma = 5, rho = 0.5, 2e5 paths) the Monte Carlo cost of the optimal policy is 5 261 618 +- 191 against the analytic 5 261 557, and the optimal-minus-TWAP paired difference is -939 +- 13 (analytic -943). The BL value-function recursion was derived here from the Bellman equation (not copied from the paper) and is verified by Gauss-Hermite Bellman consistency, the gamma = 0 closed form, and Monte Carlo.

## Avellaneda-Stoikov

Own model (paper parameters s0 = 100, T = 1, sigma = 2, k = 1.5, A = 140; 10 000 paths, seed 8001). The "paper" column is **recalled from memory and unverified**:

| gamma = 0.1 | paper (recalled) | measured |
|---|---|---|
| inventory profit (sd) | 62.94 (5.89) | 64.89 (6.49) |
| inventory sd of final q | 2.80 | 2.92 |
| symmetric profit (sd) | 67.21 (13.43) | 68.07 (13.42) |
| symmetric sd of final q | 8.66 | 8.41 |

The inventory strategy cuts the profit variance by 4.3x and the final-inventory sd by 2.9x for 3.18 +- 0.12 less mean profit. In the simulator (calibrated A = 0.184/s, k = 0.713/tick, gamma = 0.00166; 1000 paths, seed 8003, 600-s session): AS PnL 62.6 (sd 31.8), sd of final inventory 4.03; symmetric PnL 78.0 (sd 60.2), final-inventory sd 8.06; paired PnL difference (symmetric - AS) 15.4 +- 2.0. The inventory skew halves the inventory risk here but costs mean PnL.

## Tests

Release: 20/20 ctest tests pass (includes the `slow` label); ubsan: 17/17 (unit + stat); tsan: 6/6 parallel-labelled tests with no reports. The one documented expected failure is the impact concavity check (`exp 5b`, `doctest::should_fail`). See `NOTES.md` for deviations and open items.

### Impact: what was learned
The default model's impact is linear and permanent (exponent 0.997). A follow-up diagnostic (see NOTES.md, "A9 follow-up") showed this comes
from pricing background limit orders off the opposite best; an opt-in anchored rule (`fundamental_anchor`) removes it but still gives no
square-root law. This project therefore does **not** reproduce square-root impact.
