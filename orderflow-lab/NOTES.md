# orderflow-lab NOTES

## A1 - OrderBook
Done per DESIGN 3. Matching, cancel, modify, bitmap best-price scans, check_invariants. TESTPLAN 1 passes (release, ubsan).
Convention: rejected cancel/modify return `id` = the id passed (rejected submits return kNoOrder).

## A2 - ReferenceBook, differential test, bench
- ReferenceBook (std::map of deques) shares no code with OrderBook. test_book_random (unit, 200x2000, invariants every op) and
  test_book_random_slow (label slow, 50x200000) pass. Conservation identity used: entered lots = taker-filled + maker-filled
  + cancelled (incl. modify removals) + resting; the TESTPLAN wording counts only the taker side, so the maker side must be
  included (each filled lot is consumed from two orders). This is a clarification, not a spec error.
- `ofl_bench_engine` (default 5M ops, 1M warmup, 5 reps, seed 42), run on a NOT fully idle machine (load avg ~2.8, other agent building):
  CPU AMD Ryzen 7 5700U, clang 22.1.8 (zig), -O3. Verify: first 200000 ops identical to ReferenceBook.
  Throughput median 7.79 Mops/s = 128.3 ns/op (min 7.33, max 8.63 Mops/s); pool_capacity unchanged; fills/op 0.219; max live 618525; mean book levels 51.9.
  Per-op latency is dominated by the timer: the timer overhead (median back-to-back now_ns pair) measured 1396 ns on this machine
  (clock_gettime is slow here), and per-op p50 = 1397 ns equals it. Per-op latency percentiles are therefore reported but are
  timer-floor limited; the throughput ns/op is the reliable number. Overhead is not subtracted (per spec).
  latency all: p50 1397, p90 1886, p99 2305, p99.9 10966, max 498192 ns.

## A3/A4 - Hawkes simulation, likelihood, EM
- Implemented per DESIGN 4 (hawkes.{hpp,cpp}); `qc::spectral_radius_nonneg` already in qc. Extra API: `em_step` (single step, no
  preconditions; used for the hand-check test), `default_em_init`.
- Stat results (seed 1001/1002): 1-D rate 4.0121 (target 4); 2-D rates 1.0177, 0.7232 (targets 1.0128, 0.7179); KS of rescaled
  residuals far below the 1% critical value (e.g. 1-D KS 6.9e-4 vs 2.6e-3); wrong-model KS 0.065 vs 3*crit 0.0077.
- EM on 2-D model, T=1e5 (seed 1002): converged in 39 iterations, mu rel err 0.4%/0.6%, alpha rel err 0.5%, 3.8%, 0.2%, 1.3%;
  rho fitted 0.5334 vs 0.5303. ofl_hawkes_fit (seed 1001): all relative errors <= 1.33%.
- Slow 4-D EM (T=2e5, 498620 events): not converged after 500 iterations (rel_tol 1e-10; test 7 does not require convergence);
  nonzero mu within 1.5%, nonzero alpha within 2%, zero alpha all < 0.006.
- test_hawkes_em.cpp holds tests 1-6; test 7 is the label-slow executable test_hawkes_em_slow.
- ofl_hawkes_fit CSV: rel_err = |fitted-true|/|true| (no true parameter is zero in the 2-D model).

## A5 - MarketSim
- Event loop exactly per DESIGN 5.3 (tie order agent > cancel > Hawkes; RNG purposes 1..5). Additions beyond the spec API (all additive):
  `Counters::{bg_crossing_attempts, bg_live_integral}`, `bg_live_count()`, `set_fill_observer()` (used by tests, determinism check).
  Agent `on_fill` must not call back into the sim (the fill buffer is being iterated); documented.
- Stat test (seed 1003, 40 paths x 3000 s, window [600,3000]): rates 1.001, 0.9997, 0.2486, 0.2521 (targets 1,1,.25,.25);
  cancels 144622 vs nu*integral 145108 (z = -1.28); live background orders per side 75.4 / 75.8; mean spread 2.56 ticks;
  fraction of 1-s samples with both sides non-empty 1.0 (all); Var(10 s mid change) 4.14 vs 4.37 (two disjoint path sets, 5% apart).
- ofl_depth (200 paths, seed 9100, samples every 10 s on [600,1200]): mean spread 2.58 ticks, 151.2 live background orders
  (both sides); depth at the touch 4.54 (bid) / 4.56 (ask) lots; depth summed over the first 40 levels from each side's best, averaged
  over sides: 188 lots. Distance in depth_profile.csv is measured from each side's own best quote. Mid-change variance 4.2 ticks^2/10 s
  (sigma^2 ~ 0.42 ticks^2/s). Per-lot market impact of one lot at the touch is therefore small relative to depth: Q = 128 lots
  over 300 s consumes ~ the depth of the first ~25 levels of one side, replenishment aside - relevant for A9.

## A6 - Almgren-Chriss and Bertsimas-Lo (pure math)
- All TESTPLAN 6 and 7 tests pass on the first run with the DESIGN recursions; the BL recursion (a,b,c,d) was re-derived by hand
  from the Bellman equation (minimising over S gives S* = W - (theta W - rho b X)/(2a); substituting back reproduces the stated
  a_t, b_t, c_t, d_t) and verified by Gauss-Hermite Bellman consistency (to 1e-9) and Monte Carlo.
- BlParams.sigma_eps / sigma_eta are standard deviations (the test passes sqrt(0.001) for sigma_eta^2 = 0.001).
- Direct-evaluation check in test 6.4: trade k executes at S_{k-1} - eps - (eta/tau) n_k, S_k = S_{k-1} - gamma n_k + sigma sqrt(tau) xi;
  this reproduces E = gamma X^2/2 + eps sum|n| + (eta - gamma tau/2)/tau sum n^2 exactly.
- BL MC (2e5 paths, seed 3001, X0 = 0): optimal 5.26162e6 +- 191 vs V1 5.26156e6; TWAP 5.26256e6 +- 191 vs analytic 5.26250e6;
  paired diff -939.2 +- 13.3 vs -943.3 predicted. X0 = 0.1: paired diff opt - TWAP = -1653.3 +- 20.8.
- AC MC (2e5 paths, seed 2001): E and Var within 4 SE for all four lambdas (e.g. lambda 2e-6: E 1.1407e6 vs MC 1.1411e6 +- 1.0e3).

## A7 - Execution experiment
- ScheduleAgent / ProbeAgent / calibrate_ac_params / run_execution_experiment / ofl_execution implemented per DESIGN 6.3.
  (test_determinism.cpp was created here with TESTPLAN 10.3 and 10.2 (BL) - those are due in A5/A6 by the WORKPLAN table; AS and
  impact determinism are added in A8/A9.)
- Calibration (seed 7001; 400/200/500 paths): sigma^2 0.4159 ticks^2/s; epsilon 1.1347 +- 0.0314 ticks; eta 5.892 ticks per (lot/s)
  (OLS slippage slope c x tau); gamma 0.0888 +- 0.0067 ticks/lot. No throw (eta_tilde = 5.89 - 0.089*30/2 = 4.56 > 0).
  The tiny 10/10/10 calibration of test exp 1 also returns finite values (eps 1.26, eta 5.80, gamma 0.048).
- Evaluation (R = 2000, seed 7002, CRN; IS per lot in ticks): TWAP 5.458 (sd 4.089, se 0.091); AC kT=1 5.368 (sd 3.904), diff -0.090 +- 0.070;
  AC kT=3 5.112 (sd 3.166), diff -0.346 +- 0.079; BL gamma=0 identical to TWAP. unfilled_total = 0 for all.
  AC model prediction (E, sd per lot): TWAP 4.95 / 5.96; kT=1 4.97 / 5.55; kT=3 5.48 / 3.83.
  Observed: the variance reduction is real (Var ratio kT=3 / TWAP = 0.599 < 0.9) but the model's sd over-predicts and its E slightly
  under-predicts TWAP and over-predicts kT=3 (AC cost model is a linear approximation; the simulator has mean-reverting spread/depth).
  AC kT=3 also has *lower* mean IS than TWAP here (-0.346 +- 0.079, 4.4 SE), i.e. front-loading did not cost more in this simulator;
  that is a measured fact about this model (front loading takes liquidity when the book is freshly replenished), not a claim of general validity.
- BL own model (2e5 paths, seed 3001) in bl_model.csv: see A6 numbers.
