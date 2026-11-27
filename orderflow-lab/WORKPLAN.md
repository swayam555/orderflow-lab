# orderflow-lab — WORKPLAN

Prerequisite: ../PLAN.md work packages C0 (toolchain, build.sh, skeleton) and C1 (qc library) are done and green.
Each work package (WP) ends with: build `release` + `ubsan` green for its tests, a commit, and a short entry in `NOTES.md` (what was done, any deviations, measured numbers). Do the WPs in order; do not start a WP while an earlier one has a failing test. If a spec here is wrong (a test cannot pass for a mathematical reason), stop, write the evidence in NOTES.md and report — do not silently change the spec.

| WP | Content | Files | Tests that must pass | Size |
|---|---|---|---|---|
| A1 | Types + `OrderBook` (pool, levels, bitmaps, id table, matching, cancel, modify, `check_invariants`) | types.hpp, order_book.{hpp,cpp} | TESTPLAN §1 | L |
| A2 | `ReferenceBook`; randomized differential test; op-stream generator + `ofl_bench_engine` (throughput, latency, `--verify`) | reference_book.*, bench/ | §2, §11 | M |
| A3 | Hawkes params, `spectral_radius_nonneg` (in qc if not done in C1), `HawkesState`, `simulate_hawkes`, likelihood, rescaled residuals | hawkes.* | §3 | M |
| A4 | EM calibration + `ofl_hawkes_fit` app | hawkes.*, apps/ofl_hawkes_fit.cpp | §4 | M |
| A5 | `SimConfig`, `MarketSim` event loop, `Agent` interface, `ObserverAgent`, `ofl_depth` app | market_sim.*, agents.* | §5, §10.3 | L |
| A6 | Almgren-Chriss (closed form, tridiagonal reference, λ↔κT, rounding); Bertsimas-Lo (solve, policy, value, TWAP cost, MC) | almgren_chriss.*, bertsimas_lo.* | §6, §7, §10.2 (BL part) | M |
| A7 | `ScheduleAgent`, `ProbeAgent`, `calibrate_ac_params`, execution experiment, `ofl_execution` app | agents.*, experiments.*, apps/ofl_execution.cpp | §9.1, §9.2, §9.4, §9.6 | M |
| A8 | AS formulas, own-model simulator, `AsAgent` (+ symmetric), `calibrate_as_params`, `ofl_as` app | avellaneda_stoikov.*, agents.*, apps/ofl_as.cpp | §8, §9.7, §10.2 (AS part) | M |
| A9 | Impact experiment (+ profile, fit, sweep), `ofl_impact` app | experiments.*, apps/ofl_impact.cpp | §9.3, §9.5, §10.1 | M |
| A10 | Full runs of every app with defaults; README.md; tsan run; final NOTES.md | README.md | §12 | S |

## Details and ordering notes
A1
- Implement the `OrderBook` exactly per DESIGN §3; write `check_invariants()` first (it checks: FIFO links consistent both ways, level qty/count equal sums, bitmap bit set iff level non-empty, cached best = actual best, id table ↔ pool consistency, not crossed, all live remaining > 0).
- No exceptions on the hot path; rejections are return values. `check_invariants` may throw.
A2
- ReferenceBook is deliberately simple; do not share code with OrderBook (shared code would share bugs).
- The bench reports numbers; no speed assertions. Run it on an idle machine; record the CPU line and the median in NOTES.md.
A3/A4
- Use the parameter values of TESTPLAN §3–4 verbatim. EM must use the O(K·D²) recursion; no O(K²) loops anywhere except the hand-check test.
A5
- The event loop in DESIGN §5.3 is normative, including the tie order and the RNG purpose numbers. Background bookkeeping (`bg_live`) must be updated from fills (filled-out makers) and cancels.
- After A5 run `ofl_depth` and record the mean spread / depth in NOTES.md (needed to judge A9).
A6
- Pure math, no simulator. The BL recursion is the architect's derivation (DESIGN §6.2 flag); tests §7.1–7.3 validate it. If they fail, re-derive and document; do not import a remembered formula.
A7
- Calibration and evaluation use different seeds (DESIGN §6.3). If `calibrate_ac_params` throws (η̃ ≤ 0), record the estimates in NOTES.md and stop: this is a spec problem for the architect.
A8
- Soft paper-value comparisons go to NOTES.md as a table (recalled vs measured).
A9
- Run the default config first (R = 2000). Only if TESTPLAN §9.5 fails, run `--sweep=1` once and report the full table. Never change Q grid, H, R, seeds, or the fit rule.
A10 — README must contain:
- What is modelled and what is not (DESIGN §1 non-goals), with paper references.
- Benchmark numbers with the exact methodology sentence: "single thread, in-memory, pre-generated op stream, no I/O or networking, timer overhead X ns not subtracted, CPU …, clang …, -O3; not comparable to exchange tick-to-trade latency".
- Impact: measured δ̂ ± SE with the plot data path; the phrase "square-root-like" only if δ̂ ∈ [0.4, 0.6]; otherwise "concave with exponent δ̂". State that `agent_market_orders_excite = true` and that this is inspired by, not an implementation of, the Oxford LOB-Impact benchmark.
- Execution: IS table with SEs; explicitly state BL-with-information is only evaluated in its own linear model and why; BL recursion flagged as derived here and MC-verified.
- AS: own-model table vs recalled paper values (flagged as recalled), in-sim results.
- Reproduce commands (`./build.sh release`, then each app).

## Time-boxing
Expected effort: A1 L (1 day), A2–A9 M (½ day each), A10 S. If a `stat`/`slow` test fails after the implementation passes all unit tests, spend at most 2 hours debugging before writing up the evidence and reporting.
