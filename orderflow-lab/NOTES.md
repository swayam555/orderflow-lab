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
