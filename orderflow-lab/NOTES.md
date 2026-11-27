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
