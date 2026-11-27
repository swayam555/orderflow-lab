# orderflow-lab — TESTPLAN

Rules (from ../PLAN.md §5): fixed seeds; statistical assertions at ≥ 4 SE unless stated; never tune a seed or tolerance to make a test pass — record the failure in NOTES.md and escalate. Labels: `unit` (< 1 s each), `stat` (< 60 s each at -O3), `slow` (> 60 s, excluded from default `./build.sh --test`; run with `-- -L slow`), `parallel` (run under the `tsan` preset too). Every test binary also runs under the `ubsan` preset with `-L unit`.

## 1. Order book — deterministic scenarios (`test_book_basic.cpp`, label unit)
Book band [0, 1000] unless noted. After every call in every scenario: `check_invariants()`.
1. Empty book: best_bid/best_ask = nullopt; market buy qty 5 -> status Cancelled, filled 0, cancelled 5, id 1 consumed (next_id()==2).
2. Rejections: limit qty 0 -> Rejected/BadQty, id kNoOrder, next_id unchanged; price -1 and 1001 -> BadPrice; cancel(999) -> UnknownId; modify unknown -> UnknownId; modify(qty 0) and modify(price 2000) on a live order -> Rejected, order unchanged (find() identical).
3. Resting: buy 10@100 -> Resting, resting 10; best_bid 100; depth_at(Buy,100)=10.
4. Price priority: asks 5@102 (id a), 5@101 (id b); market buy 7 -> fills [b 5@101, a 2@102], maker_done = [true,false]; trade prices are maker prices.
5. Time priority: asks 3@101 (A), 4@101 (B), 5@101 (C); buy limit 6@101 -> fills A3, B3; B remaining 1; C untouched; status Filled.
6. Limit price protection: asks 5@101, 5@103; buy limit 8@102 -> fill 5@101, rest 3@102, status PartialResting; book not crossed.
7. IOC: same as 6 with IOC -> PartialCancelled, filled 5, cancelled 3, nothing rests; IOC that does not cross -> Cancelled.
8. Sell-side mirror of 4-7 (bids, sell orders).
9. Cancel: cancel middle order of a 3-order level -> FIFO of remaining is A, C; level qty updated; cancel same id again -> UnknownId; cancel last order at best -> best moves to next level (also across a gap > 64 ticks to exercise bitmap word scan).
10. Modify reduce in place: A 5@101, B 5@101 (asks); modify(A,101,2) -> Resting, cancelled 3, A still ahead of B (market buy 3 fills A2, B1).
11. Modify increase qty same price: modify(A,101,8) -> loses priority (B fills first), id unchanged.
12. Modify price into cross: bid X 5@99, asks 3@101; modify(X,101,5) -> fills 3@101 with taker = X, rests 2@101 bid, status PartialResting.
13. Fills report owners correctly and `seq` strictly increases across calls.
14. Large sweep: 100 ask levels 1 lot each, market buy 100 -> 100 fills, ask side empty, best_ask nullopt.
15. Band edges: orders at price 0 and 1000 rest/match correctly.

## 2. Order book — randomized differential (`test_book_random.cpp`)
- `unit`: 200 sequences × 2 000 ops; `slow`: 50 sequences × 200 000 ops. Seed 11, stream = sequence index.
- Op mix: 40% GTC limit, 10% IOC, 10% market, 25% cancel (70% live id, 30% random id incl. dead/never-issued), 15% modify (random field changes incl. invalid). Band [0, 300] (small to force many crossings and empty levels); price ~ U{0..300}; qty ~ U{1..10} with 2% qty ≤ 0.
- Apply each op to `OrderBook` and `ReferenceBook`; assert equal `OrderResult`, identical fill vectors (all fields), identical `levels(side, 1000)` every 50 ops and at the end, and `check_invariants()` every op in `unit`, every 1000 ops in `slow`.
- Conservation: Σ submitted qty = Σ filled (taker side) + Σ resting at end + Σ cancelled (incl. IOC/market remainders and modify reductions) — checked per sequence.
- No allocation in steady state: after a pre-sized construction (`expected_live_orders` = 10 000) `pool_capacity()` does not change in the unit runs.

## 3. Hawkes simulation (`test_hawkes_sim.cpp`)
Unit:
1. `spectral_radius_nonneg`: diag(0.3,0.7) -> 0.7; [[0,1],[0,0]] (nilpotent) -> 0 (|ρ|<1e-6); [[0.5,1],[0,0.5]] (Jordan) -> 0.5 within 1e-8; G = α/β of the 2-D model in item 7 (= [[0.4,0.15],[0.2,0.3]]) -> 0.53027756377 within 1e-9; negative entry -> throws; zero matrix -> 0.
2. `stationary_rates`, `params_from_rates` round trip (G, β, Λ) -> params -> Λ within 1e-12; `with_branching_ratio(G_shape, β, Λ, 0.8)` gives ρ(G) = 0.8 within 1e-10 and stationary rates Λ within 1e-12; validate() throws for ρ ≥ 1 and μ ≤ 0.
3. `HawkesState`: hand-computed intensities after events at t = 0.5 (type 0) and 1.0 (type 1) for the 2-D model, evaluated at t = 1.7, within 1e-14 rel; compensator over (1.0, 1.7] equals numerical integral (Simpson, 10⁴ pts) within 1e-10.
4. Log-likelihood of a 3-event hand example equals hand calculation within 1e-12.
5. Determinism: same seed/stream -> identical event vectors; different stream -> different.

Stat (seed 1001):
6. 1-D: μ = 1, α = 1.5, β = 2 (n = 0.75, Λ = 4), T = 1e5. Empirical rate N/T within 3% of 4. (Not an SE test: the variance of N/T for n = 0.75 is ≈ Λ/(T(1−n)²) = 6.4e-4, so sd ≈ 0.025 ≈ 0.6% — 3% is ~5 sd.)
7. 2-D: μ = (0.5, 0.3), α = [[0.8, 0.3], [0.4, 0.6]], β = 2 for all; ρ(G) = 0.53027756377; Λ = (1.012820513, 0.717948718); T = 1e5. Rates within 3%.
8. Time rescaling: for 6 and 7, per-dimension rescaled interarrivals -> `qc::ks_statistic` against the Exp(1) CDF is below `qc::ks_critical_001(n)` (fixed seed, so deterministic; a fail means a bug); also mean of residuals within 4 SE of 1.
9. Wrong model detection (power check, ensures the KS test is not vacuous): rescale the 1-D sample with α = 1.2 instead of 1.5 -> KS statistic > 3·`ks_critical_001(n)`.

## 4. Hawkes EM (`test_hawkes_em.cpp`)
Unit: 1. EM on the 3-event example: one iteration matches hand-computed update within 1e-12. 2. Preconditions throw (unsorted, < 10 events per dim).
Stat (seed 1002, data simulated from the 2-D model above, T = 1e5, β known):
3. Monotone likelihood: loglik_trace non-decreasing within tolerance 1e-9·|ℓ|.
4. Converged within 500 iterations.
5. Recovery: each μ_i and each α_ij within 10% relative of truth (α_ij all ≥ 0.3 so relative tolerance is meaningful); ρ(Ĝ) within 0.03 of 0.53028.
6. Fitted log-likelihood ≥ log-likelihood at the true parameters − 1e-6 (MLE property).
Slow: 7. Same as 5 on the 4-D MarketSim default params with T = 2e5 simulated by `simulate_hawkes` alone (μ within 10%, nonzero α within 15%, zero α < 0.05).

## 5. Market simulator (`test_market_sim.cpp`)
Unit:
1. Initial book matches the spec (20 levels each side, 3×3 lots, best 9999/10001, 120 background live orders).
2. Two sims with the same (seed, path) and the same sequence of `run_until` calls are identical (book levels, counters, fills), and `time()` equals t_end after each call. (Splitting one `run_until` into two calls is allowed to change the realisation, because exponentials are redrawn at the boundary; do not test equality for that.)
3. Background limit orders never cross: assert in a debug hook (counter `bg_crossing_attempts` stays 0) over 1 000 s.
4. Agent ownership: agent_cancel of a background order -> Rejected/UnknownId; on_fill called exactly for fills involving the agent with correct side/qty.
5. Agent market orders excite Hawkes: with `agent_market_orders_excite = true`, intensity of MB right after an agent market buy increases by α_MB,MB = 1.0; with false, unchanged.
Stat (seed 1003, 40 paths × 3 000 s, statistics on [600, 3000]):
6. Hawkes event counts per type / time within 5% of Λ = (1, 1, 0.25, 0.25).
7. Cancellation count ≈ ν × time-integral of live background orders: ratio within 4 SE (Poisson SE).
8. Stationarity sanity: mean live orders per side in [30, 150] (expected ≈ 75); mean spread in [1, 15] ticks; the book never has an empty side during [600, 3000] in ≥ 99% of 1-s samples (report the actual fraction).
9. Mid change variance over 10 s is finite and positive; σ² estimate across two disjoint seed sets agrees within 20%.

## 6. Almgren-Chriss (`test_almgren_chriss.cpp`, unit)
Parameters: X = 1e6, T = 5, N = 5, σ = 0.95, γ = 2.5e-7, η = 2.5e-6, ε = 0.0625, λ ∈ {0, 1e-6, 2e-6, 1e-5} (the paper's numeric example values; used only as realistic magnitudes — no table values asserted).
1. λ = 0 -> linear trajectory exactly (|x_k − X(1−k/N)| < 1e-9·X).
2. Closed form vs tridiagonal reference: max |x_k diff| ≤ 1e-10·X for all λ, and for N ∈ {1, 2, 10, 250}.
3. Σ n_k = X within 1e-9·X; x monotone decreasing for λ > 0.
4. E and Var formulas equal a direct evaluation (Σ over trades of permanent/temporary costs) within 1e-10 rel, and equal the MC mean/variance of the cost on 2e5 simulated price paths (seed 2001) within 4 SE (label stat).
5. Efficient frontier: as λ increases, Var strictly decreases and E strictly increases.
6. Optimality: for λ = 2e-6, perturbing x by ±δ e_j (j = 1..N−1, δ = 1e-3·X) increases E + λVar (quadratic, so strictly).
7. Large κT (κT = 50 and 700, via `ac_lambda_for_kappaT`): every x_k finite, in [0, X], strictly decreasing, and x_1/X equals sinh(κ(T−τ))/sinh(κT) computed in long double within 1e-12 relative (κT = 50 only; for 700, only finiteness and monotonicity).
8. `ac_lambda_for_kappaT` inverts: computed κT equals requested within 1e-12.
9. η̃ ≤ 0 throws.
10. `round_schedule`: sums exactly to X; each child within 1 of n_k; for X = 64, N = 10 TWAP (n_k = 6.4): cumulative llround(6.4k) = 6,13,19,26,32,38,45,51,58,64, so assert the exact vector 6,7,6,7,6,6,7,6,7,6.

## 7. Bertsimas-Lo (`test_bertsimas_lo.cpp`)
Unit:
1. γ = 0: S*_t = W_t/(T−t+1) for all t (with any X), and V_1 = P_0 S̄ + θ S̄²(T+1)/(2T) within 1e-12 rel, for T ∈ {1, 2, 20}.
2. T = 1: S_1 = S̄, V_1 = P_0 S̄ + θS̄² + γρX_0 S̄.
3. Bellman consistency: define Q_t(S) = E[ P'·S + V_{t+1}(P', X', W−S) ] with X' = ρX + η, P' = P + θS + γX' + ε, the expectation computed by 20-node Gauss-Hermite over η (ε enters linearly and has mean 0, so it drops out; the integrand is quadratic in η, so GH is exact). Assert Q_t(S*_t) = V_t(P,X,W) within 1e-9 relative and Q_t(S*_t ± 1) > Q_t(S*_t). Done for t = 1..T−1 at 5 states (seed 3000) with the test parameters of item 4.
Stat (seed 3001):
4. Test parameters: S̄ = 1e5, T = 20, P_0 = 50, θ = 5e-5, γ = 5, ρ = 0.5, σ_ε = 0.125, σ_η² = 0.001, X_0 = 0. 2e5 paths: MC mean of the optimal cost within 4 SE of V_1; TWAP MC mean within 4 SE of `bl_twap_expected_cost`; mean(diff) ≤ 0 and paired diff mean within 4 SE of V_1 − TWAP_expected.
5. Same with X_0 = 0.1 (information value nonzero): optimal cost < TWAP cost by > 4 SE of the paired difference. (Flag: the size of the saving depends on parameters; only the sign is asserted.)

## 8. Avellaneda-Stoikov (`test_as_model.cpp`)
Unit:
1. `as_quotes` at q = 0 symmetric around s; reservation price formula exact; spread at t = T equals (2/γ)ln(1+γ/k).
2. As γ -> 0 (γ = 1e-8) the spread -> 2/k within 1e-6 rel.
3. `as_symmetric_spread` equals the time average of the AS spread over [0,T] (Simpson) within 1e-12.
Stat (seed 8001, 10 000 paths, γ = 0.1, paper parameters):
4. Hard: Var(PnL)_symmetric / Var(PnL)_inventory > 1.5, and sd(q_T)_inventory / sd(q_T)_symmetric < 0.6.
5. Hard: mean PnL symmetric > mean PnL inventory (paired diff > 4 SE).
6. Soft (MESSAGE only, never fails; numbers written to NOTES.md): recalled paper values inventory 62.94/5.89/q sd 2.80, symmetric 67.21/13.43/q sd 8.66, flag if any deviates > 15%.
7. γ ordering: sd(q_T) decreases as γ goes 0.01 -> 0.1 -> 0.5.

## 9. Experiments (`test_experiments.cpp`)
Unit (tiny sizes; check plumbing, not statistics):
1. `calibrate_ac_params` with 10/10/10 paths returns finite values and the schedule rounds to X.
2. BL γ = 0 child schedule == TWAP child schedule (exact vector equality).
3. Impact experiment with R = 8 runs, output has Q = 0 row with I = 0 exactly (difference with itself), and CSV headers as specified.
4. ScheduleAgent with the default config, X = 64, 4 paths: unfilled == 0 and Σ agent fill qty == X on every path; IS computed by the agent equals an independent recomputation from the recorded fills.
Slow (seed per DESIGN):
5. Impact (DESIGN §8, R = 2000): I(Q) strictly increasing in Q (each step > 2 SE of the paired step difference); fitted δ̂ ∈ [0.25, 0.85]; δ̂ + 2·SE(δ̂) < 1 (concave). If 5 fails: run the documented sweep once, record in NOTES.md, mark the test as expected-fail only with written justification — do NOT change Q, H, R, seeds or the fit rule.
6. Execution (DESIGN §6.3, R = 2000): sd(IS) of AC κT = 3 < sd(IS) TWAP (F-ratio < 0.9); unfilled_total = 0 for all strategies. Means are reported, not asserted.
7. AS in sim (R = 1000): sd(q_T) AS < sd(q_T) symmetric (ratio < 0.8).

## 10. Determinism (`test_determinism.cpp`, labels unit+parallel)
1. Impact experiment with R = 16 and Q ∈ {0, 8}: CSV bytes identical for QC_THREADS = 1, 3, 8.
2. BL MC and AS own-model MC: identical RunningStats (bitwise) for threads 1 and 8.
3. Two MarketSims with same (seed, path) produce identical fill streams over 900 s; different path ids differ.
All `parallel`-labelled tests also run under `./build.sh tsan --test -- -L parallel` with zero reports.

## 11. Benchmark checks (`ofl_bench_engine --verify`, label bench, not in default ctest)
- Replay of the op stream on OrderBook and ReferenceBook (first 200 000 ops) gives identical results (proves the benchmarked engine is the correct one).
- `pool_capacity()` unchanged across timed reps.
- Report written; no assertion on speed (numbers are environment facts, not tests).

## 12. Acceptance (definition of done for the project)
`./build.sh release --test` green; `./build.sh release --test -- -L slow` green or with documented expected-fail of §9.5 only; `ubsan` unit green; `tsan` parallel green; all apps run with defaults and write their CSVs; README claims match the CSVs.
