# orderflow-lab — DESIGN

Namespace `ofl`. Depends only on `qc` (see ../PLAN.md). Units everywhere: **time in seconds, price in integer ticks, quantity in integer lots**, money = ticks x lots.

## 1. Goals and non-goals

Goals
1. A correct, fast single-instrument limit-order-book (LOB) matching engine: price-time priority, GTC limit, IOC limit, market, cancel, modify.
2. An exact multivariate Hawkes simulator with exponential kernels (Ogata thinning), branching-ratio control and EM maximum-likelihood calibration.
3. A discrete-event market simulator: Hawkes-driven background order flow + per-order cancellation hazard + pluggable agents, all on top of the engine.
4. Execution: Almgren-Chriss (AC) closed-form schedule, Bertsimas-Lo (BL) dynamic program; evaluation by implementation shortfall (IS) against TWAP.
5. Avellaneda-Stoikov (AS) market maker with inventory risk, both in its own model and inside the simulator.
6. Impact validation: measure the metaorder impact curve I(Q) produced by the simulator and test concavity (square-root-law-like), inspired by the Oxford LOB-Impact benchmark idea (arXiv:2610.27100).
7. An honest latency/throughput benchmark of the engine.

Non-goals (do not build): multiple instruments, networking/FIX/ITCH, persistence, self-trade prevention, iceberg/hidden/stop/post-only orders, auction phases, fees, multi-threaded matching, Hawkes kernels other than exponential, estimation of Hawkes decay β (β is fixed and known), reproducing the LOB-Impact benchmark protocol or claiming compatibility with it, RL agents, real market data.

## 2. Layout

```
orderflow-lab/
  CMakeLists.txt            # lib ofl (STATIC), apps, bench, tests
  DESIGN.md TESTPLAN.md WORKPLAN.md NOTES.md(created by coder) README.md(WP A10)
  include/ofl/
    types.hpp               # Tick, Qty, OrderId, OwnerId, Side, Tif, Status, Reject, Fill, OrderResult
    order_book.hpp          # OrderBook
    reference_book.hpp      # naive reference book (used only by tests/bench generator)
    hawkes.hpp              # HawkesParams, HawkesState, simulate_hawkes, likelihood, EM
    market_sim.hpp          # SimConfig, MarketSim, Agent
    agents.hpp              # ScheduleAgent (execution), ObserverAgent, ProbeAgent, AsAgent
    almgren_chriss.hpp      # AcParams, AcSchedule, almgren_chriss, tridiagonal reference, round_schedule
    bertsimas_lo.hpp        # BlParams, BlSolution, solve, policy, value, monte carlo
    avellaneda_stoikov.hpp  # AS formulas + own-model simulator
    experiments.hpp         # impact experiment, execution experiment, AS-in-sim experiment, calibration helpers
  src/ order_book.cpp reference_book.cpp hawkes.cpp market_sim.cpp agents.cpp almgren_chriss.cpp
       bertsimas_lo.cpp avellaneda_stoikov.cpp experiments.cpp
  apps/ ofl_impact.cpp ofl_execution.cpp ofl_as.cpp ofl_hawkes_fit.cpp ofl_depth.cpp
  bench/ ofl_bench_engine.cpp
  tests/ test_book_basic.cpp test_book_random.cpp test_hawkes_sim.cpp test_hawkes_em.cpp test_market_sim.cpp
         test_almgren_chriss.cpp test_bertsimas_lo.cpp test_as_model.cpp test_experiments.cpp test_determinism.cpp
  results/                  # CSV outputs (git-ignored)
```

## 3. Matching engine

### 3.1 Types (`types.hpp`)
```cpp
namespace ofl {
using Tick = std::int32_t; using Qty = std::int64_t; using OrderId = std::uint64_t; using OwnerId = std::uint32_t;
inline constexpr OrderId kNoOrder = 0;          // ids start at 1
inline constexpr OwnerId kBackground = 0;       // agents get 1,2,...
enum class Side : std::uint8_t { Buy = 0, Sell = 1 };
enum class Tif  : std::uint8_t { GTC, IOC };
enum class Status : std::uint8_t { Resting, Filled, PartialResting, PartialCancelled, Cancelled, Rejected };
enum class Reject : std::uint8_t { None, BadQty, BadPrice, UnknownId };
struct Fill { std::uint64_t seq; OrderId taker, maker; OwnerId taker_owner, maker_owner;
              Tick price; Qty qty; Side taker_side; bool maker_done; };     // price = maker's resting price
struct OrderResult { OrderId id; Status status; Reject reject; Qty filled; Qty resting; Qty cancelled; };
struct LevelView { Tick price; Qty qty; std::int32_t orders; };
struct OrderView { OrderId id; OwnerId owner; Side side; Tick price; Qty remaining; };
}
```

### 3.2 `OrderBook` API
```cpp
class OrderBook {
 public:
  // Valid prices are [min_tick, max_tick]; reserves pool for expected_live_orders and the id table for expected_total_ids.
  OrderBook(Tick min_tick, Tick max_tick, std::size_t expected_live_orders, std::size_t expected_total_ids);
  [[nodiscard]] OrderResult submit_limit(Side, Tick price, Qty qty, Tif, OwnerId);
  [[nodiscard]] OrderResult submit_market(Side, Qty qty, OwnerId);
  [[nodiscard]] OrderResult cancel(OrderId);
  [[nodiscard]] OrderResult modify(OrderId, Tick new_price, Qty new_qty);   // new_qty = new remaining quantity
  std::span<const Fill> last_fills() const noexcept;   // fills of the most recent call only; invalidated by next call
  std::optional<Tick> best_bid() const noexcept;  std::optional<Tick> best_ask() const noexcept;
  Qty depth_at(Side, Tick) const;  std::int32_t orders_at(Side, Tick) const;
  std::vector<LevelView> levels(Side, std::size_t max_levels) const;  // best first
  std::optional<OrderView> find(OrderId) const;
  std::size_t live_orders() const noexcept;
  OrderId next_id() const noexcept;                    // id the next accepted submit will receive
  void check_invariants() const;                       // O(book); throws std::logic_error with description
  std::size_t pool_capacity() const noexcept;          // for "no allocation in steady state" tests
};
```
Id assignment: every `submit_limit`/`submit_market` call that is not rejected consumes the next id (1, 2, 3, ...), including market orders and IOCs (so fills can name the taker). Rejected submits return `id = kNoOrder` and consume nothing. `modify` keeps the order's id.

### 3.3 Semantics (exact)
Validation (in this order): qty ≤ 0 -> `Rejected/BadQty`; limit price outside [min_tick,max_tick] -> `Rejected/BadPrice`. Rejections change nothing.

Matching of an incoming buy with limit price P (market buy: P = +infinity): while remaining > 0 and the ask side is non-empty and best_ask ≤ P: take the FIFO head order of level best_ask, trade q = min(remaining, maker.remaining) at price best_ask, emit `Fill`, decrement both; if maker.remaining == 0, unlink it, free its slot, mark id dead, `maker_done = true`; if level empties, clear its bitmap bit and recompute best_ask. Sell is symmetric (best_bid ≥ P).

After matching:
- GTC limit: remainder > 0 -> rest at tail of its level. Status: `Filled` (remainder 0), `PartialResting` (filled > 0, rests), `Resting` (filled 0).
- IOC limit and market: remainder is cancelled. Status: `Filled`, `PartialCancelled`, `Cancelled` (nothing filled; includes market order into an empty side).
- `cancel(id)`: unknown or dead id -> `Rejected/UnknownId`; else remove, `Cancelled`, `cancelled = remaining`.
- `modify(id, p, q)`: unknown -> `Rejected/UnknownId`; q ≤ 0 -> `BadQty`; p out of band -> `BadPrice` (order unchanged in all rejections).
  If p == price and q ≤ remaining: reduce in place, keep time priority, status `Resting`, `cancelled = old - q`.
  Otherwise: remove the order and re-submit it as a GTC limit with the same id, side, owner, price p, qty q — it goes to the back of the queue and may trade immediately (fills name it as taker). The returned result is that of the re-submission (`cancelled` = 0).
Post-condition of every call: the book is not crossed (best_bid < best_ask when both exist).

### 3.4 Data structures (decided)
- Order pool: `std::vector<Node>` with `struct Node { OrderId id; Qty remaining; Tick price; OwnerId owner; Side side; std::uint32_t prev, next; }`, free list through `next`, `kNil = UINT32_MAX`. Grows by doubling only if exhausted (allowed but should not happen in benchmarks).
- Levels: per side `std::vector<Level>` of size `max_tick-min_tick+1`, `struct Level { std::uint32_t head, tail; Qty qty; std::int32_t count; }`. Index = price - min_tick.
- Non-empty bitmap per side: `std::vector<std::uint64_t>`; best price search scans words with `std::countr_zero` (asks, upward) / `std::countl_zero` (bids, downward) starting from the previous best. Best prices cached as `std::int64_t` with sentinel (-1 = empty).
- Id table: `std::vector<std::uint32_t> slot_of_id` indexed by id (dense because ids are sequential), `kNil` when dead; grows by doubling.
- Fills buffer: `std::vector<Fill>` reserved to 1024, cleared at the start of each call.
Complexities: add O(1) + best update O(1); cancel O(1) amortised (+bitmap scan when a level empties: O(band/64) worst case); match O(#fills + levels crossed).

### 3.5 Reference book (`reference_book.hpp`)
Naive, obviously-correct model for tests: `std::map<Tick, std::deque<Order>>` per side, same API and identical id assignment rules and semantics, no performance concern. Used by the randomized differential test and by the benchmark op generator.

## 4. Hawkes process

### 4.1 Model
D-dimensional, exponential kernels (Hawkes 1971):
λ_i(t) = μ_i + Σ_j Σ_{t^j_k < t} α_ij exp(-β_ij (t - t^j_k)),  μ_i > 0, α_ij ≥ 0, β_ij > 0.
α_ij is the effect of a type-j event on type i. Branching matrix G_ij = α_ij / β_ij (= ∫ kernel). Stationary iff ρ(G) < 1; then the mean rates are Λ = (I - G)^{-1} μ.

```cpp
struct HawkesParams { int dim; std::vector<double> mu; qc::Mat alpha, beta; void validate() const; };  // throws on sizes, mu<=0, alpha<0, beta<=0, rho(G)>=1
qc::Mat branching_matrix(const HawkesParams&);
std::vector<double> stationary_rates(const HawkesParams&);
// mu = (I-G) lambda_bar, alpha = G .* beta (elementwise). Throws if rho(G)>=1 or any mu_i <= 0.
HawkesParams params_from_rates(const qc::Mat& G, const qc::Mat& beta, std::span<const double> lambda_bar);
// "Calibrate the branching ratio": G = G_shape * (target / rho(G_shape)), then params_from_rates. Requires 0 <= target < 1.
HawkesParams with_branching_ratio(const qc::Mat& G_shape, const qc::Mat& beta, std::span<const double> lambda_bar, double target);
```

### 4.2 Spectral radius (`qc::spectral_radius_nonneg`)
For a non-negative matrix, ρ(G) = lim ||G^m||^{1/m}. Algorithm: if G == 0 return 0. Let `logscale = 0`, `M = G / s` with s = max|G_ij|, `logscale = log s`. Repeat 40 times: `M = M*M; logscale *= 2; t = max|M_ij|; if t == 0 return 0; M /= t; logscale += log t`. After k squarings ||G^{2^k}||_max = exp(logscale), so return `exp(logscale / 2^40)`. (Accurate to ~1e-10 relative even for Jordan blocks of size ≤ 8 because the polynomial prefactor's 2^-40-th root is ~1.) Throw if G has negative entries.

### 4.3 State and exact simulation (Ogata 1981, "On Lewis' simulation method for point processes")
```cpp
class HawkesState {
 public:
  explicit HawkesState(const HawkesParams&);
  void decay_to(double t);        // precondition t >= time(); S_ij *= exp(-beta_ij*(t - time()))
  void add_event(int j);          // S_ij += alpha_ij for all i, at time()
  double intensity(int i) const;  // mu_i + sum_j S_ij
  double total_intensity() const;
  double time() const;
  double compensator(int i, double dt) const;   // mu_i*dt + sum_j S_ij*(1-exp(-beta_ij*dt))/beta_ij (no events in (t,t+dt])
 private: HawkesParams p_; std::vector<double> S_; double t_ = 0;   // S_ row-major D*D, starts at 0 (empty history)
};
struct HawkesEvent { double t; int type; };
std::vector<HawkesEvent> simulate_hawkes(const HawkesParams&, double T, qc::RngStream& rng);
```
Thinning (exact because intensities never increase between events):
```
state.t = 0
loop:
  lbar = state.total_intensity()             # bound valid on [t, next event)
  t_new = state.time() + rng.exponential() / lbar
  if t_new > T: break
  state.decay_to(t_new)
  lam = state.total_intensity()
  if rng.uniform() * lbar <= lam:            # accept
     u = rng.uniform() * lam; i = first index with cumsum(intensity(0..i)) >= u (i = D-1 if rounding)
     events.push_back({t_new, i}); state.add_event(i)
```
Draw order is fixed: exponential, uniform (accept), uniform (type, only if accepted).

### 4.4 Likelihood, time-rescaling residuals
log L = Σ_k log λ_{d_k}(t_k-) - Σ_i Λ_i(0,T), computed in one pass with a `HawkesState` (decay to t_k, read intensity, add event; compensator accumulated with `compensator(i, dt)` between events and up to T).
`rescaled_interarrivals(params, events, T)` -> per dimension the compensator increments between consecutive events of that dimension (first one measured from 0). Under the true model they are i.i.d. Exp(1) (time-rescaling theorem; Brown, Barbieri, Ventura, Kass, Frank 2002).

### 4.5 EM calibration with β fixed (Veen & Schoenberg 2008 style, recursive for exponential kernels)
```cpp
struct EmOptions { int max_iter = 500; double rel_tol = 1e-10; };
struct EmResult { HawkesParams params; std::vector<double> loglik_trace; int iterations; bool converged; };
EmResult fit_em(std::span<const HawkesEvent> events, double T, const qc::Mat& beta, const HawkesParams& init, const EmOptions&);
```
Preconditions: events sorted by time, all in (0,T], each dimension has ≥ 10 events (else throw).
Per iteration (all O(K·D²), K = #events):
```
A_ij = 0 (A_ij(t) = sum over past type-j events of exp(-beta_ij (t - t_l)))
num_mu[i]=0, num_a[i][j]=0
for each event k (time t_k, type i):
   decay A to t_k
   lam = mu_i + sum_j alpha_ij A_ij
   num_mu[i] += mu_i / lam                      # P(background)
   for j: num_a[i][j] += alpha_ij * A_ij / lam  # expected #children of type-j parents
   for all r: A_{r,i} += 1                      # after reading, the event becomes history
C_ij = sum over type-j events l of (1 - exp(-beta_ij (T - t_l))) / beta_ij
mu_i  <- num_mu[i] / T
alpha_ij <- num_a[i][j] / C_ij                  # if C_ij == 0 keep 0
```
Log-likelihood is evaluated before every update and appended to `loglik_trace`; stop when |ℓ_k - ℓ_{k-1}| ≤ rel_tol·|ℓ_k| (converged = true) or at max_iter. EM never decreases ℓ (tested). Default init: μ_i = 0.5·N_i/T, α_ij = 0.1·β_ij.

## 5. Market simulator

### 5.1 Background flow model (decided defaults)
Hawkes dimension D = 4: 0 = LimitBuy (LB), 1 = LimitSell (LS), 2 = MarketBuy (MB), 3 = MarketSell (MS).
Default branching matrix (rows = affected type, cols = source type), all β_ij = 2.0 s⁻¹:
```
        LB   LS   MB   MS
LB   [ 0.3  0    0    0.3 ]   # MS excites LB: replenishment of the consumed bid side
LS   [ 0    0.3  0.3  0   ]   # MB excites LS
MB   [ 0    0    0.5  0   ]   # market-order self-excitation (order splitting / herding)
MS   [ 0    0    0    0.5 ]
```
ρ(G) = 0.5 (block triangular). Target stationary rates Λ = (1.0, 1.0, 0.25, 0.25) s⁻¹, hence μ = (I-G)Λ = (0.625, 0.625, 0.125, 0.125) and α = 2G.
Cancellations are not Hawkes events: every resting background order has an independent cancellation hazard ν = 0.01 s⁻¹ (Santa-Fe/zero-intelligence style; Smith, Farmer, Gillemot, Krishnamurthy 2003), implemented as one exponential clock with rate ν·n_live and a uniformly chosen live background order.
Order attributes: size uniform on {1,2,3,4}. LB price = ref_ask − k, LS price = ref_bid + k, k uniform on {1..K}, K = 20, where ref_ask = best ask if present else the last observed best ask (symmetrically ref_bid). Background limit orders therefore never cross. Prices clamped into the band.
Initial book: band [0, 20000]; bids at 9999, 9998, ..., 9980 and asks at 10001, ..., 10020 (20 levels each), 3 orders of 3 lots per level, all owner `kBackground`, added in order best-to-worst, bids before asks.
Expected stationary background book: ~75 orders per side (inflow 1.0/s, outflow 0.25/s by market orders + ν·N), relaxation time 1/ν = 100 s. Burn-in 600 s before any measurement.

```cpp
struct SimConfig {
  qc::Mat G = default_G();  double beta = 2.0;  std::array<double,4> lambda_bar{1.0, 1.0, 0.25, 0.25};
  double cancel_hazard = 0.01;  int placement_K = 20;  Qty size_min = 1, size_max = 4;
  Tick min_tick = 0, max_tick = 20000, init_best_bid = 9999, init_best_ask = 10001;
  int init_levels = 20, init_orders_per_level = 3;  Qty init_order_size = 3;
  bool agent_market_orders_excite = true;   // agent MB/MS are injected into the Hawkes state as type 2/3 events
  void validate() const;
};
```
`agent_market_orders_excite = true` is a modelling choice (the market reacts to the metaorder: MB excites MB and LS). It stays true for all reported experiments; flagged in README.

### 5.2 Agents
```cpp
class MarketSim;
class Agent {
 public:
  virtual ~Agent() = default;
  virtual double next_wakeup() const = 0;                 // absolute time; +inf = never
  virtual void on_wakeup(MarketSim& sim) = 0;             // may place/cancel orders via sim.agent_*
  virtual void on_fill(const Fill& f, double t) { (void)f; (void)t; }   // fills where the agent is taker or maker
};
```
`MarketSim::add_agent(Agent*)` (non-owning) assigns OwnerId = 1 + index. Agents with equal wakeup times run in index order.

### 5.3 MarketSim API
```cpp
class MarketSim {
 public:
  MarketSim(const SimConfig&, std::uint64_t seed, std::uint64_t path_id);
  OwnerId add_agent(Agent*);
  void run_until(double t_end);                 // processes all events with time <= t_end, then sets time = t_end
  double time() const;  const OrderBook& book() const;
  double mid() const;                           // (bid+ask)/2 if both sides exist, else last valid mid (initial 10000.0)
  std::optional<Tick> best_bid() const, best_ask() const;
  double hawkes_intensity(int type) const;
  OrderResult agent_market(OwnerId, Side, Qty);
  OrderResult agent_limit(OwnerId, Side, Tick, Qty);  // GTC
  OrderResult agent_cancel(OwnerId, OrderId);         // rejects (UnknownId) if the order is not owned by this agent
  struct Counters { std::array<std::int64_t,4> hawkes_events{}; std::int64_t cancels = 0, agent_orders = 0, fills = 0; };
  const Counters& counters() const;
};
```
RNG streams (seed = experiment seed, stream_id = path_id): purpose 1 = Hawkes proposal exponentials, 2 = Hawkes accept/type uniforms, 3 = order attributes (k then size), 4 = cancel-clock exponentials, 5 = cancel selection; agents use purposes ≥ 16 (each agent class documents its purpose ids).

Event loop (exact; re-drawing both exponential clocks every iteration is valid by memorylessness, and the current total intensity is a valid thinning bound because intensities only decay between events):
```
run_until(t_end):
  loop:
    lbar = hawkes.total_intensity()                      # at current time
    t_h  = t + E1/lbar                                   # E1 from purpose 1
    n_c  = bg_live.size(); t_c = n_c>0 ? t + E4/(nu*n_c) : +inf   # E4 from purpose 4, drawn ONLY when n_c > 0
    t_a  = min_i agent_i.next_wakeup() (ignore < t: treat as t)
    t_n  = min(t_h, t_c, t_a)
    if t_n > t_end: hawkes.decay_to(t_end); t = t_end; return
    hawkes.decay_to(t_n); t = t_n
    if t_n == t_a (agent minimum, ties: agent first, then cancel, then hawkes): run agent wakeup(s) whose time == t_a in index order
    elif t_n == t_c: cancel bg_live[uniform_int(0,n_c-1)] (purpose 5)
    else: thinning test with U from purpose 2: accept iff U*lbar <= hawkes.total_intensity();
          if accepted: choose type (purpose 2), generate the background order (purpose 3), submit, hawkes.add_event(type)
    after every engine call: process last_fills(): remove filled-out background makers from bg_live;
          dispatch fills with an agent owner to that agent's on_fill; update last valid best bid/ask/mid
```
`bg_live` is a vector of background resting order ids plus an id->position map (`std::vector<std::uint32_t>` indexed by id) for O(1) swap-remove.
Agent market orders: after the engine call, if `agent_market_orders_excite`, `hawkes.add_event(2 or 3)`.

## 6. Execution

### 6.1 Almgren-Chriss (Almgren & Chriss 2000, "Optimal execution of portfolio transactions", J. Risk 3(2))
Discrete model, liquidate X over T in N intervals, τ = T/N, t_k = kτ, holdings x_0 = X, ..., x_N = 0, trades n_k = x_{k-1} - x_k.
Permanent impact g(v) = γv; temporary h(v) = ε·sgn(v) + η·v. Define η̃ = η − γτ/2 (must be > 0, else throw).
E[cost] = ½γX² + ε Σ|n_k| + (η̃/τ) Σ n_k²;  Var[cost] = σ² τ Σ_{k=1}^{N} x_k².
For risk aversion λ ≥ 0: κ̃² = λσ²/η̃, κ = acosh(1 + κ̃²τ²/2)/τ, and
x_k = X sinh(κ(T − t_k)) / sinh(κT)   (if κT < 1e-8 use the limit x_k = X(1 − t_k/T), i.e. TWAP).
Compute sinh ratios as `exp(-κ t_k) * (1 - exp(-2κ(T-t_k))) / (1 - exp(-2κT))` to avoid overflow for large κT.
```cpp
struct AcParams { double X, T; int N; double sigma, gamma, eta, epsilon, lambda; };
struct AcSchedule { std::vector<double> x /*N+1*/, n /*N*/; double kappa, expected_cost, variance; };
AcSchedule almgren_chriss(const AcParams&);
AcSchedule almgren_chriss_tridiagonal(const AcParams&); // reference: solve FOC -x_{j-1} + (2+κ̃²τ²)x_j - x_{j+1} = 0, x_0=X, x_N=0 (Thomas algorithm)
double ac_lambda_for_kappaT(const AcParams& base, double kappaT); // κ = kappaT/T; κ̃² = 2(cosh(κτ)-1)/τ²; λ = κ̃² η̃ / σ²
std::vector<Qty> round_schedule(std::span<const double> n, Qty X); // cumulative rounding: c_k = llround(sum_{m<=k} n_m), child_k = c_k - c_{k-1}; sum == X exactly
```
A buy program uses the same trajectory (remaining-to-buy instead of holdings).

### 6.2 Bertsimas-Lo (Bertsimas & Lo 1998, "Optimal control of execution costs", J. Financial Markets 1)
Buy S̄ shares over T periods. State before period t: (P_{t-1}, X_{t-1}, W_t), W_1 = S̄. Choose S_t (unconstrained real; may be negative), then
X_t = ρ X_{t-1} + η_t,  P_t = P_{t-1} + θ S_t + γ X_t + ε_t,  η_t ~ N(0, σ_η²), ε_t ~ N(0, σ_ε²), cost Σ P_t S_t, W_{t+1} = W_t − S_t, S_T = W_T.
Value function V_t(P,X,W) = P W + a_t W² + b_t X W + c_t X² + d_t, with
a_T = θ, b_T = γρ, c_T = 0, d_T = 0, and for t = T−1, ..., 1:
```
a_t = θ − θ² / (4 a_{t+1})
b_t = γρ + θ ρ b_{t+1} / (2 a_{t+1})
c_t = ρ² c_{t+1} − ρ² b_{t+1}² / (4 a_{t+1})
d_t = d_{t+1} + c_{t+1} σ_η²
S*_t = (1 − θ/(2a_{t+1})) W_t + (ρ b_{t+1} / (2a_{t+1})) X_{t−1}    (t < T);  S*_T = W_T
E[optimal cost] = V_1(P_0, X_0, S̄)
```
**Architect's flag:** these recursions were derived here from the model statement (they reduce to BL's well-known γ = 0 result S*_t = W_t/(T−t+1) and E[cost] = P_0 S̄ + θS̄²(1 + 1/T)/2). I do not have BL's printed coefficient formulas at hand; the tests verify them by Monte Carlo against V_1 and against the γ = 0 closed form, which is sufficient. Do not "correct" them to a remembered paper form without a failing test.
TWAP analytic cost in this model: E[P_t] = P_0 + θ S̄ t/T + γ X_0 Σ_{s=1}^{t} ρ^s, E[cost_TWAP] = (S̄/T) Σ_t E[P_t].
```cpp
struct BlParams { double S_bar; int T; double P0, theta, gamma, rho, sigma_eps, sigma_eta, X0; };
struct BlSolution { BlParams p; std::vector<double> a, b, c, d; };   // index 1..T (index 0 unused)
BlSolution bertsimas_lo_solve(const BlParams&);                       // throws if theta<=0, T<1
double bl_policy(const BlSolution&, int t, double W, double X_prev);
double bl_value(const BlSolution&, int t, double P_prev, double X_prev, double W);
double bl_twap_expected_cost(const BlParams&);
struct BlMcResult { qc::RunningStats opt, twap, diff; };              // diff = opt - twap per path (common random numbers)
BlMcResult bl_monte_carlo(const BlParams&, std::int64_t paths, std::uint64_t seed, unsigned threads);
```
MC: one path per block of 1024 paths (stream_id = block, purpose 1: per period draw η then ε, in that order); both policies use the same draws.

### 6.3 Execution inside the simulator
Program: **buy** X = 64 lots over H = 300 s starting at t0 = 600 s, N = 10 intervals (τ = 30 s). At the start of interval k (time t0 + (k−1)τ) the `ScheduleAgent` sends ONE market buy of child_k lots (+ any carried unfilled remainder). At t0 + H any remainder is sent as a final market order; still-unfilled lots are counted in `unfilled` (expected 0).
IS (ticks·lots) = Σ fills price·qty − (X − unfilled)·mid(t0). Reported per lot as IS/X.
Strategies: TWAP (κT = 0), AC κT = 1, AC κT = 3, BL with γ = 0 (must equal TWAP child list exactly — consistency check). BL with information is evaluated only in its own model (§6.2): at a 30-s interval the Hawkes intensity signal has decorrelated (kernel half-life 0.35 s), so an in-simulator BL-info strategy would be meaningless; README must say this.
AC parameters estimated from the simulator (`calibrate_ac_params`, seed 7001):
- σ²: 400 agent-free paths, mid sampled every 10 s on [600, 1200]; σ² = pooled Var(Δmid_10s)/10 (ticks²/s).
- ε, η: `ProbeAgent` on 200 paths sends market orders at t = 600 + 15m, m = 0..39, alternating buy/sell, size v_m = 1 + (m mod 8); slippage per lot = side·(VWAP − mid_before); OLS slippage = ε + c·v; η = c·τ.
- γ: 500 paths running TWAP for X = 64; γ = mean(mid(t0 + 2H) − mid(t0)) / X.
- Throw (with the estimates in the message) if η̃ = η − γτ/2 ≤ 0 or ε < 0.
Evaluation: R = 2000 paths per strategy, seed 7002, path p uses stream_id p for every strategy (common random numbers). Report mean, sd, SE of IS per lot and of the paired difference vs TWAP, plus the AC model's predicted E and sqrt(Var).

## 7. Avellaneda-Stoikov (Avellaneda & Stoikov 2008, "High-frequency trading in a limit order book", Quant. Finance 8(3))
Reservation price r = s − q γ σ² (T − t); total spread δ_a + δ_b = γσ²(T − t) + (2/γ) ln(1 + γ/k); quotes r ± spread/2. Fill intensity λ(δ) = A e^{−kδ}.
```cpp
struct AsQuote { double bid, ask; };
AsQuote as_quotes(double s, double q, double t, double T, double gamma, double sigma, double k);   // inventory strategy
double as_symmetric_spread(double gamma, double sigma, double k, double T); // (2/γ)ln(1+γ/k) + γσ²T/2  (time-average of the AS spread)
```
### 7.1 Own model (paper reproduction)
Defaults (the paper's): s0 = 100, T = 1, σ = 2, dt = 0.005 (200 steps), k = 1.5, A = 140, q0 = 0, γ ∈ {0.01, 0.1, 0.5}.
Per step (t = i·dt): compute quotes from s at step start (Symmetric strategy: s ± as_symmetric_spread/2); δ_a = ask − s, δ_b = s − bid; draw u_a, u_b, u_s (purpose 1, in that order): if u_a < min(1, A e^{−kδ_a} dt) sell 1 (q −= 1, cash += ask); if u_b < min(1, A e^{−kδ_b} dt) buy 1 (q += 1, cash −= bid); s += σ√dt·(u_s < 0.5 ? +1 : −1). PnL = cash + q·s_T. 10 000 paths, seed 8001, stream = path block (1024 paths per block), both strategies on identical draws.
**Architect's flag:** AS 2008 Table (1000 sims, γ = 0.1) is recalled as: inventory strategy profit 62.94 (sd 5.89), final q sd 2.80; symmetric 67.21 (sd 13.43), q sd 8.66; and the symmetric strategy's spread definition in the paper is recalled as "the average spread of the inventory strategy". Both are from memory; tests treat them as soft (15%) checks and make the qualitative ordering the hard check.

### 7.2 In the simulator (`AsAgent`)
Session [600, 1200] s, T = 600. Every 1 s: cancel both own quotes, compute r and spread from s = sim.mid(), q, t; bid_px = floor(r − spread/2), ask_px = ceil(r + spread/2); clamp bid_px ≤ best_ask − 1 and ask_px ≥ best_bid + 1 (when those sides exist); post 1 lot each side, except: no bid if q ≥ +20, no ask if q ≤ −20. Final PnL = cash + q·mid(1200) (ticks·lots). Symmetric agent: same spread, centred at mid.
Parameters from the simulator (`calibrate_as_params`, seed 8002): σ² as in §6.3; A, k: probe agent quoting 1 lot at fixed distance δ ∈ {1,...,6} ticks from mid on both sides (bid = floor(mid − δ), ask = ceil(mid + δ)), refreshed every 1 s, 100 paths per δ; λ(δ) = fills / (2 · 600 s); OLS ln λ = ln A − kδ. γ := 0.4 / (σ² · T) (so the reservation shift at t = 0 is 0.4 ticks per lot of inventory).
Evaluation: R = 1000 paths, seed 8003, CRN between AS and symmetric. Report mean/sd PnL, sd(q_T), mean |q| over time, fill counts.

## 8. Impact validation
Metaorder: `ScheduleAgent` buying Q lots with 1-lot market orders at times t0 + (m + 0.5)·H/Q, m = 0..Q−1, t0 = 600 s, H = 300 s. `ObserverAgent` records mid at t0 + j·H/10, j = 0..20.
Q ∈ {8, 16, 32, 64, 128} plus control Q = 0; R = 2000 paths each, seed 9001, path p uses stream_id p for every Q (CRN).
Impact I(Q) = mean_p[(mid_p,Q(t0+H) − mid_p,Q(t0)) − (mid_p,0(t0+H) − mid_p,0(t0))]; SE from the paired differences.
Fit: OLS of log I on log Q over points with I > 4·SE (need ≥ 4 such points, else the fit is reported as "insufficient"). Report δ̂ (slope), SE(δ̂), R², and Y in I = Y·σ_H·sqrt(Q/V_H), where σ_H = sd of the control's mid change over H and V_H = mean background market-order volume (both sides) in [t0, t0+H] of the control paths.
Also output the impact profile (mean paired mid change at j = 0..20) for every Q (shows decay after t0 + H).
Mechanism (why concavity is expected, not guaranteed): a depth profile that grows away from the touch (latent liquidity, Tóth et al. 2011 "Anomalous price impact and the critical nature of liquidity") and replenishment by MB→LS excitation and fresh placements. **Flag:** emergence of δ ≈ 0.5 in this zero-intelligence + Hawkes model is a hypothesis to test. Allowed sweep if the default fails the band in TESTPLAN: ν ∈ {0.005, 0.01, 0.02} × K ∈ {10, 20, 40} (9 configs, R = 1000), report the full table; the README then states the measured exponents. No other tuning.

## 9. Engine benchmark (`bench/ofl_bench_engine`)
- Build preset `release` (-O3, no -march=native). Single thread. Not pinned (state this).
- Op stream generated once (seed 42) by running the generator against a `ReferenceBook`-free fast path: use the `OrderBook` itself during generation to know live ids/best prices, then record ops; replay happens on a fresh `OrderBook` and is deterministic (same ids).
- Band [0, 65535], start mid 32768, pre-fill 2000 orders (uniform offsets 1..20 from 32768 on each side, size 1..10).
- Mix per op: 55% GTC limit (90% passive: price = opposite best ∓ uniform{1..20}; 10% aggressive: crosses by uniform{0..2} ticks), 30% cancel (uniformly random live order), 10% market (size 1..4), 5% modify (random live order; 50% reduce qty in place, 50% move price by ±uniform{1..3}). Limit sizes uniform 1..10. If no live orders exist, a cancel/modify becomes a passive limit.
- `--ops=5000000 --warmup=1000000 --reps=5`. Each rep replays warmup+ops on a fresh book constructed with capacity for all orders (no growth during replay — verified by `pool_capacity()` unchanged).
- Throughput: wall time of the measured ops (steady_clock around the loop, after warmup), median and min/max over reps; ns/op = 1e9/(ops/s).
- Latency: one extra rep timing every op individually with `now_ns()` before/after, stored in a preallocated array; report p50/p90/p99/p99.9/max overall and per op type. Also report timer overhead = median of 1e6 back-to-back `now_ns()` pairs; it is NOT subtracted.
- Also report: fills per op, mean book levels, CPU model (`/proc/cpuinfo` "model name"), compiler (`__clang_version__`), flags. Output `results/bench_engine.csv` (rep, ops, seconds, ops_per_s, ns_per_op) and `results/bench_latency.csv` (op_type, p50, p90, p99, p999, max, count).
- `DoNotOptimize`-style sink: accumulate `result.filled + result.resting` into a volatile sink.

## 10. Threading and determinism
- Engine, Hawkes simulation and MarketSim are single-threaded objects.
- Experiments run paths in parallel with `qc::parallel_blocks` (one block = one path for market experiments; 1024 paths per block for own-model MC). Each block writes into `results[block]`; summaries are computed afterwards in block order.
- Every stochastic quantity is a function of (seed, stream_id, purpose, draw index). Output files are bitwise identical for any `QC_THREADS`.

## 11. Apps (all `--key=value`, defaults = the numbers above; print summary + write CSV)
- `ofl_hawkes_fit --T=100000 --seed=1001` -> `results/hawkes_fit.csv` (param, true, fitted, rel_err) for the 2-D test model of TESTPLAN §3.
- `ofl_depth --paths=200 --seed=9100` -> `results/depth_profile.csv` (distance_ticks, mean_bid_depth, mean_ask_depth) sampled every 10 s on [600,1200].
- `ofl_impact --paths=2000 --seed=9001 [--sweep=1]` -> `results/impact.csv` (Q, I, se, n) , `results/impact_profile.csv` (Q, j, time, mean, se), `results/impact_fit.csv` (delta, se_delta, r2, Y, n_points), sweep -> `results/impact_sweep.csv`.
- `ofl_execution --paths=2000 --seed=7002` -> `results/ac_params.csv`, `results/execution.csv` (strategy, mean_is_per_lot, sd, se, mean_diff_vs_twap, se_diff, model_E, model_sd, unfilled_total), `results/bl_model.csv`.
- `ofl_as --paths=1000` -> `results/as_model.csv` (gamma, strategy, mean_pnl, sd_pnl, mean_q, sd_q, mean_spread), `results/as_sim.csv`.
