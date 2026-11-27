# orderflow-lab

### A tiny stock exchange, built from scratch in C++ — plus the trading algorithms that run on it, and an honest look at what a simulated market can and can't tell you.

```
 random buy/sell orders          the exchange                 the traders
┌─────────────────────┐   ┌──────────────────────┐   ┌─────────────────────────┐
│  Hawkes process     │──▶│  Order book          │◀──│  Execution algorithm    │
│  (orders arrive in  │   │  matches buyers and  │   │  (buy 64 lots without   │
│   bursts, like real │   │  sellers by price    │   │   moving the price)     │
│   markets)          │   │  then arrival time   │   │  Market maker           │
└─────────────────────┘   │  ~8 million orders/s │   │  (quote both sides,     │
                          └──────────────────────┘   │   manage risk)          │
                                                     └─────────────────────────┘
```

## In plain English

Every stock exchange keeps a list of people who want to buy (at the price they'll pay) and people who want to sell. When a
buyer and a seller agree on a price, a trade happens. That list is the **order book**, and the software that pairs people up is
the **matching engine**.

This project builds all of it:

1. **The exchange.** A matching engine that handles thousands of orders per millisecond, and is checked against a deliberately
   simple "obviously correct" version on randomized streams of orders.
2. **A realistic crowd.** Real markets are bursty: one trade makes the next one more likely. The simulator uses a mathematical
   model for that (a *Hawkes process*), so the fake market behaves less like coin flips and more like the real thing.
3. **The traders.** Two classic problems every trading desk faces:
   - *"I have to buy a huge amount without pushing the price up against myself."* (optimal execution)
   - *"I want to quote prices to everyone and earn the spread, without ending up stuck holding too much stock."* (market making)
4. **The reality check.** Does a simulated market push prices the way real markets do when someone trades big? This project
   *measures that* — and reports the answer even though it was not the one hoped for.

## The numbers (all measured, all reproducible)

| What | Result |
|---|---|
| Matching engine speed | **~8 million orders per second**, median **124 nanoseconds** per order, on one CPU thread |
| Is it correct? | Identical results to a naive reference book on randomized operation streams (the benchmark re-checks 200,000 operations on every run) |
| Order-flow model | Recovers its own hidden parameters to within **1–2%** |
| Smart buying vs the naive way | The Almgren–Chriss schedule narrowed the spread of outcomes by **23%** (std dev 3.17 vs 4.09 ticks) and cost slightly less: 5.11 vs 5.46 ticks per lot |
| Market maker | Inventory-aware quoting **halved** the spread of how much stock it ended up holding (std dev 4.0 vs 8.1 lots) |
| Tests | **20 tests pass**; sanitizer builds (undefined-behaviour, thread) are supported by `./build.sh ubsan --test` / `tsan` |

## The part most projects leave out: what didn't work

A good simulator should reproduce a famous real-world pattern: the **square-root law**, where doubling the size of a big order
raises the price impact by only about 40%, not 100%. I tested for it with a pre-declared pass/fail range.

**It failed.** In this model, impact grew in a straight line (exponent 0.997, not ≈ 0.5). Instead of tuning until it passed, I
tracked down why: the simulated market makers re-priced themselves off the new price after every trade, so nothing ever pulled
the price back. A fix that anchors them to a hidden "fair price" removes that artefact, but still doesn't produce the
square-root law — that needs more structure than this model has. The full investigation is in
[`orderflow-lab/NOTES.md`](orderflow-lab/NOTES.md). Reporting a negative result is the point: a simulator you can't falsify
isn't telling you anything.

Two more honest findings: the market-making strategy that manages risk *earned less* than the naive one in this simulator
(62.6 vs 78.0) — it traded profit for safety, as expected — and one of the execution methods (Bertsimas–Lo) collapses to the
simple equal-slices baseline here, which the README says rather than hiding.

## Try it

```bash
# no compiler needed: this uses zig (a C++ compiler) from PyPI
uv venv .toolchain/.venv
uv pip install --python .toolchain/.venv/bin/python ziglang cmake ninja

./build.sh release --test                         # build + run all 20 tests (QC_JOBS=4 on small machines)

build/release/orderflow-lab/ofl_bench_engine --verify=1   # speed test, with correctness check
build/release/orderflow-lab/ofl_execution                 # smart buying vs naive buying
build/release/orderflow-lab/ofl_as                        # the market maker
build/release/orderflow-lab/ofl_impact                    # the price-impact experiment
```

What `ofl_execution` prints (real output; shortfall = extra cost per lot in price ticks, lower is better):

```
strategy           IS/lot        sd   diff vs TWAP
TWAP               5.4579    4.0894         0.0000      <- naive: equal slices
AC_kappaT1         5.3675    3.9035        -0.0904
AC_kappaT3         5.1118    3.1655        -0.3461      <- smarter: front-loaded, 23% less risk
```

## What's inside

- **Order book** — integer prices, a queue per price level, a bitmap that finds the best price instantly. Supports limit, market,
  IOC, cancel and modify orders.
- **Hawkes order-flow simulator** — four order types, exact simulation (Ogata thinning), parameters fitted back from the data.
- **Execution** — Almgren–Chriss (closed form) and Bertsimas–Lo (dynamic programming, from MIT research), measured against the
  naive equal-slices schedule.
- **Market maker** — Avellaneda–Stoikov, reproduced in its original model and then run inside the simulator.
- **Reproducible by design** — same seed, same answer, on 1 thread or 16.

Built on research from Oxford (the LOB-Impact benchmark idea), MIT (Bertsimas–Lo), UCL, and the classic papers by Hawkes,
Almgren–Chriss and Avellaneda–Stoikov. Full list in [`RESEARCH.md`](RESEARCH.md).

## Go deeper

[Overview of both projects](docs/OVERVIEW.md) ·
[Full technical report](orderflow-lab/README.md) ·
[Design](orderflow-lab/DESIGN.md) ·
[Test plan](orderflow-lab/TESTPLAN.md) ·
[Measured notes & deviations](orderflow-lab/NOTES.md) ·
sibling project: [roughvol](https://github.com/swayam555/roughvol)

*A research and learning project, not a trading system. Nothing here is investment advice.*
