# orderflow-lab NOTES

## A1 - OrderBook
Done per DESIGN 3. Matching, cancel, modify, bitmap best-price scans, check_invariants. TESTPLAN 1 passes (release, ubsan).
Convention: rejected cancel/modify return `id` = the id passed (rejected submits return kNoOrder).
