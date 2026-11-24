> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Where the original quintet leaves physical slack

**Best support to price next:** negative **H3** charged rows with middle arity **2–4** and maximum tip degree **19–1998**, leaving the other tip unrestricted. Excluding founder-containing sets and the previously priced strict H3 masks, this accounts for **5.772 / 7.138 × 10⁻⁶ per core** (lower/upper) on the pooled age-20 trees and **5.727 / 7.009 × 10⁻⁶** on the separate tree. This is an empirical allocation, not a stationary lower bound or a statistical exponent estimate.

## Physical allocation

Eight preserved age-20 checkpoints, seeds `481250000`–`481250007`, contributed **896,887 cores**. The separate age-20 seed `481230001` has **235,709 cores among 1,275,795,845 vertices**. The `[parent, ordinary leaf count]` files retain the actual degrees and ancestor labels; no new path was generated. Per-tree sizes and SHA-256 hashes are in [report.log](report.log).

All entries in the table are in **10⁻⁶ normalized potential units per core**. Net credit is minus the signed outward charged action, evaluated at the original witness trial. The right columns show a *subaccount* of the credit, obtained by replacing actual uncovered donor terms with their upper boxes or their FY maxima; they must not be added to the net column.

| Family | Net lower | Net upper | Donor-gap subaccount lower | Donor-gap subaccount upper |
|---|---:|---:|---:|---:|
| f | 1.596 | 2.254 | — | — |
| g | 6.923 | 6.816 | 6.297 | 6.004 |
| H3 | 8.994 | 10.701 | 6.314 | 7.347 |
| P4 | 1.558 | 1.568 | 0.969 | 0.826 |
| S3 | 1.425 | 2.656 | 0.621 | 0.855 |
| P5 | 0.404 | 0.561 | 0.352 | 0.529 |
| B5 | 1.743 | 1.481 | 1.561 | 1.272 |
| S4 | 0.827 | 1.093 | 0.718 | 0.743 |
| **Total net** | **23.470** | **27.130** | | |

The pooled potentials per core are `1.111151 / 1.111307`. Dividing the net physical credits by these observed denominators gives **21.122 / 24.413 ppm of potential** inward of the nominal trials on these trees. Measured from the *certified*, cap-adjusted endpoints the corresponding diagnostic gaps are **21.314 / 24.540 ppm**. The original global cap budgets themselves are only **0.01745 / 0.01799 × 10⁻⁶ per core**; subtracting or repricing caps alone therefore cannot explain the empirical difference in numerators. Across the eight trees the unnormalized net per-core values span **23.153–24.189 / 26.755–27.965** in table units; held-out values are **23.709 / 27.181**. None of these ratios was fit in age or used to estimate an exponent.

## Support, direction and overlap

- **H3 is the only single higher family above 5 on both witnesses.** In the priority support, split by max tip degree: `19–37`: **1.786 / 1.925**; `38–179`: **2.803 / 3.247**; `180–1998`: **1.183 / 1.966**. Capping the other tip at 11 (and allowing middle arity up to 8) leaves only **3.547 / 3.830**. A support that excludes both moderate tips or drops the long shoulder will miss much of the available numerator. The priority cut excludes the strict nonfounder chain mask (H3 ends up to 65 with middle up to 5, terminal core arity 1 and the other two arity 2); that covers its stricter old/virgin subslice too. A separate conservative high-middle exclusion is tracked in the log; it cannot intersect the priority arity cut. The larger `max tip 12–1998, middle c≤8` cut measures **6.973 / 8.506** with the same exclusions. These cuts are nested or overlapping; choose one pricing ledger.
- **Most uncovered H3 donor slack comes from children:** **4.910 / 5.977**, versus **1.404 / 1.370** from the mother. This is for the whole H3 family, not already restricted to the priority mask. The largest individual bins have child donor degree `2–5`, receiving degree `12–37`, and receiving arity `2–8`. This is a subaccount of the full charged H3 reward; an occupation proof for the full negative row should not also claim its donor gap. Nor does this measurement certify that the existing projected row caps can be debited by that entire amount.
- **g is large but overlaps existing payments.** Total FY knot gap is **6.102 / 5.847**. Outside the named full edge-18, main FY and older-donor masks and (lower only) the older-payer extension, the measured remainder is **3.913 / 2.307**. The older-donor exclusion uses the support of the inspected `90/5/24` default cut, irrespective of whether that exact run has been booked. The largest lower remainder is when the mother donates to a host with a payer of degree ≤5 and the mother above 5 (**1.478**), then a younger child payer ≤5 with older sibling above 5 (**0.941**). Upper is more spread: mother donates with both payer and donor above 5 (**0.547**) or payer ≤5 (**0.428**). The full edge-18 cap residual measures only **0.000171 / 0.000365**. Removing FY and plain donor gaps from the complete physical edge cap slack leaves **0.626 / 0.813**; most of it lies at endpoint degrees `2–5` and `180–1998`. This last remainder is a diagnostic contrast, not a replay of every projected/tail cap. Other post-FY or exchange calculations may hit the same edge slack, and are not subtracted as if they were separate accounts.
- **Five-vertex-only is below target:** P5+B5+S4 total negative charged credit is **2.973 / 3.135**, including all degrees and arities. The held-out tree gives **3.110 / 3.261**. B5 with maximum tip degree `38–179` takes **0.463 / 0.392**, and `180–1998` takes **0.402 / 0.328**; S4 with max tip `180–1998` takes **0.239 / 0.338**. The largest B5 donor bin is a child donating at degree `2–5` to a receiver of degree `2–5` and arity `2–4` (**0.978 / 0.741**). That bin is nested inside the five-vertex total, not an extra. These observations do not give an all-age ceiling.
- **f is already the singleton account.** Its largest measured bin is `D=180–1998, c=5–8` (**0.786 / 0.978**); actual analytic states `D=2000–2011` are tiny in this census. Do not add this measured amount to a lifetime certificate for the same row. The earlier raw joined-edge and five-core seam rewards are alternative event ledgers, and have no claimed disjointness with these charged rows.

## Ledger and exact controls

For every supported unoriented embedding, take the physical gain increment, old-set promotion and the newborn extension after promotion at the updated arity, plus the original signed trial subtraction. The singleton also owns its newborn singleton. Then debit the actual centered term `(1/D_donor−r1) ΔW` from each sender and credit the unique supported extension. At exact G2 outside M3, use the original FY interpolation and move `A(1/D_donor)` between distinct paying edges instead. A physical FY gap compares `σ(x Δg−A(x))` at the actual donor reciprocal to its maximum over the **original** knots, using the actual interpolated payer slope. No new knot grid, donor projection, or independent-neighbor law enters. Uncovered plain donor gaps compare with the literal `[0,1/2]` box; their bin ranking is diagnostic, not a certified replacement for the production envelopes.

`census.cpp` uses the original wire layout on read-only mmaps and enumerates actual embeddings, with true tail/overflow marks and physical founder rates. Across all nine trees and both signs, the charged total, uncoupled raw total, and event-aggregated generator agreed to **`1.57e−12` absolute normalized units** (before dividing by core count). `exact.py` independently rebuilds literal ordinary trees and directed-arrow rates, enumerates connected sets via brute combinations, reevaluates births, and checks exact cancellation with fractions. Four targeted fixtures per side validate every family, all eight transfer channels and 24 actual FY slots per side. Exact totals are in [exact.log](exact.log). In `fixture3.tsv` the nonfounder H3 `(1,9,10)` has marks `((179,4),(8,3),(2,1))`; the lower/upper charged ticks are respectively `−13014785772891747344502119/10079227084800000000` and `−53369737938296112185593/241901450035200000`. This is a control for one occurrence; its incidence mass is not estimated by the fixture.

Replay from `.`:

```sh
bash round2/occupation-quintet/diagnostic/replay.sh
```

The replay pins `b04cc50c1e391e0420eb7b006ac781d484c68543` and the gold receipt and checks payload stat identities. Payload authentication is inherited from the prior cold audit. It uses one CPU (`CPU=1` by default), nice +19 and bounded subprocesses; the sample passes used **10.07 / 10.32 user CPU seconds**, **10.47 / 11.50 s wall**, and approximately **1.05 GiB peak RSS**. The exact replay used **4.40 user CPU seconds**, exit 0. The detailed sample bins and all new artifacts are local to this directory.
