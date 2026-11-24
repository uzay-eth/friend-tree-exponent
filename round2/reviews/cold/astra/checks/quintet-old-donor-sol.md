# Cold audit: complementary old-donor FY

**Pass on both original signs.** The matched `100/5/36; 7/5/6` cut completed all four passes below. All values are already normalized **per physical core**. The separate `170/5/70; 12/5/7` lower-witness floor also reproduced exactly; there is no matched ceiling or upper-witness result for that wider cut in this audit.

| Cut | Witness | Bound | Exact total |
|---|---|---|---|
| Matched | Lower | Floor | `49255854388977166821869/109951162777600000000000000000` |
| Matched | Lower | Ceiling | `24510597561147127159079/10995116277760000000000000000` |
| Matched | Upper | Floor | `382522191322228678519/2748779069440000000000000000` |
| Matched | Upper | Ceiling | `35825417362754219211347/54975581388800000000000000000` |
| Wide | Lower | Floor | `80878351158765429162857/109951162777600000000000000000` |

The matched intervals are approximately `[4.47979385981e−7, 2.22922586192e−6]` lower and `[1.39160762527e−7, 6.51660545605e−7]` upper. The wide lower floor is approximately `7.35584318670e−7` and replaces the matched lower floor in the sum. “Lower/upper witness” specifies the original payload/sign; “floor/ceiling” specifies the Bellman inequality on that payload.

## Original payout and support

At an exact G2 host outside M3 with at least two core neighbours, take the physical paying edge `(h,r)` and one distinct donor `t`. For `δ=Δ_h^L g(h,r)`, its outward FY contribution is `F(1/D_t)`, with `F(x)=σ(xδ−A_h(x))` and `σ=−1/+1` on lower/upper. The original bound takes `max_k F(ξ_k)` for this donor slot, so the rebate is the nonnegative difference. It is a slot in original **`g`**, not in `f` or any higher family: the separate incoming FY primal on the paying edge is unchanged. The independent exact checks interpolated the primal directly from raw ticks, matched 336 band extrema per sign, and checked that adding the selected slots back to 210 physical row cases per sign stayed under the original `g` cap. Other incidences, where unpaid, consume only their own original slots.

At actual reward time the four support tests are:

| Credit | Necessary support |
|---|---|
| D18 `g` (whole edge) | `max(D_h,D_r)≤18` |
| Stock8 FY | `max(D_h,D_r)>18`, `D_r≤5`, `D_t≤5` |
| Older-paying-neighbour FY (lower only) | `max(D_h,D_r)>18`, `D_r>5`, `D_t≤5` |
| Old-donor FY | `D_h>18`, `D_r≤5`, `D_t>5` |

Since the old-donor payer has degree at most 5, the required high endpoint of its paying edge is the host. These tests are pairwise disjoint on edge/slot ownership even when the same three physical vertices contribute in different directions or at different times. Upper degree/arity caps only further restrict support. Donor degree 5 stays unpaid; donor bins start at 6 and use extrema only over their discrete actual degrees. No transition is binned; all receives, promotions and stopped successors use actual censuses. The complementary result does not establish disjointness with a separate two-vertex exchange rebate; use max rather than sum there absent a split.

## Registration and normalization

Each path is first registered by its last promotion, paired with the unique old edge. For a mother edge the earlier birth uses the actual mother's post-promotion state and the newborn host `J`; for two children it uses the host's earlier promotion and the newborn child. The two oriented roles have separate first-edge stage caps and create one copy of each directed atom. The mother or first child may already be old when the path forms. A donor below degree 6 receives no reward but follows its ordinary continuation, so there is no extra flux at reaching degree 6. Registration and current-degree masks were checked by direct adjacency permutations of increasing labelled tree prefixes, including the 5/6 and 18/19 seams.

The three stopped test potentials have finite intensities and bounded mass per core on these finite cuts. Under the common stationary normalization the intermediate birth outputs cancel their next-stage inputs; the only exogenous term is `μV(J)`. No independent renormalization of the path or edge laws occurs. The high-discount, downward-rounded test gives `LO·V(J)` as a floor; the low-discount, upward-rounded test gives `HI·V(J)` as a ceiling. Each replay component checked as an integer newborn value after multiplying by `S·GRID` and dividing by the respective endpoint (`S=2^40`, `GRID=10^5`). D18 and stock8 receipts are in ticks/core and must each be divided by `S` once. Extended old-payer and old-donor receipts are already normalized. Consequently the exact combined edge-only floors, using the widest verified floor of each disjoint credit, are:

- Lower: `17863212515522806501670771093/10995116277760000000000000000000000` (about `1.62464971395e−6`).
- Upper: `10133550623128433496792800069/5497558138880000000000000000000000` (about `1.84328212038e−6`), with no upper-witness old-payer addition.

These are edge rebate sums, not new exponent endpoints.

## Replay and pins

```sh
taskset -c 5 nice -n 19 timeout 340s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  python3 -Bu round2/checks/quintet-old-donor-sol.py --wide
```

`PASS old-donor original payout, disjointness and normalization`; `wall=288.83 rss=68640 exit=0`. The four matched passes each visited `1621200` rows; the wide pass visited `11819808`. The script prints all exact mother-pays, mother-donates and sibling-order components. It uses the pinned owner solver for the full serial replay and separate code for selected signed reward bands, row payouts, mask checks and exact component addition. The matched prior producer receipt exists only as README decimals, which agree with this run; the wide prior exact fraction is an asserted equality. No giant verifier or H3 replay, no new payload hash scan, and no mutation to owner files.

The reviewed old-donor `README.md` is Git blob `3db9c6344669fbf4c5115243e7bd1d3424dab6bd`, SHA-256 `3c9bc4b6dd59fa59adc7842bf7b3a170253041b1cd63929f2277109a9991859a`; `replay.py` is Git blob `85f13c41926db20a41c9915a924e0414f5e0f4a5`, SHA-256 `93fad1188fd39638a83e7c904b466015e96f2dd2072ee318ef6273c1266a7ce2`. The checker also pins inherited code and the three earlier credit receipts by SHA-256. Both payloads were opened read-only in pre-authenticated mode with file identity checks.
