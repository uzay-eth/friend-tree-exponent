> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../../README.md#claim-25) for exact scope and citations.

# Gold1 FY remainder: both-high and high host arity

All values are normalized per core. Exact-G2, non-M3 hosts and the paying-edge exclusion `max(D_host,D_payer)>18` are enforced. Both-high means the payer and distinct donor have degrees at least six. This is disjoint from the three previously priced degree quadrants. Host-arity-only ceilings include previously priced quadrants, so also bound their unpriced subsets.

| Scope | Lower-witness ceiling | Upper-witness ceiling |
|---|---:|---:|
| All FY pairs at host arity 11–26 | `3.380335e−7` | `3.596892e−7` |
| All FY pairs at host arity 13–26 | `3.112114e−8` | `3.323304e−8` |
| Both-high at host arity 9–26, sum of arity 9–12 and 13–26 ceilings | `8.038226e−7` | `3.186364e−7` |
| Both-high at host arity 2–12, invariant cut 7 | `3.356814e−6` | `3.967930e−6` |

Thus high arity alone (11+) cannot supply `1e−6` on either witness. Both-high is not ruled out overall. At host arities 2–8, invariant cut 8 gives full-neighbor-degree floors **`4.67716169e−9` lower** and **`9.89331855e−8` upper**, with matching ceilings `2.745651e−6` and `3.610876e−6`. These do not justify a `1e−6` credit. The cut-11 attempt did not complete both signs; the completed cut-8 pair is the reported floor. Gold2 work superseded further widening.

## Method

`spectra.py` bounds the FY gap pointwise by the maximal signed-knot score minus the minimum of the donor score over the allowed reciprocal interval. This difference is convex in the paying slope, so its maximum occurs at the extremes of the slope interval. Exact inner, arity-overflow, tail, and infinity G2 corners enclose all physical slopes; the paying-edge gate restricts corner degrees when the host is at most 18. Multiplying the envelope by `c(c−1)` bounds every ordered pair at that host. A singleton occupation supersolution accounts for the cost of reaching physical arity: these are stationary ceilings, not pointwise ceilings on whole-tree averages. This reads only the relevant `g` table differences, not the production verifier's higher rows.

`occupation.py` keeps exact degree and stock before each neighbor crosses the invariant cut, then keeps that physical neighbor in an absorbing reciprocal band `[0,1/cut]` for all later ages. Internal transitions inside a band cancel. The physical host is exact through degree 179, the entire Gold1 G2 host degree support. Pair rewards use extrema over the physical donor reciprocal band and payer slope interval; for the slope minimum include breakpoints of the exact signed-knot upper envelope. For a floor and a host at most 18, an absorbing payer straddles the required paying-edge gate, so its reward is zero there. The ceiling instead bounds only gate-eligible slopes and can harmlessly pay that bound to the whole absorbing band. A mother can produce a registered edge from within the band; the lower value uses zero for that additional source and the upper uses the promotion cap one. No age reset or product closure is used. The exact adverse finite-state inequalities are checked at every rounded value, using the same three-stage physical birth registration as the pinned FY proof.

## Replay

From the project root with the authenticated Gold1 private copies in their original location:

```sh
R=round2/occupation-quintet/edges/fy/remainder
taskset -c 5 nice -n 19 python3 -Bu "$R/spectra.py"
taskset -c 5 nice -n 19 python3 -Bu "$R/occupation.py" --cut 8 --stock 8
taskset -c 5 nice -n 19 python3 -Bu "$R/occupation.py" --cut 8 --stock 8 --ceiling
taskset -c 5 nice -n 19 python3 -Bu "$R/occupation.py" --cut 7 --stock 12 --ceiling
```

The [spectra](spectra.log), [matched floor](invariant8.log), [matched ceiling](ceiling8.log), and [arity-12 ceiling](ceiling7.log) finished with `PASS`. No pinned source or receipt changed.
