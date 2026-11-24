# Global contraction of the signed child-occupation kernel

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 14: Claude for transfer/scalar identities; GPT for full contraction and sign transfer; inputs include Claim 12 and the cap-48 continuation certificate. See the [Astra ledger](../../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


For every `s ≥ 565658/10^6`, the upper bound proved here is

```
sup_S || Δ (s−L⁰)^−1 [p_o Σ_Y δ_Y](S) ||TV ≤ 0.8600511260.
```

Here TV is the full sum of absolute masses, with no factor 1/2. This is a bound on the unlimited kernel, uniform over shapes, not on the killed operator. Increasing `s` decreases the upper comparison integral (not necessarily the TV of the signed integral itself). In particular the bounded derivative solution is unique, with inverse norm at most `1/(1−0.8600511260)`. On bounded scores the intrinsic resolvent exists, and `U=V_s[b+C H]` gives the corresponding bounded signed resolvent. The [physical sign transfer](../critical.md) uses this contraction together with the spatial and positive-trigger results; the continuation bounds below do not resolve a numerical trial sign near criticality.

## Coupling and robust rows

Start with a single extra A leaf at a planted root. Keep the alternating parity coupling of common vertices. On a **shared** directed edge, synchronized selections use the minimum reciprocal and the excess rate on that edge belongs to the source copy with smaller degree. Separately, the higher-degree source copy can select its **private** neighbours, each at its own reciprocal. At the root the two residual selections are coupled in one joint tick as described below. Fresh leaves from a synchronized selection are paired; fresh leaves from separate selections are private. The latter includes the two separate residual choices at a root-clock tick. Common neighbours at the root have odd parity, so their reciprocal in A is at least that in B. Only A has private neighbours of the root. Retain `u,v,c,e,l,k` for common LL, common A-leaf/B-core, equal core subtree, conservatively dirty core pair, private A leaves, and private A cores, respectively. Set `B=1+u+v+c+e`, `A=B+l+k`. No core is demoted and a dirty paired child stays dirty. These are counting bounds; “dirty” need not mean that equality cannot subsequently return.

The toll is `c(1/B−1/A)+e(1/A+1/B)+v/B+k/A`, obtained by cancellation on equal core pairs, then the triangle inequality. Root gain from common leaves has rate `u` and adds a paired leaf. Private leaves add private leaves at rate `l`. Each LC sends a paired leaf at `b∈[0,1/2]` and a private leaf at `1−b`; clean CC sends a paired leaf at rate in `[0,1/2]`; dirty CC has paired/private rates `b,a−b` with `0≤b≤a≤1/2`; private cores send private leaves at a rate in `[0,1/2]`. These rates can be history dependent. Each uncertain group is linear in a scalar or simplex of rates, so the linear drift is bounded by the grouped simultaneous corners.

The common rate-one root clock selects each common neighbour or inert slot synchronously at rate `1/A`. Synchronous LL becomes a clean pair of virgin cores; LC becomes a dirty core pair. Synchronized selections on common core targets and null choices do not change the census. A residual pair of selections chooses a private neighbour of A and a common neighbour (including the inert slot) of B with joint rate `1/(AB)` per pair. The changes occur in a **single** successor: private A leaf to private core, B LL to LC, B clean CC to dirty CC. The joint choices give the correct source marginals and preserve parity because the A targets are private and the B common targets are at odd depth. An equal common child cannot first become dirty from non-root evolution: all its clocks are paired while it is equal, and its only outside input is a send by the root. Dirty descendants may evolve differently with no change to the six counts. LL/LC changes outside the root clock are only sends to the root, since the A leaf has no other neighbour. Thus no ignored physical event can increase this census toll through an unrecorded change.

### Terminals

With `j` actual child cores and degree `d`, an unlimited one-side occupation bound is `min(1/s,(j+1/s)/(1+s d))`. For the second function, degree gains can only lower it, so the promotion term bounds its generator above by `(d−1−j)/(d(1+s d))`; its discounted drift exceeds `j/d` by at least `1/(d(1+s d))`. The one-side backward table below the high boundary tightens both terminals using gain endpoints `L,L+j/2`, promotion `L/d`. After a pair exit, the sum of the two tabulated bounds is valid by the conditional jump formula and strong Markov property; the two sides need not be independent. The full toll is also at most `2/s` after integration.

For a second pair terminal write `m=v+e+l+k` and `d=B(0)`. The instantaneous toll is at most `min(2,2m/d)`: use `c(1/B−1/A)≤(l+k)/B` and `v/B+2e/B+k/A` for the rest. Under the coupling `m` gains from private leaf/core and LC/dirty CC sends at total rate at most `m`; clean sends do not increase it. At a synchronized root promotion `m` is unchanged. A residual tick increases it by at most one, at rate `(l+k)/A≤m/d`. Thus `E m_t≤m_0 exp((1+1/d)t)`, and Jensen for the concave truncated linear function bounds the integrated toll by

```
∫₀∞ e^(−s t) min(2,(2m₀/d) exp((1+1/d)t)) dt.
```

The coupon table uses a rational outward Riemann sum with step `1/64`: use the discounted factor at each left endpoint, the growing factor at each right endpoint, and pay `2/s` times the discounted factor where the grid envelope reaches 2. The exponent step is bounded by 12 positive Taylor terms and twice the first omitted term; the discount step uses `exp(−x)≤1/(1+x)`. The coupon for `m₀=1` is decreasing in `d` in the exact integral, so the certified value at degree 48 bounds all larger initial degrees regardless of the leaf count. Inside the finite domain the coupon is used at exits with the actual census, alongside both one-side continuations.

## Production and verification

`pair.cpp` uses integer arithmetic (`10^10` score scale, `__int128` products), two dense composition slices, deterministic indexing and exact integer ceilings. Same-degree moves reduce the lexicographic key `(l,u,v,c,e,k)`; all other moves increase `A`. For a fixed score the row drift maximizes each independent affine rate group. Starting at zero, increase the candidate by the exact ceiling of its maximizing policy quotient until no endpoint product has positive drift. At termination its own affine maximization tests all corners by separation. Optional `all` explicitly multiplies out the corners as a second arithmetic check. Null events cancel on both sides of the row. At `A=48` the next gain uses the analytic/table terminals, rather than killing the contribution. The solved nonnegative supersolution bounds the stopped discounted occupation by conditional drift and pays enough beyond the stop. One may localize in physical event count; the bounded score until exit and discount remove the deterministic-time terminal.

Production at cap 48 visited 22,956,300 rows. The maximum over all initial leaf counts for degrees 2–47 was `0.8600511260` at bare degree 3; the analytic high-degree tail was `0.8484503829` at degree 48. Selected bare-root bounds were `0.8495899211` (2), `0.8152448503` (4), `0.6780168587` (8), `0.5702924239` (16), `0.7959038218` (32) and `0.8446216800` (47). The initial-set maximum takes the smaller of the Bellman and coupon comparisons when both apply.

The separate `Fraction` fixture lists root selections per labelled neighbour and enumerates explicit products of uncertain endpoints, comparing exact integer outputs at cap 8. Both planted mutations were detected: splitting a residual double selection into two transitions, and dropping private-core sends. During development a conservative fixed-only initialization for a compulsory-LC row differed from the independent optimum and was replaced with initialization at zero.

### Newborn sign pilot

`newborn.py` extends the finite signed derivative by zero and bounds its full residual using intrinsic discounted first-exit mass and the contraction norm. At size 8 and `s=.565728`, the killed candidate is in `[.538120420139,.538120420140]`, but the certified radius is at most `18.189526049292`. At size 11 the true signed newborn enclosure is `[-16.283044514395,17.473459651808]`. Neither gives a sign. The global stability problem is resolved; size-cut forcing error, even after paying the exact escape transform, is too large. A useful sign needs a genuinely accurate global forcing/continuation surrogate (or a direct bound on the signed residual as seen by the newborn occupation rather than the global sup norm).

From the repository root, with `PAIR_EXE` naming a writable scratch executable:

```sh
g++ -O3 -std=c++20 round2/influence/zero/pair/pair.cpp -o "$PAIR_EXE"
(ulimit -v 6291456; "$PAIR_EXE" 48) # full production, not a short check
python3 -B round2/influence/zero/pair/fixture.py "$PAIR_EXE" 8
python3 -B round2/influence/zero/pair/fixture.py "$PAIR_EXE" 8 split-root
python3 -B round2/influence/zero/pair/fixture.py "$PAIR_EXE" 8 drop-private-core
python3 -B round2/influence/zero/pair/newborn.py 8
```
