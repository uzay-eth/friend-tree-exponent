> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# High-hub occupation from the physical degree entrance

**Result.** The physical high-hub charge law gives an exact discounted difference from the marked Poisson/inert-arm forecast after first hitting `M≥32`. A sharper Bellman solve on apex arity `k≤12` and first-arm size `≤7` gives two-sided stationary prices including **both old and new arms**. At `M=32` the total is **`[9.2621318e−7, 2.99343152e−6]` / `[3.0238689e−7, 1.79834854e−6]` per core** for the lower/upper witnesses. The `M=101` mask gives **`[5.2327439e−7, 1.51929616e−6]` / `[1.6183503e−7, 7.8516734e−7]`** and clears the known finite high-apex overlap bands. The uncapped complement is still open; no new exponent endpoint is composed here.

## 1. The entrance law and the calendar-time coupling

Start from the actual physical birth-Palm law of a nonfounder, with its living mother. Set `τ=inf{t:D(t)=M}`. If it is finite, condition on the complete history to `τ`: this includes its age, mother, all old children, and their actual geometries and ages. Subsequently write `D,k,c=k+1` for hub degree, child count and core degree, and `x` for the full planted fringe of a first child. The rate of a new child is `b=1−c/D`. The hub gains at rate `D−c+ρ`, where `ρ` is the actual core receive sum. Every old child receives a hub send at rate `1/D`; the child's opposite arrow changes only the hub. Internal branch rates use true degrees.

One exact construction splits the hub clock into a unit promotion-proposal clock thinned at `b` and separate core-send channels thinned at `1/D` (one predictable channel per core port). These are equivalent in law to the original outgoing rates; the proposal clock counts even rejected points. At each promotion proposal attach an independent inert virgin branch, started at the same calendar time and independent of the proposal process **before** conditioning on acceptances. Use its drivers in the actual branch on acceptance up to the first hub send. Old arms start from their entrance states, not virgin states; a copy for any old arm may reuse that actual geometry and omit only future hub input. No conditioning on eventual hub status is used.

For a fixed post-entrance window of length `t`, put `N_t=Pois(t)` for the proposal count and let `c_0` be the entrance core degree. Pathwise `c_t≤c_0+N_t`. Expected rejected points are at most `(c_0 t+t²/2)/M`; the expected total of sends into tracked children is at most the same (using a parent slot in the bound). Thus the young-arrival process has the independent marked Poisson local limit on fixed windows if `c_0/M→0`. At calendar time `t`, a proposal at `u` carries the physical age `t−u`. A fixed-window approximation by itself does not justify an all-future weighted calculation.

For the hub clock at *any* elapsed time, the usual [tagged-Yule construction](../../../../../astra/source/round2/fringe/hub-spectrum.md#2-a-vertexs-limiting-mass-as-seen-at-time-t) still sandwiches physical degree by the full Yule population from `M` above, and the Yule growth of `M−c_0` initial leaf tags minus a rate-one immigration Yule loss below. Start the loss families at proposal points, including unused ones; independence of the upper/lower populations and the physical marks is not needed or claimed. For the all-future error we use a sharper inverse-degree identity from the same exact Yule gain and proposal channels. Put `F_q(d)=(d−1)^−q`, `γ_q=M[1−(1−1/M)^q]`, and

```
A[q,m] = Eτ ∫₀∞ e^(−s t) (N_t)_m F_q(D_t) dt.
(s+γ_q) A[q,m] ≤ 1{m=0}F_q(M) + m A[q,m−1]
                  + γ_q[(c_0+m) A[q+1,m] + A[q+1,m+1]].
```

Indeed `d[1−(1−1/d)^q]≥γ_q`; even conditional on the whole past, physical gains are at least `D−c`, and `c/D·F_q≤(c_0+N)F_(q+1)`. `N` and `D` do not jump together. The recurrence starts with `A[q,m]≤F_q(M)m!/s^(m+1)` at large `q`, then checks four inverse orders; one can always intersect with this elementary bound. This controls all years after entrance, including the histories with unusually little degree growth. In particular it is unnecessary to delete the gamma lower-tail failures or apply a deterministic exponential replacement to those paths.

## 2. Weight the errors by changed D18 reward

The mask in the pilot keeps an occurrence when the **full first-child arm** has at most seven ordinary vertices, counting the child but omitting the exterior hub slot. Each `hV`, `P4`, or `Yt` has one such first arm; the two fork tips are counted as an unordered pair. The arm expands monotonically, so a size exit is permanent for that arm. Grade the apex as the actual nonfounder with `D≥M`, and let `R_f(D,k,x,u)` be the sum of actual signed D18 charges from these embeddings in `x`. Here `u` is the real reciprocal sum of the mother and the `k−1` other arms. Every other donor is read in the full arm. In each row the coefficient of the mother and of each other arm is the same and `0≤u≤k/2`. The original payloads are used as authenticated earlier; the replay loads only their mapped sparse rows. It does not fit a surrogate drift.

At reciprocal degree zero, `R_f^∞(x)` is independent of hub arity and exterior receive. Resolve the exact inert branch generator `L⁰` on the arm-size domain (zero at exit): `U_s=(s−L⁰)^−1 R_f^∞`. The zero-reciprocal forecast from the entrance is

```
I_s = Σ_{old actual arms x} U_s(x) + U_s(J)/s.
```

The second term is the Campbell integral of future unit-rate inert immigrants. For a finite future horizon the same construction uses the time-dependent forecast `U_(s,T−t)` and replaces `U_s(J)/s` by `∫₀^T e^(−s t)U_(s,T−t)(J)dt`. The whole-future error instead has this **signed identity** (one per family):

```
actual discounted occupation − I_s
  = Eτ ∫ e^(−s t) [ Σ_actual arms x ( R_f(D,k,x,u)−R_f^∞(x)
                                      + (U_s(x⁺)−U_s(x))/D )
                       − c U_s(J)/D ] dt .
```

Apply discounted Dynkin to the sum of `U_s` over actual arms: intrinsic evolution pays the inert resolvent, sends change only the addressed arm, and actual births add `U_s(J)` at rate `b`. The boundary expectation vanishes by bounded arm values and the linear proposal mean. This weights each hub send by the **exact change in future motif occupation**, not by a worst-case probability of any discrepancy on the hub. Rejected birth proposals cost only `U_s(J)`. The identity does not assume that arm state and hub growth are independent. For the stationary version multiply the entrance identity by the exact `e^(−μτ)1{τ<∞}` and average over the actual Palm law. Equivalently, use `π(q)=E_birth∫ μe^(−μt)q(X_t)dt`; there is no resetting of old-arm ages and no substitution of permanent mark weights for actual maternal history.

### Exact reciprocal and overflow checks

For apex `p=1/D`, original Bernstein weights are cubic in `18p`; arity controls saturate at 18. After multiplying by `(1+p)^3`, the row is a polynomial of degree at most seven: hub leaf gain has a naïve `1/p`, but the gain increment has a double zero at `p=0`. Hub promotion uses `1−(k+1)p`. Other variable-degree terms are products of a cubic with a rate at most affine in `p`; transfers onto the apex have the same cubic or gain increment form. At this arm cap all non-apex motif degree labels, including their one-step intrinsic successors relevant to the row, are inner. The donor slope has a double zero. After arity 18 the base numerator is therefore exactly affine in *actual* `k` (through the gain leaf count and promotion rate), with the same donor slope. The code reconstructs by rational interpolation, verifies separate degrees including `(101,81)` and uses Bernstein extrema on `[0,1/M]` after dividing the difference from the constant term by `p`. Thus it bounds the signed score correction by `(a + a_over (k−18)_+)/D` and `(b + b_over (k−18)_+)/D` per arm, retaining the true `k` and the sign of the send difference.

For the stationary error with no entry sample, dominate physical `D` by a Yule chain starting at 2 with gain `d−1/2`. It uses independent proposals from the rate-one hub promotion clock; the physical variables need not be independent. At the **same** `A~Exp(μ)`, `k(A)≤N(A)` and `D(A)≤Y(A)`. For `X_j=E[(N(A))_j 1{Y(A)≥M}]`, use the product of Laplace hitting factors at `LO` and tilted hitting-time exponentials plus the residual `Exp(LO)` (the common-age formula). For example, the overflow factor is bounded by `min(X₁+X₂,(X₃+2X₂)/17)` since `k(k−18)_+≤N(N)_2/17`. Do not use `π(D≥M)π(k^j)` or condition the inert marks on the final degree. Between the two discounts, differences `U_s(x⁺)−U_s(x)` are bounded using opposite endpoints for the two nonnegative terms. This yields uniform stationary errors, but at `M=32` their component widths can exceed `10^-3`. Their scale is `O(M^(−1−LO) poly(log M))` for fixed arm cap, versus a degree-tail mass scale of `M^(−μ)`. It is a convergent limit comparison, not a useful standalone remainder budget at 32.

## 3. Sharper stationary price for young arms

For a finite-threshold numerical price, keep only arms **born at hub degree ≥M** and evaluate while true hub arity `k≤17`. The whole start-at-entrance term for old arms is outside this mask. At the actual birth, evaluate a continuation at the hub's *post-promotion* arity. Reward is bracketed by exact statewise reciprocal envelopes at each `k`, before integrating the inert trajectories. For the lower solve stop on the first independent rate-`1/M` send proposal to the tracked child; for the upper solve retain both exterior send extremes (`0` and `1/M`) in the actual gain direction. Both solves retain the real arity changes via a fractional-linear Bellman test on `b∈[1−(k+1)/M,1]`. Use the **same** `b` in the denominator and continuation and stop at the true `k=18` exit. An unobserved hub gain can change these predictable rates but cannot move the compressed coordinates `(k,x)`, so the generator inequalities permit adaptive rates. The guards retain the common calendar age and are uniform after conditioning at the actual birth. All directions are exact with outward grid rounding of `10^-23`. Use the matching discount endpoint for each sign of the enclosure.

Now use birth compensation on the whole host life: from census `(2,0)` retain physical promotions at exactly `(D−k−1)/D` and bracket gain by `[D−k−1,D−k−1+(k+1)/2]` up to hitting `M`. Above that entrance, for this auxiliary birth price the only retained state change is promotion at rate between `1−(k+1)/M` and 1. Close the upper/lower Bellman inequalities with both rate endpoints. The same rate is used in each promotion reward, holding denominator and arity transition. This retains the waiting time to the threshold and subsequent birth, and applies to arbitrary live mother and sibling environments. No independence of the birth age and successful local motifs is used.

All entries below are outward seven-arm prices per core. Family order is `hV / P4 / Yt`; the main cells show the `P4` or `Yt` scales for which the method is useful. The corresponding `hV` widths show why a full D18 remainder cannot yet be booked.

| M | Witness | `hV` interval | `P4` interval | `Yt` interval |
|---:|---|---|---|---|
| 32 | lower | `[5.798e−11, 4.85805905e−6]` | `[3.7057498e−7, 1.12637983e−6]` | `[4.99953e−9, 3.796948e−8]` |
| 32 | upper | `[5.797e−11, 6.64167419e−6]` | `[4.683017e−8, 2.3514097e−7]` | `[6.713497e−8, 2.0218867e−7]` |
| 64 | lower | `[4.117e−11, 7.8597654e−7]` | `[2.8777381e−7, 6.9743289e−7]` | `[3.76583e−9, 1.273279e−8]` |
| 64 | upper | `[4.119e−11, 1.01479026e−6]` | `[3.611968e−8, 1.0270154e−7]` | `[5.008367e−8, 1.2006645e−7]` |
| 128 | lower | `[2.854e−11, 1.3444258e−7]` | `[2.0709824e−7, 4.5767377e−7]` | `[2.67805e−9, 6.52578e−9]` |
| 128 | upper | `[2.855e−11, 1.6267080e−7]` | `[2.594340e−8, 5.993138e−8]` | `[3.542935e−8, 7.776329e−8]` |

For the tighter mask also require post-birth hub arity ≥5 (so it remains ≥5). At `M=32` the lower floors sum to `1.4737502e−7` and `4.475683e−8`. This excludes complete apex crowns with at most four cores, but can still overlap previously priced high-apex supports. No credit is added to an existing ledger. Covers that retain the actual parent/sibling reciprocal sum on the dominant positive `hV` rows, and that cap only the tagged path rather than the whole arm, are plausible follow-ups. The present envelope does not settle the unbounded-arm complement or unpriced support overlap.

## Replay

From the Astra checkout, run the two witnesses separately (the paired invocation can hit the wall limit under load):

```sh
R=round2/occupation-neon/high-hub-limit
P=./round2/.venv/bin/python
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 3 nice -n 15 /usr/bin/time -p timeout 105s "$P" -Bu "$R/check.py" b3f18_lo.bin --cap 7 --degrees 32 64 128
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 15 /usr/bin/time -p timeout 105s "$P" -Bu "$R/check.py" b3f18_hi.bin --cap 7 --degrees 32 64 128
```

Both completed (67.93 s / 51.17 s wall). Each visits 84 arm states, of which 78 have at least one structural motif stencil, with extra exact polynomial checks at the finite and overflow arities. No full scan or new payload hash was run. The comparison and pilot have not had a separate proof review.

## 4. Sharper two-sided bulk price, including pre-entrance arms

[`bulk.py`](bulk.py) retains `(D,k,x)` up to a degree seam `T=101`; it pays the full supported `hV/P4/Yt` charge for `D≥M`, `k≤12`, and first-arm size `≤7`. Every actual arm promotion registers its virgin arm for the full subsequent continuation, **even when the hub has yet to hit `M`**. At a later residence, the same bounded arm may contain a mix of older slots and virgin tips. No motif-role ages are reset, and no all-nonvirgin test is imposed. The old/young split refers to the **first arm's** birth relative to the first degree hit. It includes mixed old roles inside that arm, not motifs whose tagged roles span two separate apex arms.

For `D<T` the row uses the exact reciprocal polynomial `a+b u`, with `u` the single physical mother-and-other-arm reciprocal sum. The hub gain rate is `L+1/D_x+u`. In the arm Bellman test use the **same** `u` in that gain, its successor, and the reward; the fractional-linear extrema are at `0,k/2`. Promotion changes `k` but not `D`, sends to the arm occur at `1/D` in the actual gain direction, and intrinsic arm rates are exact. For `D≥T`, `band_price` supplies the two-sided, infinite-time Bellman continuation: retain true arity transitions, interval-bound the reciprocal row at that arity, use the rate-`1/T` killing guard for the lower bound and the physical send-rate interval for the upper bound. At the hub-only boundary (pricing future births), each promotion adds both the newborn arm's value and the future-hub value. A crossing into `D≥T` carries these continuation values once. Arms born below `M` can earn past either seam until the same physical exits. Monotone size, degree and arity give a backward order in the finite part. All quotient corners and the rounding grid are rational; positive directions are converted to stationary units using the matching discount/normalization endpoints. The total and the two split rows are separately optimized; summing the component bounds is optional and generally weaker.

| `M`, witness | All first-arm births, per core | Pre-entrance births | Post-entrance births |
|---|---:|---:|---:|
| 32, lower | `[9.2621318e−7, 2.99343152e−6]` | `[4.5273610e−7, 1.67501739e−6]` | `[4.6055908e−7, 1.37032970e−6]` |
| 32, upper | `[3.0238689e−7, 1.79834854e−6]` | `[1.6097848e−7, 1.10593516e−6]` | `[1.3711452e−7, 7.2804439e−7]` |
| 101, lower | `[5.2327439e−7, 1.51929616e−6]` | `[2.9425014e−7, 8.0474448e−7]` | `[2.2532689e−7, 7.2604336e−7]` |
| 101, upper | `[1.6183503e−7, 7.8516734e−7]` | `[9.468162e−8, 3.9050589e−7]` | `[6.595770e−8, 4.0100402e−7]` |

At `M=101`, `P4` supplies `5.1602008e−7` of the lower-witness floor; `Yt` supplies `9.315526e−8` of the upper-witness floor. The upper `hV` ceilings still dominate much of the width. The first-arm size cap is a **reward** mask with permanent exit, not an upper bound on the number of vertices an apex can have. These upper bounds do not enclose the omitted arm-size or arity complement.

**Ledger.** The `M=32` total intersects [old-flow](../renewal/old-flow/README.md), the [65–100 old-hub band](../fringe/old-hub/README.md), and the published [live `hV` donor atoms](../pressure-memory/hv-lifetime/old/EAD947B.md) through degree 48; use max on overlapping objectives. `M=101` misses all three bands, the low-apex [P4 root and complementary installments](../pressure-memory/PATH.md) (including the [wider low-apex replacement](../pressure-memory/path-fast/README.md)), and the original small complete-fringe reserve since the apex's planted fringe alone has size at least `D`. The corrected `f/g`, `hC`, `Yb` and `S3` payments are different charge families. This is a support comparison to those identified releases; a newly composed endpoint would still need a full check against the exact active ledger and its normalization.

The separate [`physical_check.py`](physical_check.py) reconstructs rewards via `Local.charge` on full physical trees, with the actual mother and siblings present. It projects complete-tree arrows onto `(D,k,x)`, including aggregate exit hazard, and checks equal-sibling multiplicity, varied donor sums, and arities on both sides of the tail seam. Replayed on both pins: **1,092 exact host matches per witness**, with nonzero `hV/P4/Yt` occurrence counts `1316/854/238` each; 3.14 s wall. This is a separate evaluator check, not an external proof review.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time -p timeout 60s "$P" -Bu "$R/physical_check.py" b3f18_lo.bin b3f18_hi.bin
# Run each BLOB separately, using CPU 3 for lower and 8 for upper:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c "$CPU" nice -n 15 /usr/bin/time -p timeout 150s "$P" -Bu "$R/bulk.py" "$BLOB"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c "$CPU" nice -n 15 /usr/bin/time -p timeout 150s "$P" -Bu "$R/bulk.py" "$BLOB" --degree 32 --tail 101
```

Bulk wall times were 77.46/53.77 s at `M=101` and 105.58/55.10 s at `M=32` (lower/upper); all four exited zero.
