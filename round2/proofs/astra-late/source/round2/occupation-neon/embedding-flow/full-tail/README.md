> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Full uncharged D18 generator: effective truncation at the physical law

For `W=Φ_low+t Φ_hi`, the **entire** physical score `(L−s)W` is absolutely integrable in uniform-core units and has an effective finite truncation. At `t=1`, `|s|≤.565775`, dropping omitted rows costs less than `10⁻⁸` with `K=87,J=2^104`. Retaining the limiting weight at degrees above `J` reduces this to `K=87,J=2^43` for the same tolerance, with all arity overflow still discarded. If only the *stationary expectation* is needed, exact integrated balance gives error zero at `s=µ`, and throughout the given discount bracket the simpler potential-only limit approximation reaches `10⁻⁸` with `K=49,J=2^27`. The distinction between a termwise generator approximation and a stationary identity is essential.

## 1. Physical expansion and ownership

Work at the critical stationary probability `π` on nonfounder oriented cores. Each `v` has census `(D,k,h=1)`, leaf stock `L_v=D−k−1`, one mother, and `k` core children. There is no recharge and no extrapolated degree in a physical rate. For each occurrence of a raw D18 family `F`, write `w_F(x)` for the normalized paired potential weight. Use each physical occurrence once; labelled-role LP columns are divided by the symmetry factor (`2` for `hC,Yt`, `6` for `S3`, otherwise `1`). The pre-jump ledger is as follows:

* Discount on each extant occurrence: `−s w_F`.
* Gains at each tagged vertex: `Δ_D w_F` times its leaf rate `L_v`, plus its incoming rates `1/D_u` from each core neighbour. The latter are either internal directed edges or one actual appended exterior donor. A donor arrow increments the **receiver's** degree only. The leaf rate counts `L_v` actual leaf donors of degree one without making a new core column.
* Promotion at each tagged vertex: `b_v Δ_k w_F`, where `b_v=L_v/D_v` and `Δ_k` is at fixed degree (`L_v` decreases by one). This aggregates the actual degree-one recipient edges; no exterior core is silently used as a leaf. There are no new old-to-old edges. Newly promoted cores have exactly one old core neighbour, their promoter, so the new role in a connected motif of size `r≥2` must be a bottom tip. Delete that role, keep the connected `r−1` old roles, and charge `b_v w_F(post,virgin)` to the actual promoter `v` in the old list. The promoter is evaluated **after** promotion. The new singleton `f` is explicitly `b_v f(1,0,1)` (linearly paired), at every promoter; no predecessor motif factor applies.

There are at most `r` leaf-gain terms, `2(r−1)` directed internal gains, `r+1` exterior edge types (append a child at each role, or a mother only above the apex), `r` existing promotions and `r−1` newborn tip types for a connected `r`-role pattern. Using all these counts before dividing by symmetry is safe. When a term has several possible old-role hosts or equal siblings, the rate belongs to the specified actual incidence. No extra product of parent and donor marginals is used.

Old, post-gain and post-promotion support tests are evaluated separately; a successor can have zero weight. A degree or arity cap in the following results is a test on **pre-event** roles, not permission to zero successors. A gain taking `J` to `J+1` and a promotion taking `K` to `K+1` are still computed at their actual successor inputs. Raw motif weights and the `LF,KF` convention come from the pinned D18 reader (`round2/attacks/neon1/audit.py`, historical reference); `f2/f3`, edge center tables and routing transfers are not terms of the uncharged generator.

## 2. Pointwise term constants

Let `M_F` be the maximum absolute entry in either original table for the family, divided by `10^11`, times `1+|t|`. Read-only chunk reductions find the eight maxima in [`check.py`](check.py), with both inputs using `D0=18`, `LF=14000`, `KF=180`. For each embedded role, padded by zero on disallowed controls, the tail is a cubic in `18/D`. An adjacent control difference is bounded by `2M_F`; thus with `C_F=108 M_F`,

```
|w(D)−w(∞)| ≤ C_F/D,
|w(D+1)−w(D)| ≤ C_F/[D(D+1)]   (D≥18).
```

An inner receiver has `L≤16` (the bound `17` suffices), so `|L Δ_D w|≤34 M_F` also at `17→18`. Above the seam it is bounded by `C_F/(D+1)≤34 M_F`. The cutoff domains in the pinned `P4,Yt,Yb,S3` interior slots all end below 18; `hV,hC` can have a cubic at their capped-arity center, with the all-infinite center control zero. Thus the inequalities cover fixed, vanishing and overflow controls. For `f`, the score is a separate lookup capped at `LF` and constant in `L` for `k>KF`; its leaf-gain norm is the **finite** maximum of `L|f(L+1,k,h)−f(L,k,h)|`, scanned at `0≤L<LF, 0≤k≤KF`. It is at most `1.46044806770(1+|t|)`. In particular when `D>J≥LF+KF+2=14182`, `Δ_D f=0` and `f` and its promotion successor are already in the leaf plateau (or arity overflow). A promotion still can change the row, including the `KF→KF+1` seam. All core reciprocals are at most `1/2`; existing promotions have norm at most `2M_F` and new ones at most `M_F`. The `f` vertical gain is treated as a two-role pattern for **both** directions.

## 3. Re-rooted tail bounds

Each summand type uses a connected pre-event pattern of `p≤5` cores. If all arities are at most `K`, re-root at a labelled role with `D>J`; upward steps are forced and each downward step costs at most `K` choices. The expected number of omitted incidences is bounded by `p(K+1)^(p−1) π(D>J)` per type. If any arity exceeds `K`, assign the whole incidence to one of its **maximal-arity** vertices, with a deterministic tie rule. The expected count of such incidences is at most `p π[(k+1)^(p−1); k>K]`. Since `p≤r+1` for an original `r`-body family, the weight on this vertex is at most `(k+1)^r`. Do not apply a root-aritiy tail to an apex that is not the offending vertex. These are physical incidence reversals: expand along true parents and actual children on the oriented core tree. The physical law satisfies them by mass transport, including when one donor is above the apex.

By the [physical stable-age law](../../../../../../astra/source/round2/degree-law.md#birth-times), at common exponential age a newborn starts at degree two and arity zero. The gain process is dominated by the pure birth chain of rate `d−1/2` at `d`, and the arity by its own unit-rate Poisson clock. The mechanisms can be coupled by the disjoint incoming-gain and outgoing-ring proposal channels, with the same elapsed age, as in [pressure-memory](../../pressure-memory/README.md#an-effective-joint-degree-and-child-count-tail); these are not independent physical censuses. At rate `s0=.565673≤µ`, `N(Exp(s0))` is geometric with tail `(1+s0)^(-K−1)` and, conditional on `N>K`, is distributed as `K+1+N'` for a fresh geometric `N'`. Its factorial moments are `m!/s0^m`. This bounds the full arity polynomial loss by finite rational arithmetic.

For degree use the hitting transform `P_J=∏_(d=2)^J (d−1/2)/(d−1/2+s0)`. There is no product loop up to the huge cutoff: `log(1+x)≥x/(1+x)` and `d−1/2+s0≤d+1` imply `P_J≤(3/(J+2))^s0`, by an integral bound on the harmonic sum. For `J=2^j`, use the exact, looser dyadic enclosure `2^(-floor(9(j−2)/16))`, since `s0>9/16`. Multiply each incidence bound by the absolute norm from §2; the arity and degree events are disjoint in this accounting. This controls the **sum of the absolute omitted terms**, not only their signed sum. For `η=π/µ` multiply the resulting per-core error by at most `1/s0`.

### Keeping the leading high-degree score

On `k≤K`, replace the census of any pre-event `D>J` by a formal infinity of the **same arity** only in the approximant. Keep topology and exact low censuses. Use the tail limit of each motif weight, zero reciprocal for a high core, and `b=1` for a high promoter. For the old gain at a high receiver use zero, including for the leaf term; this is the limit of `L Δ_D w`, not the formal product `∞·0`. A low receiver's physical gain successor is still evaluated at its true degree, even when that successor crosses the cap. `f` retains its plateau indexed by its actual arity. This produces a finite degree alphabet with an overflow symbol. For a high tagged weight the error is bounded by `C_F/J` per changed argument (telescope the tensor); a gain/promotion difference pays for both evaluations. For a high **donor** the reciprocal error is at most `1/J`; a high promoter has `1-b≤(K+1)/J`. The leaf gain error for a high receiver is `C_F/J`; otherwise `L|Δ_D w−Δ_D w_∞|≤34(p−1) C_F/J` when its low degree is inner. For a low **tail** receiver, differentiate in its own tail coordinate: the change in any of its four controls from the other high tags is at most `(p−1)C_F/J`. Cubic differentiation then gives `L|Δ_D w−Δ_D w_∞|≤108(p−1)C_F/[J(D+1)]`, within the same bound. All the remaining constants are listed as `frozen error*J` in [`terms`](check.py), including outside-donor patterns and newborn `f`. The same re-rooted degree probability now multiplies an extra `1/J`.

The overflow symbol needs no new infinite hierarchy row. For fixed `K,J`, sum over the relevant high-role censuses using projection and tightness, or approximate each overflow cylinder by finite sub-slices up to a further `J'`: its bounded omitted mass is controlled by the same re-rooting bounds at `J'`. If an implementation only has the fixed-column identification constraints, the rigorous evaluation prescription is a *fixed finite* approximation with the physical remainder attached. To let the objective support grow along a cofinal LP sequence, impose the same nonnegative arity/degree slice budgets on the corresponding labelled term counts (including order five for exterior donors), not just on the total mass missing at the final cap. They are necessary physical constraints by §3; they give uniform errors for feasible relaxations by comparison to any fixed intermediate cap. Fixed-column convergence alone does not control objective values computed with no such cuts. Neither procedure requires a cutoff on number of ordinary leaves in a realized tree.

### Matched subtraction using the generator balances

An alternate sharp finite objective uses formal subtraction first, without an overflow symbol. For each free D18 slot subtract its limiting weight on the region `D≥18`; more precisely use the product of commuting operators `I−1_{D_i≥18} E_i^∞` on each family weight, including all combinations when several tags are large. For `f`, subtract its leaf plateau **globally** as a function of arity. This leaves exactly zero `f` at any high degree and zero anchored `f` gain/promotion there; the newborn anchored `f` term is integrated explicitly as `µ f_anchored(virgin)` (retain this constant outside the cap). The subtraction seam at `17→18` belongs to the derivative of the anchored potential. The anchored table norm is at most `2^m M_F` for `m` free slots; anchored high-slot weights and their successors are `O(1/D)` and high-slot gains `O(D⁻²)`. If a high role is only an exterior donor, use `1/D` on that term. These yield the `anchor` line with the same `P_J/J` loss. Promotion cannot reintroduce a high capped-interior role. The anchored higher family's newborn still uses its matched post-promotion weight.

Let `A` be this anchored potential and `B=W−A` the removed potential. This is a **matched replacement**, not subtraction of an unpaired running reward: restore `(µ−s)π B` to the truncated anchored generator. Justification is the integrated identity below (also valid for `A,B` by the same moment bounds). At `s=µ` the restoration is exactly zero. At other trials one must include it; the anchor rows in the table budget the residual generator only. The stationary-limit line shows a direct alternative at arbitrary trial in the bracket.

## 4. Stationary integration and the completed hierarchy

For an integrable labelled occurrence weight, summed over its copies, the physical count at time `u` per starting cell has expected value `e^(µu)` times its stationary per-core mean. Equivalently this is the eigenmeasure relation for `Q_u`, applied to nonnegative and negative parts of the occurrence statistic. When writing the jump increment as the pre-event ledger of §1, some births produce motifs whose apex is an older core, and some incoming directed gain clocks belong to other cells. Transport the event payment to its prescribed apex; changing this owner leaves the stationary sum invariant by oriented mass transport. Over a bounded time interval, the integral of the absolute jump decompositions is finite by the same stationary majorants from §§2–3, multiplied by `e^(µu)`. One may first enumerate all finite-census bounded occurrence cylinders, use compensation and add their event ledger, then pass with the summable absolute tails. This gives the full balance

```
π L W = µ π W,                 π (L−s) W = (µ−s) π W.
```

For a direct finite-row justification, cap **the potential's pre/post census tests** at fixed degree and arity and sum its finite labelled balances from the completed hierarchy. When removing these caps, the additional degree seam is at `D=J`, proportional in magnitude to its gain intensity. For each affected labelled role re-root at that receiver with `D=J`: the multiplicity with other arities bounded by `K` is at most `r(K+1)^(r−1)`. Its integrated gain intensity is `µ π(D>J)` by the [physical degree crossing balance](../../../../../../astra/source/round2/degree-law.md#point-probabilities), so the seam vanishes as `J→∞` at fixed `K`, even though `πD=∞` need not be ruled out. Arity seams have bounded promotion rate and the same arity moment bound (use `k≥K` instead of `k>K`). Newborn and outside donor contributions pass by absolute convergence. Then remove the arity cap. This also explains why pointwise weight approximation alone cannot justify passing a hard cutoff through the generator.

The stationary score at `s=µ` is exactly zero; this does not say the pointwise score or its absolute discarded tail is zero. Away from µ, integrate the full gain and birth ledger by the identity and approximate only `πW`, scaled by `|µ−s|≤.000102` over the original trial bracket. Degree and arity caps use the discount-only line with this small gap in place of `|s|`. Both forms apply to the physical signed pair for any rational `t` simply by changing `1+|t|`; cancellation within blobs is not needed for the stated guarantees.

## Results and replay

The grid chooses the least `K` meeting half the budget in the arity term, then the least dyadic degree cap meeting the whole budget at that `K`. These are effective caps, not claims of minimal physical cutoffs or of a practical order-five LP at these sizes. Values below are rounded upwards.

| Scheme | At `10⁻⁶`: `(K, log₂ J)`, error | At `10⁻⁸`: `(K, log₂ J)`, error |
|---|---|---|
| Full generator, omit | `(76,91)`, `7.49748e−7` | `(87,104)`, `9.17706e−9` |
| Full generator, retain limiting rows | `(76,38)`, `9.44524e−7` | `(87,43)`, `8.36419e−9` |
| Matched anchored generator (restore correction exactly) | `(81,40)`, `9.28274e−7` | `(92,45)`, `7.80221e−9` |
| Integrated score, omit potential tail | `(37,45)`, `7.87112e−7` | `(49,59)`, `7.33218e−9` |
| Integrated score, retain limiting potential | `(37,22)`, `7.17238e−7` | `(49,27)`, `6.72092e−9` |

From the repository root:

```sh
taskset -c 8 nice -n 18 python3 -Bu round2/occupation-neon/embedding-flow/full-tail/check.py
taskset -c 8 nice -n 18 env OPENBLAS_NUM_THREADS=1 round2/.work/venv/bin/python -Bu \
  round2/occupation-neon/embedding-flow/full-tail/check.py --scan b3f18_lo.bin b3f18_hi.bin
taskset -c 8 nice -n 18 env OPENBLAS_NUM_THREADS=1 round2/.work/venv/bin/python -Bu \
  round2/occupation-neon/embedding-flow/full-tail/check.py --sparse b3f18_lo.bin b3f18_hi.bin
```

The ordinary budget replay uses no payload. The optional single read-only scan obtains potential maxima and the finite `f` leaf-gain norm; it retains no payload or extracted tables. Blob authentication is inherited from the original D18 verification, not rerun here. The exact finite fixture compares arity moments in two bases, checks cubic bounds on both ordinary and extreme degrees, and checks complete jump differences against old-role and newborn ownership on physical directed arrows. The sparse replay passed 28 whole-potential gain/promotion comparisons, with a nonzero newborn term in every one of the eight production families. No bracket endpoint is recomputed.
