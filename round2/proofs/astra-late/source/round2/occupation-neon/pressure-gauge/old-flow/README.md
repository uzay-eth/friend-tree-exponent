> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Terminal-degree gauge on older charged rows

**The terminal antiderivative cancels the selected receive coefficient exactly. It gives zero improvement to either optimized integrated bound with the same flow constraints.** On ranked `P4`/`Yt` faces the full weight also cancels all unresolved gains, and leaves a simple known-rate running term plus signed birth and seam terms. `hV` has a genuine routing curl. This comparison books no charge or endpoint movement.

## Exact accounting

Write the row (one signed D18 witness, or their contrast) as `c=a+Σ r_i e_i`. At fixed other tags and terminal `(k,h)`, let `z(L)=-Σ_{u=L_min}^{L-1} e_j(u)`, with `L_min=max(0,2−k−h)`. For all controls, `cʹ=c+(A−μ)z` therefore has terminal coefficient zero. Keep the actual gain rate in `A` even when the extractor used a center in `a`; this adds a pressure coefficient back on centered faces if `Δz` there is nonzero. `A` also includes known leaves, tagged receives, all tracked promotions and their exits. A promotion changes arity and hence changes the antiderivative's anchor.

For a full tag intensity `M` and its extension/birth forcing `F`, `M[(μ−A)z]=F(z)`, so the stationary observable is `M(cʹ)+F(z)`. The terminal antiderivative vanishes at the newest tag `(1,0,1)` and has zero full-tag birth correction in the unmasked hierarchy. A full-weight gauge need not vanish there. At the preceding prefix action the missing term is `+b z(x_new)` with the *post-promotion* prefix in `x_new`. In a stopped region, also retain the terminal payoff `−z(x_exit)`, or equivalently use a zero extension of the mask and keep the resulting jump on every crossing, including entry from the complement in a stationary calculation. These conventions are alternatives for the same flux.

In the stopped Bellman problem with discount `s`, use `cʹ=c+(A−s)z`. For every policy, `Vʹ=V−z` with boundary `Vʹ_exit=−z_exit` when the raw boundary is zero. The same follows locally from `cʹ+(A−s)(V−z)=c+(A−s)V`. Both extrema and the residual of every fixed continuation are unchanged. This applies at the true `μ` and throughout the bracket; the numerical pilot uses the common fixed `s=.565775` and does not use discount monotonicity on a signed contrast. It is also an obstruction for the corresponding exact occupation LP, where the correction is a linear combination of flow constraints. A pointwise interval formed by separately boxing the new running row and the transported flux need not be invariant and need not improve.

## Routing curl and a pressure-free representation

For two active directions `i,j`, take the circulation `e_i(x)+e_j(x+i)−e_j(x)−e_i(x+j)`. A scalar correction leaves it invariant. Only compare two directions with actual free pressures on the four edges: a zero column for an absent surplus donor is not a gradient constraint. If the apex parent and child pressures both occur, their coefficients must agree to cancel them with one degree difference. On the sampled faces their coefficients agree.

* `hV`: the apex and terminal child square has nonzero circulation on **50/50** tested squares for each sign and the contrast. On the maximal contrast square, `x=((12,5,1),(1,1,1),(0,1,1))`, the terminal rate was routed into higher rows, while the apex remains free. The contrast circulation has absolute value `41886438996851/171475000000000000`; some corrected edge has coefficient at least `6.1067851e−5` in magnitude for any scalar gauge (at least half the circulation on one of the two nonterminal edges if both terminal edges are forced to zero). This is a coefficient obstruction, not a stationary mass or endpoint bound. The terminal-only correction leaves the running receive-box width unchanged on average here; the continuation and exits still have to transform.
* `P4`, `Yt`: there is no higher routing. On live gain faces all free coefficients are differences of the scaled signed motif weight `W`, including the parent pressure. Thus `z=−W` cancels all of them simultaneously. The running term after full cancellation consists of the known transfer from supported smaller motifs plus the per-witness `(s−λ)W` (apply the contrast after this operation). This decomposition was independently replayed on every ranked face; do not transfer it across a support seam without the signed seam term. To exploit it beyond the invariant comparison, bound the **joint** weighted known transfer and its birth/seam counterpart using stronger physical occupation constraints. A separate bound on the pressure-free running term would not suffice: on the paired sample its ranges were `[5.856e−5,4.724e−4]` for `P4` and `[−2.402e−3,−4.126e−4]` for `Yt`, much larger than the corresponding integrated values.

## Exact ranked replay

Nonfounder descriptors `(L,k,h)` with apex `D∈{18,24,32,64,65}`, `k=5`; the other slots are listed in `seeds` in [replay.py](replay.py). This includes diagnostic shells at terminal degrees 4–8 and at apex degree 65. These are **comparison supports only**, overlapping booked old-flow and other work; no additional ownership is claimed. The following widths are unweighted face means of the *running* receive box, for `c=c_lo−1.2415 c_hi`; they are not stationary estimates. Each stopped interval is from the face ranked by maximal local width drop, with per-tag caps one degree and one arity above the initial state, and exit on the first cap crossing.

| Family | Faces | Static width, raw → terminal gauge | Exact matched stopped interval (rounded) |
|---|---:|---:|---:|
| `hV` | 50 | `2.29409e−5 → 2.29409e−5` | `[-2.08175e−8, 1.69208e−8]` |
| `P4` | 20 | `1.63477e−5 → 1.91859e−6` | `[-5.81759e−7, 4.71528e−7]` |
| `Yt` | 20 | `4.05164e−4 → 6.27320e−5` | `[-8.55153e−6, 7.99454e−6]` |

The replay asserts rational equality for both extrema of each witness and the contrast, including the full-weight gauge for `P4`/`Yt`. Deliberately dropping the exit terms shifts the apparent paired `Yt` interval by roughly `+1.7e−3`; this is a bookkeeping error, not a recovery. Exact gain from the coordinate change under matched constraints is **0**, so it cannot meet either complementary contrast target in [ratio-transfer/averaged](../../ratio-transfer/averaged/README.md) by itself.

```sh
taskset -c 8 nice -n 17 timeout 20s env OPENBLAS_NUM_THREADS=1 ./round2/.venv/bin/python -Bu round2/occupation-neon/pressure-gauge/old-flow/replay.py b3f18_lo.bin b3f18_hi.bin
```
