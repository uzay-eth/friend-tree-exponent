> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Claim 21 extension: size 20 and finite high-parent occupation

The physical occupation replay gives the outward bracket **`[0.565676634491, 0.565771874464]`** when added to the [cumulative size-20 ledger](../renewal/next/README.md). It adds stationary signed reward of at least `5.520306761e−7` on the lower witness and `2.799995635e−7` on the upper. The fixed denominators give inward endpoint moves of at least `4.068523300e−7` and `2.048023029e−7` before endpoint rounding. The mixed edge is the main contributor; sibling pairs are positive but do not account for a large fraction of the residual budget.

## Exact composition

| Term | Lower witness | Upper witness |
|---|---:|---:|
| Original trial `λ` | `0.565673` | `0.565775` |
| First-body floor | `2.512846393e−6` | `2.025752646e−6` |
| Size-12 higher floor | `1.276902849680e−6` | `1.501029458674e−6` |
| Joint size-20 suffix | `58961819710386489/10^23` | `46635843154908967/10^23` |
| Joint high-parent addition | `55203067615510485475692014463/10^35` | `139999781774155931356847859909/(5·10^35)` |
| All-degree ceiling | `1.356833021311` | `1.367169995210` |
| Exact endpoint | `76752873706217324196974475692014463/135683302131100000000000000000000000` | `386753165449918989299233643152140091/683584997605000000000000000000000000` |

Sum the four rewards and divide by the single potential ceiling for that witness. Add the quotient to the original lower trial and subtract it from the original upper trial, then round outward at 12 places. The reviewed extension conservatively rounds the first two rewards down and the ceiling up relative to [Claim 21](../release/README.md). It does not use an earlier outward endpoint as a trial.

## Finite-degree floors

A virgin child has one leaf and no core children, hence degree two. Each line is an **exact** degree event, not a threshold or an asymptotic fit. Columns are stationary expectations per physical core with at most ten children. The pair count requires at least four core children in total and counts each unordered pair of virgins once. The reward columns price the sum of the mixed edge and sibling-pair terms in normalized signed witness units. All table entries are rounded down.

| Parent degree | Hub mass | Virgin edges | Virgin pairs | Lower reward | Upper reward |
|---:|---:|---:|---:|---:|---:|
| 18 | .0040065499839020 | .0018337532593034 | .0001470372798880 | `8.27960616e−8` | `1.5893277e−9` |
| 32 | .0015422659941728 | .0008014182916395 | .0001112680873731 | `8.6374088e−9` | `8.2739333e−9` |
| 64 | .0004993082834405 | .0002863716702961 | .0000583778279790 | `4.125305e−10` | `1.0268948e−9` |
| 128 | .0001642405677524 | .0000999828225428 | .0000251719748303 | `1.82219e−11` | `8.77708e−11` |

The main replay visits all integer degrees 18–191 on each side. The four exact-degree results are slices of that same support, not extra allowances.

| Reward over degrees 18–191 | Lower witness | Upper witness |
|---|---:|---:|
| Mixed `g` alone | `5.432244944e−7` | `2.694676445e−7` |
| Virgin-pair `hC` alone | `8.7345295e−9` | `1.02830012e−8` |
| Joint | `5.520306761e−7` | `2.799995635e−7` |

## Joint occupation without a product law

Follow a physical newborn with a parent. The projected state `(d,k,j)` is its degree, number of core children, and number of virgin core children; it starts at `(2,0,0)`. Put `l=d−k−1`. The only changes in the projection, with hazards conditional on the full physical filtration, are:

* Own gain → `(d+1,k,j)`, with rate in `[l+j/2, l+(k+1)/2]`. Every virgin contributes exactly `1/2`; each other neighboring core contributes at most `1/2`.
* Own promotion → `(d,k+1,j+1)` at `l/d`. The new arm is virgin, and the parent's degree stays fixed.
* Loss of one virgin → `(d,k,j−1)` at `j(3/2+1/d)`. For each virgin the loss comes from one leaf sending a gain at rate one, its own promotion at rate `1/2`, or a send from the tracked vertex at rate `1/d`. Virgin status cannot recur: degree never decreases, and a degree-two core which promoted its only leaf cannot lose that core child. Arrows out of the virgin toward the hub are already included in the hub-gain rate.

All other physical changes are allowed. In particular, child and parent arrows remain live, and the unobserved part of the gain rate may depend on all the past. This calculation tracks arms together; it needs neither a fixed-age coupling loss nor an independence assumption conditional on future success. Both `d` and `k` are monotone. Kill only at the upper degree and arity cuts, set continuation to zero there, and process `d,k` downward and `j` upward. On a row with reward `r`, promotion rate `b`, virgin loss rate `a`, and gain rate `g`, check

```
(HI + g + b + a) V(d,k,j)
  <= r + g V(d+1,k,j) + b V(d,k+1,j+1) + a V(d,k,j−1).
```

The inequality is affine in `g` and is tested at both endpoints using exact rationals; the value is rounded down on the `10⁻²⁴` grid. This is a bounded discounted subsolution for every adapted exterior. Multiply the newborn value by `LO≤μ` and use the physical stable-age law; `HI≥μ` and the reward is nonnegative. For an exact-degree column the reward is zero before the stated degree and the process is stopped on the next gain. This accounts for final-degree restriction directly.

## Physical charge and disjointness

Use §3 of [Neon’s pinned GMP charging decomposition](../../../../../../reviews/cold/astra/occupation-neon/check/neon-proof-pin.md#motif-definitions-and-transfers), source commit `cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`, path `round2/proof-b3f-motif4.md`, blob `3bbda54c3cbdbd235a6cf4af8ecf80ca75d16ae7`. The tracked verifier (`round2/attacks/neon1/round2/cert/cert_b3f_verify_gmp.c`, historical reference) and GMP receipt archive (`round2/attacks/neon1/outputs/neon1-production-rehearsal.tar.gz`, historical reference) preserve the code and the two pinned production runs. Use exact physical Bernstein weights. A virgin edge keeps its two included directed reciprocals. It receives the `f` transfers from both vertices. At the tail parent it sends the prescribed portion `T = 3t²(1−t)(F₀−F₁)/18 + t³Yref`, `t=18/d`, of each external gain derivative to a higher row. Thus the external reciprocal sum multiplies `Δg−T`, while the exterior centers multiply `T`. Use the full `[0,k/2]` box for that sum. Both possible promotion terms include the newly created higher motif allocated to this edge. No first-body charge is claimed.

For a pair of virgins at a tail hub, `hC` includes both directed contacts with each virgin and the two transfers from the incident `g` rows. There is no outgoing Level-4 transfer or Level-4 newborn on this face. The untracked exterior reciprocal sum at the hub lies in `[0,(k−1)/2]` and multiplies its physical `hC` gain derivative. Check the signed row floors before multiplying by `j` or `binom(j,2)`. The joint Bellman solve can share a gain-rate choice across distinct rewards; only the better of joint and summed separate bounds is used. Sixteen star fixtures (eight per witness, each at both exterior endpoints) recompute the drift through full potential differences and check exact agreement with the sum of these rows, first-body terms, and omitted exterior transfers.

| Booked reward | Current charged-row support |
|---|---|
| First body | Nonfounder `f`, exact rows through degree 64 and all-stock continuation in five bands. |
| Size-12 higher | Top-owned `hV,hC,P4,Yt,Yb,S3` on complete planted fringes of size ≤12 in the interior domain. |
| Size-20 joint | The same six higher families on complete fringes of size 13–20, with 3–4 cores and every degree ≤16; plus the minimum charged `g` per actual inner→inner edge from the root (both degrees <18), throughout the retained size ≤20, cores ≤4 solve. The `g` floor can pay before size 13. |
| High-parent joint | `g` from a nonfounder of degree 18–191 with 1–10 core children to a virgin `(L,k,h)=(1,0,1)`; `hC` at the same kind of host with 4–10 core children, one row per unordered pair of virgins. |

The inner-edge floor never pays a high-parent mixed edge. A planted degree-18 owner cannot fit in the size-12 class, and a high-parent pair has at least five cores at its owner and is beyond the size-20 degree seam. The exact-degree slices and separate edge/pair prices are checks on the same high-parent support. The [separate audit](../../../../../../reviews/cold/astra/checks/neon-high-parent-sol.md) records the competing support and why it is excluded.

## Replay

From the repository root, set `PYTHON`, `D18_LOWER` and `D18_UPPER` as in the [Claim 21 instructions](../release/README.md). The [extension manifest](inputs.txt) pins the added numerical sources, the separate audit script, and the base verifier and its manifest.

```sh
OPENBLAS_NUM_THREADS=1 "$PYTHON" -Bu round2/occupation-neon/high-parent/verify_release.py \
  --lower "$D18_LOWER" --upper "$D18_UPPER"
OPENBLAS_NUM_THREADS=1 "$PYTHON" -Bu round2/checks/neon_high_parent_sol.py "$D18_LOWER" "$D18_UPPER"
```

The release verifier authenticates both original payloads, re-extracts the first-body and inner-edge rows, reruns the size-12 and size-20 solves and all-degree ceiling, and prices the joint high-parent contribution on 1,740 row descriptors per sign. The separate audit independently checks the charged rows and high-parent recurrence and redoes the exact final arithmetic. The exhaustive GMP runs and the prior are inherited.
