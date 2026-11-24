# Cold audit: residual pressure memory after signed EDGE

**CORRECT in scope.** The exact joint price strengthens `02d0bd87c64b53a2bda2004a4af65de9142dda13` to **`[0.565678056632067, 0.565770469400619]`** outward. This inherits the original D18 GMP B verifications, the physical mean/stable-age identities and incoming interval `[0.565673570322, 0.565775]`. The frozen EDGE baseline is `[0.565677980686114, 0.565770513279953]`. The inward moves from its **exact** endpoints are approximately `7.5956e−8` and `4.3879e−8`. No board entry was made.

## What is left after EDGE

Use signed charges at the original `λ`, with `ε=+1` lower, `−1` upper, `δ=ε(s−λ)>0`. Let `Q_x` be the recharged pointwise upper allocation used by EDGE. Suppressing zero continuations, the two remainders are

```text
R_f(x) = F(x) + δ(M−Q_x) + (A−s)z(x) + b V(post(x),J)
R_g(e) = G(e) − G_res(e) + (A_edge−s)V(e)
```

The `f` continuation sees the new edge exactly once through the parent's promotion; its endpoint is **post**-promotion. Each old edge sees the two gain successors and two promotion successors of that surviving edge. Newborn vertex `z−z(J)` is zero, and founder edge births are excluded. Outside the `V` cut both the value and its old-edge generator vanish by coordinate monotonicity; the clamped census tail has no retained new edge. Ignoring the unpriced remainders is valid because their lower bounds are nonnegative. The leftover from the original potential allocation inequality is not priced separately.

The paid masks are 170 `f` states through `D=18`, 54 `g` states with both `D≤4`. `f` and `g` are distinct charged rows. The `g` mask shares **actual edges** with the frozen size-20 inner reward: the uniform inner debit must be subtracted on **every** paid edge, including those for which the size-20 trajectory might never book a reward. It has empty intersection with the high-parent mixed-edge mask. All six higher-body charge families, including `hC`, and the frozen high-parent payment remain intact. A different argument for a row's occupation is not a different reward. The old first-body-only residual or another EDGE `f` price cannot be added to the joint price.

The companion (`edge-memory-sol.py`) rebuilt all 2,079 finite `f` rows and the clamped arity tail from the authenticated signed witness, and independently contracted the 54 paid physical `g` rows and their shared-rate corners. Omitting the newborn edge in the paid `f` residual makes a corner negative in 108 lower / 98 upper states. Adding newborn continuations to the old-edge equation would add spurious value on 45 paid states each side. Omitting the same-edge **debit** restores an already-reserved positive amount on all 54 states. Omitting the same-edge **objective** instead is safe: verified memory `f`-only prices `1.016956395e−7` / `5.988262812e−8` alone give the outward control `[0.565678055600371, 0.565770469479848]`; the joint prices below are used in place of those floors.

## Why the physical law is feasible in the dual

Use two measures: the uniform-core vertex law and the parent→child incidence law per core, each of mass one. The latter has endpoint counts `k` at the sending cell and `1` at the receiving cell. The reciprocal-weighted counts are `kp=1−b−p` and `q`, with `p=1/D_parent`, `q=1/D_child`, `b=L/D` for the vertex in question. Reversal gives the downward and upward input marginals, `πk=1`, `πρ=1−μ`, and birth compensation gives `πb=μ`. The census balance is discounted Dynkin at `μ` from `J`; one common interpolation of `μ` in its incoming bracket lifts it exactly. Founder contributions vanish in the normalized mean; the bounded tests and the finite core-child moments in the flow proof (`occupation-neon/flow/README.md`) justify the weighted reversals. No independence or equality between laws of adjacent descriptors is used.

On a paid edge the full receives at the parent and child are `u+q` and `p+v`. Outgoing ports sum to the parent's `kρ`, incoming ports to the child's `ρ`. The *same* vertex `ρ` is used in the price and census gain. On exact child bins, the gain weight `w=(L+p+v)/D_child` uses that same edge input; `q` is fixed so the memory source is exactly `w q/(1+q)`. The other edge bins can independently mix endpoint reciprocals and child gain corners: the cell rectangle covers `w` by `b≤w≤(1+b)/2`. For the memory source at a coarse child bin, the chord/tangent bracket for the concave `q/(1+q)` is affine in `q`; interpolate between the two at the actual `q` and lift to reciprocal corners with that same interpolation weight. Independently mixing the actual `w` within its cell bounds then preserves the product. These are marginal constraints, not a prescribed law of the two receives.

For a finite parent, track `B=k/2−Σ_child q`. Own gain leaves `B` unchanged; own promotion adds one newborn with `q=1/2` and increases `k/2` by `1/2`, also leaving `B` unchanged before the census test changes. A child's gain increases `B` by `q²/(1+q)` at rate `L_child+ρ_child=w/q`; the product gives the source just described. The finite-cell discounted balance for `B·1_cell` has zero newborn value and includes exactly this source. At own gain its coefficient is `(L+p_parent+k/2)B−B²`. The physical `(B,B²)` on each interval `[a,b]` is the quadratic Bézier combination of `(a,a²)`, `((a+b)/2,ab)`, `(b,b²)` with weights `(1−t)²,2t(1−t),t²`. This also preserves linear receives in the objective. Both LPs keep the identical hull and columns; removing the 20 `B` rows is the only change.

### Unbounded columns and rational repair

- The degree tail is not killed: its reciprocal range includes zero. For fixed overflow `k<12` keep `b=1−(k+1)p`, keep promotion into the next arity, and make gain invisible to the census test. For `k≥12`, use the four corners of `0≤p≤1/19, 0≤b≤1−13p`; gain and promotion are both invisible. Keep `kp=1−b−p` exactly, and use the excess-stub column for `k−12` in the unweighted counts. The saturated point still carries the downward receive in the occupation column; the separate `tailup` column carries the entire incoming **child** sum in both the reversal and global receive rows. Its weight is at most `πρ<1`; the stub weight is at most `πk=1`. The tail `w` lift uses the endpoints `b,(1+b)/2` and a common interpolation fraction across the physical `(b,p)` corners, preserving the `child_w` marginal. No finite cap is imposed on actual arity or receive.
- The fixture rebuilt **37 vertex and 35 edge** local contributions as exact convex combinations, including nonzero saturated and unsaturated tail traffic, a degree `10^10` vertex, unequal donors, and nontrivial memory sources. For each proposal, the separate verifier rounds the equality dual, computes the maximum positive violation over every column of each class, and checks every repaired inequality in rational arithmetic. Class totals are at most one, so subtracting those four defect maxima from `y·rhs` is conservative. The checks cover all **9,674** columns in each model: `3,528` occupation, `6,144` edge, one stub, one tailup, including all zero-price and overflow columns. HiGHS only proposes the multipliers.

## Exact composition

Use the unchanged EDGE trial and negative local margin, and add the frozen reserve and the **one** joint lower price in signed numerator units. The old all-degree denominator ceiling **plus** the one-sided census correction remains necessary; on the lower sign also add `sup V` since there is at most one counted edge per core. On the upper sign `−ΣV≤0`. Memory is a dual certificate for occupation of the remainder; it does not change the potential or its denominator. Both net numerators are positive, so division by the full ceiling has the requisite outward direction.

| | Lower | Upper |
|---|---|---|
| Flow-only joint price | `3333116239993075574694023/36450000000000000000000000000000` | `231029505304032388502591/4500000000000000000000000000000` |
| Memory joint price used | `18402664634233592780094779/178500000000000000000000000000000` | `101984612445685206436931/1700000000000000000000000000000` |
| Full denominator | `1357493815758092773/1000000000000000000` | `273436004214812983/200000000000000000` |
| Exact endpoint | `2586244278310822397274177587149/4571936719109803018500000000000` | `6574835699364880591816419305857/11621030179129551777500000000000` |

## Replay

From `[local path]`, the following passed. The producer's combined manifest (`occupation-neon/interpolate/compose/residual/edge-projection/combined/inputs.sha256`) returned six `OK` entries, and the EDGE manifest returned 24. The companion uses the original full-payload SHA-256s and checks the EDGE, projection, canonical root memory and transitive flow sources. It reuses the frozen Bellman tables' solver and the model column builder; its physical lifts, residual contractions, proposal repair and endpoint arithmetic are separate. Keep `pressure-memory/marginal/replay.py` at the repository root and include its `round2/occupation-neon/flow/check.py` dependency in any clean archive.

```sh
sha256sum -c round2/occupation-neon/interpolate/edge/inputs.sha256
sha256sum -c round2/occupation-neon/interpolate/compose/residual/edge-projection/combined/inputs.sha256
ulimit -v 10485760
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
P=round2/.work/venv/bin/python
A=round2/checks/edge-memory-sol.py
nice -n 15 taskset -c 2 timeout 320s "$P" -Bu "$A" b3f18_lo.bin
nice -n 15 taskset -c 4 timeout 320s "$P" -Bu "$A" b3f18_hi.bin
```

Each companion run printed the exact fractions above and the corresponding outward bound, with rational repair of the memory joint numerator `43576407219905221/178500000000000000000000000000000` lower and `242314793563069/1700000000000000000000000000000` upper. The producer source reviewed with 2f74 was `combined/replay.py` SHA-256 `7d536fab426f568f338af51f6d2f6bd104683ddfb9ff17f12315a5561cd14239`; its manifest hash was `6abcfc4bc9dc9811992f1afc425aeecacc3bbb1ea7049b4be55dabe9223b577e`.
