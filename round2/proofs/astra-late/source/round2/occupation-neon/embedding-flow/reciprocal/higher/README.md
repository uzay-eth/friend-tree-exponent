> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Original D18 hV on genuine reciprocal cells

**Yes.** The original allocated nonfounder `hV` charge, and its lower-minus-`θ`-upper contrast, are exact rational functions on *every fixed-arity cell* of the oriented chain `apex → middle → tip`. This includes two or three simultaneous high roles. For each high role use `y_i=18/D_i` and `a_i=1+y_i/18`. A common denominator is `Δ=∏_high a_i³`. No other pole is needed for this family. The common numerator has degree at most **7 per high coordinate** on the triple-mass block and **6 per high coordinate** on each live pressure block. Finite coordinates are exact constants. This is a cost reader for the physical measure; it does not solve or extend the shared triple LP.

## Original controls and allocation

[`hv.py`](hv.py) uses the pinned `Cert` layout and the [`Extract.weight`](../../../higher/extract.py) table lookup at a control. It reconstructs the original physical cubic, not the GMP verifier's extremal bounds for a control:

* For `D<18`, the old exact state is `(D−k−1,k,1)`. For `D≥18`, the control labels are `INF,T1,T2,T3` with `κ=min(k,18)`, and coefficients `(1−y)³,3y(1−y)²,3y²(1−y),y³`. The tensor product gives the weight of a chain. `INF` has no middle `hV` control and is padded with zero there. The entire allocated row vanishes for middle `k>17`. This matches the original `FULL hV` loop, including the absence of a lower-edge transfer at overflow middle arity.
* Promotion increments **actual** `k` at unchanged `D`; recompute capped labels only for table lookup. A promotion at middle `k=17` kills the old `hV` weight. Gain increments `D` only. A finite degree-17 receiver enters its exact `T3` weight at degree 18; it does not access a high-moment coordinate or an arbitrary control corner. Physical inverse degrees are always `1/D`, including for a capped or overflow label. For fixed high arity, the continuous control interval is `0≤y≤min(1,18/(k+1))`; a discrete rate evaluation requires the true integer `D≥max(18,k+1)`.
* Put `−λ hV` into the base for each witness, using that witness's original trial `λ`. Each gain retains the leaf-stock contribution and **each true internal sender**. The only external ports are the apex mother, and the off-path children at each of the three roles. If the original route to `P4`, `Yb`, or `Yt` is available on the pre-event census, retain only `count × original center × gain difference` here; the `(actual pressure − count × center)` difference was already allocated to the higher family. Otherwise the slope is the original `hV` gain difference, multiplying the aggregate of the actual live exterior reciprocals. No clipping, hindsight extrema, or subsequent credit allocation is folded into these slopes.
* Promotion keeps the existing difference times `(D−k−1)/D`. The `hV` row also owns three possible newborn level-four payments: `Yb` when the apex promotes, `Yt` when the middle promotes, `P4` when the tip promotes. Evaluate the newborn term on the **post**-promotion eligibility of the promoter, at its true degree, and on the virgin `(1,0,1)`. If a required finite slot is actually high, the term is zero even when the free coordinates carry nonzero tail controls. The newborn `hV` itself belongs to its old edge row; it is not added again here.

### Incoming `g → hV` transfers

Write the chain `(p,u,w)`. There are exactly two lower-edge allocations: from the apex donor, on the old edge `(u,w)`, and from the tip donor, on the old edge `(p,u)`. Each equals `(1/D_donor − original center_at_u) × transfer_coefficient`. The center is the matching `rpar/rch` for inner `u`, and `trpar/trch` for high `u` (`h=1` index `k_u`).

For an inner middle, the coefficient is the exact gain difference of `g` at that middle. This includes a `17→18` seam evaluated at true degree 18. For a high middle the original transfer is **not** that full gain difference. It is

```
B2(y_u) * (g(INF, w) − g(T1(u), w))/18 + B3(y_u) * Yp_ref(u,w)
B2(y_u) * (g(p, INF) − g(p, T1(u)))/18 + B3(y_u) * Yw_ref(p,u)
```

respectively. `Bj` are the original Bernstein basis functions, not the new storage basis. Apply the physical control mixture to the other old-edge endpoint in each expression. The terms are already polynomial, with zero contribution at `INF` and `T1` in the **receiver** coordinate. The singleton-to-edge `c2_f/c3_f` allocations belong to `g`, not `hV`; there is no direct `f→hV` term. This implements the signed physical identities underlying the original GMP `hV` loop (`round2/attacks/neon1/outputs/production-final/clone/round2/cert/cert_b3f_verify_gmp.c#L1320`, historical reference) and original proof (`round2/attacks/neon1/outputs/production-final/clone/round2/proof-b3f-motif4.md#31-exact-zero-sum-telescoping-transfers`, historical reference), before it bounds anything by Bernstein hulls or port corners.

`Reader.row` returns the raw physical witness drift divided by `S=10¹¹`. If `C_- = +raw_lower`, `C_+ = −raw_upper`, then `contrast(lower,upper,roles,θ)` returns `C_-−θ C_+ = raw_lower+θ raw_upper`. Both witnesses see the **same** live ports. In particular the upper mode's sign is applied exactly once. A signed cost on a supported cell is not a positive unpaid payment or an exponent sign.

## Denominators and moment controls

For `D→D+1`, `y→y/a`. Write the cubic in the receiver coordinate as `P=∑_{j=0}³ c_j y^j`, where `c_j` includes any other high axes. Its cleared gain difference is

```
U = ∑_j c_j y^j (a^(3−j) − a³)
L Δ_D P = [18 (U/y) − (k+1)U] / a³ .
```

`U` has degree at most six and the `j=0` term vanishes. Every `j≥1` term is divisible by `y^(j+1)`, so dividing by `y` introduces **no pole** (the code checks divisibility by `y²`). The formula also stays bounded at `y=0`, without subtracting two separately unbounded leaf rates. Products of simultaneous high *roles* use the tensor denominator; there is still only one receiver gain per event. Other terms have at most a single additional inverse degree or promotion factor in a given coordinate. Thus the original tensor cubic plus clearing (degree 6) grows to at most 7 in a mass coordinate, while a port slope is at most 6 because its rate is in the port measure.

The cubic pole can be genuine: on `((tail,2),(4,2),(4,1))`, the `θ=7/5` apex-mother slope numerator at `y0=−18` is `−2025082539/31250000000`. No inner weight is evaluated there; this substitution only proves algebraic noncancellation. This pole lies outside every physical cell. Saturating the arity and crossing the finite seam create *separate control labels*, not rational factors. The extra transfer references are stored constants times Bernstein controls, not evaluations at a shifted reciprocal. There is **no** `D−1`, neighbor-shift, leaf-stock cap, or rescaled tail variable in these original `hV` row terms. This calculation does not establish the same denominator for a different family, inverse transition test, or altered score.

For the owner's storage convention, `W=∏_high a_i^m`, call `moment_controls(row,m,bounds)` on one physical reciprocal slab. It multiplies the numerator by the exact polynomial `W/Δ` and converts to tensor Bernstein coefficients over that slab; full physical intervals are the default. This requires `m≥3` in general. The required axis degrees are at most **`m+4` for mass and `m+3` for ports**; the existing degree `2m+2` accommodates both. At `m=3` the existing storage degree is 8 (nine coordinates per high axis) although an independent degree-7 mass net suffices. There is no degree in the pressure variable beyond the affine factor: apply that net to the triple mass divided by `W`, and each port net to that **same triple's** true reciprocal-flux measure divided by `W`. A cell with all exact finite roles uses one control, and the oriented chain carries symmetry factor one. This operation alone does not identify missing mass or flux outside a represented cell.

## Bounded replay

```sh
OPENBLAS_NUM_THREADS=1 taskset -c 11 nice -n 19 /usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M' \
  timeout 30s ./round2/.venv/bin/python -Bu \
  round2/occupation-neon/embedding-flow/reciprocal/higher/check.py \
  b3f18_lo.bin b3f18_hi.bin
```

Replay output: `exact fixtures {'rows': 80, 'terms': 1122, 'ports': 176, 'extract': 6, 'controls': 1700, 'zeros': 8, 'moment': 8} nonzero branches 19 genuine cubic pole numerator at y0=-18 -2025082539/31250000000`. Actual final `wall=23.84 cpu=14.18+0.95 rss=84612`. The numeric fixture uses a separate `Cert.embed` path with exact physical `D,D+1` values, split by charge term. All four centered routes and all three newborn branches are nonzero in at least one fixture. It also checks inner `Extract.row` parity, `17→18`, every one of the three high roles, simultaneous tails, middle `17→18` **arity** exit, middle overflow zero, and exact evaluation past degree `10⁹`. It checks both Bernstein conversion and stored-moment contraction. Only sparse reads and header checks were used; full blob authentication is inherited from the release replay (optional `--authenticate` is available, not run here). No LP, big scan, new booking, or shared triple file edit was performed.
