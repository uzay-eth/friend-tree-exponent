> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Signed high-core block at the analytic seam

At `s=141431/250000`, under the physical uniform-core law `π`, the original gold pair `Φ_lower+Φ_upper` gives the following **per-core** enclosure for the block defined below:

```
−5.829193e−7 ≤ π H ≤ −1.116081e−7       (width < 4.714e−7).
```

This same enclosure holds for the corresponding polynomial-retained block `H̃`. Both have support at all `D≥2000`, with the same pre-event arity cutoff `K=106`. There is no degree cap in either target. This result evaluates the high core's own singleton drift and its entire newborn contribution. A core of this degree cannot be an interior vertex of an extant higher motif; it can still be a tip or a donor. Those other contributions are tracked in the signed complement, not estimated here.

## Same support, original raw coefficients

Let `F(D,c)` be the sum of the two analytic singleton weights and `J=(2,1)`. Set `N(D,c)=F(J)+g_lower((D,c+1),J)+g_upper((D,c+1),J)`. A promotion at a high core creates exactly this singleton and edge among nonzero newly created weights: any larger new motif would have that core as an interior vertex, violating the pinned support masks. For a high receiver let `ρ_K=Σ_{u∼core o, c_u≤106} 1/D_u`. The exact target is

```
H(o) = 1{D≥2000, c≤106} [−s F + (D−c+ρ_K) Δ_D F
                         + (1−c/D)(Δ_c F + N)].
```

The indicators are on **pre-event** roles. In particular a promotion at `c=106` uses the real successor `c=107`; a gain at `D=1999` is in the complement, even when its successor has degree 2000. The sum of `H` on a tree equals the selected part of the parent ledger even though the ledger splits gain receipts among roles and this definition places them at the high receiver. Physical incidence reversal gives equality of their integrated expectations, with no rescaling. `H̃` replaces only `Δ_D F` with the parent ledger's second-order reciprocal expansion on exactly the same support.

For `x=1/D`, write `F=a+u(c)x+v(c)x²` with `u(c)=−kc−h` and `v(c)=−jc²+hc`. The constants come from the two original headers. `N=n+z_c x`, with `z_c` depending only on the clamped post-promotion arity. Two affine reads at the knot and twice the knot recover the actual `INF` value `n`; each `z_c` then needs only one knot entry in the column for `J`. The needed joint moments are thus degree and arity moments plus one *linear* incoming sum. No joint ages or higher-dimensional donor-state positions are priced for this block. In particular this calculation does not multiply separate physical marginals for degree, arity, or neighbors.

## Two-sided physical pricing

For the solve only, use the full physical incoming sum `ρ` and split the reward into disjoint bands. Compute the signed finite band `2000≤D<10000, c≤55` using the polynomial gain. From birth `(2,1)` to the band and inside it, let `G` be the actual gain and `b=1−c/D` the actual promotion rate. The root only moves to `(D+1,c)` or `(D,c+1)`. The future discount value satisfies the robust inequalities obtained by taking the minimum/maximum of

```
[ h(D,c,G) + G V(D+1,c) + b V(D,c+1) ] / (σ + G + b),
```

with continuation zero at `D=10000` or `c=56` and reward zero below 2000. Here `h` is the complete signed bracket defining the finite-band surrogate. Check all four corners `G∈{D−c,D−c/2}` and `σ∈{.565702,.565749}` after rational rounding; the cross-multiplied constraint is affine in each rate and their roles are shared between reward and transition. This comparison holds for the actual possibly adaptive exterior filtration. Discounted compensation from the actual newborn gives bounds for its expected signed lifetime, then the stable-age identity multiplies by `μ`. The multiplier uses the minimum/maximum of both exact gold endpoints times the corresponding bound, so a negative value is not priced with the wrong endpoint. We keep `s=.565724` fixed in the reward throughout; `σ` ranges only to cover the physical age discount `μ`.

For the unbounded suffix `D≥10000, c≤55`, use the exact reward with full receive. Extract `r₀+r₁(c)/D`, where

```
r₀ = n−s a
r₁(c) = z_c−k+(1+s)h + c[(1+s)k−n].
```

The error is bounded by `(A+B c+H₂ c²)/D²`, computed from `F,N` as follows. Write `Z=max|z_c−k|`. For the gain after extracting `−u(c)/D` the error numerator is at most `(kc+h)(1+c)+3(jc²+hc)`, because `D−c≤G≤D`. For the promotion add `c Z+|−j(2c+1)+h|`, using `0≤b≤1`; for discount add `s(jc²+hc)`. Sum the resulting nonnegative coefficients to obtain the envelope in the code. Negative linear coefficients are bounded on this same degree suffix; positive ones only occur at higher arity and use the intersection bound `min(π(D≥10000) upper, π(c≥min positive class) upper)`. Do not multiply these probabilities. The same quadratic envelope and a geometric arity tail bound pay `D≥2000, 56≤c≤106`, with no overlap with either band.

To integrate the quadratic envelope, use the common-age proposal domination from the parent proof. The focal gain proposal chain has rate `d−1/2`; an independent rate-one Poisson proposal clock dominates promotions. They are independent **drivers**, not independent physical marks, and are evaluated at the same independent exponential age at rate `.565702` for an upper bound. The degree hitting transform is a product of `g/(g+.565702)`; conditional on the hit before observation, the time through degree `d` is tilted exponential with rate `g+.565702`, and the residual age is independently exponential at `.565702`. Their complete homogeneous moments, times `i!`, give the first three factorial moments of the proposal count on the hit event. This bounds `1,c,c²` jointly with the degree mask. Geometric moments pay the high-arity band without assuming any degree dependence there.

Finally compare the two **physical expectations on identical retained support**:

* `|π H−π H̃| ≤ 1.947e−13` per core. Use the two derivative remainders at each high receiver, total gain at most `D`, and the degree-hit upper mass. This bound applies unchanged to a restriction to either band, or to full receive. The solve combines a polynomial finite band with exact suffix bounds; the common interval allows twice this error to cover both target variants.
* Replacing full receive by `ρ_K` changes the exact expectation by at most `8.277e−26`. Re-root the omitted sends at the actual `c>106` donor; its total rate to cores is at most `c/D≤1`. The supremum of the high receiver's `|Δ_D F|` times the donor's geometric tail pays the difference. When working with the parent ledger this send belongs to arity overflow, not to the retained complement.

| Physical units, per core | Rational replay (display rounded) |
|---|---:|
| Signed finite band (polynomial, full receive) | `[-4.119973435081e−7, -1.897165947688e−7]` |
| Signed degree suffix (exact, full receive) | `[-1.709214090652e−7, +7.810792947128e−8]` |
| Absolute arity-band error | `8.427050801076e−14` |
| Common raw/retained width | `4.713110343555e−7` |

## Signed complement and replay

The parent retained score is `π R̃=π H̃+π C̃`, with `C̃` defined by deleting precisely the singleton discount, singleton gain, singleton promotion change, and two newborn weights at the selected high action; all other terms keep their **original sign and ownership up to incidence reversal**. Thus it includes the low-action singleton block, old `g` and higher-motif discount/gain/promotion, newborns not in `H̃`, and high-donor sends into other receivers when not already in `H̃`. In terms of the raw score, define `C` the same way with the exact gain and keep the dropped arity terms separate. The parent `3.569e−7` raw-to-retained bound applies to `R̃` as a whole. There is no need to add the tiny comparison error separately if substituting the certified interval for `π H̃` into that whole-score formula. The sign of `π C̃` is still unknown. In particular, this signed block is not an extra credit against the original `f/g/FY` bounds at their different trials.

```sh
taskset -c 5 nice -n 19 python3 -IBu round2/occupation-quintet/true-drift/moments/check.py
```

The replay checks both Bellman inequalities on the exact integer grid, exact reward arithmetic at rate corners, and signed complements of directed ordinary births using independent connected-subset enumeration at `D=1999,2000,2001`. It also checks a high receiver with an overflow donor, and displays fixtures at `4096` and `65536` with the signed complement kept separate. It checks the original detached source pin and reads 87 sparse words per payload including fixtures, without scanning any whole table.
