# Cold audit: reciprocal oriented chains and original hV tail price

**Verdict: the stationary shared receive-field projection is sound in both extension directions, and the CLI embeds the original unbounded-degree `hV` charge exactly on its selected support.** There are two reproducible direct-call footguns (unaligned cutoff and missing base edge ports), plus the usual finite founder boundary. This audits the current `reciprocal/{chains,price_chains,check}.py` and the local `reciprocal/higher/hv.py`; `8d8f` did not resolve as a local Git object or through the configured GitHub repository. The verifier pins the actual files by SHA-256. The future `bands/` scheme is outside scope.

## 1. Pre-existing field: before versus after the edge

For a physical edge `(u,v)` and stored edge coordinate `j`, write `q=B_j/W_edge` for its nonnegative tensor selector. A prepend has one eligible slot (mother of `u`); an append has `k_v` eligible slots, one for each actual core child of `v`. Write `c` for that multiplicity and `n` for the number resolved by the chain mask. The disintegrated terms integrate **the same physical edge** once per resolved extension. Multiplying a stored chain column by the new axis's `W` cancels exactly its extra denominator, leaving `q`; this works in either direction, including at overlapping slab boundaries with the same atom assignment.

For each old exterior receive field `F_l` (a sum of actual reciprocal sender degrees), the essential reconstruction is:

| Extension | Old selector field reconstructed from each resolved chain | All other old receive fields |
|---|---|---|
| `(a,u,v)` prepended | The old mother flux into `u` is `1/D_a`; middle has no exterior mother port. | Old roles shift by one and the old off-edge child sets are unchanged. |
| `(u,v,z)` appended | The old child flux into `v` is `1/D_z` **plus** the middle's spare-child field. | Old roles keep their indices and their exterior sender sets. |

Thus `given_l` in `connect` is precisely the covered integral of `q n F_l`, not just the flux from the chosen donor and not the uncorrected chain port. With `V=∫q(c−n)` and `E_l=∫q F_l`, the row pair is

```
0 ≤ c E_l − given_l = ∫q(c−n)F_l ≤ (slots_l/2) V.
```

The cap uses the **old** field's slot count, not the new spare count. In a fully covered append with `k_v>1`, each extension carries the *same reconstructed old field*, so the sum is `k_v E_l`. Using `E_l` instead, or restoring only the selected child, would break that identity. The selector has an additional pair of rows: the known selected reciprocals sum only the identified donor fluxes, and the unresolved part is bounded by the maximum of `1/H` and the reciprocals of every omitted explicit donor degree. Omission can be caused by **any** of the role, arity or wide-count restrictions; `1/H` alone would be unsound for an omitted low donor. Counting only the mass of one chain per census would also be unsound when several actual placements have identical labels.

**Finite founder effect.** An ordinary founder omitted as a planted apex is an unresolved prepend donor. If her reciprocal exceeds the computed missing-donor cap `u`, the selector's upper inequality needs the finite excess `q(1/D_founder−u)_+` on that edge. This is the only boundary residual in the independent checks. Since `u≥1/H`, a positive excess requires `D_founder<H`, with at most `H−1` direct children then. The excess vanishes under stationary normalization. The pre-existing-field transport pair remains valid on the finite image: the missing founder slot is still counted as missing and the old port retains the real send. The supplied `check.py` models a planted exterior donor of degree 23, so it is useful coverage but not a test of the low-degree founder exception.

## 2. Masks, one-jump balance and storage bounds

- The first two roles require at least one core child. The tip may be virgin or bare. Each role has its own arity cap; a high role must also pass the base high-arity cap `K`, whereas an explicit role has the full feasible `k< D` before the chain cap. The wide count includes **explicit** degrees `>inner` as well as high cells, and only listed `wide_roles` can be wide. A masked-out occurrence is zero in the tested chain image; this includes the case of a fixed arity becoming high overflow solely on a degree gain.
- Degree and arity never decrease, so an already existing chain does not re-enter a disallowed arity or wide-count mask. The balance assembly correctly keeps exits and does not confuse an existing-chain change with a new chain. On promotion of `v` on edge `(u,v)`, the new chain `(u,v,new)` enters with *postpromotion* `k_v`, same `D_v`, and a virgin tip. Its pre-event edge lies in the base census whenever the target chain is retained; in particular a high retained `k_v+1` implies pre-event `k_v` is represented. Each existing chain's loss/gain uses that chain's literal receive fields and both true neighbors of a middle gain.
- High-axis closure uses global `H/D` and the source-axis weight; the within-tail difference cancels the unbounded leaf stock and also covers changes across reciprocal slabs. At `D=H−1→H` the changed coordinate contributes raw test `1` with **no** high weight at the source. Nonuniform slabs just restrict the Bernstein controls to their interval. The supplied lifts put an exact cut point in the lower reciprocal interval; in particular `D=18` belongs to the selected `D≥18` band.
- One ordered ancestor pair per tip gives total chain mass at most `1/μ`. The receive groups inherit the same bound with the per-occurrence `1/2` or `max(K,H−2)/2` factors, since the chain caps only restrict a base cell. Tensor selectors sum to `1/W≤1` and the mass lanes sum to the unsplit atom. Thus inherited rational group repair applies to the chain columns on the physical image as well.

**Direct-call option pitfall:** `ChainModel(..., project=True, edge_balances=False)` with default sharing constructs base edge mass but no base edge port blocks, and `connect` raises `KeyError` on the base mother port. The price CLI's `--no-chain-balances` disables only the chain rows, so it does not hit this failure. Reject the incompatible constructor flags explicitly if exposing them to callers, or allocate base edge ports independently of edge balance rows.

## 3. Exact physical D18 hV embedding

The CLI selects one high role and uses the other roles' **exact** degrees at most `inner`; high degree means selected threshold slab(s), not a separate tail-price label. `minimum≥18`, `H_model≤18` and `m≥3` are enforced for the physical charge. The mapping to the reader is `y=18/D=(18/H_model)t`; the reciprocal slab endpoints are rescaled in exactly that way, while the atom still has `W=(1+1/D)^m`. The cubic old table weight is evaluated at the physical `D` and the physical shifted `D+1`; clearing its denominator contributes `(1+y/18)^3`. The original inner degree-17 successor is exactly the degree-18 `T3` value. The high middle's two incoming `g→hV` terms use their original `T2` divided difference and `T3` reference arrays, not a replacement by the full edge gain difference. Promotions use postpromotion ownership for `P4/Yt/Yb`; middle arity `>17` is an original zero row. Each mode uses the original sign: lower is raw lower, upper is minus raw upper, difference is raw lower plus `θ` raw upper, and sum is raw lower minus `θ` raw upper.

`moment_controls` multiplies each base or slope numerator by `W/Δ` before local Bernstein conversion. Degrees `m+4` (base) and `m+3` (slope) fit `2m+2` for `m≥3`. Each slope is attached to its own **chain** receive block, not an edge projection or marginal rate. Both witnesses contract against the same atom, including the same ports; no symmetry factor is needed for the ordered path. The first slab extends to `t=0`, so this is an exact unbounded-degree price on the requested one-high-role, arity-masked support, without an upper truncation or a claim on the complement. The authenticated column check includes all three roles, all signed modes, the degree-18 boundary and degree `10^9+7`.

**Direct-call cutoff pitfall:** `high_charge` filters whole slabs by their upper endpoint. It does not insert or require the cutoff `H_model/minimum` as an existing boundary. A direct call with cuts at `D=22` and `D=18`, requesting `D≥19`, succeeds but omits a nonzero `D=19` charge while still pricing `D=22`. The CLI inserts the exact cutoff, and its mass/reciprocal comparisons use that first slab. A direct caller asking for the entire band should align the boundary and place equality on the selected side, or the API should explicitly reject unaligned cutoffs. For custom very large arity caps, respect the reader's additional physical edge `D≥k+1`: `moment_controls` rejects a passed interval extending past `y=min(1,18/(k+1))`; the default CLI caps are far below that restriction.

## Evidence

All runs used one BLAS/OMP thread and `nice`, without an LP solve or a production copy. `sympy` was available through a temporary install outside the checkout (`PYTHONPATH=[local path]`). From the project root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice timeout 60 round2/.work/venv/bin/python -B round2/occupation-neon/embedding-flow/reciprocal/check.py
PYTHONPATH=[local path] OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice timeout 95 round2/.work/venv/bin/python -B round2/occupation-neon/embedding-flow/reciprocal/higher/check.py b3f18_lo.bin b3f18_hi.bin --authenticate
PYTHONPATH=[local path] OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice timeout 110 round2/.work/venv/bin/python -B round2/checks/embedding-chain-sol.py b3f18_lo.bin b3f18_hi.bin
```

- Supplied chain checker: `PASS 4896 chain balances, 464 nonzero, 133 distinct entrance tests`; `PASS 46720 reciprocal and persistent receive projections including degree 18 and out-of-census mothers`.
- Original reader replay: `rows=80, terms=1122, ports=176, extract=6, controls=1700, zeros=8, moment=8`; 19 nonzero branches; cubic pole numerator `−2025082539/31250000000` at `y0=−18`. This compares against the separate numeric `Cert.embed` path and authenticates both payloads.
- Independent literal-arrow audit: 15,092 balance rows (983 nonzero), 462,088 projection row values, 55 positive finite founder residuals, 27 multiple-child and 9 partially covered append incidences, 36 appends with surviving spare-child flux, 368 restored selected arrows, and 369 covered incidences with nonzero other old fields. It reconstructs inequalities from literal physical arrows, never using `check.py`'s synthetic planted lift or the chain projection's fields for the expected row values.
- Independent price-layer contractions: 144 exact signed objective images, 40 nonzero at each selected role, and both direct-call failures reproduced. This compares the assembled columns with direct physical evaluations of the authenticated reader; the separate original reader replay checks its coefficients against numeric D18 tables.
