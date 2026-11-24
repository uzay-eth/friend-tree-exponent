> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Pinned Neon D18 charged first-body slack — Astra integration

For the two original D18 payloads accepted in Astra-on-Neon-1, let `Φ` be the published normalized potential (each stored table entry divided by `S=10^11`). The sign is `+1` for `λ=565673/1000000` (lower) and `−1` for `λ=22631/40000` (upper). If `N_{L,k,h}` counts physical core vertices of state `(L leaves,k core children,h parent flag)`, then for either certificate, on **every finite physical tree** with a core,

```math
 \mathop{\mathrm{sign}}\nolimits (\mathcal L\Phi-\lambda\Phi)
 \ \ge\ \sum_{L\ge0,\ 0\le k\le4} a_k N_{L,k,1}.
```

The following `a_k` are exact *all-L* floors. `k=0` requires `L≥1`; `k≥1` includes `L=0`. The indicated worst legal `L` attains the computed first-body row minimum, not necessarily the global drift minimum.

| k | worst L | lower signed floor | upper signed floor |
|---:|---:|---:|---:|
| 0 | 1 | `256918177771749/100000000000000000000` | `20697884432099/10000000000000000000` |
| 1 | 0 | `4272535757583/4000000000000000000` | `861071198001/1000000000000000000` |
| 2 | 0 | `6713770119873/5000000000000000000` | `5412445357533/5000000000000000000` |
| 3 | 0 | `932799090941/625000000000000000` | `6016226911043/5000000000000000000` |
| 4 | 0 | `78795139786299/50000000000000000000` | `31764287097597/25000000000000000000` |

For integrating life histories, a few tighter per-state rows (minutes or other clocks do not enter the coefficients):

| `(L,k,h)` | lower | upper |
|---|---:|---:|
| `(1,0,1)`, newborn state | `256918177771749/100000000000000000000` | `20697884432099/10000000000000000000` |
| `(2,0,1)` | `508265145797027/150000000000000000000` | `16379780018491/6000000000000000000` |
| `(3,0,1)` | `36914485349857/10000000000000000000` | `2973168588449/1000000000000000000` |
| `(1,1,1)` | `25500867560131/10000000000000000000` | `123288571208369/60000000000000000000` |
| `(1,2,1)` | `12286814155743/5000000000000000000` | `990210891789/500000000000000000` |

Use a state-specific coefficient **in place of** its all-L band floor, or add only the nonnegative difference. Each vertex is counted once in the first-body sum. The edge charges below are separate summands of the pinned decomposition.

## Local charge and global validity

Write `d=L+k+h`, `B=(L−1,k+1,h)` when `L>0`, `b=(1,0,1)`, `Δ=f(L+1,k,h)−f(L,k,h)`, `ρ=Σ 1/D_donor` over core neighbors, and `R=h r_parent+k r_child` using the published raw integer centers divided by `Sr`. Write the formula with **unscaled** table entries, then divide by `S` at the end. The common first-body term is

```
A = L Δ + (L/d) [f(B)+f(b)+g(embed(B),b)−f(L,k,h)] − λ f(L,k,h).
```

The bracketed term is absent at `L=0`. The newborn *vertex and edge* are assigned here, once, to the parent that promotes a leaf; `d` does not change when it buds. Inner charged row is `A+R Δ`. Its missing donor credit is precisely `(ρ−R)Δ`, which the checked `g` rows receive split by the actual neighbor and role, with the exact reciprocal obtained from that neighbor's Bernstein embedding. For finite tail `L<LF`, `k≤KT`, the charged row is `A+ρ(Δ−ψ)+Rψ` with `ψ=c2/d²+c3/d³`. The missing `(ρ−R)ψ` goes into the T2 and T3 edge charges, with coefficients `c2/(3D0²)` and `c2/D0²+c3/D0³`. This is why neither dropping the donor term nor adding the edge recharge here is correct. For `k>KT` the transfer is zero. Exact affine extrema are at `ρ∈{0,(k+h)/2}`. For `L=LF` the table gain is zero and the capped transfer term uses the convex Bernstein controls of `ψ`; for `L≥LF+1` use the same valid transfer hull, the restricted newborn-edge Bernstein controls, and the two bud rate endpoints. This independently rederives the pinned GMP first-body formulas at `cert_b3f_verify_gmp.c:1004–1167`; the complementary edge source terms are at `1218–1255`. The accepted GMP replay proves the remaining signed charged rows are nonnegative. Therefore ignoring them is safe in the same sum of charged terms. The `¼ Σ(k−h)=0` recharge belongs only to the `c0` comparison, unchanged here.

## Consequence without an occupation estimate

For either sign let `a=a0`, `b=a1`. The exact values give `a>2b>0`. The affine floor `a−(a−b)k` is nonpositive for `k≥2`, so the signed first-body charge of **every nonroot** dominates that floor, including unscanned large `k` by the accepted GMP nonnegativity. Since `Σ_nonroot k≤C−1`, the global signed drift is at least `b(C−1)`. This alone shifts the endpoints by `b/C0`. With signed side inequalities and `Φ≤C0 C`, the effective rates are `η_−=λ_−+b_−/C0_−` and `η_+=λ_+−b_+/C0_+`. For lower, `Ψ=Φ−b_−/η_−` satisfies the lower drift at `η_−`; start iteration at a sufficiently large deterministic size with `EΨ>0`, which exists by the original positive lower product for `EΦ`. For upper use `Ψ=Φ+b_+/η_+`, which satisfies the upper drift at `η_+` and still dominates `c0 C`. The new exact enclosure from this floor alone is

```
[1388727487581492461/2454997964200000000, 369301187979920833/652735583940000000]
≈ [0.5656735801139580, 0.5657745602756464].
```

The inward shifts are about `5.80e−7` and `4.40e−7` (1% of the prior production width in total). These are conservative algebraic gains, not rounding effects: the state slack is order `10^−6`, versus the `10^−11` table storage quantum. No new occupation lower bound is claimed. A census/regenerative calculation can use the additional exact differences above the affine floor.

## Inner-edge extension and unconditional allowance

The directed physical core edge `p→w` is oriented by ancestry. For two inner endpoints (`2≤D_p,D_w≤17`) the exact `g` charge comprises:

- `−λg(p,w)` and both gain terms with physical **included** donor credit `1/D_w` or `1/D_p` and *centered external* donor credit; positive-degree external slots not included here transfer their displaced rate into the appropriate higher motif;
- both promotion terms, with `L/D` of the promoting endpoint, plus the new `hC` (when `p` promotes) or `hV` (when `w` promotes) assigned **to this g charge** exactly once;
- inward first-body transfers `(1/D_w−r_child(p)) Δf(p)` and `(1/D_p−r_parent(w)) Δf(w)`. Both are needed: omitting them yields 3,671 lower and 40,068 upper wrong-sign rows in the exact negative control.

This is the sum of the pinned `g` source and exact-inner transition lines 1211–1277. `check_edges.py` independently contracts all 43,776 legal ordered pairs per sign, includes root-parent states, asserts full range coverage (parent `k≥1`; child `h=1`), and leaves mixed/tail edges at the original GMP allowance of zero. Thus for the entire tree, with `E_II` the number of both-inner core edges, one can add `m E_II` to the first-body inequality, where

| scope | lower signed floor | upper signed floor |
|---|---:|---:|
| all inner→inner (`m`) | `25227936004759/1700000000000000000000` | `5046841434259/340000000000000000000` |
| nonroot inner→newborn state | `25354426150427/1700000000000000000000` | `745983027279/50000000000000000000` |
| nonroot degree ≤4 → newborn state | `374615525873/25000000000000000000` | `59975051779/4000000000000000000` |
| `(0,1,1)→(1,0,1)` | `374615525873/25000000000000000000` | `749900066669/50000000000000000000` |
| `(1,1,1)→(1,0,1)` | `562427823161/37500000000000000000` | `59975051779/4000000000000000000` |
| `(1,2,1)→(1,0,1)` | `1501397469751/100000000000000000000` | `150088975619/10000000000000000000` |

The listed edges are physically attainable: form a rooted core path through `p,w` with an extra core child when `k_p=2`; give each terminal core one leaf, then add leaves to reach the specified counts. Such a finite tree is attainable from `T3` by attaching physical leaves and later promoting them (each prescribed recipient has a positive attachment rate). Add only the difference above `m` if an edge already receives the common credit. The rates remain near `1.5e−8` per edge, well above the storage quantum but much smaller than first-body gains.

One can realize the common credit **globally**, without a lower bound on `E_II`, using unused first-body surplus. Write `A(k)=a−(a−b)k`. The exact sum over nonroots is `b(C−1)+(a−b)k_root`. For each nonroot tail vertex, the signed `f` charge minus `A(k)` covers `m(k+1)` (the number of its incident core edges): for `k≤4`, `check.py` verifies this using the exact tail-only minima and capped rows; for `k≥5`, nonnegativity gives surplus at least `(a−b)k−a≥m(k+1)`, checked at `k=5` and by the slope. At a tail root, the unused root term `(a−b)k_root` covers all its incidences. Each edge not counted in `E_II` has at least one tail endpoint, so these reimburse at least `m(C−1−E_II)`. No `f` charge or edge is double booked (a two-tail edge needs only one credit). Combining gives `(b+m)(C−1)` for either sign, and hence

```
[885313785947169446267/1565061202177500000000,
 376687206692677815401/665790295618800000000]
≈ [0.5656735881736863, 0.5657745526954197].
```

The additional inward shifts over `f` alone are `8.06e−9` and `7.58e−9`. For integration above this combined floor retain only excess edge credit above `m`, and subtract the charged tail-incidence budget from the `f` surplus (or explicitly assign each missing edge only once). The `hV/hC` charges themselves have not been spent.

## Physical bare-core corrector (no new cross-scan)

Let `s` be the sign, `B(T)` count nonroot bare cores (`L=0`), and `ε=7/10^7`. For the **new Astra potential** `Φ*=Φ−sε B`, the exact additional signed drift at the old `λ` is `ε(λB−ℒB)`. A nonroot becomes bare only on its last-leaf promotion at rate `1/D`, and loses bare status only when a core donor sends it a leaf, with total rate `ρ`. A newly promoted core has one leaf. Thus the extra signed row is `ελ+ερ` for `L=0`, `−ε/D` for `L=1`, and zero otherwise. Dropping `ερ≥0` is safe. No independent occupation hypothesis enters.

Set `a*=a−ε/2`, `b*=b+ελ`, and use the new affine floor `A*(k)=a*−(a*−b*)k`. The corrector check verified every nonroot inner descriptor (152 per sign) against `A*` and reused the exact tail minima printed by `check.py` for `k=0..4` to check tail reimbursement `m(k+1)`, even with a worst-case `ε/18` loss. For `k≥5`, accepted nonnegativity plus the affine slope covers this loss and `m(k+1)`. The root's unused `(a*−b*)k_root` still covers all its tail incidences. Consequently the signed drift of `Φ*` is at least `(b*+m)(C−1)`. For lower, `Φ*≤C0 C` and `Φ*≥(c0−ε)C` with positive `c0−ε`; for upper, the corresponding bounds are `Φ*≤(C0+ε)C` and `Φ*≥c0 C`. The same constant-offset argument yields

```
[885314122522604446267/1565061202177500000000,
 376687206692677815401/665790533618800000000]
≈ [0.5656738032294518, 0.5657743504481111].
```

This buys another `2.15e−7` and `2.02e−7` beyond the unmodified `f+g` bound. The gain comes from reweighting the physically attainable bare class, whose `k=1,L=0` row limited the old floor. The earlier `f+g` bound is still the one for the original potentials. When stacking refinements, use the corrected local rows and comparison constants, and subtract the new affine and tail budgets only once. The exploratory arithmetic check exercised the bare count drift on 76 small trees (59 entries into the bare class and 19 exits). The bare corrector is outside the Claim 21 replay.

## Portable replay

`pinned-fg-slabs.json.gz`, SHA-256 `63aab8b10d07b43abd80e5852e60b72f5513bc23618d83ce6b0fd3ba3d5ab321`, contains first-body raw table slices and inputs for the inner–inner edge calculations. Recheck on one core, no NumPy or full certificate files needed:

```
taskset -c 3 nice -n 15 /usr/bin/time -f 'wall=%e s RSS=%M KiB' python3 -B check.py
taskset -c 3 nice -n 15 python3 -B check_edges.py
```

`check.py` prints separately reproducible canonical JSON SHA-256 roots for each sign and slab (`json.dumps(...,sort_keys=True,separators=(',',':')).encode()`) and checks all slab ranges before the inequalities. `check_edges.py` prints the exact edge extrema and negative control. `charged-f.json` gives signed corners and operands for the 2,079 nonroot descriptors through degree 64 per sign (`extract.py --buffer`). The [release verifier](../release/README.md) re-extracts the first-body data directly from the original payloads.

To reproduce the byte-to-slab link **from the original files without copying them**, set `PYTHON`, `D18_LOWER`, and `D18_UPPER` as in the [release instructions](../release/README.md) before running from this directory:

```sh
set -o pipefail
OPENBLAS_NUM_THREADS=1 taskset -c 3 nice -n 15 "$PYTHON" -B extract.py --portable --edges --verify-sha "$D18_LOWER" "$D18_UPPER" | gzip -n > replayed-fg-slabs.json.gz
sha256sum replayed-fg-slabs.json.gz
```

The extraction streams each original once for SHA, checks exact header/data sizes and no change in inode/size/mtime through extraction, then reads selected mmap pages. Each original is 3,495,091,968 bytes.
