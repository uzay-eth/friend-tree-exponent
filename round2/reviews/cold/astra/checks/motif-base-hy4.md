CORRECT

# Independent referee review: finite all-degree four-vertex raw drift theorem

Scope. I reviewed the compact base certificate and its two proofs:
`proof.md` (`motif-cert/proof.md`), `higher-rows.md` (`motif-high/higher-rows.md`), and the frozen/verifier sources
`snapshot.hpp` (`motif-cert/snapshot.hpp`), `export.c` (`motif-cert/export.c`),
`low.cpp` (`motif-cert/low.cpp`), `envelopes.cpp` (`motif-cert/envelopes.cpp`),
`tail.py` (`motif-cert/tail.py`), `comparison.cpp` (`motif-cert/comparison.cpp`),
`low_oracle.py` (`motif-cert/low_oracle.py`), `tree_check.py` (`motif-cert/tree_check.py`),
`mutations.py` (`motif-cert/mutations.py`), `replay.py` (`motif-cert/replay.py`),
`oracle.py` (`motif-high/oracle.py`), `physical_check.py` (`motif-high/physical_check.py`),
`check.cpp` (`motif-high/check.cpp`), `negative.py` (`motif-high/negative.py`) and
`certificate.json` (`motif-cert/certificate.json`). The raw tables
`mc-upper.bin` (`outputs/mc-upper.bin`) and `mc-lower.bin` (`outputs/mc-lower.bin`) were opened read-only; their digests match the certificate. I changed no target file and made no commit.

My own check is `motif_base_hy4.py`. It is a separate implementation: I derived the
analytic-tail algebra, the envelope inequality, the generator/FY/newborn bookings and the comparison arithmetic
from the model in the proofs, wrote my own small-tree generator and my own power→tensor-Bernstein conversion, and
queried the author's compiled binaries only as a read-only table oracle (exact GMP row values, frozen weights).

## 1. Analytic f tail (independent re-derivation)

From `w=a+bq+cq^2+eqx`, `q=L/D`, `x=1/D`, the gain successor `(q+x)/(1+x), x/(1+x)` and the promotion successor
`(q-x,x)`, I recomputed `(1+x)^2 Δ⁺w/x = b(1-q)(1+x)+c(1-q)(2q+(1+q)x)+e[(1-2q)x-qx^2] = N` and
`Δ^bw = -bx+c(-2qx+x^2)-ex^2 = U`. Hence the signed residual is `σ{ρN+(1+x)^2[q(U+C+dx)-λw]}` with
`ρ∈{q,(1+q)/2}`, `C=f(J)+g(∞,J)`, `d=d_2[g(T(min(k+1,K_2+1)),J)-g(∞,J)]`. This is exactly proof.md:52-65
and tail.py:100-105. I checked the identity against a direct evaluation of the residual at `D∈{D_F,D_F+1,3D_F+17}`
and `L∈{0,1,⌊D/3⌋,D-1}` for both `ρ` endpoints and both signs — exact equality.

I then bounded `p` by my own tensor-Bernstein conversion over 16 `q`-panels and `x∈[0,1/D_F]`, applying the
division by `(1+1/D_F)^2` only when the maximum is negative (correct because the residual is `p/(1+x)^2`, whose
maximum for negative `p` occurs at `x=x_max`). Finally I semantically sampled the direct residual on a rational grid
including `q=0` and small `x` and required it to be ≤ the cap. Recomputed caps:

```
upper  cap=-2532   (cap_rational -79531159920317/31404816000, sampled worst -1.47e7 <= cap)
lower  cap=95165   (cap_rational 951645367447/10000000,        sampled worst 0     <= cap)
```

These match the frozen `upper-analytic.json` (`outputs/mc-package-verification/upper-analytic.json`) /
`lower-analytic.json` (`outputs/mc-package-verification/lower-analytic.json`). The seam argument at `D_F-1` and the
`D|Δ⁺f|=O(1)` conclusion follow from the displayed `N`; no unproved floating input survives.

## 2. Affine envelope

For a fine label with projection weights `(n/den, (den-n)/den)` the required enclosure is
`|den·δ − n·m₁ − (den−n)·m₂| ≤ n·r₁ + (den−n)·r₂`, which is exactly the convexity of the expectation over the
product interpolation: if each coarse corner lies in `[m−r,m+r]`, so does the physical mixture. For two coordinates
this becomes the four-corner weighted check (envelopes.cpp:68-75), with one-tail radii established before
two-tail cells. I verified the projection masses sum to one and sampled 960 one-dimensional enclosures per side
directly from the stored weights (recomputing `δ` from `weight("E", …)` rather than from the builder's difference
function). All passed. The full 1D=2,249,983 and 2D=7,418,600 inequalities are author-run production sweeps;
I did not rerun them.

## 3. FY transfers and newborn bookings

I rebuilt the generator identity on my own trees. For each physical event I computed the direct drift
`Σ_x r_x(Φ⁺−Φ)+b_x(Φᵇ−Φ)` with `Φ` from the frozen tables, then decomposed it into per-row increments, moving each
donor contribution `(i_y−z)Δ` between adjacent rows exactly once (centering `r_par,r_ch,x_par,x_ch,x_cc` for exact
states, signed FY knots `A_pc,A_cp,A_cc` otherwise) and booking each newborn edge/triple/cherry on its promoter's
post-promotion state. I asserted:
- the transfer identity `Σ(rows) = σ(direct−λΦ)` holds on every tree (no drift lost or double-booked),
- each FY donor is an active edge,
- each credited row is ≤ its exact GMP corner bound.

Both sides passed on 53 trees / 416 events / 208 newborns / 208 higher rows / 246 corner queries, including the
degree seams `D₂=72`, `D₃=28`, `D₄=14`, clamped high child counts, three simultaneous tails and degree 117.
The FY array ordering (proof.md:17: `A_pc,A_cp,A_cc`; low.cpp `fy[1]` at the top with count 1 = `A_cp[s_p](ι_w)`,
`fy[0]` at the bottom with count `k_w` = `A_pc[s_w](ι_p)`, `fy[2]` with count `k_p−1` = `A_cc`; and the conjugate
counts `1,k,k−1`) is internally consistent across low.cpp, low_oracle.py, tree_check.py and `bracket-motifs.md`:212.
I initially suspected a label swap and confirmed there is none.

## 4. Multiplicity over all D and signs

The finite comparison covers every `(D,k)` with `2≤D<D_F`, i.e. `Σ_{D=2}^{1399}D = 979,299 =
f_inner+f_table` (certificate.json:11), which I recomputed. The binomial multiplicities
`C(k,2)` cherries/Y, `C(k,3)` tridents and `k_s−1` forks were checked for all `k<D_F`; the fork count uses the
active apex with `k_s≥2`. Both signs are the two witnesses: upper `σ=+1`, lower `σ=−1`; the high-degree singleton
uses the worst endpoint of the neighbour input `L+(k+1)/2`, the lower tail uses
`q(B₁−B₂q−C/D) ≥ q(B₁−B₂−⌈C/D_F⌉/D_F)` with the checked positivity `coeff≥0` and
`q_min=max(0,(D_F−k−1)/D_F)` (0 for the overflow class), and `Σ(k−1)=−1`. No degree or sign is exempted.

## 5. Retained disjoint Q4 and Q1

`Q_4` is the sum of `J(σΔ_xh₄,i_y)` over four-core external donors at an exact contact; `Q_1` is the sum of
`J(σΔf,i_y)` over high-degree singleton neighbours. They come from different additive summands of the generator
(the drift of `h_P,h_F,h_Y,h_T` versus the drift of `f`), and in my implementation occupy disjoint keys of one
nonnegative reserve: a key indexed by a four-set support for Q4 and by a singleton core for Q1. I confirmed each
reserve term is ≥0 (`J(a,p)=max(a,0)/2−pa≥0` for `p≤1/2`) and that the whole-tree inequality
`actual + Q_4 + Q_1 ≤ e·C_core` holds simultaneously. The planted double-credit (doubling the Q4 reserve) is
rejected by the corner bounds; the "duplicate claiming" control in `negative.py` (`motif-high/negative.py`) and
the `double-gap` plant in `tree_check.py` (`motif-cert/tree_check.py`) cover the author's side.

## 6. Comparison and bracket

I recomputed the per-core cap sum and the outward bracket with exact `Fraction` arithmetic from certificate.json:

```
upper  λ=56581165/10^8   err=155608   min_weight=258643133503   base=0.5658122516320553
lower  λ=5656453/10^7    err=115075   min_weight=314288655137   base=0.5656449338556988
```

These are the exact endpoints behind the outward bracket `[0.5656449338, 0.5658122517]` in
`README.md` (`motif-cert/README.md`); the displayed decimals are rounded outward, as required for an enclosure.
`comp["lower"]` is the positive lower bound on `SΦ/C` used as denominator in both directions, so the signed
negative upper analytic cap (`−2532`) still helps rather than hurts.

## 7. My independent tests, commands and planted errors

From the repository root, one CPU:

```
taskset -c 0 python3 -B round2/checks/motif_base_hy4.py upper
taskset -c 0 python3 -B round2/checks/motif_base_hy4.py lower
```

Observed (both sides), ending `INDEPENDENT REFEREE TEST: PASS`:

```
analytic tail identity + sampling PASS cap=-2532 / 95165
affine envelope sampling PASS n=960
multiplicity over all D and signs PASS finite_rows=979299
own small-tree generator/FY/newborn/credit PASS {'trees': 53, 'events': 416, 'newborns': 208} rows=208 corners=246
planted FY/birth/credit mutations rejected
comparison arithmetic PASS base=0.565812251632 / 0.565644933856
```

Planted errors that my harness demonstrably rejects: (i) dropping the child-to-other-child FY transfer breaks the
transfer identity; (ii) dropping a `F` newborn breaks the transfer identity; (iii) doubling the Q4 credit makes a
T-row reward exceed its GMP corner bound. Each is a test of the test, showing the checked steps are load-bearing.

## 8. What I did not recompute

I did not rerun the 74,226,396 finite rows, the 2,249,983+7,418,600 projection inequalities, the 66,763,695
four-core corners, the high-degree scans or the full `replay.py`/UBSan replay. Those were pinned and run by the
author; I relied on the frozen digests, the recorded outputs under `outputs/` and certificate.json
for the aggregate row counts and per-family caps, and verified all downstream arithmetic exactly. The founder and
occupation lemmas and the separate reciprocal-degree tail are outside this review. I found no gap in the material
under review.

OPEN: none.
