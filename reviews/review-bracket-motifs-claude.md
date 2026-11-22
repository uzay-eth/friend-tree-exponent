> Packaging note: one internal directory name on line 16 was redacted for public release; all mathematical content is unchanged.

VERDICT: CORRECT

# Cold referee audit: Strogatz 17/30 refutation, bracket [113/200, 1133/2000]

- **Referee:** Claude Opus 5.5 (`claude-opus-5-5`), Claude Code harness. Independent and cold.
- **Date:** 2026-11-22.
- **Object audited:** `attempts/strogatz-exact/PROOF.md` at commit `07919c1`.
- **Certificates:**
  - `build/upper_5665.bin`, sha256 `3a4aef888a5b567a8ada05f9546f98c979d361da9f9b33ba6e794eb91fb30abe`
  - `build/lower_5650.bin`, sha256 `9e966fba9aee6eae3da659b693d3bb9c481c1cea66635b5218d0887f182a60d9`

  I recomputed both hashes before starting and they match.
- **Files read:** only `PROOF.md`, `cert/cert_format.h`, the two `.bin` files, and `route-exact/exact_rationals_n21.txt`.
- **Files not read:** `cert_verify_gmp.c`, `verify_python_cross.py`, `cert_build.c`, `run_all.sh`, `RESULT.md`, `route-edge/*`, `checks/*` (other than this file, which I wrote) and the shared coordination board.
- **Repository changes:** none. All scratch code and logs are in `build/ref-claude/` (gitignored). Everything ran under `nice -n 10` on at most 3 cores.

## 0. Summary

| Item | Upper (`1133/2000`) | Lower (`113/200`) |
|---|---|---|
| Exact checks `f / g / hV / hC` | 484,138 / 160,545 / 13,519,800 / 3,164,625 | same |
| Total exact checks | **17,329,108**, 0 failures | **17,329,108**, 0 failures |
| `g(tail,tail)=0` (full 29×29 tail block, 841 entries) | all zero | all zero |
| Header, nesting and knot structure | OK | OK |
| Global constant (Lemma 5), over 241,900 states | `c0 = 2989767809/S`, at (L,k,h)=(0,1,1) | `C0 = 17527873930/S`, at (668,0,1) |
| `Phi(T3) = f(2,0,0)` | `14606753031/S` | `14608704616/S` |
| Prefactor, recomputed by me | **A_up = 14606753031/2989767809** (≈4.8856) | **A_lo = 7304352308/8763936965** (≈0.8335) |
| Minimum exact slack in S-units: `f` | 24776800469/250000 | 0 (equality is allowed) |
| Minimum exact slack in S-units: `g` | 2507610629/1312500 | 32318417927/16500000 |
| Minimum exact slack in S-units: `hV` | 2621009681/5750000 | 18244173547/40250000 |
| Minimum exact slack in S-units: `hC` | 5253400669/11500000 | 2515030267/5500000 |
| Exact tree-level tests (§3) | 1,698 trees, 0 failures | 1,698 trees, 0 failures |

- The prefactors match PROOF.md Theorem 1 digit for digit.
- The check counts match PROOF.md §6 exactly. I derived the counts myself and did not take them from PROOF.md.
- Theorem 1's finite-n bounds hold against exact `E[C_n]` for every n = 3..21 (§4).
- Every planted false control was rejected (§5).
- I found no fatal flaw. There are three minor notational imprecisions (§6).

## 1. Analytic audit of PROOF.md (adversarial)

### §2 Canonical rooting
- **T_3 and the root.** T_2 is a single edge, and every move attaches the new leaf to an endpoint. So T_3 is always K_{1,2} and its centre r is unique.
- **The core is a connected subtree containing r.**
  - Degrees never decrease, so r stays in the core.
  - In a tree with n ≥ 3, the internal vertices form a connected subtree: every interior vertex of a path is internal.
  - Every leaf's neighbour has degree ≥ 2 once n ≥ 3, so every leaf hangs off a unique core vertex.
  - With h, L and k defined as in PROOF.md, D_u = L_u + k_u + h_u holds.
- **Generator.**
  - The probability of attaching to w is (1/n) Σ_{v~w} 1/D_v.
  - For a core vertex u this gives (L_u + ρ_u)/n, because each leaf neighbour contributes 1.
  - For a leaf child of u it gives 1/(n D_u).
  - Budding sends u to (L−1, k+1, h) with D unchanged, and creates a newborn vertex (1,0,1).
  - The formula for (𝓛Φ)(T) is correct. I confirmed it exactly against a fully explicit enumeration over every pair (V, W), weighted 1/deg V, on 27 trees per certificate (§3).
- **The bound ρ_u ≤ (k+h)/2** holds because every core D_v ≥ 2.

### Forward invariance of H_u and C_u, and creation events
- Under leaf gain, D goes up by 1 and k is fixed. Under budding, k goes up by 1 and D is fixed. h never changes. So a vertex that leaves H_u or C_u never comes back.
- Entry into H_u happens only when (L,0,1) with D < D_hu buds. Entry into C_u happens only when (L,1,h) with D < D_cu buds. Leaf gain can never cause entry.
- **Wording point.** A newborn path (p→u′→nb) is created at *every* budding where u′ ∈ H_u, not only on entry. Likewise, newborn cherries {w, nb} are created at every budding where u′ ∈ C_u. The formulas in §4.5 handle the general case correctly: they use the indicators `1_{w^bud∈H_u}` and `1_{p^bud∈C_u}`. Only the remark in Lemma 4.3 ("k_u=0→1 and k_u=1→2") reads more narrowly; see §6.
- No other structure is created: nb has no children.
- I attacked each boundary with dedicated trees, and all passed (§3):
  - u = (L,0,1) budding into H_u;
  - u = (L,1,h) budding into C_u;
  - leaving H_u at D = D_hu − 1 or at k = K_hu;
  - leaving C_u at k = K_cu;
  - D_0 − 1 → D_0;
  - k = K_T → K_T + 1;
  - L = L_F and L > L_F;
  - k = K_F and k > K_F.

### Lemma 1: PL Fenchel–Young
- The proof is a pure convex-combination identity on the containing subinterval. No convexity of A is needed.
- A is exactly affine on [0, 1/D_0] because ξ_1 = 10/240 = 1/24 = 1/D_0, which I checked against the header.
- max_j is convex and min_j is concave. So A*(affine in t) is bounded by interpolating between its corner values, which is the right direction for both UPPER and LOWER.
- The coefficients k_w, k_p − 1 and 1 in front of A* are all non-negative, which the argument needs.

### Lemma 2: 1/D tail
- **Leaf gain.** Δt = −t/(D+1). With F(μ_p) affine in t, this gives (L+ρ)ΔF = t·γ·(F[INF] − F[T1]).
- **Upper bound on γ.** L + ρ ≤ D − (k+h)/2 < D + 1, so γ < 1. This holds for every finite tree.
- **Lower bound on γ.** γ ≥ L/(D+1) = 1 − (k+h+1)/(D+1).
  - For κ ≤ K_T we have k = κ exactly, so this is non-decreasing in D. Its minimum is at D = D_0, which gives γ_min = (D_0 − κ − h)/(D_0 + 1).
  - For κ = K_T + 1, γ_min = 0 is trivially valid.
- **Budding.** b = 1 − (k+h)/D ≥ (D_0 − κ − h)/D_0 for κ ≤ K_T, and b < 1.
  - For κ ≤ K_T at a tail vertex, L ≥ 11 > 0, so b_min > 0 is consistent.
  - For κ = K_T + 1, budding leaves κ unchanged, so ΔF = 0.
- Every rate enters linearly, so checking the two endpoints of each rate interval suffices.

### Lemma 3: quadratic and bi-affine credit
- The bounds t·((1−t)d_0 + t·d_1) ≤ t·max(d_0,d_1) and the bi-affine split are correct.
- I checked the corner assignment against §4.5:
  - corner (T1, INF) receives max(d00, d10);
  - corner (T1, T1) receives max(d01, d11) from p and max(d10, d11) from w;
  - corner (INF, T1) receives max(d00, d01).
- At the INF corners the actual tail credit is 0, so no term is lost.
- **Load-bearing check:** removing these terms makes both certificates fail (control 12/13 in §5).

### Lemma 4 and §4.5: term-by-term bookkeeping
I enumerated every term of (𝓛Φ)(T) myself.

**Leaf gain at an inner u:**
- Δf: split by r_par and r_ch, with the residuals going to the edges.
- Δg(p→u): L_u + 1/D_p stay on the edge. Each child term 1/D_x goes to x_ch plus an h_V source if u ∈ H_u. Otherwise it is split by Fenchel–Young with A_cp: k_u·A*_cp goes to (p→u), and A_cp(1/D_x) goes to (u→x).
- Δg(u→w): L_u + 1/D_w stay on the edge. The parent credit goes to x_par plus an h_V source, or else to Fenchel–Young with A_pc: k_u·A_pc goes to (p→u). The sibling credits go to x_cc plus a cherry source, or else to Fenchel–Young with A_cc: (k_u − 1)·A_cc(1/D_w) plus (k_u − 1)·A*_cc go to (u→w).
- h_V with u as top, middle or bottom, and h_C with u as apex or child: full rate, bounded by the stated intervals.

**Leaf gain at a tail u:** Lemma 2 covers the full rate on every structure incident to u.

**Budding:**
- f(u′) − f(u) + f(nb) + g(u′→nb) goes to the vertex.
- The new state's changes on every incident structure are counted.
- Newborn h_V(p→u′→nb) goes to the edge (p→u).
- Newborn cherries go to the edges (u→w).

Everything matches §4.5. Signs and factors (k_u, k_u − 1, h_u, the (h+k−1)/2-type ranges, L/D) all agree.

**Exact decomposition test.** I tested the decomposition exactly rather than only by reading:
- For every tree in §3, I computed each piece's actually assigned value. These are: vertex, edge, vertical 2-path and cherry pieces, at actual rates and actual t; each Fenchel–Young product xy is split as A(x) + (xy − A(x)).
- Σ(pieces) equals (𝓛Φ)(T) **exactly** on all 3,396 trees (both certificates).
- Every single piece satisfies piece_LHS ≤ λ·piece (UPPER), or ≥ λ·piece (LOWER).
- This identity test catches deliberately planted bookkeeping bugs (§5).

### Lemma 5: global comparison
- Each non-root core w receives exactly one g term (its incoming edge) and at most one h_V term (when its parent is in H_u). Each u ∈ C_u receives C(k_u, 2) cherries. Every g, h_V and h_C term is assigned exactly once.
- With ext(0, ·), any h_V or h_C entry that is 0 or absent is covered.
- For tail w, the term is a convex combination of the INF and T1 corners, so the bound ext(g_ext(INF), g_ext(T1(κ,1))) is valid.
- For L > L_F, f clamps to L_F, and g_eff depends only on κ. So the finite box [0..800]×[0..150]×{0,1} plus VSTAR_0 and VSTAR_1 covers every D, L and k.
- I recomputed c_0 and C_0 myself; both are > 0. Φ(T_3) > 0 in both certificates.
- Lemma 5 (Φ ≥ c_0·C or Φ ≤ C_0·C) also holds on every test tree.

### Induction, Theorem 1 and Corollary 2
- E[Φ(T_{n+1})] = E[Φ(T_n)] + E[(𝓛Φ)(T_n)]/n. The pointwise drift bound plus Lemma 5 then gives the stated bounds, with A = Φ(T_3)/c_0 (or /C_0).
- The Gamma-function asymptotic ∏_{j=3}^{n−1}(1 + α/j) = 2Γ(n+α)/(Γ(3+α)Γ(n)) is correct.
- **Corollary 2** is valid:
  - Markov's inequality along n_k = 2^k gives a summable sequence, and Borel–Cantelli applies.
  - C_n is non-decreasing, so C_n ≤ C_{2n_k} ≤ (2n)^{a+ε}.
  - Intersecting over rational ε finishes the argument.
- The refutation of E[C_n] = n^{17/30+o(1)} follows from Theorem 1 alone. The margin is 17/30 − 1133/2000 = 1/6000.
- The Remark identifying μ ∈ [0.565, 0.5665] depends on an outside existence theorem by the "Astra lead", which I did not audit. It is not needed for Theorem 1 or Corollary 2.

## 2. Independent exact verifier: `build/ref-claude/refv.c` (GMP mpq, no floating point)

- **Written from scratch.** I wrote it from PROOF.md §4.5 and §5 only. I implemented the state ordering, `pars`, `chs`, `hu_list`, `cu_list` and `vidx` myself rather than including `cert_format.h`.
- **Arithmetic.** Every inequality is homogeneous in the table entries, so I work in S-units. Centering constants are R/S_r and knot abscissae are knot_num/240. Every comparison is an exact `mpq_cmp`.
- **What it implements:**
  - every vertex case: VSTAR_h with b ∈ {0,1} and both g_nb; exact L = L_F; tail L ≥ L_F + 1 with both b and both g_nb; inner box; outer box with r ∈ {L, L+(h+k)/2};
  - the edge, h_V and h_C inequalities, with all A_tot, Fenchel–Young conjugate, γ/b endpoint and Lemma 3 quadratic terms;
  - the tail×tail zero test on all 841 tail pairs;
  - header and nesting checks: 2 ≤ D_cu ≤ D_0 ≤ L_F; D_hu ≤ D_0; K_T + 1 ≤ K_F; K_hu ≥ 1; K_cu ≥ 2; K_F + 1 ≥ D_0;
  - knot checks: 0 = ξ_0 < ξ_1 = 1/D_0 < … < ξ_8 = 1/2, strictly increasing;
  - payload length, with exactly 0 trailing bytes;
  - the Lemma 5 global comparison.
- **Overflow guard.** All table values have |x| < 1.7·10^10 S-units, so no int64 overflow occurs before conversion to mpq.

Output, condensed from `build/ref-claude/full_{upper,lower}.log`:
```
upper_5665: f checks=484138 fails=0 | g 160545/0 | hV 13519800/0 | hC 3164625/0
            TOTAL checks=17329108 fails=0 header_ok=1 knots_ok=1 global_ok=1  RESULT PASS exit=0
            GLOBAL c0 = 2989767809/S (L=0,k=1,h=1); A_up = 14606753031/2989767809
lower_5650: f 484138/0 | g 160545/0 | hV 13519800/0 | hC 3164625/0
            TOTAL checks=17329108 fails=0 ... RESULT PASS exit=0
            GLOBAL C0 = 17527873930/S (L=668,k=0,h=1); A_lo = 7304352308/8763936965
```

## 3. Exact tree-level tests from the model definition (`trees.py`, Fraction only)

**What is computed:**
- Φ(T) directly from §3: f(π_f), multilinear μ-evaluation of g, h_V and h_C.
- (𝓛Φ)(T) = Σ_u [(L_u + ρ_u)·(Φ(T + leaf at u) − Φ) + (L_u/D_u)·(Φ(bud at u) − Φ)], with Φ recomputed from scratch on every successor tree.
- For trees with n ≤ 45, additionally the explicit generator Σ_V Σ_{W~V} (1/deg V)·ΔΦ over the full graph, re-rooted at r. It agrees exactly with the grouped version on 27 of 27 trees per certificate.
- Per tree, four checks:
  1. drift: (𝓛Φ) ≤ λΦ for UPPER, or ≥ λΦ for LOWER;
  2. the Lemma 4 identity;
  3. per-piece domination;
  4. Lemma 5.

**Tree sets (each run on both certificates):**
- Small trees: T_3, K_{1,3}, P_4, and P_5 rooted at its centre and at an end.
- 12 random friend trees with n = 12..45 (explicit generator).
- 49 handcrafted extreme trees, covering:
  - D ∈ {D_0−2, D_0−1, D_0, D_0+1}, at the root and at a child;
  - tail→tail edges in both orientations;
  - tail→H_u→tail paths with k_u = 1 and 2;
  - H_u with k = K_hu = 8 and D_u = 15;
  - C_u apexes, h ∈ {0,1}, with two tail children and with mixed children;
  - entry into H_u under a tail parent and under an inner parent;
  - entry into C_u for a child and for the root;
  - L = 0 chains;
  - k ∈ {11, 12, 13, 14} with D ∈ {23, 24, 25}, rooted and non-rooted;
  - L ∈ {799, 800, 801, 803};
  - k ∈ {150, 151, 152} (VSTAR_0 and VSTAR_1);
  - root (L_F, K_F);
  - root (L_F+2, K_F+1);
  - the κ clamp at D_0−1 → D_0 with k = 22, 23;
  - C_u at k = K_cu exiting by budding.
- 1,560 random adversarial trees, with up to 30 core vertices and n in the thousands. Degrees were drawn from {2, 3, 11, 12, 13, 15, 16, 17, 22..26, 30, 48, 100, k+h+{0,1,2}, k+h+{799, 800, 801, 807}}.
- 46 random friend trees with n up to 500.

**Result:** 1,698 trees per certificate. For every tree: drift OK, identity exact, 0 per-piece failures, and Lemma 5 OK. That is 0 failures on 3,396 tree evaluations.

## 4. Theorem 1 against exact E[C_n]
- I wrote my own exact enumeration (`ecn.py`) over rooted shapes. It reproduces column 3 of `exact_rationals_n21.txt` exactly for n = 3..17, with 399,466 states at n = 17.
- Using the file's values for n = 3..21, A_lo·∏(1 + 113/(200j)) ≤ E[C_n] ≤ A_up·∏(1 + 1133/(2000j)) holds for **all n = 3..21**. For example, at n = 21: E/upper ≈ 0.2548 and E/lower ≈ 1.498. At n = 3: E[C_3] = 1, against upper 4.8856 and lower 0.8335.

## 5. Planted false controls (all rejected)

The exact verifier `refv` exits with code 2 and REJECT on each of these:

| # | Control | Outcome |
|---|---|---|
| 1 | upper with λ = 1132/2000 | first fail VSTAR, slack −548580466653/75500 |
| 2 | lower with λ = 114/200 | first fail inner f(2,0,0) |
| 3 | upper with λ = 17/30 | g fails |
| 4–5 | f[(2,0,0)] ∓ 10^10 (upper and lower) | 1 fail each |
| 6 | g[3, 0] − 10^10 | 1 fail |
| 7 | planted g[INF, T1] = 1 (tail-tail nonzero) | structural reject |
| 8–9 | hV[5000] ∓ 10^9 (upper and lower) | 1 fail each |
| 10–11 | hC[777] ∓ 10^9 (upper and lower) | 1 fail each |
| 12–13 | broken Lemma 3: tail-credit quadratic terms omitted (upper and lower) | many g failures |

The `nofy` and `outer_rate` variants only *loosen* the verifier, by dropping or narrowing upper-side terms. They pass, as a loosening should. I list them for transparency; they are not counted as controls.

Tree-level controls, run on both certificates, all rejected:
- f(2,0,0) ∓ 1.0 on T_3;
- g[(1,1,0),(1,0,1)] on P_4;
- hV[T1(1,0),(3,1,1),T1(0,1)];
- hC[(0,2,0); T1(0,1), T1(0,1)];
- λ = 1132/2000 or 114/200 on an RFT with n = 200. This caught 35 or 36 piece failures, even where the whole-tree drift still passed.

Decomposition mutations, all detected by the exact identity test:
- dropping newborn h_V;
- dropping newborn h_C;
- dropping the h_V x_par source;
- dropping the A_pc Fenchel–Young residual. This one needs a tree with an inner, non-root vertex outside H_u, so I added two targeted trees; it is detected on both certificates.

## 6. Minor non-fatal imprecisions
1. **§4.5 definition of v⁺.** It writes v⁺ = T1(k_v, h_v) when D_v = D_0 − 1. It should be T1(κ(k_v), h_v) with κ = min(k, K_T+1). Inner states with D = 23 have k up to 23 > K_T + 1 = 13. My verifier uses the clamp, and the tree tests at k = 14, 22, 23 confirm that this reading is correct and passes.
2. **Lemma 4.3, last bullet.** The phrase "(k_u = 0→1 and k_u = 1→2)" suggests newborn h_V and h_C appear only on entry. In fact a newborn path or cherry is created at every budding with u′ ∈ H_u or C_u. The inequalities in §4.5 already include the general case.
3. **Small wording points:**
   - Lemma 2 says γ_p ∈ [γ_min, 1]; in fact γ_p < 1 strictly. This is harmless.
   - The VSTAR newborn bound uses D = K_F + 1 + h, but a VSTAR state that can bud has D ≥ K_F + 2 + h. This is conservative and valid.

I did not verify the historical attribution to Strogatz, the claim of "sharpest known" bounds, or the outside existence theorem cited in the Remark. None of them affects the mathematics of Theorem 1 or Corollary 2.

## 7. Reproduce

All commands run inside `build/ref-claude/`:
```
gcc -O2 -o refv refv.c -lgmp
nice -n 10 ./refv ../upper_5665.bin ; nice -n 10 ./refv ../lower_5650.bin          # full replay
nice -n 10 ./refv ../upper_5665.bin --lam 1132 2000                              # control (exit 2)
python3 -I run_trees.py ../upper_5665.bin all 60 ; python3 -I tree_more.py ../upper_5665.bin 1500
python3 -I ecn.py 17 > ecn.log ; python3 -I bounds.py
```
