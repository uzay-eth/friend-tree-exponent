UPPER VERDICT: HOLDS
LOWER VERDICT: HOLDS

# Adversarial referee report on Neon's "strogatz-exact" claim (Opus 5.5, attack-neon1-opus)

Claim checked: for every n >= 3,
`(7304352308/8763936965) prod_{j=3}^{n-1}(1+113/(200j)) <= E[C_n] <= (14606753031/2989767809) prod_{j=3}^{n-1}(1+1133/(2000j))`,
from `[neon-branch]/attempts/strogatz-exact/PROOF.md` and the certificates
`checks/attack-neon/{upper_5665,lower_5650}.bin`.
The copies of `cert_build.c, cert_format.h, cert_verify_gmp.c, run_all.sh` in `attack-neon/` are byte-identical to the ones
in `[neon-branch]/.../cert/` (`cmp` printed `same` for all four).
Certificate hashes (sha256): upper `3a4aef88...0abe`, lower `9e966fba...60d9`.

**Bottom line.** I could not break either bound. I re-derived Lemma 4 from the definition of Phi and the two events, and
every term of §4.5 matches my derivation: no missing term, no wrong sign, coefficient, endpoint or multiplicity. The
verifier checks exactly the §4.5 list, in exact arithmetic, for every table entry. An independent verifier I wrote from
the text alone agrees: 0 failures, same counts, same prefactors. On 34 776 trees, Lemma 4 holds as an exact
rational identity (Σ_τ E[τ] = (LΦ)(T)). Every per-term local inequality and both drift and floor inequalities hold on
all of them. I found only nits (§5). None affects correctness, and the certificates as checked prove both bounds.

---

## 1. Logic audit (task 1)

### 1.1 Exact drift
For n >= 3 every leaf has exactly one (core) neighbour, so the new leaf goes to core u with probability
(L_u + ρ_u)/n and to a given leaf child of u with probability 1/(n D_u). This gives
(LΦ)(T) = Σ_u [(L_u+ρ_u) Δ_{L,u}Φ + (L_u/D_u) Δ_{bud,u}Φ]. A vertex's state (L,k,h) changes only under leaf gain at it
(L+1) or budding at it (L-1,k+1, D unchanged), so each term of Φ changes only under events at its own 1, 2 or 3
vertices, plus newborn terms. Correct.

### 1.2 My derivation of Lemma 4, compared with §4.5
For each term τ and each vertex x of τ, the rate L_x + ρ_x splits as follows: L_x, plus credit from neighbours inside
τ (kept exactly), plus credit from neighbours outside τ (transferred, Fenchel–Young-decoupled, or bounded by a range).

* **f(u), u inner**: ρ_u Δ_L f = (h r_par + k r_ch)Δf + h(1/D_p − r_par)Δf [to edge p→u] + Σ_w (1/D_w − r_ch)Δf [to each edge u→w].
  The pieces sum back to ρ_u Δf for every k and h. At the root (h = 0) the r_par piece is absent; the verifier
  multiplies r_par by h. **Tail u**: no transfer, range ρ ∈ [0,(h+k)/2]. **L ≥ L_F**: Δ_L f = 0 through the clamp π_f, and
  b ∈ [(L_F+1)/(L_F+1+k+h), 1] with t ∈ [0, D0/(L_F+1+k+h)] is a valid rectangle. **VSTAR**: Δ_L f = Δ_bud f = 0, and
  only f(nb) + g(u'→nb) remains, with b ∈ [0,1]. **Newborn**: f(nb) + g(u'→nb), evaluated at the exact t of the
  budding vertex. All match §4.5 item 1.
* **g(p→w), p inner**: ρ_p = 1/D_w [internal] + h_p/D_{pp} [x_par split into hV(pp→p→w) if p ∈ H_u, else FY with A_pc]
  + Σ_{c≠w} 1/D_c [x_cc split into the cherries {w,c} if p ∈ C_u, else FY with A_cc, multiplicity k_p − 1].
  **w inner**: ρ_w = 1/D_p [internal] + Σ_{c∈ch(w)} 1/D_c [x_ch split into hV(p→w→c) if w ∈ H_u, else FY with A_cp,
  multiplicity k_w].
  The FY halves go to the right edges: the A*_pc[p], (k_p−1)A*_cc[p] and k_w A*_cp[w] halves stay on (p→w); edge (p→w)
  receives A_cp[p](1/D_w), (k_p−1)A_cc[p](1/D_w) and k_w A_pc[w](1/D_p) from its neighbouring edges. Every product
  (1/D)·Δg is decoupled exactly once. This is exactly A_tot plus the three conjugate terms of §4.5 item 2, with the
  same conditions (p ∉ H_u, h_p = 1 / p ∉ C_u, k_p ≥ 2 / w ∉ H_u, k_w ≥ 1).
  **Tail p or tail w**: Lemma 2 takes the whole rate, (L+ρ)Δ_L F(μ) = t·γ·(F[INF] − F[T1]) with γ ∈ [γ_min, 1). That
  identity and the γ_min, b_min bounds are correct: L/(D+1) is increasing in D for fixed k ≤ K_T. Budding at a tail
  vertex moves κ only when k ≤ K_T.
  **Newborns at budding**: the hC(p'; w, nb) newborn goes to edge (p→w), and the hV(p→w'→nb) newborn goes to edge
  (p→w), both with indicators p' ∈ C_u and w' ∈ H_u. That is everything a budding event creates: nb has k = 0.
* **Quadratic terms**: in each term, the only products of a tail parameter with a quantity depending on that same
  parameter are:
  * (1/D_p)Δ_w g and (1/D_w)Δ_p g in g;
  * (1/D_p)Δ_u hV and (1/D_w)Δ_u hV in hV;
  * (1/D_{w_i})Δ_u hC in hC.

  Each is bounded by Lemma 3 at the T1 corner of that parameter. The verifier uses (d0 at INF, d1 at T1) evaluated at
  the other parameter's corner, which is the right bilinear corner bound. The source terms (1/D_p)·Δ_u g(u→w) etc.
  are bilinear and exact, and are evaluated with inv(T1) = 1/D0, inv(INF) = 0. A*(affine in t) is convex (UPPER: max)
  or concave (LOWER: min), so the corner interpolation bounds it in the right direction. A(t/D0) is exactly affine
  because ξ_1 = 1/D0 (Lemma 1(1)). After these bounds every LHS is multi-affine in the tail parameters, with
  non-negative weights μ, so checking the corners suffices.
* **hV(p→u→w), u ∈ H_u**:
  * rates: p in [L_p+1/D_u, L_p+1/D_u+(h_p+k_p−1)/2], w in [L_w+1/D_u, L_w+1/D_u+k_w/2], u in
    [L_u+1/D_p+1/D_w, …+(k_u−1)/2];
  * tail p and tail w: Lemma 2;
  * source (1/D_p − x_par)Δ_u g(u→w) + (1/D_w − x_ch)Δ_u g(p→u).

  The **drop to zero** when u leaves H_u (D_u+1 = D_hu by leaf gain, or k_u = K_hu → K_hu+1 by budding) is in the
  check (hu_up = 0, hu_bud = 0). Same for **hC**, with u-range (h_u+k_u−2)/2 and cu_up/cu_bud = 0 on exit from C_u.
* **Forward invariance (load-bearing, holds)**: a 2-path or cherry whose middle vertex is outside H_u or C_u has value
  0 and gets no inequality. This is sound only because such a term stays 0 after every event except the budding
  entries k: 0→1 (into H_u) and 1→2 (into C_u). Those create only the newborn terms already assigned to edges. D and k
  are nondecreasing and h is fixed, so this holds.
* **Root (h = 0), VSTAR (k > 150), L > 800, clamp κ = min(k, 13)**: all handled. T1(κ,0) and T1(κ,1) states with
  κ ≥ 1 are in P. A tail vertex with k ≥ 13 has γ_min = 0 and no κ change. The VSTAR budding target goes through
  cl_vidx → VSTAR.

### 1.3 Lemma 5 and the induction
* **Lemma 5**: each non-root core w receives one g term and at most one hV term, and each u ∈ C_u receives C(k_u,2)
  cherries. The corner extrema bound the interpolated values. For tail u the extremum over {INF, T1(κ_u,1)} is used,
  and the L ≥ L_F tail states are covered by the L = L_F entry, which is itself tail. So Φ_up ≥ c0·C(T) and
  Φ_lo ≤ C0·C(T) with c0 = W_min, C0 = W_max. The minimum and maximum are attained at (L=0,k=1,h=1) for UPPER and at
  (L=668,k=0,h=1) for LOWER (my verifier, §2).
* **Induction**: E[Φ(T_{n+1})] = E[Φ(T_n)] + E[LΦ(T_n)]/n gives (1 ± λ/n)-recursions from Φ(T3) = f(2,0,0). Both
  directions are valid without any sign assumption on Φ. Dividing by c0 or C0 gives exactly the stated prefactors.
* **Asymptotics**: Γ-asymptotics of the product, plus Borel–Cantelli for Corollary 2; correct.

**Conclusion of task 1: no false step and no missing term found.**

## 2. Verifier vs §4.5 (task 2)

I compared `cert_verify_gmp.c` line by line with §4.5. It checks exactly the listed inequalities, with the same terms,
signs and extremes (`choose_rate` picks r_hi or r_lo by the sign of Δ for each mode). It loops over every table entry
and uses exact mpq/mpz arithmetic.
* The entry counts match full enumeration: g = 555·290 − 27·15 (tail–tail) = 160 545, hV = 555·84·290 = 13 519 800,
  hC = 75·290·291/2 = 3 164 625.
* The integer divisions D_all/D0 are exact: they occur only when one endpoint is T1, so D_all contains D0.
* The tail–tail zero condition is checked.

Nothing in the text is checked more weakly, or omitted.

Neon's verifier, run by me (`taskset -c 0,1`):
```
$ ./cert_verify_gmp ../attack-neon/upper_5665.bin
  [1/5] f checks:      484138  fails: 0  (PASS)
  [2/5] g checks:      160545  fails: 0  (PASS)
  [3/5] hV checks:   13519800  fails: 0  (PASS)
  [4/5] hC checks:    3164625  fails: 0  (PASS)
  [5/5] Global comparison: Phi(T3)=14606753031/10000000000, min W(u)=2989767809/10000000000, exact prefactor A = 14606753031/2989767809 (4.88558108)  (PASS)
  [6/6] Exact brute-force trees: 37 tested, fails: 0 (PASS)
=== SUMMARY: 17329108 exact GMP inequality checks, 0 failures => VERIFIED OK ===
$ ./cert_verify_gmp ../attack-neon/lower_5650.bin
  ... f 484138 / g 160545 / hV 13519800 / hC 3164625, fails 0 each
  [5/5] Global comparison: Phi(T3)=14608704616/10000000000, max W(u)=17527873930/10000000000, exact prefactor A = 7304352308/8763936965 (0.83345560)  (PASS)
=== SUMMARY: 17329108 exact GMP inequality checks, 0 failures => VERIFIED OK ===
```

**Independent verifier** (`attack-neon1-opus/verify_indep.c`, written from §§3–5 only; it does not include Neon's header
or code). It uses GMP rationals for f and g and exact `__int128` with common denominator
M = 2·D0·(D0+1)·S_r·dp·du·dw for hV and hC. It also checks the structural header conditions, g(tail,tail) = 0, and
Lemma 5:
```
$ ./verify_indep ../attack-neon/upper_5665.bin 2
structural header conditions: OK
layout: INNER_N=572 NG=601 npars=555 nchs=290 NHU_ACT=84 NCU_ACT=75 NV=241904
g(tail,tail)==0: OK (0 nonzero)
[f ] checks 484138 fails 0 min slack 9.910720e-06
[g ] checks 160545 fails 0 min slack 1.910560e-07
[hV] checks 13519800 fails 0 min slack 4.558278e-08
[hC] checks 3164625 fails 0 min slack 4.568174e-08
[L5] Phi(T3)=14606753031/S  min W(u)=2989767809/S at (L=0,k=1,h=1)  => prefactor A = Phi(T3)/W = 14606753031/2989767809  (positive OK)
TOTAL failures: 0  => ALL LISTED INEQUALITIES HOLD
$ ./verify_indep ../attack-neon/lower_5650.bin 2
[f ] checks 484138 fails 0 min slack -0.000000e+00      (equality case, see N4)
[g ] checks 160545 fails 0 min slack 1.958692e-07
[hV] checks 13519800 fails 0 min slack 4.532714e-08
[hC] checks 3164625 fails 0 min slack 4.572782e-08
[L5] Phi(T3)=14608704616/S  max W(u)=17527873930/S at (L=668,k=0,h=1)  => prefactor A = Phi(T3)/W = 7304352308/8763936965  (positive OK)
TOTAL failures: 0  => ALL LISTED INEQUALITIES HOLD
```
(Slacks are in natural units, i.e. table value / S.)

**Negative controls for my verifier** (`indep_negative_controls.log`):
* Stricter λ: UPPER 1132/2000 → 1 236 682 failures; UPPER **11329/20000 = 0.56645 → 281 912 failures**;
  LOWER 114/200 → 2 519 701; LOWER 1131/2000 → 982 934. Both certificates are therefore tight to about 5·10⁻⁵ in the
  exponent.
* Single-entry perturbations of entries that are actually used are all rejected:

  | perturbation | failures |
  |---|---|
  | hV[0] ∓ 10⁷ | 1 |
  | hC[0] ∓ 10⁷ | 1 |
  | g[(1,1,0)→nb] −10⁷ | 1 |
  | f[802] ∓ 10⁷ | 810 |
  | f[1000] −10⁵ | 2 |
  | xpar[(1,1,1)] +0.1 | 2506 |
  | xch[(1,1,1)] −0.1 | 141 587 |
  | xcc[(2,2,1)] +0.1 | 1054 |
  | rpar[(1,1,1)] +0.1 | 1 |
  | rch[(0,2,0)] −0.1 | 290 |
  | kpc[(3,9,1)] knot 3 −10⁸ | 289 |
  | kcc[(0,7,0)] knot 2 −10⁸ | 228 |

  Perturbing unused entries (for example x-constants of states outside H_u, or g rows with an h = 0 child) correctly
  changed nothing.

## 3. Numerical attack (task 3)

The tools are `decomp.py` (my own Φ evaluator and Lemma-4 decomposer, exact `Fraction`) and `run_decomp.py`. For every
tree and both certificates:

1. **(LΦ)(T)** is computed exactly by applying every event.
2. **E[τ]** for every term τ (f, g, hV, hC) is computed with the true rates, exact credit transfers and exact FY
   products. The **identity Σ_τ E[τ] == (LΦ)(T)** is checked exactly; this is Lemma 4 with no slack anywhere.
3. **E'[τ]** is the same with each FY product replaced by A(x) on the receiving edge plus A*(y) on the edge itself,
   i.e. the §4.5 left-hand side at the tree's true rates and true t = D0/D. **Every single term** is checked:
   E'[τ] ≤ λτ (UPPER), ≥ λτ (LOWER). This separates term groups, so a violated local inequality cannot hide behind
   slack elsewhere.
4. The tree drift inequality and the floor inequality are checked.

**Controls.**
* My Φ equals `neon_brute.phi` on 1 884 (tree, certificate) pairs: `compared 1884 (tree, cert) pairs, mismatches 0`.
* Planted decomposition errors make the identity fail:

  | planted error | identity failures |
  |---|---|
  | drop r_par from the vertex term | 216 (n ≤ 8) |
  | drop the hV newborn | 216 |
  | drop the hC newborn | 170 |
  | drop the parent-credit FY product | 2 on a tree with a k = 9 vertex |

  The last one gives 0 on n ≤ 8, as expected: such small trees contain no vertex with a parent and a child that is
  outside H_u.

**Results** (no identity, local, tree or floor failure anywhere). Minimal local slack is absolute, in units of S;
"rel" is the relative slack (λτ − E')/|λτ| of that same term. Tree-drift slack is (λΦ − LΦ)/(λΦ).

| family | trees | identity fails | min local slack f / g / hV / hC (UPPER) | min local slack f / g / hV / hC (LOWER) | min tree-drift rel slack UP / LO |
|---|---|---|---|---|---|
| all rooted trees n ≤ 14, root deg ≥ 2 (`small14.log`) | 32 972 | 0 | 9.9e-6 / 1.9e-7 / 4.9e-8 / 4.9e-8 | 7.1e-6 / 2.0e-7 / 4.9e-8 / 4.9e-8 | 7.56e-5 / 5.43e-5 |
| hub families part 1 (`hubs.log`, indices 0–799: H_u-type vertex repeated 3/12 times under the root or under a non-root hub, kids of degree 2…31, chains) | 800 | 0 | 9.9e-6 / 1.9e-7 / 5.0e-8 / 5.0e-8 | 7.1e-6 / 2.0e-7 / 4.9e-8 / 5.4e-6 | 1.04e-4 / 6.50e-5 |
| hub families part 3 (`hubs_part3.log`, indices 1200–1563: repeated cherries with inner / D0−1 / tail children, tail hubs with 13–40 copies of a child type, D = 22…25 vertices with k up to 14, repeated L = 0 chains, k > K_T, VSTAR root k = 152, VSTAR child k = 151, child L = 801, child L = 800 with k = 2) | 364 | 0 | 9.9e-6 / 0* / 5.0e-8 / 4.7e-8 | 0** / 0* / 5.0e-8 / 4.6e-8 | 6.88e-5 / 5.93e-5 |
| random histories, 8 mixed policies (natural, hub, budding, budding small vertices, deepest-vertex budding, steering degrees to 11/12/15/21–23); 40 runs to n = 160 (`policy.log`) | 320 | 0 | 9.9e-6 / 0* / 4.9e-8 / 5.0e-8 | 7.1e-6 / 0* / 4.9e-8 / 5.0e-8 | 1.10e-4 / 6.21e-5 |
| same, 30 runs to n = 400 (`policy400.log`) | 240 | 0 | 9.9e-6 / 0* / 4.9e-8 / 5.0e-8 | 7.1e-6 / 0* / 5.0e-8 / 5.0e-8 | 1.13e-4 / 7.24e-5 |
| same, 60 runs to n = 700 (`policy700.log`) | 420 | 0 | 9.9e-6 / 0* / 4.9e-8 / 5.0e-8 | 7.1e-6 / 0* / 4.9e-8 / 5.0e-8 | 1.85e-4 / 9.97e-5 |

\* 0 is a tail–tail edge (g ≡ 0 and E' ≡ 0): equality, as designed.
\** A VSTAR vertex in the LOWER certificate, where f(VSTAR) = 0 (equality case of the VSTAR b = 0 check; see N4).

Smallest nonzero relative local slack seen: 2.3e-6 (LOWER hV, `Hrep reps=12 Lu=0 ku=1 kid=14`). Smallest floor slack:
LOWER 4.6e-2 (policy), UPPER 0.98 (chain repetition, n = 361).

`hubs_part2b.log` reruns hub indices 800–1199 (repetition 30, capped at n ≤ 1200) plus 1200+. It was still running
when I finished this report. Its last line then was
`--- hubs 900/1564 (50s): trees=62 failures(identity/local/tree/floor)=0`, with 0 violations. It stops at its 420 s
time limit; any later lines in the log are genuine output that this report does not cover. Most hub trees in index
range 800–1199 have n > 1200 and were skipped as too slow for the exact brute-force drift. The original hub run was
killed by the service restart at index 800; its output up to 800 is in `hubs.log`.

## 4. Verdicts (task 4)

* **UPPER (1133/2000): HOLDS.**
  * Lemma 4's decomposition is correct term by term, and it is confirmed as an exact rational identity on 34 776
    trees.
  * The §4.5 inequalities are exactly those the verifier checks, and they pass in Neon's GMP verifier and in my
    independent int128/GMP verifier.
  * Lemma 5 and the induction are correct, and the prefactor recomputes to 14606753031/2989767809.
  * The certificate is tight: λ = 0.56645 already fails. 1133/2000 < 17/30 by exactly 1/6000.
* **LOWER (113/200): HOLDS**, with the same evidence. The prefactor recomputes to 7304352308/8763936965.

## 5. Nits (no effect on correctness)

* N1. §4.5 defines v⁺ := T1(k_v, h_v) for D_v = D0 − 1. It should read T1(κ(k_v), h_v): k_v can be up to 23 > K_T + 1.
  The code clamps (`cl_t1_idx`, and `t1()` in mine).
* N2. Lemma 4 is titled "Exact … Decomposition", but it is an inequality because of FY and the rate ranges. The text
  itself uses "≶" correctly.
* N3. The never-re-enter property of H_u and C_u is presented as a remark, but it is needed (see §1.2). It is true.
* N4. The LOWER certificate has f(VSTAR_0) = f(VSTAR_1) = 0. The VSTAR inequality with b = 0 (0 ≥ λ f(VSTAR)) then
  holds with equality, which "≥" allows.
* N5. For VSTAR, the text's bound D ≥ K_F + 1 + h is loose: budding needs L ≥ 1, so D ≥ K_F + 2 + h. The check covers a
  superset, which is fine.
* N6. Neon's verifier hard-codes `num_threads(3)`. I ran it under `taskset -c 0,1`.
* Out of scope: RESULT.md's remark combining this with Astra's existence theorem was not checked here.

## 6. Files and commands (all in `checks/attack-neon1-opus/`)

* `neon_verify_{upper_5665,lower_5650}.log`: Neon's verifier output.
* `verify_indep.c`, `indep_{upper_5665,lower_5650}.log`, `indep_negative_controls.log`:
  `gcc -O2 -fopenmp verify_indep.c -lgmp -o verify_indep; taskset -c 0,1 ./verify_indep CERT 2 [TABLE INDEX DELTA | lam P Q]`.
* `decomp.py` (the evaluator; set `DECOMP_PLANT=rpar|hvnb|hcnb|pc` for the planted controls), `run_decomp.py`
  (modes `small N`, `hubs [START]` with `HUB_NMAX`, `policy RUNS NMAX SEED EVERY`, `tree PARENTARRAY`), and
  `cmp_neon_phi.py`.
* Logs: `small9.log`, `small14.log`, `hubs.log`, `hubs_part3.log`, `hubs_part2b.log`, `policy.log`, `policy400.log`,
  `policy700.log`.

Everything used at most two cores (`taskset -c 0`, `-c 1`, `-c 0,1`). Nothing outside this folder and this report was
written.
