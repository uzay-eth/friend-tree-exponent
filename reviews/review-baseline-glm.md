> Packaging note: one local filesystem path in the Strogatz quote row was redacted for public release; all mathematical content is unchanged.

VERDICT: CORRECT
BOUNDS: [0.543, 0.586] — for every n ≥ 3, 0.50232622·Γ(n+0.543)/Γ(n) ≤ E[C_n] ≤ 2.16187785·Γ(n+0.586)/Γ(n) (so 0.502·Γ(n+0.543)/Γ(n) ≤ E[C_n] ≤ 2.17·Γ(n+0.586)/Γ(n) as claimed), and a.s. limsup log C_n/log n ≤ 0.586 — both final certificate tables verified exactly by my own verifier in checks/strogatz-glm/

Independent check by a different model family (GLM). The author's verifier code (`certs*.py`) and the
author's `independent/` folder were **not read**. The certificate tables were regenerated with the
author's `dump_tables.py` (found via `run_all.sh`) and read from the dumped text files only.

All source URLs below were read on 2026-11-21 (UTC).

## 1. Model and literature (REPORT §1–2)

| Claim in REPORT | My check | Result |
|---|---|---|
| Model: start from a single edge, V uniform, W uniform among V's neighbours | KR abstract "complete redirection: a new node randomly selects an existing target node, but attaches to a random neighbor" ([arXiv:1706.01514](https://arxiv.org/abs/1706.01514)); ABBDDKL abstract "sequential manner one by one by selecting an existing target vertex and connecting to one of its neighbours (or friends), chosen uniformly at random" ([arXiv:2403.20185](https://arxiv.org/abs/2403.20185)) | same model in all three sources ✓ |
| KR call it isotropic redirection (IR), r = 1 | TeX source (https://arxiv.org/e-print/1706.01514), "We define this process as isotropic redirection (IR)", complete redirection | ✓ |
| KR quote, §5.1 eq. (32): 𝒩 ∼ N^μ with μ ≈ 0.566, slopes 0.5652–0.5673, "accurate to within 0.001" | verbatim at lines 637–651 of paper.tex; 31 numbered display equations precede `\label{core}`, so it is indeed eq. (32); the growth/sum equations {growth:core},{nu:sum} are then (34),(35) | ✓ exact quote, ✓ eq. (32), ✓ "sum rule (34)–(35)" |
| "17/30" absent from papers | grep of KR TeX, of the ABBDDKL arXiv-v1 PDF text and of the published RSA PDF | ✓ nowhere (nor in the 2017 paper's TeX) |
| ABBDDKL, Random Structures & Algorithms 69(1) e70082, online 2026-07-30, arXiv v1 only (2024-03-29) | CrossRef record for DOI 10.1002/rsa.70082 (Random Structures & Algorithms 69(1), published 2026-07-30); arXiv submission history (v1 only, 29 Mar 2024); published PDF at https://luc.devroye.org/Addario+Briend+Devroye+Donderwinkel+Kerriou+Lugosi-Random_friendship_trees-2026-.pdf | ✓ all three details exact |
| Thm 4.7 "0.1 < δ ≤ λ < 0.9" (arXiv numbering) | rendered text of https://arxiv.org/pdf/2403.20185v1 | ✓ verbatim |
| "their proof gives 0.9142" (Prop 6.22: δ = 1 − γ/2, γ = 3 − 2√2, "omitting the details of the coupling") | Prop 6.22 verbatim: "Let 𝛿 = 1 − 𝛾/2, where 𝛾 = 3 − 2√2. Then, for any 𝛿 < 𝜆 < 1, n^{−𝜆}X_n^{≥2} → 0 a.s."; proof says "Omitting the details of the coupling, we conclude that for 𝜆 ∈ ((1+𝛽)/2, 1)…"; 1 − (3−2√2)/2 = √2 − 1/2 ≈ 0.91421 | ✓ arithmetic ✓, and the REPORT's reading — that Thm 4.7's own λ < 0.9 (and the abstract's n^{9/10}) is stronger than what the written proof delivers (λ > 0.9142) — is a fair and, as far as I can tell, correct reading |
| δ < (√13−3)/6 ≈ 0.1009 (Prop 6.20) | Prop 6.20 verbatim: α = (√13−3)/2, any 0 < δ < α/3 ≈ 0.101; (√13−3)/6 = 0.100925… | ✓ |
| Strogatz post quote | archived local copy of post `1572287647715934210`: "they discovered a weird exponent of 0.566 in it. It's an open problem to prove that 17/30 is really the exponent" | ✓ verbatim (x.com itself blocks fetching; replies not retrievable, agreeing with the REPORT) |

Not read, with reason: the Briend thesis (HAL bot-blocked, per REPORT; not attempted) and the Hartle
thesis (HTTP 403, per REPORT; not attempted); the IOP published PDF of the 2017 paper (IOP page
returned a bot-block HTML; content verified from the arXiv TeX source instead); the context-only
citations of REPORT §2 rows 4–8 (Cannings–Jordan 2013, KR 2024 review, Hartle et al. 2026, Levens
et al. 2022, Zhao et al. 2026 etc.) were not re-fetched, since none of them enters the proof.

## 2. Proof by hand (REPORT §3.1–3.5)

**Lemma 1 (one-step dynamics).** P(W_n = w | T_n) = (1/n)Σ_{v~w} 1/D_v follows directly from V
uniform, W uniform among V's neighbours. (i) for u ∈ K: leaf neighbours contribute 1 each and core
neighbours give ρ_u. (ii) a leaf ℓ of u has its single neighbour u, so P(W = ℓ) = 1/(nD_u), and
u → (L−1, c+1) while ℓ becomes core with (L,c) = (1,1). "Every leaf of T_n (n ≥ 3) has a non-leaf
neighbour": otherwise the leaf-leaf edge would be an isolated component, forcing n = 2. For n ≥ 3,
C_{n+1} − C_n = 1{W_n is a leaf}, hence E[C_{n+1} − C_n | T_n] = Σ_u L_u/(n D_u) = S_n/n (★).
I verified (★) exactly by enumerating all 362,880 trees at n = 10 (below), and Lemma 1's invariants
(degrees never decrease, core stays core, T_3 = path a.s., T_2 the initial edge) are all elementary
and correct.

**Lemma 2 (event list / local states).** The move list is exhaustive: a move attaches the new vertex
to W, which is either core (E1 for W; E3/E4/E5/E6/E7-type class changes for W's neighbours; label
changes) or a leaf (E2 for the leaf's neighbour; the newborn core vertex is the converted leaf).
Checked the classification algebra by hand: pendant class ℓ means (L,c)=(ℓ,1), class B means (0,2);
every "other" neighbour has degree ≥ 3 (a degree-2 core vertex is (1,1), i.e. a class-1 pendant, or
(0,2)); c never decreases so pendant→B/"other" transitions are one-way, and "other" stays "other"
for ever; a remembered creator is "other" at birth, hence for ever. The rate formulas are exact
where claimed: E3 = m_ℓ(ℓ + 1/D) because a class-ℓ pendant's credit seen from w is exactly 1/D_u;
E4 = m_ℓ ℓ/(ℓ+1) from (ii) with D_w = ℓ+1; E5's per-neighbour rate 1/D_u + 1/D_x ∈ [1/D, 1/D + 1/2];
the credit intervals follow from 1/D_v ≤ 1/2 for core v (pendant/B credits exact, others ≤ 1/3, the
static label P = min(D_birth, 20) satisfying 1/D_current ≤ 1/P since degrees only grow). All of this
is *independently confirmed by brute force* below (exact equality of the (★★) drift with true rates,
per move, on random trees). I found **no error** in the final-bound-relevant event list.

**Lemma 3 (comparison) and the certificates→bounds step.** 𝒜f is affine (multilinear) in each
uncertain rate and the rates enter separate terms, so min/max over the product of intervals is the
sum of per-term endpoint choices — correct, no monotonicity of f needed. True rates lie in the
intervals (Lemma 2, verified), so min 𝒜f ≥ δf bounds the true drift. T_3 is a.s. the path (out of
T_2, V ∈ {1,2} and W is forced to be the other vertex), so Φ_3 = f(centre) with label/no-claim is
deterministic. Induction over n ≥ 3 gives E[Φ_n] ≥ Φ_3 Π and ≤ Φ_3 Π; with f ≥ 0, sup f for the
lower bound (C_n ≥ Φ_n/sup f pointwise) and with f ≥ 1 (upper, C_n ≤ Φ_n) yields the theorem. The
Γ conversion Π_{m=3}^{n−1}(1+a/m) = 2Γ(n+a)/(Γ(n)Γ(3+a)) is an exact rational identity for integer
n (telescoping Γ(z+k)/Γ(z)); the constants reproduce (float lgamma, cross-checked against an
independent Stirling evaluation): 0.87551816·2/Γ(3.543) = 0.50232622 and 3.95462373·2/Γ(3.586) =
2.16187785 — matching the REPORT to every printed digit. Base case n = 3: product empty, so
0.87551816 ≤ 1 = E[C_3] ≤ 3.95462373. The rounded theorem constants 0.502 < 0.50232622 and 2.17 >
2.16187785 are sound roundings (NOTE: REPORT §1 says "0.50·", §3.1 "0.502·" — both ≤ 0.50232622).

**Almost-sure corollary.** C_n is nondecreasing (core vertices keep degree ≥ 2). Markov:
P(C_{2^j} > 2^{j(λ+ε)}) ≤ E[C_{2^j}] 2^{−j(λ+ε)} = O(2^{−jε}) via the mean upper bound; summable;
Borel–Cantelli (no independence needed) gives C_{2^j} ≤ 2^{j(λ+ε)} eventually, and monotonicity gives
C_n ≤ C_{2^j} ≤ (2n)^{λ+ε} for 2^{j−1} < n ≤ 2^j. Hence limsup log C_n / log n ≤ λ + ε for all ε,
i.e. ≤ 0.586. Correct.

**The "infinitely many states outside the finite table" argument (most important).** The state
space is covered by: box (L ≤ 60, c ≤ 10), tail L > 60 (flat extension kills E1 and makes the
newborn type constant; the rate extrema L/D ∈ {61/(61+c), 1}, 1/D ∈ {0, 1/(61+c)} are the achieved
endpoints over the whole tail and enter separate affine terms — sound), and c > 10 (lower: f = 0
makes the requirement automatic once f ≥ 0 everywhere; upper: f = G and Af = (L/D) f(nb) ≤ max f(nb)
≤ λG — verified exactly: max over all newborn-type values = 2.91903 ≤ λG = 2.92147, which turns out
to be the binding constraint of the whole upper certificate, matching its minimum slack
2.4452e-03 = λG − max f(nb)). I could not break this argument.

**Minor notes, none of which affects the result** (by section and line):
* §3.5(c): "max f(nb) ≤ λG" is correct only with "nb" read as the *newborn states* (a newborn
  always has c = 1 and L = 1); as luck would have it, that stronger-looking claim also holds for
  this table. A careless reader could mistake it for the max over all c = 1 states, which is
  4.5005 > λG (those states, e.g. (60,1,1,0,7), exist in the box but are not newborn states).
  My verifier checks both readings and reports both pass.
* §3.5: the 90 upper-table states with L + c < 2 (e.g. (0,0) and (0,1) entries, where f = 0) can
  never occur as the state of a core vertex (degree ≥ 2 by definition), so the drift inequality is
  correctly not required there; the certificate's λf = 0 would be violated by a positive drift at
  those unreachable states, so they must be (and my verifier does) skip them explicitly.
  The REPORT does not spell this out.
* §3.5(b): the interval endpoints are infima in places (ρ lower endpoint, 1/D → 0, L/D → 1); since
  the drift is affine and the slack inequality is one-sided, using infima is conservative and sound.
* Upper certificate check counts: my verifier's 313,938 box + 5,148 tail = 319,086 exact checks
  coincides exactly with the author's "319 086 exact checks" — strong evidence the two verifiers
  interpret the state space identically. (The lower counts differ: mine 1,186,328 + 8,008 =
  1,194,336 vs the author's 1,205,767; presumably a bookkeeping difference over degenerate/region-c
  states; every box state except nothing is covered by mine, so this is immaterial.)

## 3. Independent verifier (my code, from §3.1–3.5 only): checks/strogatz-glm/

`verify_glm.py` (exact rationals; table entries parsed as dyadic rationals via `float.fromhex` →
integer multiples of 2^−62, exactness asserted; all rates `Fraction`; the certificate inequality as
an exact rational comparison). Run via `run_verify_all.sh` (log: `verify_glm.log`):

* **Final lower table** (lower9b; state (L,c,m1,m2,b,n3,n4,n5), K0 = 2, Dcap = 6, δ = 543/1000):
  box 1,186,328 checks + tail 8,008 checks, **0 failures**. Min slack 0 at dead states (L = 0, f = 0);
  min slack among f > 0 states **1.0834e-06** (the author's independent check reported 1.1e-06).
  sup f = 0.63906883, Φ_3 = f(2,0,0,0,0,0,0,0) = 0.55951636, Φ_3/sup f = **0.87551816** ⇒
  E[C_n] ≥ 0.87551816·Π = **0.50232622**·Γ(n+0.543)/Γ(n). **PASS**.
* **Final upper table** (state (L,c,m1,b,P), K0 = 1, Pmax = 20, λ = 293/500): f ≥ 1 on every
  degree-≥2 box state; G = 4.985445 ≥ 1; box 313,938 + tail 5,148 checks (90 degree<2 states skipped
  as unreachable), **0 failures**; min slack **2.4452e-03** (author: 0.00245); region c > 10 holds
  in both readings. Φ_3 = f(2,0,0,0,3) = 3.9546237 ⇒ E[C_n] ≤ 3.95462373·Π =
  **2.16187785**·Γ(n+0.586)/Γ(n). **PASS**.
* **Negative controls**: (i) *perturbed table must fail*: f raised by ~2×(min slack)/δ at the
  minimum-slack state — the lower re-check fails at 1 state, the upper at 3 states: **detected**;
  (ii) final lower table at δ = 0.55: fails (1,186,310 box failures), as expected; (iii) final
  upper table at λ = 0.58: fails (308,808 box failures), as expected. (iv) As a warm-up I also
  verified the 50×9 lower table (δ = 27/50): PASS, with constants 0.87910243/0.50607536 matching
  the REPORT's exactly.
* **Per-step expected change by exact enumeration** (`check_drift_glm.py`, Part A): **all** trees
  reachable by the process up to n = 10 ((n−1)! = 362,880 trees at n = 10, each with its exact
  rational probability; total probability mass exactly 1). For every 3 ≤ n ≤ 9,
  E[C_{n+1}] − E[C_n] and E[S_n]/n are *exactly equal* rationals (Lemma 1 (★)), and the E[C_n]
  agree exactly with the values recorded in `outputs/run_all.log`: E[C_3] = 1, 4/3, 19/12, 647/360,
  17189/8640, 4717423/2177280, 101568829/43545600, 200684675729/80621568000 = 2.489218216...
* **Drift identity (★★) per move** (`check_drift_glm.py`, Part B): for 12 random trees (sizes up to
  45) and a pseudo-random rational test function, the event-list decomposition with *true* rates
  (computed exactly from the tree, including the static-creator-label bookkeeping for the upper
  scheme and the tracked-degree state for the lower scheme) reproduces exactly, in `Fraction`
  arithmetic, the conditional expectation of the change of Φ computed by brute force over every
  (V, W) move. Deliberately wrong rates (E3's 1/D term set to 0; E5's per-neighbour rate set to 1/2)
  are detected — the test is discriminating. **PASS**.

## 4. Simulation (small size only)

Compiled the author's `sim.c` (`gcc -O2`) and ran it fresh: `./sim 10000 400 7001` (400 runs to
n = 10^4, seed not used by the author). Analysed with my own `sim_identity_check.py`
(`sim_identity.log`):
* At every checkpoint n, the local log-slope of E[C_n] from finite differences over neighbouring
  checkpoints vs the ratio of run means E[S_n]/E[C_n]: agree within 0.018 across n = 13 … 6310,
  consistent with 400-run sampling noise (the identity must hold in expectation; the exact
  enumeration above already proves it for all states).
* Independent per-step check with a tiny Python simulator of my own (different RNG and code path):
  mean of the indicator 1{W_n is a leaf} (= C_{n+1} − C_n) vs mean of S_n/n, paired, 4000 runs at
  n ∈ {30, 61, 137, 316, 717, 1622}: max |z| = **1.29**. **PASS** — the reported identity
  E[C_{n+1}] − E[C_n] = E[S_n]/n is confirmed numerically at small size, and exactly at all sizes.

## 5. What is left unchecked

* The long simulation runs (10^7, 10^8) and everything derived from them (μ ≈ 0.5657 plateau,
  the case against 17/30 = 0.56667) were **not** rerun, per instructions; the "17/30 is the
  exponent?" question remains open, exactly as the REPORT says.
* The secondary certificates (certs6 lower δ = 0.531, certs7 dynamic label δ = 0.535, certs8
  recorded degrees, the simple (L,c) pair, and the sensitivity experiments) were not re-verified in
  full; the *final* pair (certs9 tracked degrees δ = 0.543; static label λ = 0.586) — the only two
  tables the published constants rest on — are verified here, plus the 50×9 table (same state as
  the final lower, δ = 0.54) as a warm-up.
* The tail states were checked at the analytically derived rate extrema rather than at every L;
  this covers all L > 60 simultaneously as argued in §3.5(b) and re-derived by hand here.
* The two theses (Briend, Hartle) and the X replies could not be fetched (see "not read" list).

## 6. Files in this directory

| file | role |
|---|---|
| `verify_glm.py` | my independent exact verifier for the two final tables (from REPORT §3.3–3.5 only) |
| `run_verify_all.sh` | runs all checks + negative controls; regenerates the tables via the author's `dump_tables.py` if absent (tables not committed: ~89 MB) |
| `check_drift_glm.py` | exact enumeration of all trees to n = 10; per-move brute-force of the (★★) drift identity for both final state schemes; broken-rate negative controls |
| `sim_identity_check.py` | small-size simulator identity check (author's `sim.c` output + my own Python simulator) |
| `verify_glm.log`, `check_drift_full.log`, `sim_identity.log` | outputs of the runs described above |
