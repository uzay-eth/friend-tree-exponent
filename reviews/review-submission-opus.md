# Referee report: "Growth of the non-leaf population in the random friend tree"

**R1 · Claude Opus 5.5 · maximum reasoning.** Independent cold report on §§1–5 at `526f350`.

**Manuscript.** `paper/friend-tree-exponent.tex` (one author; dated 23 November 2026), read at repository snapshot `526f350`.

**Standard applied.** Annals of Applied Probability / Electronic Journal of Probability.

**Scope of this report.** Sections 1–5 were close-read in full, including every proof behind the three headline claims and the whole motif drift argument of Section 3. Sections 6–7 were checked only for consistency with the data files and for exposition. The referee did not consult any earlier review of this work.

**How the audit was done.** Every logical step of Sections 1–5 was re-derived by hand. All computer-assisted parts were rebuilt from source and rerun. In addition, the referee wrote four independent programs that share no code with the supplement:

- an exact enumeration of the chain on unlabelled trees;
- an implementation of the motif potential taken from the paper's text alone, with a brute-force drift test;
- an implementation of the banded label-update rules taken from the paper's text alone, with a brute-force drift test on long histories;
- a small PDF text extractor, used to check the cited sources.

All arithmetic in these checks is exact (Python `fractions` or the supplement's GMP and 128-bit integer code).

---

## Summary for the editor

The paper proves that the number $`C_n`$ of non-leaves of the random friend tree satisfies $`\mathbb E C_n=n^{\mu+o(1)}`$ and $`C_n=n^{\mu+o(1)}`$ almost surely, for one deterministic $`\mu`$. It also proves the explicit bracket $`113/200\le\mu\le1133/2000`$, so $`\mu<17/30`$. A second, independent certificate gives $`\mu\le 5657/9983`$.

**No proof step was found to be false, and none was found that cannot be justified.** A few steps need a short repair, written out below; most of these are missing one-line computations or missing statements of a coupling fact. Every number in Sections 1–5 that the referee checked agrees with the proofs, the certificate files and the checker output. That covers all constants, counts, digests and rate formulas. The computer-assisted parts reproduce digit for digit. The referee's independent direct drift tests also pass in the large-degree and band-crossing regimes, which the paper's own direct tests never reach.

Status of the three headline claims:

| Claim | Status |
|---|---|
| The exponent is below $`17/30`$ | **Established.** Two independent certificates (motif, $`1133/2000`$; banded, $`5657/9983`$) were reproduced and tested directly beyond the paper's test range. Note that $`\limsup_n \log\mathbb E C_n/\log n<17/30`$ and the almost-sure $`\limsup`$ do not need Theorem 1.1. |
| Existence (Theorem 1.1, mean and almost sure) | **Established.** No false or unjustifiable step. The argument is a dense sketch, and three places need written repairs (E7, E22, and stating E2 as a lemma). |
| Bracket (Theorem 1.2) | **Established, computer-assisted.** The reduction (Lemmas 3.1–3.2, Section 4) is sound and every finite check was reproduced. The constants match. The paper must identify the motif arrays it certifies (see Issue 2). |

**Recommendation: major revision.** The mathematics appears correct. The revision concerns verifiability and presentation at journal standard: Section 2 compresses a long coupling argument into a few paragraphs; the certified objects are identified only partly in the paper; and several citations need correction. The three main issues are listed in §8.

---

## 1. The claims, in the referee's words

**Model.** $`T_2`$ is one edge. Given $`T_n`$, pick a uniform vertex $`V_n`$ and a uniform neighbour $`W_n`$ of $`V_n`$, and attach a new leaf to $`W_n`$. This is the 1-step random-walk attachment of Saramäki–Kaski, the "isotropic redirection" of Krapivsky–Redner, and the "random friend tree" of Addario-Berry et al. $`C_n`$ is the number of vertices of degree at least 2 (the core).

**Theorem 1.1.** There is a deterministic $`\mu\in[0,1]`$ with $`\log\mathbb E C_n/\log n\to\mu`$ and $`\log C_n/\log n\to\mu`$ almost surely.

*Method.* Embed the chain in continuous time with rate-one clocks. In $`k+1`$ time units, any eligible core vertex can plant a prescribed $`k`$-vertex rooted seed with probability at least $`p(M,k)>0`$. The planted branch is then compared, through an acceptance-weighted coupling across its cut edge, with a free copy of the seed. This gives the approximately super-multiplicative inequality

```math
F(t+h+s)\ \ge\ F(t)\,\frac{p}{4k}\,e^{-\eta s}\bigl(F(s)-k-O_k(e^{-s})\bigr),\qquad F(t)=\mathbb E\,C(t),
```

with loss $`\eta\to0`$. It yields a common growth rate for $`F`$, which transfers to $`\mathbb E C_n`$. For the almost-sure statement, at least $`C/6`$ sites with disjoint clocks run independent planting experiments. Chebyshev's inequality and a martingale/Kronecker argument then give a geometric lower growth on almost every path. Markov's inequality and Borel–Cantelli give the upper limit.

**Theorem 1.2.** For every $`n\ge3`$, with $`P_a(n)=\prod_{j=3}^{n-1}(1+a/j)`$:

```math
\frac{7304352308}{8763936965}P_{113/200}(n)\le\mathbb E C_n\le\frac{14606753031}{2989767809}P_{1133/2000}(n),
\qquad
\mathbb E C_n\le\frac{1539738114735}{3908512529}P_{5657/9983}(n).
```

*Method.* The proof uses potentials $`\Phi`$ that are sums of local weights, with $`\mathcal L\Phi\le\lambda\Phi`$ (or $`\ge`$) for the scaled generator $`\mathcal L\Phi(T)=n\,\mathbb E[\Phi(T_{n+1})-\Phi(T_n)\mid T_n=T]`$. Two such potentials are built:

- **Motif potential (Section 3).** Weights sit on vertices, directed core edges, vertical 2-paths and cherries. Degrees $`D\ge 24`$ are represented by the exact two-point interpolation $`t=24/D`$. The drift inequality on all trees reduces to about $`1.7\times10^7`$ rational inequalities per direction, which are checked in GMP.
- **Banded potential (Section 4).** Each core vertex carries a weight that depends on a history-dependent label of its neighbourhood. Two exact transfers redistribute the drift. The weights come from a 226,545,617-row integer recursion, closed by trial rewards for newborns and an affine tail for $`D\ge2400`$.

**Corollaries stated in Sections 1–4.** $`0.5650\le\mu\le0.5665`$, with gaps $`1/6000`$ and $`1/299490`$ below $`17/30`$. The mean of each exact-degree count $`N_{=k}(n)\le C_n`$ obeys both upper bounds. Almost surely $`\limsup_n\log C_n/\log n\le 1133/2000`$, and $`\le 5657/9983`$.

---

## 2. What was rerun and what was computed independently

All runs used the snapshot's sources on x86-64 Linux with GCC 13.3 and GMP 6.3. Relative paths are relative to the repository root. The independent referee scripts are described at the end of this report.

| Check | Outcome |
|---|---|
| Motif arrays rebuilt by `verify/motif-bracket/cert_build.c` for $`(1,1133,2000)`$ and $`(2,113,200)`$ | SHA-256 `3a4aef88…0abe` (upper) and `9e966fba…60d9` (lower). These equal the digests in the supplement (`proofs/bracket-motifs.md`, §8) and are bit-identical for 1, 5 and 16 OpenMP threads. |
| `cert_verify_gmp.c`, each direction | 484,138 / 160,545 / 13,519,800 / 3,164,625 checks, 0 failures. Comparison constants $`2989767809/10^{10}`$ and $`17527873930/10^{10}`$; $`10^{10}f_\pm(2,0,0)=14606753031,\ 14608704616`$; prefactors exactly as in Theorem 1.2. The 37 built-in trees pass. |
| Motif negative controls | 10 of 10 rejected with exit code 2: exponent override in both directions, and single-entry perturbation of $`f,g,h_V,h_C`$ in both directions. |
| `verify/banded-upper/check.cpp`, built with UBSan | 226,545,617 rows, minimum slack 3, newborn reserve 4,340,525, weights in $`[3908512529,\,1764716280715]`$, centre weight 1539738114735. The SHA-256 of `recipe.txt.gz` equals the digest printed in §4.4. |
| `verify_banded.cpp` (independent verifier), built with UBSan | Same rows and band counts; slack 0; reserve 4,341,454; centre 1539738113358; minimum 3908512428. It chooses different (also admissible) weights. |
| Self-tests of the independent banded verifier (`test_banded_audit.py`) | The smoke certificate (exponent $`3/5`$) passes; 436 exhaustive histories ($`n=3`$–$`7`$) and 50 multi-band trees pass; all 14 negative controls (3 tree-level and 11 command-line) are caught, as §5 states. |
| Tail inequalities of §4.3 | $`5657A-9983(B+a_*)=6235`$ and $`5657(A+B)-9983s=4755`$, as stated. |
| Supplement's `verify_python_cross.py --max-tree-size 14` | Brute-force drift on all 32,972 rooted trees ($`3\le n\le 14`$), both directions: pass. Theorem 1.2 at $`3\le n\le 21`$: pass. |
| Supplement's `verify/tree-drift/motif_identity.py` (term-by-term check of Lemma 3.1) | 172 sampled or constructed trees (up to $`n=254`$) plus 5 boundary trees (up to $`n=1104`$). The identity, every local inequality and the comparison with $`\lvert K\rvert`$ all hold: 0 failures. |
| **Referee:** exact $`\mathbb E C_n`$ by own enumeration (`exact_mean.py`) | Identical to the supplement's file for $`3\le n\le17`$; e.g. $`\mathbb E C_{13}=863016923161025341999/295667991134208000000`$. |
| **Referee:** Lemma 2.1, site count of §2.4, promotion rate | Over every tree with at most 17 vertices: eligible/core $`\ge4/11`$; sites/core $`\ge3/10`$; $`\sum_uL_u/D_u\ge1`$. |
| **Referee:** motif $`\Phi`$ from the paper's text plus brute-force drift (`motif_direct.py`) | The local update agrees with global recomputation on validation trees. 135 trees per direction: friend trees up to $`n=3500`$, plus adversarial trees with tail endpoints, $`D\in\{23,24\}`$, $`\kappa=K_T+1`$, $`L\approx L_F`$, $`k>K_F`$. Altogether 646 (upper) and 566 (lower) tail vertices. Results: upper $`\max\mathcal L\Phi/\Phi=0.56642885<0.5665`$; lower $`\min=0.56503801>0.565`$; the comparison with $`\lvert K\rvert`$ holds; 0 failures. |
| **Referee:** banded labels from the text of Lemma 4.1 plus brute-force drift (`banded_direct.py`) | Two batches, 16 histories in all: 4 to size 30, and 12 to size 220 (6 friend-tree histories and 6 adversarial ones that favour hubs and conversions). They reach degree 207 and cross the band boundaries 16, 32, 80 and 200. The invariant of Lemma 4.1 was asserted after every step. 2,712 history times; 1,016 + 34,004 label-state queries, answered by the production verifier through its read-only observer. Results: $`\max n\mathbb E[\Delta\Phi]/\Phi=0.566661576<5657/9983=0.566663328`$; 0 failures. |

The referee's direct tests are the only ones that exercise the tail interpolation, the leaf and child overflow states, band crossings and pooled-mark re-encoding on actual trees. The paper's own direct tests (§5) stop at 14 vertices and at size 9, where none of this logic is active: a tail state needs at least 25 vertices, and a band crossing needs a vertex of degree 16.

---

## 3. Step-by-step audit of Sections 1–5

Verdict codes: **H** = holds; **H+R** = holds with a repair, given in full; **NJ** = not justified as written; **F** = false. Line numbers refer to `paper/friend-tree-exponent.tex`.

### 3.1 Section 1 (lines 26–61)

| ID | Location | Step | Verdict | Notes |
|---|---|---|---|---|
| I1 | l.27 | Model, $`C_n`$, notation $`N_{\ge k},N_{=k}`$ | H | Agrees with [AB, §1]. |
| I2 | l.29–31 | Literature statements | H+R | Factual content checked (§5 of this report). Repair: the Strogatz post states a conjectured value; it does not ask a question (see §5). Cite [AB] in its published version. |
| I3 | l.34–40 | Theorem 1.1 | H | Proved in Section 2 (§3.2 below). |
| I4 | l.41–51 | Theorem 1.2 | H | Proved in Sections 3–4; constants reproduced exactly. |
| I5 | l.52 | $`\log P_a(n)=a\log n+O(1)`$; $`0.5650\le\mu\le0.5665`$; gaps $`1/6000`$, $`1/299490`$; $`N_{=k}\le C_n`$ | H+R | All arithmetic checked. Repair: say that the bracket on $`\mu`$ uses Theorem 1.1, while $`\limsup\log\mathbb E C_n/\log n\le1133/2000`$ does not. |

### 3.2 Section 2: existence and the sample-path exponent (Theorem 1.1)

| ID | Location | Step | Verdict | Notes |
|---|---|---|---|---|
| E1 | l.64 | Continuous-time embedding; jump chain independent of holding times; non-explosion | H | Waiting time at size $`j`$ is $`\mathrm{Exp}(j)`$; $`\sum_j(E_j-1)/j`$ converges in $`L^2`$. |
| E2 | l.66–76 | Coupling principle: fresh exponentials and marks at rings or deterministic times, induction over jump indices, then non-explosion | H+R | Correct but informal. Repair: state it as a lemma. "Any construction that, at each jump, gives the correct conditional law of (holding time, source, target) given the revealed past produces the friend-tree process; this applies at stopping times by the strong Markov property of the clocks." The split choice "first $`d`$ with acceptance $`d/(d+1)`$" is a special case. |
| E3 | l.78 | Yule decomposition; $`\mathbb E N=ke^t`$, $`\mathbb E N^2=k(k+1)e^{2t}-ke^t`$ | H | From $`\mathrm{Var}(\mathrm{Geom}(e^{-t}))=e^{2t}-e^t`$. |
| E4 | l.81–84 | Lemma 2.1: at least $`C/4`$ eligible core vertices | H | The argument gives $`E\ge(C+1)/4`$. Exhaustive check up to 17 vertices: $`E/C\ge4/11`$. The constant is asymptotically sharp: take a cubic core tree with every core edge subdivided by a bare vertex. |
| E5 | l.90 | Lemma 2.2, stage 1: probability $`\ge e^{-2}(1-e^{-1})/2`$ | H | In the bare case $`w`$ keeps degree 2 until its first ring because $`u,z`$ are silent. The leaf case has probability $`e^{-1}`$, which is larger. |
| E6 | l.90 | Stage 2: $`e^{-1}\Pr(\mathrm{Pois}(1)\ge M)`$; $`L_u\ge M+1`$ and $`c_u\le2`$ are preserved | H | $`c_u`$ changes only at a ring of $`u`$. |
| E7 | l.92 | Stage 3: "a mark below $`M/(M+2)`$ forces a promotion"; probability $`e^{-3}M/(M+2)`$ | H+R | The inequality is not shown. Repair: exterior arrivals only add leaves at $`u`$, and $`c_u\le2`$ until $`u`$ rings, so at the ring $`L_u/D_u\ge(M+1)/(M+3)>M/(M+2)`$, because $`(M+1)(M+2)-M(M+3)=2`$. With leaves ordered first, a uniform mark below $`M/(M+2)`$ therefore selects a leaf. |
| E8 | l.92–96 | Stage 4: $`(e^{-(k+1)}/(k+1))^{k-2}`$ | H | The source has at most $`k-1`$ neighbours, so the stated bound $`k+1`$ is loose but valid. Silence requirements fall on disjoint clocks. |
| E9 | l.97 | $`D_u\ge M`$; the branch contains no old core | H | $`D_u\ge M+1`$ in fact. |
| E10 | l.99 | Branches planted at distinct old cores stay disjoint | H | Each branch is a component of $`T\setminus K_t`$ hanging from a distinct old core. |
| E11 | l.102–111 | Cut-edge comparison and inequality (2.1) | H | The construction gives uniform actual choices. Marks of $`u`$ in $`[0,1/M]`$ form a rate-$`1/M`$ Poisson process independent of the reference. The two tests imply agreement, and agreement implies that reference cores are actual cores. The bound conditions on $`\mathcal F_{t+h}`$. |
| E12 | l.113 | $`b>\epsilon`$; $`b\log(b/\epsilon)-b+\epsilon\ge6`$; $`A\log A-A+1>5`$ | H | $`b\log b-b\ge-1`$; $`8\log8-7=9.64`$. |
| E13 | l.116 | Good reference histories; exceptional probability $`O_k(e^{-4s})`$; contribution $`O_k(e^{-s})`$ | H | A mark below $`1/M`$ at an original vertex of degree $`<M`$ raises an original degree; marks counted up to $`\epsilon s`$ are Poisson with mean $`k\epsilon s/M\ge6s`$; Cauchy–Schwarz with the Yule second moment. |
| E14 | l.118 | Comparison weight $`\ge e^{-\eta s}`$ on good histories assigned to the root; pigeonhole choice of root | H | Factors are $`\ge1/2`$ for at most $`bs`$ rings and $`\ge e^{-1/M}`$ afterwards. |
| E15 | l.121 | $`F(s)\le k+\max_{\lvert T\rvert=k}\mathbb E C_T(s)`$ | H | Strong Markov at the hitting time of size $`k`$, plus monotonicity. |
| E16 | l.121 | Path functionals $`X_u`$; disjointness; $`\sum_uX_u\le C(t+h+s)`$; summing bounds from site-specific couplings | H | The key point is correct: each $`X_u`$ is a functional of the tree path, so $`\mathbb E[X_u\mid\mathcal F_t]`$ may be bounded with any valid coupling, and the bounds add. |
| E17 | l.122–124 | Block inequality for $`F`$ | H | Factors $`1/4`$ (Lemma 2.1), $`p`$ (Lemma 2.2) and $`1/k`$ (pigeonhole). |
| E18 | l.126 | Promotion rate $`\ge1`$ after the first birth; $`F(t)\to\infty`$ | H | Checked on all trees with at most 17 vertices. |
| E19 | l.126 | Iteration gives $`\liminf\ge\limsup-\eta`$; letting $`\eta\downarrow0`$ gives a common limit in $`[0,1]`$ | H | The order of limits ($`s\to\infty`$ for fixed $`M,\epsilon,k`$, then $`\eta\to0`$) is stated and is right. |
| E20 | l.128 | Transfer to $`a_n=\mathbb E C_n`$ | H | $`a_n\le eF(\log n)`$; $`\Pr(N>n)\le2(1-e^{-t})^{\lfloor n/2\rfloor}`$; $`\mathbb E[N;N>n]=o(1)`$ via the second moment. |
| E21 | l.131 | At least $`C/6`$ sites with disjoint clock sets | H | $`\#\{c\ge3\}<\#\text{ends}\le A_0`$; at most $`2A_0-1`$ bare paths; $`C\le6A_0+3J-2`$. Exhaustive check up to 17 vertices: $`\ge3/10`$. |
| E22 | l.133 | Independent retained flags with a common probability $`p`$ from site-specific inputs, valid with a shared exterior | H+R | Correct but compressed. Repair, to be written out: (i) at sites with a leaf, thin stage 1 by an independent coin of probability $`(1-e^{-1})/(2e)`$, so both site types have stage-1 probability $`e^{-2}(1-e^{-1})/2`$. (ii) At the ring of $`u`$ in stage 3 (a stopping time), redefine the future clocks of the selected leaf and of its newborn from the site's own randomness. This is allowed because no other site reserves these clocks: the leaf is adjacent only to $`u`$. Then thin by $`e^{-2\tau}`$. (iii) On a retained flag the script succeeds whatever the exterior does, because exterior moves can only add leaves at $`u`$ (see E7) or at $`z`$. Adjacent sites and triples were checked case by case, and no interference is possible. |
| E23 | l.135–137 | The $`Y`$ variables are independent given the block-start past, dominated by the endpoint core count, and satisfy $`\mathbb E Y\ge(p/k)e^{-\eta s}(F(s)-k-O_k(e^{-s}))`$ | H | Each $`Y`$ uses only fresh site randomness; the branches are disjoint (E10). |
| E24 | l.139 | Truncation at $`R`$ and Chebyshev: $`\Pr(\cdot)\le16R^2/(m(\mathbb E Y)^2)\le96R^2/(C(\mathbb E Y)^2)`$ | H | Arithmetic re-derived. |
| E25 | l.141–145 | Infinitely many promotions, so $`C(t)\to\infty`$ a.s. | H | $`\prod_{n=n_0}^{n_1-1}(1-1/n)=(n_0-1)/(n_1-1)`$. |
| E26 | l.147–157 | Block indicators; $`L^2`$-bounded martingale; Kronecker; $`\liminf\log C(t)/t\ge\log g/(h+s)`$ | H | Failed blocks have nonnegative increments. |
| E27 | l.158 | Let $`s\to\infty`$ through the integers, then use countably many parameter choices with $`\eta\downarrow0`$; the case $`\mu=0`$ | H | |
| E28 | l.160 | Upper limit from Markov and Borel–Cantelli at integer times; time change $`\tau_n/\log n\to1`$ | H | |

**Verdict on Section 2.** The proof of Theorem 1.1 is correct. No step is false, and every gap closes with the repairs above. At AoAP/EJP standard, though, the section is a sketch. E2, E16 and E22 each deserve a formal statement and proof, and Lemma 2.2 should be stated with explicit filtrations: "conditional probability at least $`p`$ given $`\mathcal F_t`$ on the event that $`u`$ is eligible at $`t`$".

### 3.3 Section 3: motif bounds and the full motif drift lemma

| ID | Location | Step | Verdict | Notes |
|---|---|---|---|---|
| M1 | l.175–179 | Rooting at the centre of $`T_3`$; non-root core vertices have core parents; $`s_u=(L,k,h)`$ | H | A leaf neighbour of a core vertex is always its child. |
| M2 | l.180–195 | Rates $`r_u=L_u+\rho_u`$ and $`b_u=L_u/D_u`$; the two events are exhaustive; $`\rho_u\le c_u/2`$ | H | |
| M3 | l.199–258 | Projection $`\pi`$, measures $`\zeta`$, index sets, domains $`H,J`$, knot set $`\Xi`$ | H+R | Counts re-derived: 572 inner states, 601 indices, 555 parent indices, 290 child indices, 84 path centres, 75 cherry centres, $`2\cdot801\cdot151+2`$ vertex entries; $`\sum_a\zeta_s(a)i(a)=1/D`$. Repair: "$`\pi^{\rm end}`$" (l.236) is undefined and should read $`\zeta`$. |
| M4 | l.262–278 | Newborn terms; $`H`$ is entered only from $`k=0`$ and $`J`$ only from $`k=1`$; exits via the zero extension | H | |
| M5 | l.280–295 | Vertex redistribution $`\theta_u`$ and $`U_e`$ | H | Inner vertices keep $`h\,r_{\rm par}+k\,r_{\rm ch}`$; the differences go to edges. |
| M6 | l.297–316 | Splitting the credits $`i_zd`$ by orientation (centring value inside the domain, knot function outside) | H | Each triple $`(x,y,z)`$ is split exactly once, and each sent $`A(i_z)`$ is received by the edge $`xz`$. |
| M7 | l.317–354 | Multiplicities $`n_1,n_2,n_3`$; $`\sigma_p,\sigma_w`$; $`I_e`$, $`R_e`$, $`E_e`$ | H | Sender/receiver bookkeeping re-derived for all three orientations, including the root ($`h=0`$, so $`n_1=0`$ and the root is not in $`H`$). |
| M8 | l.355–367 | $`P_{puw}`$ and $`B_{u12}`$ | H | Each path receives one residual from each of its edges; each cherry receives one per ordered child. |
| M9 | l.368–387 | **Lemma 3.1** (exact decomposition) | H | Follows from M4–M8. Only one old state changes per event, plus at most one newborn. The supplement's term-by-term identity test, rerun by the referee, passes on 177 trees, including trees with tail and overflow states. |
| M10 | l.391–395 | Rate-interval bound $`rd\preceq M_I(d)`$ | H | |
| M11 | l.397–471 | **Lemma 3.2**, rules (1)–(4) | H | Rule (1) is multilinear interpolation. Rule (2): the extrema of a piecewise-linear $`z\mapsto zd-A(z)`$ lie at knots, and Jensen applies to the convex max (concave min). Rule (3): the gain moves $`t\mapsto t-t/(D+1)`$, and $`r_v/(D+1)\in[(D_0-\kappa-h)/(D_0+1),1]`$, with $`0`$ when clamped; budding gives $`b_v\in[(D_0-\kappa-h)/D_0,1]`$. Rule (4): $`i_v\Delta_xF=(t/D_0)((1-t)d_0+td_1)`$. |
| M12 | l.473–503 | Five vertex families | H | They match the checker exactly. The count reconciles: $`8+572+2\cdot241{,}024+5\cdot302=484{,}138`$. |
| M13 | l.505–538 | Edge corner expressions | H | $`555\cdot290-27\cdot15=160{,}545`$. A tail–tail edge gives $`0=0`$, since no newborn motif can be centred in $`H`$ or $`J`$. |
| M14 | l.540–586 | Path and cherry corner expressions | H | The products of a credit and an edge jump are separately affine. The rate tables (top, bottom, middle, child, apex) were re-derived. Counts: $`555\cdot84\cdot290`$ and $`75\cdot290\cdot291/2`$. |
| M15 | l.575–586 | Assembly: $`X(T)\preceq\sum_c\omega_c\widehat X(c)\preceq\lambda\tau(T)`$, hence $`\mathcal L\Phi\preceq\lambda\Phi`$ | H | Confirmed directly on 135 large trees per direction (§2). |
| M16 | l.590–618 | Comparison of $`\Phi`$ with $`\lvert K\rvert`$ through $`g_m,v_m,c_m`$ and a finite enumeration | H | Including $`0`$ in each extremum is conservative. Using both tail indices also at $`t=1`$ is valid. |
| M17 | l.620–671 | **Lemma 3.3**: the finite verification and its constants | H+R | Reproduced exactly. Repair: the lemma refers to "rational tables" without saying which. Print the two SHA-256 digests (they are in the supplement) and the build command, and state that a rebuild with a different digest needs a fresh verifier run, because the builder uses floating point. |
| M18 | l.653–658 | Clearing denominators by $`Q\,10^{10}\cdot2D_1D_uD_2(D_0+1)\,10^6`$ | H | The code uses $`D_{\rm eff}=1`$ at $`\infty`$ and $`D_0`$ at $`q`$, as stated. |
| M19 | l.673–693 | Iteration and $`\Phi_\pm(T_3)=f_\pm(2,0,0)`$ give the two product bounds | H | The lower direction needs no positivity of $`\Phi_-`$. |
| M20 | l.695–701 | Almost-sure upper bound $`1133/2000`$ | H | Dyadic Borel–Cantelli and monotonicity of $`C_n`$. |

**Verdict on Section 3.** The motif drift lemma, Lemmas 3.1–3.3 together with the four inequality families, is correct as a reduction. The reduction is faithfully implemented by `cert_verify_gmp.c`: the referee checked each block of the vertex, edge, path, cherry and comparison code against the text. It is confirmed by an independent implementation of $`\Phi`$ on trees containing every tail and overflow case.

### 3.4 Section 4: the banded upper bound $`5657/9983`$

| ID | Location | Step | Verdict | Notes |
|---|---|---|---|---|
| B1 | l.705–713 | $`\lambda=5657/9983=17/30-1/299490`$ | H | |
| B2 | l.717–729 | Gain and conversion rates; the newborn has type $`(2,1)`$; the rates at $`T_3`$ are $`2`$ and $`1`$ | H | |
| B3 | l.731–770 | Alphabet: band caps $`h_c(D)`$, $`K_c`$, witnesses, pooled marks, the largest-mark rule, ordering | H | Caps are nonincreasing in $`D`$ and in $`c`$; $`K_c`$ is nonincreasing; $`h_1=620`$ in every band. The per-band state counts in the table equal the verifier's counts. |
| B4 | l.772–813 | **Lemma 4.1** (history invariant) | H | The referee implemented the update rules from this text alone and asserted the invariant after every step of 16 histories, covering band crossings, pooling, saturation and the creation, update and dropping of witnesses. The newborn witness choice (first entry of the creator's old multiset) is valid because the caps are at least 6 and $`K_c\ge5`$ for $`c\le3`$. |
| B5 | l.817–835 | Generator identity summed over $`I_u=\{u\}\cup S_u\cup J_u`$ | H | Events outside $`I_u`$ cannot change $`z_u`$; the sets are disjoint by the tree property. |
| B6 | l.837–849 | Signed transfer $`g`$: exact cancellation; readable from both ends | H | The checker enforces the four observability conditions on every key. |
| B7 | l.851–866 | Nonnegative transfer $`a`$ and the amended row $`R_u`$ | H | |
| B8 | l.870–881 | $`Q_z`$ bound; dependent rates may be bounded termwise | H | |
| B9 | l.883–895 | Central targets and rate intervals, including pooled marks with $`b<T`$ | H | Using no subtraction when $`b<T`$ is conservative because $`a\ge0`$. |
| B10 | l.897–925 | Neighbour and witness rate intervals, including saturation | H | Each interval re-derived and matched to `check.cpp`, lines 304–321. |
| B11 | l.926–940 | Sender charge $`F_u`$; row inequality (4.1) | H | |
| B12 | l.942–964 | Tail $`D\ge H`$: affine bound $`(B+a_*)(1-q)+sq`$ and the two endpoint inequalities | H | $`r_u(1-q)/(D+1)-q/D\le1-q`$, using $`r_u\le D`$. Slacks verified. |
| B13 | l.966–976 | Recipe file: counts and digest | H | 12,697 rewards, 63,309 signed entries, 18 nonnegative entries; the digest matches. |
| B14 | l.978–1038 | **Lemma 4.2**: order of computation, superset enumeration, integer acceptance test, results | H+R | Reproduced exactly. The order claim is enforced at run time, since the checker aborts on any uncomputed successor. Repair: the weights depend on a floating-point *proposal* (`check.cpp`, lines 324–337). Any accepted integer is valid, but the constants $`m`$ and $`z_{(2,0)}`$ in Theorem 1.2 belong to this particular admissible system. The integer-only "least admissible" rule of the second verifier gives $`1539738113358/3908512428`$. Say this explicitly, or state the theorem with the integer-only rule. |
| B15 | l.1040–1060 | Proposition (unlabelled): the product bound | H | |
| B16 | l.1062–1071 | Almost-sure upper bound $`5657/9983`$ | H | |

**Verdict on Section 4.** Correct as a computer-assisted proof. The text is specific enough that an independent label tracker, written from Lemma 4.1 alone, agrees state by state with the production verifier's state space on long histories. The drift inequality holds on all 2,712 history times tested, with a smallest relative margin of $`3\times10^{-6}`$.

### 3.5 Section 5: computer-assisted verification

| ID | Location | Claim | Verdict | Notes |
|---|---|---|---|---|
| V1 | l.1075–1080 | Two-step structure (identity, then arithmetic) | H | |
| V2 | l.1082–1095 | Check counts; 128-bit integers; UBSan | H | Reproduced. |
| V3 | l.1097–1110 | Direct tests: 32,972 rooted shapes; means for $`3\le n\le21`$; 23,116 histories | H+R | The counts are right: $`32{,}972=A(14)-A(2)`$ for the rooted-tree numbers $`A`$, and $`23{,}116=1+3+12+60+360+2520+20160`$. The 32,972-tree test was rerun and passes. Repair: say that these sizes contain no tail state, no overflow state and no band crossing, and add tests that do (the referee's tests in §2 show that such tests are feasible). |
| V4 | l.1112–1116 | Negative controls | H | 10 of 10 rejected in the referee's runs. |
| V5 | l.1117–1126 | Commands and requirements | H+R | `run_all.sh` installs pinned Python packages, which needs network access. The two headline suites need only a C/C++ compiler and GMP; document how to run them offline. |

---

## 4. Which finite statements imply each theorem, and whether the paper's description suffices

**Theorem 1.1 (existence, mean and almost sure).** No computation is involved. The proof rests on Lemma 2.1, Lemma 2.2, inequality (2.1), the good-history estimate (E13–E14) and the block inequality (E17) for the mean. For the sample path it adds the site selection (E21), the independent flags (E22), Chebyshev (E24) and the martingale argument (E26). The description suffices once the repairs to E2, E7 and E22 are written out.

**Theorem 1.2, motif pair.** For each direction $`(\lambda,\preceq)\in\{(1133/2000,\le),(113/200,\ge)\}`$, the proof needs exactly these finite facts:

1. $`g[a,b]=0`$ for the $`27\times15`$ tail–tail index pairs;
2. the vertex inequalities of the five families (484,138);
3. the edge inequalities for $`(a,b)\in\mathcal P\times\mathcal W`$ except tail–tail (160,545);
4. the path inequalities on $`\mathcal P\times H\times\mathcal W`$ (13,519,800);
5. the cherry inequalities on $`J\times\binom{\mathcal W}{2}_{\rm rep}`$ (3,164,625);
6. the extremum of the per-vertex lower (upper) bound for $`W_u`$ over $`0\le L\le800`$, $`0\le k\le150`$, $`h\in\{0,1\}`$ and the two overflow states. Its values are $`2989767809/10^{10}`$ and $`17527873930/10^{10}`$;
7. the values $`f_\pm(2,0,0)`$.

Lemma 3.1 and Lemma 3.2 turn items 1–5 into $`\mathcal L\Phi\preceq\lambda\Phi`$ on every tree; items 6–7 give the prefactors. *Sufficiency of the description:* yes, for a reader who has the arrays. The referee implemented $`\Phi`$ from the text alone, read the shipped arrays, and recovered the drift inequality on large trees. The paper must, however, identify the arrays (Issue 2) and fix the notation slips (M3).

**Theorem 1.2, banded bound.** The proof needs:

1. the 226,545,617 row inequalities (4.1) on the enumerated superset of finite states, with targets defined by the recursion;
2. the newborn closure $`\min_t\{s(t)-z(2,1;t)\}\ge0`$ over the 12,697 newborn labels, including $`t=H`$;
3. the two tail inequalities;
4. the minimum weight $`m`$ (finite and tail) and the centre weight $`z_{(2,0)}`$.

Lemma 4.1, the generator identity and the two exact transfers turn items 1–3 into $`R_u\le\lambda z_u`$ for every core vertex; item 4 gives the prefactor. *Sufficiency:* yes. The geometry, caps, constants and the hashed recipe determine everything. The one weak point is the weight-selection rule (B14). A second verifier, written from the specification, accepts the same input with different weights, and the referee's own label tracker agrees with the verifier's state space.

**The claim "below 17/30".** It follows from item set 1–5 of the motif upper bound alone, or from the banded rows alone: either gives $`\limsup_n\log\mathbb E C_n/\log n<17/30`$, and Borel–Cantelli gives the almost-sure $`\limsup`$. Theorem 1.1 is needed only to speak of a single exponent and to make the lower bound $`113/200`$ hold almost surely.

---

## 5. Bibliography and attributed claims

Each cited source was retrieved and read at the stated location, except where noted.

| Ref. | Bibliographic data | Attributed claim | Where in the source | Verdict |
|---|---|---|---|---|
| [SK] | Physica A **341** (2004) 80–86; arXiv:cond-mat/0404088 (single version). Correct. | The rule belongs to the random-walk attachment family of SK. | Model section | Correct. SK grow networks by random walkers, with walk length as a parameter. |
| [CJ] | Electron. Commun. Probab. **18** (2013), paper 77, 1–5, doi:10.1214/ECP.v18-2518. Correct. | For walks of length one the fraction of leaves tends to one. | **Theorem 2.1** (almost-sure convergence) | Correct; give the theorem number. |
| [KR] | J. Stat. Mech. (2017) 073405; arXiv:1706.01514 v2. Correct (journal DOI 10.1088/1742-5468/aa7a3f not given). | "Isotropic redirection"; nucleus exponent 0.566; fitted slopes 0.5652–0.5673; scaling law for $`N_k`$ (exact degree). | §5.1, eqs. (32)–(33); intro eq. (1) | Correct. |
| [Strogatz] | Post of 20 September 2022 (created 18:12 UTC). The X page needs a login; the text was read through the public syndication endpoint. | "Strogatz asked whether 17/30 could be the exponent"; "in that question … $`N_k`$ counts vertices of degree exactly $`k`$". | Post text and attached image | **Inaccurate.** The post reports a "weird exponent of 0.566" found by Redner and Krapivsky and says that "It's an open problem to prove that 17/30 is really the exponent." So it presents 17/30 as the conjectured value, not as a question. The post's text does not mention $`N_k`$; its attached photograph shows a screen with a log–log plot labelled $`N_2,\dots,N_{100}`$ in the style of [KR, Fig. 8]. Quote the post accurately, and cite an archival copy if the journal allows social-media sources. |
| [AB] | arXiv:2403.20185 v1 (29 Mar 2024). **Now published:** *Random Friend Trees*, Random Structures & Algorithms **69** (2026), e70082, doi:10.1002/rsa.70082. | Almost-sure bounds $`(\sqrt{13}-3)/6\approx0.1009`$ and $`\sqrt2-1/2\approx0.9142`$ (Props. 6.20, 6.22); comparability of $`X^{\ge k+1}`$ and $`X^{\ge k}`$. | Props. 6.20 and 6.22 in both versions; the comparison theorem is **Thm 4.8 in arXiv v1 but Thm 4.6 in the published version** | Content correct. Cite the published version and give its numbering (l.31, l.1188). AB's own Theorem 4.7 (4.5 in RSA) states $`0.1<\delta<\lambda<0.9`$, which conflicts with Prop. 6.22. The manuscript rightly cites the proposition. |
| [HKRZ] | Phys. Rev. E **114** (2026) 034314; arXiv:2604.01540 v1. Correct (DOI 10.1103/9lmp-pydp not given). | Modified rules with tractable leaf-degree recurrences; 0.569 by nonlinear least squares over $`N=10^2`$–$`10^6`$ with 2,400 realisations; weighting unspecified. | Fig. 6 and caption | Correct. |
| [supplement] | GitHub revision `c4653036b993` | Proofs and verifiers | — | That revision is an ancestor of the reviewed snapshot. `proofs/` and `verify/` are identical between the two; only READMEs, the paper and the review files differ. |

**Relevant omissions.**
- T. S. Evans and J. Saramäki, *Scale-free networks from self-organization*, Phys. Rev. E **72** (2005) 026138. It is the companion random-walk model, and its preferential-attachment claim is the one [CJ] corrected (see [AB, §1]).
- P. L. Krapivsky and S. Redner, *Organization of growing random networks*, Phys. Rev. E **63** (2001) 066123. This is the origin of the redirection mechanism.
- [AB, §7, Open Question (1)] asks whether the gap between the polynomial bounds can be closed and whether $`n^{-\mu}X_n^{\ge2}`$ converges. The manuscript answers the first part on the logarithmic scale and leaves the second open; it should say so.
- The argument in §2.3 is an approximate super-multiplicativity (Fekete-type) argument; a pointer to the standard lemmas would help readers.

---

## 6. Exposition

1. **Density.** Section 2 condenses what would normally be 10–15 pages into four. The informal phrases "in every initial tree", "the marginal continuation preserves the tree law" and "by thinning conditional successes in each stage" each stand for a lemma.
2. **Overloaded notation.** Many symbols carry two or more meanings:
   - $`M`$: degree threshold in §2 and max/min operator in §3;
   - $`H`$: path domain in §3 and $`2400`$ in §4; $`h`$: block preparation length, parent indicator and band caps $`h_c`$;
   - $`A`$: $`8`$ in §2, the knot functions $`A_{pc}`$ and the tail constant in §4;
   - $`B`$: $`B(u)`$, $`B_s`$, $`B_{u12}`$ and the tail constant;
   - $`J`$: the cherry domain and the number of triples; $`T`$: trees and the threshold $`T=8`$; $`p`$: probability, parent vertex, and the numerator in the code;
   - one event is called "promotion", "budding" and "conversion".

   Please use distinct symbols per section.
3. **Section 3.** Define $`\zeta`$ consistently (fix $`\pi^{\rm end}`$). One fully worked edge inequality and a figure of the four motif types would help readers. Explain why the knots and the domains $`H,J`$ were chosen.
4. **Section 4.** A table of label types and their transitions, plus one short worked history, would make Lemma 4.1 readable.
5. **Computer-assisted lemmas.** Give a short paragraph on what is trusted (compiler, GMP, the 128-bit narrowing checks) and list the digests of every input that determines a constant.
6. **Section 6.** Outside the close-read scope, but it states theorem-like results whose proofs exist only in the supplement. A journal will want them proved in the paper, moved to a separate paper, or turned into remarks.
7. **"How this work was produced."** The disclosure is welcome. The sentence "No human referee has yet read the proofs" belongs in a cover letter rather than in a submitted paper. The editor should check the paper against the journal's policy on AI-assisted research and author responsibility.
8. **Finite-size bounds.** The prefactors are loose. At $`n=21`$, $`\mathbb E C_n/P_{1133/2000}(n)=1.245`$ against the certified $`4.886`$, and the banded prefactor is $`393.9`$. One sentence saying the constants are not optimised would avoid misreading.

## 7. Minor corrections

1. l.29: describe the Strogatz post accurately (see §5).
2. l.31, l.1188, l.1321–1323: cite the published [AB] and use its theorem numbering (4.8 → 4.6).
3. l.52: "together with Theorem 1.1" before "the first pair gives $`0.5650\le\mu\le0.5665`$".
4. l.86–97: add the inequality of E7. Write the source-degree bound in stage 4 as $`k-1`$.
5. l.133: write out the thinning and clock-redefinition details of E22.
6. l.236: $`\pi^{\rm end}`$ should read $`\zeta`$.
7. l.391–395: rename the extremum operator $`M`$.
8. l.620–671: put the digests and the build command into Lemma 3.3.
9. l.1032–1037: explain the dependence of the constants on the weight-selection rule (B14).
10. l.1040: give the Proposition a number and label.
11. l.1102–1110: state the regimes that the direct tests do not reach.
12. Bibliography: add DOIs for [SK], [KR] and [HKRZ], and the theorem number for [CJ].

---

## 8. Recommendation and top three issues

**Recommendation: major revision.** The referee found no false step and no step that cannot be justified. The three headline claims ($`\mu<17/30`$, existence of $`\mu`$, and the explicit bracket) are, in the referee's judgement, correct. The finite checks reproduce exactly, and independent direct tests in the regimes the paper never tests show no failures. The revision is about making the proofs verifiable from the paper, as AoAP/EJP require.

**Top three issues.**

1. **Section 2 must become a full proof.** Write the coupling principle (E2), the summation over site-specific couplings (E16) and the construction of independent site flags (E22) as stated lemmas with proofs. Add the missing inequality in Lemma 2.2 (E7). Give Lemma 2.2 explicit conditioning. Every step was found correct, but several are asserted rather than shown.
2. **Pin down the certified objects and test them where they matter.**
   - Print the SHA-256 digests of the two motif arrays in Lemma 3.3. They exist only in the supplement, and the arrays come from a floating-point builder.
   - State that the banded constants depend on the weight-selection rule, which includes a floating-point proposal, or switch to the integer-only rule.
   - Extend the reported direct tests beyond 14 vertices and size 9, so that they reach the tail interpolation, the overflow states and band crossings.
3. **Correct the literature context.**
   - The Strogatz post states 17/30 as the conjectured exponent, not as a question, and does not itself define $`N_k`$.
   - [AB] is published in Random Structures & Algorithms, where the cited Theorem 4.8 is Theorem 4.6.
   - Add Evans–Saramäki (2005), Krapivsky–Redner (2001) and AB's Open Question (1).
   - Decide the status of the supplement-only results in Section 6.

---

### Appendix: referee scripts (not part of the supplement)

The referee used the following independent audit scripts. The packaged public checkers are described in [the verification guide](../verify/README.md).

- `exact_mean.py`: exact law of the chain on unlabelled trees, with checks of Lemma 2.1, the site count and the promotion rate.
- `motif_direct.py`: the motif potential implemented from §3.2, and a brute-force drift test with local and global cross-validation.
- `banded_direct.py`: banded labels from Lemma 4.1, invariant assertions, and brute-force drift using weights from the production verifier's read-only observer (`verify/tree-drift/read_weights.cpp`).
- `pdftext.py`: a minimal PDF text extractor used to read [AB], [KR] and [CJ].
