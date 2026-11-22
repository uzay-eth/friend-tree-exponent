CLAIM 1 VERDICT: HOLDS
CLAIM 2 VERDICT: HOLDS

# Adversarial referee report on strogatz-astra `existence.md` (Claim 1) and `sample-exponent.md` (Claim 2)

Referee: Claude Opus 5.5 (`claude-opus-5-5`, claude-code harness), 2026-11-22.

Scope and ground rules I followed:
- I judged only the text of `[astra-branch]/attempts/strogatz-astra/existence.md` and `.../sample-exponent.md`. I read nothing else in `[astra-branch]` (I listed the directory names once and opened none of the files) and wrote nothing there.
- All scripts are mine and live in `checks/attack-astra-opus/`. Every run was pinned to core 0 (`taskset -c 0`), so total CPU use stayed at one core. Interpreter: `PY=attempts/strogatz-fable/.venv/bin/python` (numpy 2.1.3, scipy 1.14.1).
- The outputs pasted in section 3 are copied verbatim from the `out_*.txt` files. `build_report.py` pasted them into this report mechanically.

## 0. Bottom line

I could not break either proof. Every sentence I checked is correct as written, or needs only a one-line rewording (listed in section 4). None of these rewordings changes a constant or the logic. Every numerical test of a sub-lemma passes, and every deliberately planted error was caught.

These are the places I pushed hardest, and why each survives:

1. **Eligible >= C/4.** The argument is right and actually gives E >= (C+5)/4 for C >= 2. The constant 1/4 is sharp asymptotically: edge-subdivided cubic caterpillars give E = (C+5)/4 exactly. No tree with up to 17 vertices violates the bound (81 136 trees).
2. **Planting script.** No outside clock can undo a stage, for four reasons:
   - a leaf of u can be promoted only by u's own ring;
   - c_u can grow only through u's own ring;
   - deg(w) can change only through rings of u or of w's other neighbour x;
   - the branch touches the rest of the tree only through u, which is silent from its single stage-3 ring until t+h.

   The only probability that depends on the exterior is the stage-3 leaf choice, which is >= (M+1)/(M+3) >= M/(M+2), and outside arrivals only raise it. A forced-script simulation in a hostile exterior confirms all the conclusions in 24 000 of 24 000 runs.
3. **Coupling.** The marks reproduce the exact law of the actual process. E_u is a function of u's own clock and marks, so it is independent of the reference. On E_u ∩ Acc the internal trees coincide. All three facts hold numerically, and planted errors in the ordering or the acceptance probability are caught.
4. **Degree-M condition in the good histories.** It is correct, but only asymptotically. At eps = 0.01, M = 50, k = 6M/eps = 30000, the claimed bound P(Pois(6s) <= kM) <= e^{-4s} fails at s = 20 and s = 50: the probability there is 1 to double precision. It starts to hold only at s0 ≈ 2.47·10^6 ≈ 9.9 M²/eps. The text says "for large s" and fixes all parameters before s → ∞, so this is not a gap. Anyone quoting the constants should know about it, though.
5. **Claim 2 independence.** The promoted leaf's identity does depend on outside arrivals. The proof's stage-3 flag does not depend on that identity: it is (u's ring count, u's own mark below M/(M+2) with leaves ranked first, and silence of the freshly re-clocked v and v'). Numerically, with two adjacent sites, the proof's flags are independent (z = +1.09 and -0.87 at 4·10^6 runs). The "natural" flags ("u's ring picks a leaf") are clearly dependent (z = -21.9). So the mark device is necessary, and the proof uses it.

---

## 1. `existence.md`, sentence by sentence

Notation used below: n1 and n3 are the numbers of core vertices with core degree 1 and >= 3, and E is the number of eligible vertices. S_u is the planting event at u, and X_u is the number of core vertices in the branch at u at time t+h+s.

### 1.0 Set-up
| # | Sentence (abridged) | Verdict |
|---|---|---|
| 0.1 | Continuous-time embedding with rate-one clocks; ring → uniform neighbour → new leaf there. | VALID. The total rate is N, and the jump chain picks V uniformly and then W uniformly among V's neighbours, which is the discrete model. |
| 0.2 | N_T is a Yule process (rate j at size j). | VALID. |
| 0.3 | E[N_T(t)^2] = k(k+1)e^{2t} - k e^t. | VALID. N_T - k ~ NegBin(k, e^{-t}). Checked to 6 decimals by exact summation (`poisson_tails.py`, last block). |

### 1.1 Finding places
| # | Sentence | Verdict |
|---|---|---|
| 1.1 | Core vertices induce a tree when C > 0. | VALID. Deleting all leaves of a tree leaves a subtree; for C = 1 it is a point. |
| 1.2 | Definitions: c, L, type B (c=2, L=0), eligible = (c<=2 and L>=1) or (B with a B neighbour). | – |
| 1.3 | "There are at least C/4 eligible vertices." | VALID; see the chain below. |
| 1.4 | "For C=1 the centre is eligible." | VALID: c = 0 and L >= 2. For C = 0 (the initial edge) the bound reads 0 >= 0. |
| 1.5 | B vertices with no B neighbour form an independent set; their 2I core edges are distinct; I <= (C-1)/2. | VALID. The core tree has C-1 edges. |
| 1.6 | n3 <= n1 by summing c-2. | VALID. In fact Σ(c-2) = -2 gives n3 <= n1 - 2. |
| 1.7 | Each core-degree-one vertex has a leaf and is eligible. | VALID. Its degree is >= 2 with only one core neighbour, so L >= 1 and c = 1 <= 2. |
| 1.8 | Every ineligible vertex is among the I or has c >= 3. | VALID. c = 0 is impossible for C >= 2, c = 1 is eligible, c = 2 with L >= 1 is eligible, and c = 2 with L = 0 is type B. |
| 1.9 | "The claimed fraction follows." | VALID. C <= E + (C-1)/2 + (E-2), hence E >= (C+5)/4 >= C/4. |

The count is tight in two senses:
- **Smallest cases.** Among the trees with up to 11 vertices, 6 trees with 11 vertices meet E = (C+5)/4 exactly; they are listed in output (i-a). The smallest ratio E/C is 3/7, reached at C = 7 by the subdivided K_{1,3} with a leaf on each end.
- **Asymptotically.** `families.py` builds edge-subdivided cubic caterpillars with a leaf hung on each skeleton leaf. They give E = (C+5)/4 exactly, so E/C → 1/4 (0.2516 at C = 803). The constant 1/4 cannot be improved.

### 1.2 The planting script
For every requirement below, the probability counts only clocks and marks that the script names. Their increments on the stated interval are independent of the past, by the Markov property at deterministic times or the strong Markov property at the ring times.

| # | Requirement | Probability | Can the exterior destroy it? |
|---|---|---|---|
| S1a | (L_u >= 1) u silent on unit 1 | e^{-1} | No. A leaf of u can be promoted only by u, and c_u can grow only when u's own ring promotes a leaf. VALID. |
| S1b | (u of type B with B neighbour w; N(w) = {u, x}) u and x silent, and the first ring of w falls in unit 1 and chooses u | e^{-2}(1-e^{-1})/2 | No. Only rings of u or x can add a vertex to w, so deg w = 2 at its first ring. "Before that ring w has degree two" is VALID; the forced runs show deg w = 2 in 12 000 of 12 000 runs. A planted error (x not silent) gives deg w up to 7. |
| S1 end | "u now has a leaf, and still has at most two core neighbours" | – | VALID. |
| S2 | u silent on unit 2, and a leaf ℓ of u (chosen by a fixed rule at the start of unit 2) rings >= M times | e^{-1} P(Pois(1) >= M) | No. ℓ's only neighbour is u, so each ring adds a leaf to u, and nothing else can touch ℓ. "Stays a leaf" is VALID. |
| S3 | exactly one ring of u in unit 3; it picks a leaf; v and v' then silent until the unit ends | >= e^{-1}·M/(M+2)·e^{-2} | No. At the ring L_u >= M+1 and c_u <= 2, so P(leaf) >= (M+1)/(M+3) >= M/(M+2). Arrivals only increase L_u. In the forced runs the smallest L/D seen equals (M+1)/(M+3) exactly. VALID. |
| S4 | u silent; per added vertex: the prescribed internal source rings once and picks the prescribed target, all other internal clocks are silent, and the newborn is silent until the unit ends | >= e^{-1}·e^{-1}·(1/k)·e^{-(k-2)}·e^{-1} per unit | No. During S4 the branch's only outside neighbour is u, which is silent. The source's degree is <= k; the text's "k+1" is a harmless overestimate. Parent-before-child order always has an internal neighbour available, because the partial branch is connected and has >= 2 vertices. k = 2 needs no additions. VALID. |

Consequences:
- **h and p.** h = k+1 and p(M,k) >= [e^{-2}(1-e^{-1})/2]·[e^{-1}P(Pois(1)>=M)]·[e^{-3}M/(M+2)]·[e^{-(k+1)}/k]^{k-2} > 0, uniformly in the tree. This gives the lemma as stated (VALID).
- **D_u >= M at t+h and forever after.** VALID. After S3, L_u >= M, plus the edge to v, so D_u >= M+1, and degrees never decrease.
- **The branch is exactly (T,v) at t+h.** VALID. The forced runs show 0 isomorphism failures in 24 000 runs; a planted error (u not silent in S4) gives 4–17 failures per 300 runs for every seed with k >= 3 (K2 has no stage 4, so it shows 0).
- **"Extra arrivals at u ... can only help stages 1–3".** VALID.
- **Disjointness of branches.** VALID. At any later time, the far side of u_i v_i consists of v_i, vertices that were leaves at t, and vertices born later; the tree grows only by leaves. So it never contains an old core vertex. Suppose v_j were on the far side of u_i v_i. Then u_j, which is old core, differs from u_i and from v_i, and is joined to v_j by an edge other than u_i v_i, would be on that far side as well. That is a contradiction.
- **No independence needed.** VALID. Σ_u 1_{S_u} X_u <= C(t+h+s) holds pointwise for whatever subset of the events S_u occurs. Scripts at neighbouring vertices may be mutually exclusive (for example, u_1 may be the x of u_2's script), and that costs nothing under linearity.

### 1.3 Comparing a planted branch with a free tree
| # | Sentence | Verdict |
|---|---|---|
| 3.1 | Inside the branch, couple clocks, choices and newborn identities up to the first discrepancy. | VALID. |
| 3.2 | At v: keep the reference's internal choice with probability D_v/(D_v+1), otherwise choose u; acceptance marks are independent. | VALID. In the actual tree, each internal neighbour gets (D/(D+1))(1/D) = 1/(D+1) and u gets 1/(D+1), which is the exact law. Test (P3): the actual mean from the coupled construction matches a direct simulation (abs(z) <= 1.28 in 5 configs). A planted acceptance probability of (D+1)/(D+2) is caught (z up to 16.4). |
| 3.3 | At u: put v first; with no mark in [0,1/M], u never chooses v because D_u >= M. | VALID: {v chosen} = {U <= 1/D_u} ⊂ {U <= 1/M}. |
| 3.4 | "This event has probability e^{-s/M} and is independent of the reference history." | VALID. The count of u-rings with mark <= 1/M is Poisson(s/M), and it is a function of u's clock and marks, which are not inputs of the reference. Test (P2): P(E_u) matches e^{-s/M}, and E[C_ref 1(E_u∩Acc)] matches e^{-s/M}E[C_ref W] (abs(z) <= 0.92). |
| 3.5 | On E_u and acceptance at every ring of v, the internal trees agree throughout. | VALID. Only u can push a vertex into the branch from outside, and every internal ring has the same effect in both processes. Test (P1): 0 violations in 200 075 runs where E_u∩Acc held. A planted error (v last in u's ordering) produces 6 144–11 176 violations per configuration. |
| 3.6 | Reference core ⇒ actual core. | VALID. Degrees agree except at v, which has one more in the actual tree. |
| 3.7 | Weight e^{-s/M}·Π D_v/(D_v+1). | VALID, read as E[X_u \| F_{t+h}] >= e^{-s/M} E[C_T(s) Π_v D_v/(D_v+1)] on S_u. Given the reference history, the acceptance probability is exactly that product. |

Direct Monte Carlo of the inequality (`comparison_mc.py`, 42 configurations: 7 rooted seeds, M ∈ {3, 6, ∞}, s ∈ {1, 2}, 10^5 runs per side):
- **The claimed inequality holds everywhere**, very slackly; the smallest z is 47.9.
- **Planted weight e^{+s/M}·Π (D+1)/D:** rejected in 42 of 42 configurations (z <= -54).
- **Planted weight e^{+s/M}·Π D/(D+1) (sign error only):** rejected in 9 of 28 finite-M configurations, namely the high-degree roots.
- **Direction of the effect without the weight.** When v has degree >= 2 in T, E[actual] < E[C_ref]: for example 3.924 vs 4.220 for K_{1,5} rooted at the centre with M = 3, s = 2, and 3.395 vs 3.868 for P3 rooted at the middle. The same holds for K2 at s = 2. When v is a leaf of T, E[actual] > E[C_ref], because v counts as an actual core vertex but not as a reference one. So an unweighted comparison would be false, and the weight is genuinely needed.

### 1.4 A high-degree vertex among the starting vertices
| # | Sentence | Verdict |
|---|---|---|
| 4.1 | A log A - A + 1 > 5 for A = 8. | VALID (9.6355). |
| 4.2 | b > eps with b log(b/eps) - b + eps > 5, b → 0; b = 7/log(1/eps) works for small eps. | VALID. In fact it works for every eps in (0,1), since the rate is 7 + b log b - b + eps >= 6 + eps; see the grid in output (iii). |
| 4.3 | k >= 6M/eps; everything fixed before s → ∞. | VALID. |
| 4.4 | Bad set has probability O(e^{-4s}) for large s: ring counts via the Chernoff bound. | VALID. At s = 20 and 50 the exact log-tails are -128.18 and -314.43 (the Chernoff bounds are -122.53 and -306.32) for (eps, b), and -196.04 and -485.56 for A = 8. All are <= -5s, and the union over k clocks costs only a factor k. |
| 4.5 | Degree condition: until some original degree reaches M, every original mark in [0,1/M] raises an original degree; at most kM such raises; the mark count is Poisson(k eps s/M) with mean >= 6s; its lower tail at kM is <= e^{-4s} for large s. | VALID but asymptotic only. The logic is right: while all original degrees are below M, 1/M < 1/D_y, so the fixed original neighbour is chosen; there are at most k(M-1) - 2(k-1) < kM such raises. The threshold, however, is enormous. At eps = 0.01, M = 50, k = 30000 we have log P(Pois(120) <= 1.5·10^6) = 0 at s = 20 and log P(Pois(300) <= 1.5·10^6) = 0 at s = 50, so the bound fails there. It first holds at s0 ≈ 2.467·10^6. Because the text says "for large s", this is not a gap. |
| 4.6 | Separate marks for the internal choice and for acceptance. | VALID. |
| 4.7 | Cauchy–Schwarz: the bad part contributes O(e^{-s}). | VALID: E[N^2]^{1/2} P(bad)^{1/2} <= (k(k+1))^{1/2} e^{s}·O(e^{-2s}). |
| 4.8 | For histories assigned to v the weight is >= e^{-ηs}, η = b log 2 + (A+1)/M. | VALID. There are <= bs rings of v before eps·s, each with factor >= 1/2 since D_v >= 1. After that there are <= As rings, each with factor >= M/(M+1) >= e^{-1/M} since D_v >= M. The u-factor is e^{-s/M}. |
| 4.9 | Pigeonhole: some v accounts for >= (E C_T(s) - O(e^{-s}))/k. | VALID. |
| 4.10 | "The choice of T,v is made from the laws, before planting." | VALID. p and h are uniform over all rooted k-vertex trees. |

### 1.5 Iteration and return to discrete size
| # | Sentence | Verdict |
|---|---|---|
| 5.1 | If size k is not reached by s, then C(s) <= k. | VALID. |
| 5.2 | F(s) <= k + max_T E C_T(s). | VALID. On {τ_k <= s}, C(s) <= C(τ_k+s) by monotonicity; then apply the strong Markov property at τ_k. |
| 5.3 | F(t+h+s) >= F(t)[p/(4k)] e^{-ηs}[F(s)-k-O(e^{-s})]. | VALID. E[C(t+h+s) \| F_t] >= p·#elig(t)·e^{-ηs}(E C_T(s) - O(e^{-s}))/k, with #elig >= C(t)/4 and E C_T(s) >= F(s) - k. |
| 5.4 | F(s) → ∞ because the promotion rate Σ L_u/D_u >= 1. | VALID. Checked on all 81 135 trees with 3 to 17 vertices, with 0 violations. |
| 5.5 | Blocks of length h+s give liminf t^{-1} log F >= (h+s)^{-1} log B_s. | VALID for every s with B_s > 0. |
| 5.6 | A sequence of s realising the limsup gives liminf >= limsup - η. | VALID, since log(F(s)-k-O(e^{-s})) = log F(s) + o(1) and h, p are fixed. |
| 5.7 | M → ∞, eps ↓ 0. | VALID. η = 7 log 2/log(1/eps) + 9/M → 0, and the order of limits is right: s → ∞ first, then the parameters. |
| 5.8 | a_n = E C_n: F(t) = E a_{N(t)}; N from 2 is a sum of two Geometric(e^{-t}); a_n = O(F(log n)); F((1-α)log n) <= a_n + o(1). | VALID. The jump chain is independent of the holding times because the total rate equals the size. P(N(log n) >= n) → P(Gamma(2,1) >= 1) = 2/e. a_n is nondecreasing, and a_m <= m together with the geometric tail exp(-n^α) gives the o(1). So log a_n/log n → μ. |

---

## 2. `sample-exponent.md`, sentence by sentence

### 2.1 Disjoint sets of clocks
| # | Sentence | Verdict |
|---|---|---|
| A.1 | A counts core vertices with c <= 2 and L >= 1; the c >= 3 vertices number at most A. | VALID (n3 <= n1 - 2 <= A - 2). |
| A.2 | Each maximal B-path lies between two non-B vertices. | VALID. B vertices have core degree 2, and core leaves are never of type B. |
| A.3 | Suppressing the paths gives a tree on the non-B vertices, so there are <= 2A-1 paths. | VALID. There are at most 2A-2 non-B vertices, so at most 2A-3 paths. |
| A.4 | C <= 6A + 3J. | VALID. C <= A + (A-2) + 3J + 2(2A-3). |
| A.5 | At least C/6 sites. | VALID: A + J >= (6A+3J)/6. The constant is sharp: doubly subdivided cubic caterpillars give sites/C → 1/6 (0.1678 at C = 1204). |
| A.6 | A one-core tree supplies its centre. | VALID. |
| A.7 | Reserved sets ({u, one leaf} or {u, w, z}) are disjoint. | VALID. A-vertices are not of type B, the triples are disjoint, and a leaf has a single neighbour. In a triple, N(w) = {u, z} exactly; this was asserted on all trees with up to 17 vertices, and the disjointness check passed with 0 violations. |

### 2.2 The preparation
| # | Bullet | Verdict |
|---|---|---|
| B.1 | Unit 1. A-sites: u silent. Triples: u and z silent, and w's first ring picks u (u first in the ordering, mark < 1/2). Flags are independent; thin to a common probability. | VALID. deg w = 2 at its first ring (same argument as S1b), so "mark < 1/2" is exactly "chooses u". The flags are functions of disjoint reserved clocks and marks. |
| B.2 | Unit 2: u silent and >= M rings of the tagged leaf; the tagged leaves are distinct and their future clocks are independent standard clocks. | VALID. At a triple the tagged leaf is the newborn of w's first ring, which has a fresh clock. At an A-site the reserved leaf is still a leaf because u was silent. Increments on unit 2 are independent of everything before. |
| B.3 | Unit 3: exactly one ring of u, mark < M/(M+2) with leaves ranked first; "always promotes some leaf"; the identity may change but promotion cannot fail; v and v' silent, then thinned to e^{-2}. | VALID; this is the point the brief asked about. With leaves first, L_u >= M and c_u <= 2, U < M/(M+2) implies ceil(U·D_u) <= L_u whatever the arrivals, so the first two requirements depend only on u's own clock and mark. The third is evaluated on the clocks of v and v' after τ. v is F_τ-measurable and τ is a stopping time, so re-sampling those clocks at τ preserves the law; this is the strong Markov property, and the "For clarity" paragraph says exactly this. The flag is therefore a function of independent site-specific inputs, and the product law holds. Thinning by e^{-2τ} is VALID. Test: see output (iv). |
| B.4 | The next k-2 units: seed growth with u silent; disjoint interiors; common p by thinning. | VALID. Only fresh branch clocks and u's reserved clock are involved. |
| B.5 | "For clarity ... sampled at selection ... thinned ... product Bernoulli law ... only the edge to u ..." | VALID. |

### 2.3 Independent contributions in a block
| # | Sentence | Verdict |
|---|---|---|
| C.1 | At each prepared site: a fresh reference, acceptance marks at v, the exterior test on u's own clock; Y = reference count on the event (good ∧ assigned to v ∧ exterior test ∧ all acceptances), else 0. | VALID. Y_i = prep_i·Z_i, where Z_i depends only on that site's branch clocks and marks after preparation, its acceptance marks, and u_i's clock and marks on the comparison window. The Y_i are i.i.d. given the tree at the block start. |
| C.2 | E Y_s >= (p/k) e^{-ηs}[F(s)-k-O(e^{-s})]. | VALID; the same computation as in Claim 1. |
| C.3 | The sum of the contributions is <= the actual core count at the block end. | VALID. The branches are disjoint, and on Y_i > 0 the branch equals the reference tree. |
| C.4 | Truncate at R_s with E min(Y,R) >= EY/2; Chebyshev gives P(< C·EY/24) → 0 as C → ∞. | VALID. With m >= C/6 sites, P(S < C·EY/24) <= P(S < m·EY/4) <= 16R²/(m(EY)²) <= 96R²/(C(EY)²), uniformly in the tree. |
| C.5 | "holds in their own filtration" | VALID. The event lies in the natural σ-field and the bound is a function of C(t) alone, so it passes to the coarser filtration by the tower property. |
| C.6 | C → ∞ a.s. (conditional Borel–Cantelli with step probability >= 1/n; Yule never explodes). | VALID. |
| C.7 | The fraction of failed blocks → 0 a.s. (bounded martingale differences, weights 1/j, Kronecker). | VALID. |
| C.8 | liminf log C(t)/t >= log(EY_s/24)/(h+s) a.s. when EY_s/24 > 1. | VALID. Failed blocks lose nothing, by monotonicity. |
| C.9 | "The existence and positivity of μ give the lower bound μ. (If μ = 0 ... automatic.)" | NEEDS REPAIR (wording only). Positivity of μ is never proved, but the parenthetical μ = 0 case makes it unnecessary. Read "existence" for "existence and positivity". The actual argument (μ > η, s → ∞ along a sequence, then η ↓ 0 over countably many parameter choices) is VALID. |
| C.10 | Upper bound by Markov at integer times, Borel–Cantelli and monotonicity. | VALID: P(C(j) > e^{(μ+δ)j}) <= e^{(-δ+o(1))j}, which is summable. |
| C.11 | The hitting time of size n divided by log n → 1 a.s., which gives the discrete statement. | VALID. The variances Σ 1/j² are finite (Kolmogorov), and C_n = C(τ_n) pathwise. |

---

## 3. Numerical tests (commands and real output)

### (i) Exhaustive enumeration of all trees (= all reachable trees)
Every unlabelled tree with >= 2 vertices is reachable: a leaf can be attached to any W by picking a neighbour of W as V. Tree counts per n match OEIS A000055 (1, 1, 2, 3, 6, 11, 23, 47, 106, 235, ..., 48629).

(i-a) n <= 11, with the tight cases E = (C+5)/4 listed:
```
$ taskset -c 0 $PY enum_trees.py 11 --show-tight
n  #trees  minC  maxC  min E/C (C>=1)   #trees with E=ceil(C/4)  min sites/C   #E=(C+5)/4
 2       1  (only C=0 trees)
 3       1     1    1         1 (C=1,E=1)        1                     1 (C=1,sites=1)   0
 4       2     1    2         1 (C=2,E=2)        1                     1 (C=2,sites=2)   0
 5       3     1    3       2/3 (C=3,E=2)        1                   2/3 (C=3,sites=2)   1
 6       6     1    4       2/3 (C=3,E=2)        1                   1/2 (C=4,sites=2)   1
 7      11     1    5       2/3 (C=3,E=2)        1                   1/2 (C=4,sites=2)   2
 8      23     1    6       3/5 (C=5,E=3)        1                   1/2 (C=6,sites=3)   2
 9      47     1    7       1/2 (C=6,E=3)        1                   3/7 (C=7,sites=3)   3
10     106     1    8       3/7 (C=7,E=3)        1                   3/7 (C=7,sites=3)   4
   tight E=(C+5)/4: (((()()))((()))((()))) C=7 E=3 I=3 n1=3 n3=1
   tight E=(C+5)/4: (((()))((()))((()))()) C=7 E=3 I=3 n1=3 n3=1
   tight E=(C+5)/4: ((()()()()()()())(())) C=3 E=2 I=1 n1=2 n3=0
   tight E=(C+5)/4: ((()()()()()())(()())) C=3 E=2 I=1 n1=2 n3=0
   tight E=(C+5)/4: ((()()()()())(()()())) C=3 E=2 I=1 n1=2 n3=0
   tight E=(C+5)/4: ((()()()())(()()()())) C=3 E=2 I=1 n1=2 n3=0
11     235     1    9       3/7 (C=7,E=3)        1                   3/8 (C=8,sites=3)   6
total trees checked: 435
violations:
  E>=C/4                   0
  E>=(C+5)/4 (C>=2)        0
  I<=(C-1)/2               0
  n3<=n1-2                 0
  sites>=C/6               0
  C<=6A+3J                 0
  paths<=2A-1              0
  n3<=A                    0
  reserved disjoint        0
  stage-1 preconditions    0
  sum L/D >= 1             0
```
(i-b) n <= 17 (81 136 trees):
```
$ time taskset -c 0 $PY enum_trees.py 17
n  #trees  minC  maxC  min E/C (C>=1)   #trees with E=ceil(C/4)  min sites/C   #E=(C+5)/4
 2       1  (only C=0 trees)
 3       1     1    1         1 (C=1,E=1)        1                     1 (C=1,sites=1)   0
 4       2     1    2         1 (C=2,E=2)        1                     1 (C=2,sites=2)   0
 5       3     1    3       2/3 (C=3,E=2)        1                   2/3 (C=3,sites=2)   1
 6       6     1    4       2/3 (C=3,E=2)        1                   1/2 (C=4,sites=2)   1
 7      11     1    5       2/3 (C=3,E=2)        1                   1/2 (C=4,sites=2)   2
 8      23     1    6       3/5 (C=5,E=3)        1                   1/2 (C=6,sites=3)   2
 9      47     1    7       1/2 (C=6,E=3)        1                   3/7 (C=7,sites=3)   3
10     106     1    8       3/7 (C=7,E=3)        1                   3/7 (C=7,sites=3)   4
11     235     1    9       3/7 (C=7,E=3)        1                   3/8 (C=8,sites=3)   6
12     551     1   10       3/7 (C=7,E=3)        1                   1/3 (C=9,sites=3)   8
13    1301     1   11       3/7 (C=7,E=3)        1                  3/10 (C=10,sites=3)   12
14    3159     1   12       2/5 (C=10,E=4)        1                  3/10 (C=10,sites=3)   16
15    7741     1   13      4/11 (C=11,E=4)        1                  3/10 (C=10,sites=3)   23
16   19320     1   14      4/11 (C=11,E=4)        1                  3/10 (C=10,sites=3)   31
17   48629     1   15      4/11 (C=11,E=4)        1                  3/10 (C=10,sites=3)   45
total trees checked: 81136
violations:
  E>=C/4                   0
  E>=(C+5)/4 (C>=2)        0
  I<=(C-1)/2               0
  n3<=n1-2                 0
  sites>=C/6               0
  C<=6A+3J                 0
  paths<=2A-1              0
  n3<=A                    0
  reserved disjoint        0
  stage-1 preconditions    0
  sum L/D >= 1             0

real	0m22.577s
user	0m21.181s
sys	0m0.016s
```
(i-c) Planted errors. `--plant-eligible` drops "B with a B neighbour" from eligibility; `--plant-I` lets I count all B vertices; `--plant-sites` drops the triple sites. All three are caught:
```
$ for f in --plant-eligible --plant-I; do taskset -c 0 $PY enum_trees.py 11 $f; done; taskset -c 0 $PY enum_trees.py 16 --plant-sites
=== --plant-eligible (n<=11)
violations:
  E>=C/4                   1
     first counterexample: (((((()))))((((())))))  {'n': 11, 'C': 9, 'E': 2, 'I': 0, 'n1': 2, 'n3': 0, 'A': 2, 'J': 2, 'sites': 4, 'npaths': 1, 'disjoint': True, 'pre_ok': True, 'promo': Fraction(1, 1)}
  E>=(C+5)/4 (C>=2)        41
     first counterexample: (((()))(()))  {'n': 6, 'C': 4, 'E': 2, 'I': 0, 'n1': 2, 'n3': 0, 'A': 2, 'J': 0, 'sites': 2, 'npaths': 1, 'disjoint': True, 'pre_ok': True, 'promo': Fraction(1, 1)}
  I<=(C-1)/2               0
  n3<=n1-2                 0
  sites>=C/6               0
  C<=6A+3J                 0
  paths<=2A-1              0
  n3<=A                    0
  reserved disjoint        0
  stage-1 preconditions    0
  sum L/D >= 1             0
=== --plant-I (n<=11)
violations:
  E>=C/4                   0
  E>=(C+5)/4 (C>=2)        0
  I<=(C-1)/2               72
     first counterexample: (((()))(()))  {'n': 6, 'C': 4, 'E': 4, 'I': 2, 'n1': 2, 'n3': 0, 'A': 2, 'J': 0, 'sites': 2, 'npaths': 1, 'disjoint': True, 'pre_ok': True, 'promo': Fraction(1, 1)}
  n3<=n1-2                 0
  sites>=C/6               0
  C<=6A+3J                 0
  paths<=2A-1              0
  n3<=A                    0
  reserved disjoint        0
  stage-1 preconditions    0
  sum L/D >= 1             0
=== --plant-sites (n<=16)
violations:
  E>=C/4                   0
  E>=(C+5)/4 (C>=2)        0
  I<=(C-1)/2               0
  n3<=n1-2                 0
  sites>=C/6               3
     first counterexample: (((((((()))))))((((((())))))))  {'n': 15, 'C': 13, 'E': 13, 'I': 0, 'n1': 2, 'n3': 0, 'A': 2, 'J': 3, 'sites': 2, 'npaths': 1, 'disjoint': True, 'pre_ok': True, 'promo': Fraction(1, 1)}
  C<=6A+3J                 0
  paths<=2A-1              0
  n3<=A                    0
  reserved disjoint        0
  stage-1 preconditions    0
  sum L/D >= 1             0
```
(i-d) Asymptotic sharpness of 1/4 and 1/6:
```
$ taskset -c 0 $PY families.py
d = 1
  r=   1 n=   10 C=    7 eligible=   3 E/C=0.4286 (C+5)/4=3.00  sites=   3 sites/C=0.4286
  r=   2 n=   15 C=   11 eligible=   4 E/C=0.3636 (C+5)/4=4.00  sites=   4 sites/C=0.3636
  r=   5 n=   30 C=   23 eligible=   7 E/C=0.3043 (C+5)/4=7.00  sites=   7 sites/C=0.3043
  r=  10 n=   55 C=   43 eligible=  12 E/C=0.2791 (C+5)/4=12.00  sites=  12 sites/C=0.2791
  r=  50 n=  255 C=  203 eligible=  52 E/C=0.2562 (C+5)/4=52.00  sites=  52 sites/C=0.2562
  r= 200 n= 1005 C=  803 eligible= 202 E/C=0.2516 (C+5)/4=202.00  sites= 202 sites/C=0.2516
d = 2
  r=   1 n=   13 C=   10 eligible=   9 E/C=0.9000 (C+5)/4=3.75  sites=   3 sites/C=0.3000
  r=   2 n=   20 C=   16 eligible=  14 E/C=0.8750 (C+5)/4=5.25  sites=   4 sites/C=0.2500
  r=   5 n=   41 C=   34 eligible=  29 E/C=0.8529 (C+5)/4=9.75  sites=   7 sites/C=0.2059
  r=  10 n=   76 C=   64 eligible=  54 E/C=0.8438 (C+5)/4=17.25  sites=  12 sites/C=0.1875
  r=  50 n=  356 C=  304 eligible= 254 E/C=0.8355 (C+5)/4=77.25  sites=  52 sites/C=0.1711
  r= 200 n= 1406 C= 1204 eligible=1004 E/C=0.8339 (C+5)/4=302.25  sites= 202 sites/C=0.1678
```

### (ii) The comparison inequality of step (c)
Actual run: the seed (T,v) hangs from u, with D_u = M counting v. u's other neighbours are a core vertex x carrying 2 leaves, plus M-2 leaves, and the whole tree evolves. The M = inf rows mean u never chooses v. This is the most stringent case, because the exterior then gives the branch no help. Columns: E act, the unweighted E C_ref, the claimed RHS e^{-s/M}E[C_ref W], and the planted RHS e^{+s/M}E[C_ref Π(D+1)/D].
```
$ time taskset -c 0 $PY comparison_mc.py 100000 20261122
NRUNS=100000 seed=20261122
seed                  M   s |          E act |  E Cref (unw.) |  exp(-s/M)E[CW] | z=(act-RHS)/se |   PLANTED RHS    z_pl
------------------------------------------------------------------------------------------------------------------------
K2 (v leaf)           3   1 |   1.655+-0.003 |   1.470+-0.003 |    0.582+-0.001 |          345.2 |   5.518+-0.03  -122.6
K2 (v leaf)           3   2 |   2.806+-0.006 |   2.960+-0.006 |    0.657+-0.002 |          364.1 |  31.156+-0.36   -78.7
K2 (v leaf)           6   1 |   1.639+-0.003 |   1.468+-0.003 |    0.688+-0.002 |          299.2 |   4.649+-0.03  -112.9
K2 (v leaf)           6   2 |   2.750+-0.006 |   2.950+-0.006 |    0.915+-0.002 |          305.6 |  21.521+-0.23   -81.7
K2 (v leaf)         inf   1 |   1.629+-0.003 |   1.472+-0.003 |    0.813+-0.002 |          247.7 |   3.994+-0.02  -104.5
K2 (v leaf)         inf   2 |   2.706+-0.005 |   2.950+-0.006 |    1.277+-0.003 |          228.5 |  15.584+-0.17   -76.6
P3 v=end (leaf)       3   1 |   2.804+-0.003 |   2.120+-0.004 |    0.960+-0.002 |          474.6 |   7.333+-0.04  -106.0
P3 v=end (leaf)       3   2 |   4.457+-0.007 |   3.848+-0.007 |    0.861+-0.003 |          494.2 |  42.694+-0.42   -91.0
P3 v=end (leaf)       6   1 |   2.777+-0.003 |   2.114+-0.004 |    1.127+-0.003 |          404.1 |   6.169+-0.03  -113.3
P3 v=end (leaf)       6   2 |   4.346+-0.007 |   3.872+-0.007 |    1.202+-0.004 |          416.7 |  31.107+-0.31   -86.3
P3 v=end (leaf)     inf   1 |   2.754+-0.003 |   2.116+-0.004 |    1.335+-0.003 |          323.0 |   5.182+-0.02   -99.6
P3 v=end (leaf)     inf   2 |   4.262+-0.007 |   3.864+-0.007 |    1.688+-0.005 |          312.6 |  21.998+-0.23   -78.4
P3 v=middle           3   1 |   1.869+-0.003 |   2.121+-0.004 |    1.028+-0.001 |          235.5 |   5.009+-0.02  -188.8
P3 v=middle           3   2 |   3.395+-0.007 |   3.868+-0.007 |    1.151+-0.002 |          316.6 |  16.123+-0.06  -200.9
P3 v=middle           6   1 |   1.867+-0.003 |   2.112+-0.004 |    1.213+-0.001 |          179.0 |   4.203+-0.01  -168.5
P3 v=middle           6   2 |   3.372+-0.007 |   3.853+-0.007 |    1.606+-0.003 |          243.3 |  11.515+-0.05  -173.1
P3 v=middle         inf   1 |   1.857+-0.003 |   2.112+-0.004 |    1.434+-0.002 |          112.9 |   3.559+-0.01  -141.9
P3 v=middle         inf   2 |   3.342+-0.007 |   3.863+-0.007 |    2.242+-0.004 |          143.3 |   8.302+-0.03  -147.0
K1,3 v=leaf           3   1 |   2.994+-0.004 |   2.172+-0.004 |    0.973+-0.002 |          461.8 |   7.661+-0.04  -119.3
K1,3 v=leaf           3   2 |   4.875+-0.008 |   4.066+-0.008 |    0.870+-0.003 |          493.2 |  48.327+-0.46   -94.5
K1,3 v=leaf           6   1 |   2.959+-0.004 |   2.171+-0.004 |    1.149+-0.003 |          393.2 |   6.459+-0.03  -110.2
K1,3 v=leaf           6   2 |   4.759+-0.007 |   4.063+-0.008 |    1.208+-0.004 |          422.0 |  34.593+-0.29  -104.1
K1,3 v=leaf         inf   1 |   2.932+-0.004 |   2.168+-0.004 |    1.353+-0.003 |          322.7 |   5.474+-0.03   -92.5
K1,3 v=leaf         inf   2 |   4.666+-0.007 |   4.056+-0.008 |    1.694+-0.006 |          325.3 |  24.535+-0.21   -95.4
K1,3 v=centre         3   1 |   1.983+-0.004 |   2.171+-0.004 |    1.174+-0.002 |          203.3 |   4.275+-0.01  -186.1
K1,3 v=centre         3   2 |   3.680+-0.007 |   4.057+-0.008 |    1.417+-0.002 |          291.7 |  12.689+-0.04  -232.4
K1,3 v=centre         6   1 |   1.978+-0.004 |   2.181+-0.004 |    1.388+-0.002 |          144.3 |   3.652+-0.01  -155.9
K1,3 v=centre         6   2 |   3.657+-0.007 |   4.061+-0.008 |    1.977+-0.003 |          208.8 |   9.124+-0.03  -190.5
K1,3 v=centre       inf   1 |   1.968+-0.004 |   2.174+-0.004 |    1.636+-0.002 |           78.8 |   3.080+-0.01  -119.8
K1,3 v=centre       inf   2 |   3.651+-0.007 |   4.058+-0.008 |    2.756+-0.004 |          104.3 |   6.545+-0.02  -135.4
K1,5 v=centre         3   1 |   2.085+-0.004 |   2.225+-0.004 |    1.328+-0.002 |          171.4 |   3.815+-0.01  -173.1
K1,5 v=centre         3   2 |   3.924+-0.008 |   4.220+-0.008 |    1.702+-0.003 |          262.2 |  10.772+-0.03  -240.6
K1,5 v=centre         6   1 |   2.089+-0.004 |   2.225+-0.004 |    1.568+-0.002 |          113.4 |   3.232+-0.01  -130.0
K1,5 v=centre         6   2 |   3.938+-0.008 |   4.225+-0.008 |    2.375+-0.004 |          174.7 |   7.738+-0.02  -178.4
K1,5 v=centre       inf   1 |   2.081+-0.004 |   2.223+-0.004 |    1.851+-0.003 |           47.9 |   2.731+-0.01   -85.0
K1,5 v=centre       inf   2 |   3.920+-0.008 |   4.215+-0.008 |    3.313+-0.005 |           62.7 |   5.518+-0.01   -99.4
P5 v=end              3   1 |   5.137+-0.004 |   4.493+-0.004 |    2.025+-0.004 |          587.9 |  15.723+-0.07  -155.9
P5 v=end              3   2 |   7.847+-0.009 |   7.409+-0.009 |    1.645+-0.004 |          650.8 |  85.522+-0.91   -85.6
P5 v=end              6   1 |   5.107+-0.004 |   4.499+-0.004 |    2.395+-0.004 |          471.8 |  13.283+-0.06  -141.1
P5 v=end              6   2 |   7.773+-0.009 |   7.414+-0.009 |    2.309+-0.006 |          531.0 |  61.580+-0.68   -79.3
P5 v=end            inf   1 |   5.081+-0.004 |   4.492+-0.004 |    2.829+-0.005 |          356.7 |  11.150+-0.05  -122.7
P5 v=end            inf   2 |   7.663+-0.008 |   7.424+-0.009 |    3.216+-0.008 |          383.1 |  44.847+-0.69   -54.0
smallest z for the claimed inequality: (47.938155630250705, 'K1,5 v=centre', 'inf', 1.0)

real	5m35.419s
user	2m45.371s
sys	0m0.112s
```
The sign-error-only planted weight e^{+s/M}·Π D/(D+1), evaluated from the same runs:
```
$ taskset -c 0 $PY planted_expplus.py out_comparison_mc.txt
K2 (v leaf)        M=3 s=1: E act=1.655  exp(+s/M)E[CW]=1.134  z=+145.8
K2 (v leaf)        M=3 s=2: E act=2.806  exp(+s/M)E[CW]=2.492  z=+32.4
K2 (v leaf)        M=6 s=1: E act=1.639  exp(+s/M)E[CW]=0.960  z=+165.7
K2 (v leaf)        M=6 s=2: E act=2.750  exp(+s/M)E[CW]=1.782  z=+135.3
P3 v=end (leaf)    M=3 s=1: E act=2.804  exp(+s/M)E[CW]=1.870  z=+190.0
P3 v=end (leaf)    M=3 s=2: E act=4.457  exp(+s/M)E[CW]=3.266  z=+89.1
P3 v=end (leaf)    M=6 s=1: E act=2.777  exp(+s/M)E[CW]=1.573  z=+233.8
P3 v=end (leaf)    M=6 s=2: E act=4.346  exp(+s/M)E[CW]=2.341  z=+191.4
P3 v=middle        M=3 s=1: E act=1.869  exp(+s/M)E[CW]=2.002  z=-37.3
P3 v=middle        M=3 s=2: E act=3.395  exp(+s/M)E[CW]=4.367  z=-94.1
P3 v=middle        M=6 s=1: E act=1.867  exp(+s/M)E[CW]=1.693  z=+52.6
P3 v=middle        M=6 s=2: E act=3.372  exp(+s/M)E[CW]=3.128  z=+26.8
K1,3 v=leaf        M=3 s=1: E act=2.994  exp(+s/M)E[CW]=1.895  z=+196.8
K1,3 v=leaf        M=3 s=2: E act=4.875  exp(+s/M)E[CW]=3.300  z=+113.2
K1,3 v=leaf        M=6 s=1: E act=2.959  exp(+s/M)E[CW]=1.604  z=+234.1
K1,3 v=leaf        M=6 s=2: E act=4.759  exp(+s/M)E[CW]=2.353  z=+229.7
K1,3 v=centre      M=3 s=1: E act=1.983  exp(+s/M)E[CW]=2.287  z=-54.4
K1,3 v=centre      M=3 s=2: E act=3.680  exp(+s/M)E[CW]=5.376  z=-164.3
K1,3 v=centre      M=6 s=1: E act=1.978  exp(+s/M)E[CW]=1.937  z=+8.4
K1,3 v=centre      M=6 s=2: E act=3.657  exp(+s/M)E[CW]=3.851  z=-21.2
K1,5 v=centre      M=3 s=1: E act=2.085  exp(+s/M)E[CW]=2.587  z=-89.8
K1,5 v=centre      M=3 s=2: E act=3.924  exp(+s/M)E[CW]=6.457  z=-182.1
K1,5 v=centre      M=6 s=1: E act=2.089  exp(+s/M)E[CW]=2.188  z=-20.4
K1,5 v=centre      M=6 s=2: E act=3.938  exp(+s/M)E[CW]=4.626  z=-61.6
P5 v=end           M=3 s=1: E act=5.137  exp(+s/M)E[CW]=3.944  z=+136.2
P5 v=end           M=3 s=2: E act=7.847  exp(+s/M)E[CW]=6.241  z=+91.1
P5 v=end           M=6 s=1: E act=5.107  exp(+s/M)E[CW]=3.342  z=+256.9
P5 v=end           M=6 s=2: E act=7.773  exp(+s/M)E[CW]=4.497  z=+222.1
configs: 28; planted inequality violated (z<-3) in 9
```
(ii-b) The coupling of step (c) implemented literally, with pathwise and law checks:
```
$ time taskset -c 0 $PY coupling_pathwise.py 100000 4242
NRUNS=100000 seed=4242 plant_order=False plant_accept=False

K2 (v leaf)  M=3 s=2
  (P1) runs with E_u & Acc: 25524; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.5120 +- 0.0016   exp(-s/M) = 0.5134
       E[C_ref 1(E_u&Acc)] = 0.6536 +- 0.0046   exp(-s/M)E[C_ref W] = 0.6545 +- 0.0016   z = -0.18
  (P3) E[act] coupled = 2.8166 +- 0.0057   direct = 2.8063 +- 0.0057   z = 1.28

P3 v=end (leaf)  M=3 s=2
  (P1) runs with E_u & Acc: 22236; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.5117 +- 0.0016   exp(-s/M) = 0.5134
       E[C_ref 1(E_u&Acc)] = 0.8646 +- 0.0061   exp(-s/M)E[C_ref W] = 0.8675 +- 0.0026   z = -0.43
  (P3) E[act] coupled = 4.4527 +- 0.0068   direct = 4.4544 +- 0.0068   z = -0.17

K1,3 v=centre  M=6 s=2
  (P1) runs with E_u & Acc: 54162; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 6
  (P2) P(E_u) = 0.7164 +- 0.0014   exp(-s/M) = 0.7165
       E[C_ref 1(E_u&Acc)] = 1.9678 +- 0.0079   exp(-s/M)E[C_ref W] = 1.9755 +- 0.0031   z = -0.92
  (P3) E[act] coupled = 3.6485 +- 0.0074   direct = 3.6559 +- 0.0074   z = -0.70

K1,3 v=leaf  M=6 s=1
  (P1) runs with E_u & Acc: 52751; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 6
  (P2) P(E_u) = 0.8469 +- 0.0011   exp(-s/M) = 0.8465
       E[C_ref 1(E_u&Acc)] = 1.1500 +- 0.0045   exp(-s/M)E[C_ref W] = 1.1467 +- 0.0029   z = 0.62
  (P3) E[act] coupled = 2.9600 +- 0.0036   direct = 2.9604 +- 0.0036   z = -0.07

P5 v=end  M=3 s=1
  (P1) runs with E_u & Acc: 45402; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.7176 +- 0.0014   exp(-s/M) = 0.7165
       E[C_ref 1(E_u&Acc)] = 2.0278 +- 0.0076   exp(-s/M)E[C_ref W] = 2.0227 +- 0.0037   z = 0.61
  (P3) E[act] coupled = 5.1291 +- 0.0038   direct = 5.1327 +- 0.0038   z = -0.67

real	1m16.510s
user	0m40.920s
sys	0m0.050s
```
Planted errors: `--plant-order` puts v last in u's ordering, which breaks P1; `--plant-accept` uses acceptance probability (D+1)/(D+2), which breaks the law, P2 and P3.
```
$ for f in --plant-order --plant-accept; do taskset -c 0 $PY coupling_pathwise.py 100000 4243 $f; done
NRUNS=100000 seed=4243 plant_order=True plant_accept=False

K2 (v leaf)  M=3 s=2
  (P1) runs with E_u & Acc: 25485; of these with a discrepancy or act<ref: 9467;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.5114 +- 0.0016   exp(-s/M) = 0.5134
       E[C_ref 1(E_u&Acc)] = 0.6507 +- 0.0045   exp(-s/M)E[C_ref W] = 0.6568 +- 0.0016   z = -1.27
  (P3) E[act] coupled = 2.8177 +- 0.0057   direct = 2.8173 +- 0.0057   z = 0.05

P3 v=end (leaf)  M=3 s=2
  (P1) runs with E_u & Acc: 22377; of these with a discrepancy or act<ref: 8356;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.5151 +- 0.0016   exp(-s/M) = 0.5134
       E[C_ref 1(E_u&Acc)] = 0.8632 +- 0.0061   exp(-s/M)E[C_ref W] = 0.8614 +- 0.0026   z = 0.27
  (P3) E[act] coupled = 4.4449 +- 0.0068   direct = 4.4478 +- 0.0068   z = -0.31

K1,3 v=centre  M=6 s=2
  (P1) runs with E_u & Acc: 54467; of these with a discrepancy or act<ref: 8869;  min D_u seen at u-rings: 6
  (P2) P(E_u) = 0.7192 +- 0.0014   exp(-s/M) = 0.7165
       E[C_ref 1(E_u&Acc)] = 1.9963 +- 0.0080   exp(-s/M)E[C_ref W] = 1.9801 +- 0.0031   z = 1.89
  (P3) E[act] coupled = 3.6623 +- 0.0074   direct = 3.6507 +- 0.0074   z = 1.11

K1,3 v=leaf  M=6 s=1
  (P1) runs with E_u & Acc: 52726; of these with a discrepancy or act<ref: 6144;  min D_u seen at u-rings: 6
  (P2) P(E_u) = 0.8476 +- 0.0011   exp(-s/M) = 0.8465
       E[C_ref 1(E_u&Acc)] = 1.1567 +- 0.0045   exp(-s/M)E[C_ref W] = 1.1534 +- 0.0029   z = 0.62
  (P3) E[act] coupled = 2.9667 +- 0.0036   direct = 2.9631 +- 0.0036   z = 0.70

P5 v=end  M=3 s=1
  (P1) runs with E_u & Acc: 45237; of these with a discrepancy or act<ref: 11176;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.7165 +- 0.0014   exp(-s/M) = 0.7165
       E[C_ref 1(E_u&Acc)] = 2.0167 +- 0.0076   exp(-s/M)E[C_ref W] = 2.0231 +- 0.0037   z = -0.77
  (P3) E[act] coupled = 5.1332 +- 0.0038   direct = 5.1328 +- 0.0038   z = 0.08
NRUNS=100000 seed=4243 plant_order=False plant_accept=True

K2 (v leaf)  M=3 s=2
  (P1) runs with E_u & Acc: 30891; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.5142 +- 0.0016   exp(-s/M) = 0.5134
       E[C_ref 1(E_u&Acc)] = 0.8210 +- 0.0050   exp(-s/M)E[C_ref W] = 0.6526 +- 0.0016   z = 31.85
  (P3) E[act] coupled = 2.9514 +- 0.0058   direct = 2.8173 +- 0.0057   z = 16.44

P3 v=end (leaf)  M=3 s=2
  (P1) runs with E_u & Acc: 28398; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.5125 +- 0.0016   exp(-s/M) = 0.5134
       E[C_ref 1(E_u&Acc)] = 1.0964 +- 0.0067   exp(-s/M)E[C_ref W] = 0.8600 +- 0.0026   z = 33.02
  (P3) E[act] coupled = 4.4957 +- 0.0068   direct = 4.4469 +- 0.0068   z = 5.07

K1,3 v=centre  M=6 s=2
  (P1) runs with E_u & Acc: 56507; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 6
  (P2) P(E_u) = 0.7171 +- 0.0014   exp(-s/M) = 0.7165
       E[C_ref 1(E_u&Acc)] = 2.1025 +- 0.0081   exp(-s/M)E[C_ref W] = 1.9796 +- 0.0031   z = 14.13
  (P3) E[act] coupled = 3.7113 +- 0.0075   direct = 3.6486 +- 0.0074   z = 5.96

K1,3 v=leaf  M=6 s=1
  (P1) runs with E_u & Acc: 61456; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 6
  (P2) P(E_u) = 0.8451 +- 0.0011   exp(-s/M) = 0.8465
       E[C_ref 1(E_u&Acc)] = 1.3364 +- 0.0046   exp(-s/M)E[C_ref W] = 1.1498 +- 0.0029   z = 34.60
  (P3) E[act] coupled = 2.9684 +- 0.0036   direct = 2.9647 +- 0.0036   z = 0.72

P5 v=end  M=3 s=1
  (P1) runs with E_u & Acc: 52164; of these with a discrepancy or act<ref: 0;  min D_u seen at u-rings: 3
  (P2) P(E_u) = 0.7157 +- 0.0014   exp(-s/M) = 0.7165
       E[C_ref 1(E_u&Acc)] = 2.3251 +- 0.0077   exp(-s/M)E[C_ref W] = 2.0191 +- 0.0036   z = 36.08
  (P3) E[act] coupled = 5.1568 +- 0.0039   direct = 5.1347 +- 0.0038   z = 4.07
```

### (iii) Poisson tail inequalities of the good-history section
```
$ time taskset -c 0 $PY poisson_tails.py
eps=0.01 M=50 A=8.0 b=7/log(1/eps)=1.520031 k=30000  plants: A=False chernoff=False k=False
A log A - A + 1 = 9.6355   (need > 5)
b log(b/eps) - b + eps = 6.1265  (need > 5);  b > eps: True
b(eps) check over a grid of eps:
   eps=0.5      b=10.0989 rate=20.7540 b>eps=True
   eps=0.2      b=4.3493 rate=9.2443 b>eps=True
   eps=0.1      b=3.0401 rate=7.4401 b>eps=True
   eps=0.05     b=2.3367 rate=6.6965 b>eps=True
   eps=0.01     b=1.5200 rate=6.1265 b>eps=True
   eps=0.001    b=1.0134 rate=6.0011 b>eps=True
   eps=1e-06    b=0.5067 rate=6.1488 b>eps=True
   eps=1e-12    b=0.2533 rate=6.3988 b>eps=True

--- s = 20 ---
 [cross-check scipy logsf: -128.18, -196.04]
 log P(Pois(eps s=0.2) >= bs=30.40) = -128.18; Chernoff exponent -s*rate = -122.53; -5s = -100
 log P(Pois(s) >= As=160) = -196.04; Chernoff exponent = -192.71; -5s = -100
 degree condition: marks in [0,1/M] up to eps s ~ Pois(k eps s/M = 120);  threshold kM = 1500000
   log P(Pois(120) <= 1500000) = -0;  -4s = -80  -> holds: False
 ring-count bounds hold: True

--- s = 50 ---
 [cross-check scipy logsf: -314.43, -485.56]
 log P(Pois(eps s=0.5) >= bs=76.00) = -314.43; Chernoff exponent -s*rate = -306.32; -5s = -250
 log P(Pois(s) >= As=400) = -485.56; Chernoff exponent = -481.78; -5s = -250
 degree condition: marks in [0,1/M] up to eps s ~ Pois(k eps s/M = 300);  threshold kM = 1500000
   log P(Pois(300) <= 1500000) = -0;  -4s = -200  -> holds: False
 ring-count bounds hold: True

threshold s0 where log P(Pois(k eps s/M) <= kM) <= -4s:
   s0 ~ 2.46696e+06  (kM / (k eps / M) = M^2/eps = 250000)
   at s = 1.01*s0 = 2.49162e+06: log P(bad) <= -1.00009e+07, -4s = -9.9665e+06; holds: True

Yule second moment E[N^2] = k(k+1)e^{2t} - k e^t (exact sum):
   k= 1 t=1.0: sum=12.059830 formula=12.059830
   k= 3 t=2.0: sum=633.010632 formula=633.010632
   k= 6 t=1.5: sum=816.702416 formula=816.702416
   k=10 t=3.0: sum=44176.311915 formula=44176.311915

ring-count bound failures: 0
```
A note on method: my first version used `scipy.stats.poisson.logcdf`. It underflows to -inf far in the tail, which gave a spurious s0 ≈ 2.58·10^5. The script now sums the log-pmf directly (`gammaln` + `logsumexp`), and the result s0 ≈ 2.467·10^6 agrees with the Chernoff estimate y/3 = 1 + log y, y = 6s/(kM) ≈ 9.87. The tails at s = 20 and 50 are cross-checked against scipy, where scipy is safe.

Planted errors: A = 4 (rate 2.55 < 5), a Chernoff exponent that is too strong, and k = 3M/eps. All three are caught:
```
$ for f in --plant-A --plant-chernoff --plant-k; do taskset -c 0 $PY poisson_tails.py $f | grep -E "plants|A log A|log P|holds|s0|none found|failures"; done
===== --plant-A
eps=0.01 M=50 A=4.0 b=7/log(1/eps)=1.520031 k=30000  plants: A=True chernoff=False k=False
A log A - A + 1 = 2.5452   (need > 5)
 log P(Pois(eps s=0.2) >= bs=30.40) = -128.18; Chernoff exponent -s*rate = -122.53; -5s = -100
 log P(Pois(s) >= As=80) = -53.73; Chernoff exponent = -50.90; -5s = -100
   log P(Pois(120) <= 1500000) = -0;  -4s = -80  -> holds: False
 log P(Pois(eps s=0.5) >= bs=76.00) = -314.43; Chernoff exponent -s*rate = -306.32; -5s = -250
 log P(Pois(s) >= As=200) = -130.54; Chernoff exponent = -127.26; -5s = -250
   log P(Pois(300) <= 1500000) = -0;  -4s = -200  -> holds: False
threshold s0 where log P(Pois(k eps s/M) <= kM) <= -4s:
   s0 ~ 2.46696e+06  (kM / (k eps / M) = M^2/eps = 250000)
   at s = 1.01*s0 = 2.49162e+06: log P(bad) <= -6.34163e+06, -4s = -9.9665e+06; holds: False
ring-count bound failures: 2
===== --plant-chernoff
eps=0.01 M=50 A=8.0 b=7/log(1/eps)=1.520031 k=30000  plants: A=False chernoff=True k=False
A log A - A + 1 = 9.6355   (need > 5)
 log P(Pois(eps s=0.2) >= bs=30.40) = -128.18; Chernoff exponent -s*rate = -152.93; -5s = -100
 log P(Pois(s) >= As=160) = -196.04; Chernoff exponent = -352.71; -5s = -100
   log P(Pois(120) <= 1500000) = -0;  -4s = -80  -> holds: False
 log P(Pois(eps s=0.5) >= bs=76.00) = -314.43; Chernoff exponent -s*rate = -382.32; -5s = -250
 log P(Pois(s) >= As=400) = -485.56; Chernoff exponent = -881.78; -5s = -250
   log P(Pois(300) <= 1500000) = -0;  -4s = -200  -> holds: False
threshold s0 where log P(Pois(k eps s/M) <= kM) <= -4s:
   s0 ~ 2.46696e+06  (kM / (k eps / M) = M^2/eps = 250000)
   at s = 1.01*s0 = 2.49162e+06: log P(bad) <= -1.00009e+07, -4s = -9.9665e+06; holds: True
ring-count bound failures: 4
===== --plant-k
eps=0.01 M=50 A=8.0 b=7/log(1/eps)=1.520031 k=15000  plants: A=False chernoff=False k=True
A log A - A + 1 = 9.6355   (need > 5)
 log P(Pois(eps s=0.2) >= bs=30.40) = -128.18; Chernoff exponent -s*rate = -122.53; -5s = -100
 log P(Pois(s) >= As=160) = -196.04; Chernoff exponent = -192.71; -5s = -100
   log P(Pois(60) <= 750000) = -0;  -4s = -80  -> holds: False
 log P(Pois(eps s=0.5) >= bs=76.00) = -314.43; Chernoff exponent -s*rate = -306.32; -5s = -250
 log P(Pois(s) >= As=400) = -485.56; Chernoff exponent = -481.78; -5s = -250
   log P(Pois(150) <= 750000) = -0;  -4s = -200  -> holds: False
threshold s0 where log P(Pois(k eps s/M) <= kM) <= -4s:
   none found up to s = 1.09951e+12 (rate k eps/M = 3 <= 4: lower tail >= exp(-3 s))
ring-count bound failures: 0
```

### (iv) Claim 2, product law of the third-stage flags with outside arrivals
Two adjacent sites u1–u2, each with L leaves, run for one unit with every clock free. The proof's flag (u's mark < thr with leaves first, plus silence of v and v') is compared against the "natural" flag (u's ring picks a leaf), which serves as the planted error.
```
$ time (taskset -c 0 $PY stage3_independence.py 4000000 1 0.5 11; taskset -c 0 $PY stage3_independence.py 4000000 3 0.5 12)
NRUNS=4000000 L=1 thr=0.5 seed=11
  proof flags (mark<thr): P(F1)=0.07960 P(F2)=0.07976 P(F1F2)=0.006393 P(F1)P(F2)=0.006349  cov=+4.34e-05 (se~4.0e-05, z=+1.09)
     theory P(F_i) = 0.07952; z1=+0.58 z2=+1.78
  natural flags (chooses a leaf): P(F1)=0.10231 P(F2)=0.10252 P(F1F2)=0.009430 P(F1)P(F2)=0.010488  cov=-1.06e-03 (se~4.8e-05, z=-21.90)
NRUNS=4000000 L=3 thr=0.5 seed=12
  proof flags (mark<thr): P(F1)=0.07947 P(F2)=0.07941 P(F1F2)=0.006276 P(F1)P(F2)=0.006311  cov=-3.45e-05 (se~3.9e-05, z=-0.87)
     theory P(F_i) = 0.07952; z1=-0.36 z2=-0.87
  natural flags (chooses a leaf): P(F1)=0.13351 P(F2)=0.13319 P(F1F2)=0.017590 P(F1)P(F2)=0.017782  cov=-1.91e-04 (se~6.6e-05, z=-2.91)

real	3m39.714s
user	2m6.461s
sys	0m1.209s
```

### (v) Planting lemma, forced-script simulation in a hostile exterior
The script's clocks are forced and every other clock runs free. Checks: (a) the branch is isomorphic to (T,v); (b) D_u >= M; (c) the branch contains no old core vertex. Also recorded: L/D at the stage-3 ring, deg w at its forced ring, and the stage-4 source degrees.
```
$ time taskset -c 0 $PY planting_forced.py 1000 77
NRUNS=1000 seed=77 plant_x=False plant_u4=False
case L>=1 M=3 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=0 (k=2)  
case L>=1 M=3 P3 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=3)  
case L>=1 M=3 P3 v=mid       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=3)  
case L>=1 M=3 K1,3 v=centre  fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=4)  
case L>=1 M=3 K1,3 v=leaf    fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=4)  
case L>=1 M=3 P4 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=4)  
case L>=1 M=6 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=0 (k=2)  
case L>=1 M=6 P3 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=3)  
case L>=1 M=6 P3 v=mid       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=1 (k=3)  
case L>=1 M=6 K1,3 v=centre  fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=1 (k=4)  
case L>=1 M=6 K1,3 v=leaf    fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=4)  
case L>=1 M=6 P4 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=4)  
case B  M=3 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=0 (k=2)  deg w at ring: min 2 max 2
case B  M=3 P3 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=3)  deg w at ring: min 2 max 2
case B  M=3 P3 v=mid       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=3)  deg w at ring: min 2 max 2
case B  M=3 K1,3 v=centre  fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=4)  deg w at ring: min 2 max 2
case B  M=3 K1,3 v=leaf    fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=4)  deg w at ring: min 2 max 2
case B  M=3 P4 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=4)  deg w at ring: min 2 max 2
case B  M=6 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=0 (k=2)  deg w at ring: min 2 max 2
case B  M=6 P3 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=3)  deg w at ring: min 2 max 2
case B  M=6 P3 v=mid       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=1 (k=3)  deg w at ring: min 2 max 2
case B  M=6 K1,3 v=centre  fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=1 (k=4)  deg w at ring: min 2 max 2
case B  M=6 K1,3 v=leaf    fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=4)  deg w at ring: min 2 max 2
case B  M=6 P4 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=4)  deg w at ring: min 2 max 2

real	1m32.691s
user	1m32.396s
sys	0m0.031s
```
Planted errors: `--plant-x` (x not silent in stage 1) is caught because deg w at the ring exceeds 2; `--plant-u4` (u not silent in stage 4) is caught because the isomorphism check fails:
```
$ for f in --plant-x --plant-u4; do taskset -c 0 $PY planting_forced.py 300 78 $f; done
NRUNS=300 seed=78 plant_x=True plant_u4=False
case L>=1 M=3 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=0 (k=2)  
case L>=1 M=3 P3 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=3)  
case L>=1 M=3 P3 v=mid       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=3)  
case L>=1 M=3 K1,3 v=centre  fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=4)  
case L>=1 M=3 K1,3 v=leaf    fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=4)  
case L>=1 M=3 P4 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.714 (M/(M+2)=0.600)  max srcdeg=2 (k=4)  
case L>=1 M=6 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=0 (k=2)  
case L>=1 M=6 P3 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=3)  
case L>=1 M=6 P3 v=mid       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=1 (k=3)  
case L>=1 M=6 K1,3 v=centre  fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=1 (k=4)  
case L>=1 M=6 K1,3 v=leaf    fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=4)  
case L>=1 M=6 P4 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=2 (k=4)  
case B  M=3 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=0 (k=2)  deg w at ring: min 2 max 4
case B  M=3 P3 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=3)  deg w at ring: min 2 max 3
case B  M=3 P3 v=mid       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=3)  deg w at ring: min 2 max 3
case B  M=3 K1,3 v=centre  fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=4)  deg w at ring: min 2 max 6
case B  M=3 K1,3 v=leaf    fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=4)  deg w at ring: min 2 max 3
case B  M=3 P4 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=2 (k=4)  deg w at ring: min 2 max 6
case B  M=6 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=0 (k=2)  deg w at ring: min 2 max 5
case B  M=6 P3 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.818 (M/(M+2)=0.750)  max srcdeg=2 (k=3)  deg w at ring: min 2 max 5
case B  M=6 P3 v=mid       fail(a,b,c)=[0, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=1 (k=3)  deg w at ring: min 2 max 5
case B  M=6 K1,3 v=centre  fail(a,b,c)=[0, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=1 (k=4)  deg w at ring: min 2 max 7
case B  M=6 K1,3 v=leaf    fail(a,b,c)=[0, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=2 (k=4)  deg w at ring: min 2 max 7
case B  M=6 P4 v=end       fail(a,b,c)=[0, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=2 (k=4)  deg w at ring: min 2 max 5
NRUNS=300 seed=78 plant_x=False plant_u4=True
case L>=1 M=3 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=0 (k=2)  
case L>=1 M=3 P3 v=end       fail(a,b,c)=[10, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=3 (k=3)  
case L>=1 M=3 P3 v=mid       fail(a,b,c)=[7, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=3)  
case L>=1 M=3 K1,3 v=centre  fail(a,b,c)=[13, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=4)  
case L>=1 M=3 K1,3 v=leaf    fail(a,b,c)=[12, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=3 (k=4)  
case L>=1 M=3 P4 v=end       fail(a,b,c)=[11, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=3 (k=4)  
case L>=1 M=6 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.818 (M/(M+2)=0.750)  max srcdeg=0 (k=2)  
case L>=1 M=6 P3 v=end       fail(a,b,c)=[6, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=3 (k=3)  
case L>=1 M=6 P3 v=mid       fail(a,b,c)=[6, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=1 (k=3)  
case L>=1 M=6 K1,3 v=centre  fail(a,b,c)=[5, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=1 (k=4)  
case L>=1 M=6 K1,3 v=leaf    fail(a,b,c)=[5, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=4 (k=4)  
case L>=1 M=6 P4 v=end       fail(a,b,c)=[8, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=3 (k=4)  
case B  M=3 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=0 (k=2)  deg w at ring: min 2 max 2
case B  M=3 P3 v=end       fail(a,b,c)=[6, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=3 (k=3)  deg w at ring: min 2 max 2
case B  M=3 P3 v=mid       fail(a,b,c)=[6, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=3)  deg w at ring: min 2 max 2
case B  M=3 K1,3 v=centre  fail(a,b,c)=[13, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=1 (k=4)  deg w at ring: min 2 max 2
case B  M=3 K1,3 v=leaf    fail(a,b,c)=[12, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=4 (k=4)  deg w at ring: min 2 max 2
case B  M=3 P4 v=end       fail(a,b,c)=[17, 0, 0]  min L/D=0.667 (M/(M+2)=0.600)  max srcdeg=4 (k=4)  deg w at ring: min 2 max 2
case B  M=6 K2             fail(a,b,c)=[0, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=0 (k=2)  deg w at ring: min 2 max 2
case B  M=6 P3 v=end       fail(a,b,c)=[5, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=3 (k=3)  deg w at ring: min 2 max 2
case B  M=6 P3 v=mid       fail(a,b,c)=[4, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=1 (k=3)  deg w at ring: min 2 max 2
case B  M=6 K1,3 v=centre  fail(a,b,c)=[5, 0, 0]  min L/D=0.800 (M/(M+2)=0.750)  max srcdeg=1 (k=4)  deg w at ring: min 2 max 2
case B  M=6 K1,3 v=leaf    fail(a,b,c)=[7, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=3 (k=4)  deg w at ring: min 2 max 2
case B  M=6 P4 v=end       fail(a,b,c)=[10, 0, 0]  min L/D=0.778 (M/(M+2)=0.750)  max srcdeg=3 (k=4)  deg w at ring: min 2 max 2
```

---

## 4. Verdicts

**Claim 1 (existence.md): HOLDS.** I found no false step and no step that is unjustified as written. Optional repairs, none of which changes the argument:
1. In the good-history section, "lower tail at kM is at most e^{-4s} for large s" is true, but "large" means s >= s0(k,M,eps) ≈ 10·M²/eps when k = 6M/eps (≈ 2.47·10^6 at eps = 0.01, M = 50). At s = 20 and 50 the bound is vacuous. It would be worth stating s >= C M²/eps explicitly.
2. In stage 4, "degree at most k+1" can be sharpened to k. This is harmless.
3. The weight statement in §3 is best written as the conditional inequality E[X_u | F_{t+h}] >= e^{-s/M} E[C_T(s) Π D_v/(D_v+1)] on S_u, which is what is actually used.

**Claim 2 (sample-exponent.md): HOLDS.** I found no false step and no step that is unjustified as written. Optional repairs:
1. "existence and positivity of μ": positivity is not proved and is not needed, because the μ = 0 case is handled. Drop the word.
2. The re-clocking at the stage-3 stopping time is stated briefly ("can be sampled at selection, with independent increments"). It is a correct use of the strong Markov property of the clock family, and it is essential. My test shows that with the naive event "u picks a leaf" the flags at adjacent sites really are dependent (cov = -1.06·10^-3, z = -21.9). The proof does not use that naive event.

Remaining limits of this check:
- The numerical tests cover small cases only: trees with up to 17 vertices, seeds with up to 6 vertices, and s <= 2. The verdicts rest on the line-by-line verification in sections 1–2, and the numerics corroborate them.
- Neither proof identifies the value of μ, and neither claims to.
