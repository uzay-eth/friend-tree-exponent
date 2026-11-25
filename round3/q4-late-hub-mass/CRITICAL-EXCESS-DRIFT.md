> **Public checkpoint export note:** Edited export of `late-mass-astra/CRITICAL-EXCESS-DRIFT.md` (frozen original SHA-256 `4c27a043e261f9c2f3e9db1dd529dfce70cce0c8f548e875039fb85e262fa4e8`). Changes are limited to local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

# Critical conditional bias for the observed excess-degree tail

Astra addendum following local snapshot b404108. Fix a deterministic cutoff sequence 2≤m=m(n)<n with

    m^(1+μ) (log n)²/n → λ ∈ (0,∞).

For the original edge-seeded friend tree let q_n=Σ_{m<v≤n}(d_v(n)−1)/(n−2), R_m=Σ_{v>m}Z_v, C_m be the number of cores at size m and X=lim C_m/m^μ. Write c for the reviewed late-mass multiplier and M_0 for the total positive birth-mark intensity. All logarithms are natural. Let J be the archived inert virgin: its root has one ordinary leaf and one exterior parent slot. The slot never sends, is counted in the root degree, and a root choice of the slot has no internal effect. Put

    b = E_J[1/D_root(∞)],       1/∞ = 0,
    δ(X) = (1−b) M_0/(1+μ) · √(λ X/c).

Then almost surely the conditional distribution given the observed size-n history satisfies

    sup_z |P((R_m−q_n)/√(q_n/n) ≤ z | ℱ_n) − Φ(z−δ(X))| → 0.       (1)

Here Φ is the standard normal CDF. In observable units the conditional bias is E[R_m|ℱ_n]−q_n=(1−b)t(m,n)/n+o_ω(√(q_n/n)), where t(m,n) counts the observed crossing core edges. In particular (E[R_m|ℱ_n]−q_n)/√(q_n/n)→δ(X) almost surely. The Gaussian noise has conditional variance one in the limit law. Its positive shift retains the random core-growth amplitude. These statements use q_n as their center throughout.

Let t=t(m,n) be the observed number of core edges crossing the label cutoff, and L=log(n/m). For the integrated drift V_n defined below, the two separate assertions that identify the shift are

    t/(C_m L) → h=M_0/μ                                      a.s., (2)
    E[|n V_n/t − (1−b)| | ℱ_n] → 0                           a.s.  (3)

Ratios at the finitely many zero denominators may be defined arbitrarily. In fact (1) allows the fully tree-dependent shift (1−b)t/√(nq_n) in place of δ(X). We prove (2)–(3) directly, rather than extracting a coefficient from the earlier two-sided inequalities.

## Inputs and scales

The reviewed notes give C_s/s^μ→X>0, R_m~cX m^(μ−1), q_n~R_m for diverging sublinear cutoffs, and H_m/C_m→h for H_m=#{x≤m:Z_x>0}. The marked rectangle for labels through m gives #{x≤m:Z_x>r/m}/C_m→h(r), with h(r)=∫₀¹ v^(μ−1)ν((rv,∞))dv increasing to h as r↓0. Convergence on a countable sequence of positive r suffices. The limiting uniform-core degree law is proper. `PROOF.md` (reviewed Q4 baseline note, not exported), `POSITIVE-LABELS.md` (reviewed Q4 note, not exported), `FINITE-TREE.md` (reviewed Q4 note, not exported).

The number K_m of labels through m that **ever** become cores, including current cores, is O_ω(C_m). The independent old-leaf-promotion bounds of SMALL-POWERS.md bound K_m itself on dyadic starts; interpolate by monotonicity of K_m and regular variation of C_m. In continuous time every vertex has a rate-one source clock and a uniform neighbour mark. At size epoch T_s, uniformly for v≥s→∞, T_v−T_s=log(v/s)+o(1) a.s. The physical degree dichotomy says a vertex of unbounded degree has positive limiting share. It holds simultaneously for the countably many physical labels. The inert virgin is used solely as a comparison process with its parent slot retained. `SMALL-POWERS.md` (reviewed Q4 note, not exported), [archived degree dichotomy](../../round2/proofs/astra-late/source/round2/open/age-order/late/README.md), [archived inert dynamics](../../round2/proofs/astra-late/source/round2/open/forest-age/inert/README.md).

The critical hypothesis gives log m/log n→1/(1+μ), L/log n→μ/(1+μ), and C_m≥m^(μ/2) eventually. We use countably many deterministic gates K_m≤B C_m, B integer. The gate need not be known in advance: for martingale estimates stop at the first time the **current** count of early cores exceeds B C_m. On each typical path some B suffices forever for all large starts. For discrete hit counts with total predictable mean ≤B C_m O(L), the difference from that mean is o(C_m L) a.s. by bounded-jump Bernstein and Borel–Cantelli. Even an error C_m is valid since exp(−a_B C_m/L) is summable on the size gate. The same estimate for a continuous counting process follows from the exponential bound for a count minus its compensator, obtained by exponentiating the count at parameter u and subtracting (exp(u)−1) times the integrated rate. Conditional at T_m, C_m is a known number. These estimates need no independence across n.

## Upper bound for the crossing count

A new crossing occurs exactly when an early core x selects a young leaf attached to x, promoting it. Put d=max(3,ceil(√L)). Register an early label at the first **degree gain after size m** that takes its degree to at least d. At that gain it has a fresh young leaf. Before registration any one early label creates at most d crossing edges. Indeed before any post-m gain it has no young neighbours, and until a qualifying gain its degree is <d. Each distinct crossing consumes a different young neighbour. This includes early labels promoted after size m. In total these early crossings cost at most d K_m.

At registration start a private linear birth-death process from one with rates j and j/d. Its population can be kept below x's physical leaf count. Use its births for true leaf-source births and supplement at the remaining physical leaf rate. At an x ring take a virtual death for a uniform mark ≤min(1,j/d) and add deaths at rate (j/d−1)_+. Order physical leaves first; when physical and virtual leaf counts coincide, the probability of losing a physical leaf is j/D_x≤j/d. When they differ, a loss of one cannot reverse their order. Zero is absorbing. The virtual drivers consist of the host's fresh source tape and private births and extra deaths. They can be assigned independently in registration order. Physical gains from other cores add leaves, distinct hosts have disjoint leaf stocks, and an old leaf promoted to a registered host later takes fresh clocks. The construction preserves all physical arrow intensities. It is the simultaneous version of the archived fresh-leaf degree-hitting test. [Prior single-tag comparison](../../round2/proofs/astra-late/source/round2/open/age-order/adjacent/README.md).

The virtual embedded walk has downward/upward probability ratio 1/d, hence extinction probability 1/d by gambler's ruin. On survival it is unbounded. Every registered label with a surviving comparison belongs to H_m by the physical dichotomy. Pad the registration tapes out to ceil(B C_m). Their number of extinctions is binomial with probability 1/d, regardless of which of these tapes are ultimately used. A conditional Chernoff bound shows it is ≤2 ceil(B C_m)/d eventually a.s. on the size gate. Therefore on K_m≤B C_m the number of registered vertices at every elapsed time is ≤H_m+O_B(C_m/d). Registration is adapted. Its criterion does not use survival; survival was used only to bound the count afterwards.

Post-registration crossings are at most all the rings of registered sources. The compensated count up to elapsed time L+1 has error O(C_m) on the gates, as above. Its integrated intensity is at most (H_m+O_B(C_m/d))(L+1) on the good event. The actual interval T_n−T_m lies below L+1 eventually. Adding the pre-registration cost dK_m gives limsup t/(C_m L)≤h. If comparison tapes are sampled separately for different starts, they may be joined conditionally on the physical path for this countable set of experiments. Equivalently apply each summable bound in its own coupling and transfer the resulting physical inequality. No physical hub decisions are treated as independent.

## Old cores that account for almost all the crossings

Fix 0<ε<1, r>0 and an initial degree cap D. Put j_0=floor(m exp(εL)) and let G be the set of early labels of degree ≤D at m and degree ≥r j_0/(4m) at j_0. This is measurable at j_0. Almost surely, eventually every early label with Z_x>r/m and starting degree ≤D is in G. Indeed restarting the continuous-time tree at j_0, the total terminal amplitude S has law Gamma(j_0,1), and S Z_x is dominated by a unit-rate Yule terminal Gamma(d_x(j_0),1) via the gain-rate bound r_x≤d_x. This is a marginal bound for each x; the union bound needs no independence. The probability of S<j_0/2 is exponentially small in j_0. For any initial degree below r j_0/(4m), the probability of the dominating gamma exceeding r j_0/(2m) is ≤exp(−a r j_0/m) by exponential Markov at 1/2. The union over m is summable in n since j_0/m grows polynomially. Thus liminf |G|/C_m≥h(r)−π_core(d>D).

Every x in G has at most D originally old neighbours. Its additional non-young-leaf neighbours can only be young vertices promoted at a ring of x. Almost surely the maximum number of rings of any one label through size n² is O(log n): for each label this is a sum of independent Bernoulli source choices with mean ≤1+2 log n; exponential Markov and a union bound over n² labels give a summable bound with a sufficiently large deterministic constant. Hence throughout sizes j_0,…,n the fraction of neighbours of x that are young leaves is 1−o(1), uniformly on G. At each birth the conditional chance of a crossing from G is the sum of those fractions divided by the current size. The counting estimate then gives, for t_G consisting of crossings from G **after j_0** through size n,

    liminf t_G/(C_m L) ≥ (1−ε)[h(r)−π_core(d>D)].                  (4)

This uses predictable fractions for the observable set G, including any nonhubs in it. First send D→∞ and r↓0, then ε↓0 on fixed countable parameter sets. Together with the upper bound, (4) proves (2). It also shows that crossings outside the selected cohort have arbitrarily small limiting fraction in this order of limits. This is stronger than a lower bound on the number of productive old hosts.

## Reciprocal degrees on the cut

At each crossing, the promoted endpoint y has degree two with its early parent x and one newborn leaf. Give its branch an independent inert virgin tape in global crossing order. Run its internal arrows in both copies up to the first later arrow x→y. A root send y→x creates a leaf at x in the physical tree and makes no internal change in the inert copy. All new sources inside the young component can use the same tapes. The cut tree has one early boundary edge per young core component, so different registered components have disjoint interiors for this cutoff. After an incoming x→y arrow the physical continuation may use fresh clocks. The unused part of the inert tape still defines its terminal mark. The tapes can be generated independently before launch and used progressively. In particular their marks need not be independent of the random number or birth times of crossings.

Let D_i^0(s) be the inert root degree at age s of tape i and U_i=1/D_i^0(∞). The U_i are iid terminal marks in [0,1/2], of mean b. We will show, on the physical law, that

    sup_{n≤j≤n²} | (1/t) Σ_{y crossing at n} 1/d_y(j) − b | → 0 a.s. (5)

First bound incoming arrows into the selected components of (4) up to size n². The host degree is ≥r j_0/(4m). Each eligible old host has O(log n) possible crossing core recipients by the source-count bound. At step j≥j_0 their total conditional incoming-arrow probability is at most

    B C_m · O(log n) · 4m/(r j_0 j)

on the gates. Stop if the observed early-core or individual-source-count gate fails. Summing to n² bounds the compensator by O_B(C_m (log n)² m/(r j_0))=o(C_m). A bounded-jump exponential bound makes the actual number O(C_m) eventually a.s. Each such arrow can invalidate at most one branch. This step allows repeated sends and covers the full future interval used in (5).

Fix an integer age S_0>0. The number of crossings created after floor(n exp(−2S_0)) and before n is O_{ω,S_0}(C_m) by the same gate and the bounded log-size interval. Every earlier crossing is at least S_0 old by T_n for large n. For an unspoiled selected crossing at least this old, at every j in (5) the difference between its reciprocal and U_i is between zero and V_i(S_0)=1/D_i^0(S_0)−U_i. This bounded variable decreases to zero as S_0→∞, so E V_i(S_0)→0 without any time-rate assumption.

Sum over the full prefix of t tapes, including the discarded indices when upper-bounding the sum of V_i. The prefix means of U_i and V_i(S_0) tend to b and E V_i(S_0) a.s. in these triangular couplings: for a known initial C_m, Hoeffding and a union over all integer prefixes ≥C_m give O(exp(−a C_m)) tails for each fixed accuracy. By (2), t≥C_m eventually. All excluded crossings together cost at most half their count in the reciprocal sum. Equation (4), the incoming-arrow bound and the age bound show that the limsup error in (5) is at most E V_i(S_0) plus half the discarded fraction from (4). Send the fixed cohort errors and then the age error to zero. The conclusion involves only the original path. In particular this argument does not assume that a bounded-age approximation is accurate at a prescribed deterministic rate.

## Integrate the crossing-edge process

At size j write A_j=Σ_{m<v≤j}(d_v(j)−1), Y_j=A_j/j, t_j=t(m,j), k_j for the number of early cores, and a_j=Σ_{early cores x} (# late neighbours of x)/d_x. Exactly,

    R_m−q_n = M_n + V_n − 2q_n/n,
    V_n=Σ_{j≥n} [t_j − Σ_{crossing at j}1/d_y(j) + a_j]/[j(j+1)].

Here M_n=Σ_{j≥n}(I_j−p_j)/(j+1), where I_j is the indicator of a late recipient and p_j its conditional probability. The preceding post-review proof gives E[M_n|ℱ_n]=0 and conditional N(0,1) convergence after division by √(q_n/n). It also gives a finite bound on the new crossings and early-source term:

    E[ Σ_{j≥n} (t_j−t+k_j)/[j(j+1)] | ℱ_n] ≤ 24 k_n/n.

This uses E[k_j|ℱ_n]≤12 k_n and new-crossing probability ≤k_j/j. Each crossing already present contributes its own 1−1/d_y(j), and each newer crossing contributes between zero and one. Also 0≤a_j≤k_j. `ADDENDUM-CUTOFF-THRESHOLD.md` eqs. (2), (3), (6), (8) (audited in [`REVIEW-CRITICAL-EXCESS-DRIFT.md`](REVIEW-CRITICAL-EXCESS-DRIFT.md), not separately exported).

Let B_n denote the bounded supremum error in (5), set to zero at t=0. Its conditional mean tends to zero a.s.: for n≥N use the bounded full-path majorant sup_{r≥N} B_r, martingale convergence of its conditional expectations, and then N→∞. The normalized tail of the summation weights past n² is ≤1/n. Since nΣ_{j≥n}1/[j(j+1)]=1, for t>0 this gives the finite comparison

    E[|nV_n/t−(1−b)| | ℱ_n] ≤ E[B_n|ℱ_n] + 1/n + 24k_n/t.        (6)

Since k_n≤K_m=O_ω(C_m), (2) now proves (3). Finally t/√(nq_n)→h [μ/(1+μ)]√(λ X/c). Use conditional L¹ comparison with the observable shift (1−b)t/√(nq_n), the conditional noise limit, and continuity of the Gaussian CDF. The subtraction 2q_n/√(nq_n) vanishes. This proves (1) and its conditional mean statement. There is no independence assumption between the realized drift and the martingale.

## Scope

The number-of-hubs input, the inert virgin and the one-tag high-degree survival bound have prior provenance as indicated above. The asymptotic intensity of crossings from a growing early cohort, the uniform reciprocal-degree replacement (5), and the integrated drift limit at a positive finite critical ratio are the contributions of this addendum. The proof provides an expectation defining b and does not require a closed form for it. A rational enclosure gives 0.8833499≤1−b≤0.9347984. The finite scheme converges to this coefficient at unrestricted precision. `CRITICAL-EXCESS-EVIDENCE.md` (audited in [`REVIEW-CRITICAL-EXCESS-DRIFT.md`](REVIEW-CRITICAL-EXCESS-DRIFT.md), not separately exported). All conclusions here concern the observed excess-degree center. No successor claim about a different center is needed.
