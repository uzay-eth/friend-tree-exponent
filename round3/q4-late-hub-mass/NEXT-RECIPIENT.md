> **Public checkpoint export note:** Edited export of `late-mass-astra/NEXT-RECIPIENT.md` (frozen original SHA-256 `0cd50c23647af388b9e84e8663d5e70d1a3b71a0023743a1d4aea1f7a7f3c87c`). Changes are limited to local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.** Corresponds to audited revision SHA-256 `0cd50c23647af388b9e84e8663d5e70d1a3b71a0023743a1d4aea1f7a7f3c87c` (distinct from earlier unreviewed pin `646ca34a...`). See [`REVIEW-NEXT-RECIPIENT.md`](REVIEW-NEXT-RECIPIENT.md).

# The next-recipient center has a different Gaussian threshold

Astra addendum to local snapshot 16f0c12. Start with the original edge seed. At size n and label cutoff 2≤m≤n let p=p(m,n) be the conditional probability that the next attachment goes to a label >m. Write R_m for the terminal share of those labels, C_m for the core count at size m, X=lim C_m/m^μ, and c for the late-share multiplier. Set h=M_0/μ, the limiting number of positive-share labels through m divided by C_m. Let b=E_J[1/D_root(∞)] for the inert virgin J (a degree-two root, one ordinary leaf, one silent exterior slot), and put κ=1−b.

For deterministic m→∞, m=o(n), the conditional law of (R_m−p)/√(p/n) converges almost surely to standard normal **if and only if**

    m^(1+μ)/n → 0.                                                   (1)

At a finite positive limit λ of the ratio in (1), its conditional law tends almost surely, in Kolmogorov distance, to

    N(κ h √(λ X/c), 1).

Its conditional mean tends to the displayed shift. When the ratio tends to infinity the standardized error escapes to plus infinity in conditional probability almost surely. Thus the log-squared critical window of the observed excess-degree center is still in the centered regime here.

More precisely, for every deterministic m≥n^ρ with any fixed ρ>0 and n/m→∞,

    E[R_m | ℱ_n] − p = (C_m/n)(κ h + o_ω(1)).                      (2)

The leading conditional bias separates into three contributions at scale C_m/n:

| Recipient of the gain | Coefficient |
| --- | ---: |
| An early label, including a later promotion of an early leaf | 0 |
| A young leaf promoted by its early host | h/2 |
| The young core endpoint of an existing crossing edge | h(1/2−b) |

Here “existing” in the last row refers to the time of that gain. The last row also splits according to whether the edge already crossed at the observation. The two coefficients are given below. The positive-label count H_m/C_m→h and C_m/m^μ→X are as recorded in the `POSITIVE-LABELS.md` (reviewed Q4 note, not exported) and `PROOF.md` (reviewed Q4 note, not exported). The positive value κ was bounded between 0.8833499 and 0.9347984 in the separate inert calculation. Only the inert constant and its definition are needed here.

## Exact finite increment

At size j≥m denote physical degrees by d_w, incoming rate numerators by r_w=Σ_{z~w}1/d_z, and the number of neighbours of w with labels ≤m by e_w. Let w_j be the recipient of birth j+1 and I_j its indicator of being late. The unnormalized late rate is H_j=j p(m,j)=Σ_w (d_w−e_w)/d_w. If w receives, the change of H is

    H_(j+1)−H_j = I_j + g_j(w),       g_j(w)=e_w/[d_w(d_w+1)].      (3)

The old late neighbours of w lose incoming rate, the newborn contributes its parent's reciprocal new degree, and w receives an extra unit if late. In particular p(m,j) is a bounded submartingale. Its limit is R_m for fixed cutoff by the finite-tree convergence. Subtracting p(m,n) and summing gives the physical identity

    R_m−p(m,n) = M_(m,n) + B_(m,n),
    M_(m,n)=Σ_{j≥n} (I_j−p(m,j))/(j+1),
    B_(m,n)=Σ_{j≥n} g_j(w_j)/(j+1) ≥ 0.                            (4)

The centered martingale series converges in L² and has conditional mean zero. B_(m,n) is pathwise nondecreasing in the cutoff m at fixed n. This monotonicity concerns g_j on one physical path. It does not require monotonicity of the martingale. The finite one-step mean of g is

    j E_j g_j(w_j)
      = (1/2) Σ_{early cores x} L_x/d_x
        + Σ_{early cores x} d_x(m) r_x/[d_x(d_x+1)]
        + Σ_{crossing late cores y} r_y/[d_y(d_y+1)],               (5)

where L_x is the number of present leaves of x. All leaves with an early neighbour occur in the first sum, including early leaves. Formula (2)'s age split instead follows the actual recipient labels.

For an early recipient x, e_x=d_x(m), since no new early neighbour can appear. Its unweighted charges after any observation telescope as d_x(m) times the decrease in 1/d_x, hence total at most one. A late recipient has e=0 unless its older parent x is early. Its first gain gives exactly 1/2; thereafter its charges telescope as decreases in its own reciprocal degree. Charge both kinds at that young endpoint to its early parent x. Early leaves that never promote have no charge. These relations include all cross-edge births and repeated sends.

## Averaging near the observation

Write T_m for the size epoch of the rate-one source clock model, L=log(n/m), and measure elapsed time s from T_m. We use the following claim, proved in the [averaging argument](NEXT-RECIPIENT-AVERAGING.md). For any fixed real v, sum all the charges g in (3) at elapsed times s>L+v with weight exp(L−s). Almost surely the sum divided by C_m tends to

    exp(−v) h [ 0 + 1/2 + (1/2−b) ].                                (6)

The three terms converge separately. In the proof one compares productive old hosts with rate-one sources continually promoting new inert virgins. This comparison controls late gains on each side of the cutoff, including old hosts that freeze and young branches that ever receive another exterior send. It does not rely on the previous asymptotics of the excess-degree drift.

Uniformly in j≥n, T_(j+1)−T_m−log((j+1)/m)→0 a.s. The first charge in (4) occurs strictly after T_n. Sandwich its cutoff time between L−ε and L+ε and its weights n/(j+1) between exp(L−s)exp(±ε). Apply (6) and send ε↓0. Thus nB_(m,n)/C_m→κh a.s., separately for the components. The averaging argument also gives a deterministic K with

    P(n B_(m,n)/C_m>K, C_m≥m^(μ/2)) = O_a(n^(−a))                 (7)

for every a>0 and m≥n^ρ. To pass to conditional L¹, the size gate belongs to ℱ_n and eventually holds a.s. By Markov and Borel–Cantelli, the conditional probability of the exceptional event in (7) is eventually ≤n^(−4) a.s. Also E_n[(nB/C_m)²]≤4n² from (4), since 0≤p,R≤1 and E_n M²≤1/n. Cauchy–Schwarz controls the conditional exceptional tail. For the capped variables use the bounded full-path tail supremum of their absolute error from κh: conditioning it on ℱ_n and sending its initial index to infinity proves almost-sure conditional L¹ convergence. The same argument applies to each nonnegative part. This proves (2).

More explicitly, let Γ(du)=E_J[−d(1/D_root(u))] denote the expected positive Stieltjes measure of **post-promotion** reciprocal drops in the inert virgin; its mass is 1/2−b. A virgin begins at promotion of the young leaf. At unit stationary immigration rate, a root gain of age u after observation has promotion time either before or after the observation. Integrating exp(−s) over observation-forward times s≥0 divides its mass into factors 1−exp(−u) for s<u and exp(−u) for s≥u respectively. Accordingly the part from crossing edges already present at n has limiting coefficient h∫(1−exp(−u))Γ(du); the part from edges promoted after n has coefficient h∫exp(−u)Γ(du), in addition to their initial promotions h/2. The same stationary averaging proof applies with these indicators. To replace the deterministic dividing epoch L by T_n, the only extra discrepancy comes from edges born in [L−ε,L+ε]; at productive hosts this costs at most one per such birth, with mean 2ε per host and exponential concentration. The discarded hosts obey the same square bound as before. Per productive host the future post-promotion contribution from new crossings belongs to [0.1191026,0.1310137], while that from crossings already present belongs to [0.2523362,0.3156958]. These are expansions in conditional L¹ for the actual correction B, so no independence with M is needed. `NEXT-RECIPIENT-EVIDENCE.md` (author replay note, not exported; intervals independently certified in [`REVIEW-NEXT-RECIPIENT.md`](REVIEW-NEXT-RECIPIENT.md)).

## The fluctuation scale and its converse

For every diverging sublinear cutoff the finite-tree estimates give p~q~R_m~cX m^(μ−1), where q is the observed excess-degree proportion. The previously proved **compensated-noise** theorem gives the conditional N(0,1) limit for M/√(q/n), hence also for M/√(p/n). This input is solely about the identical centered sum of hit indicators in (4). See `ADDENDUM-CUTOFF-THRESHOLD.md` eq. (6) (audited in [`REVIEW-CRITICAL-EXCESS-DRIFT.md`](REVIEW-CRITICAL-EXCESS-DRIFT.md), not separately exported) and `FINITE-TREE.md` (reviewed Q4 note, not exported).

When m^(1+μ)/n→λ∈(0,∞), C_m/√(np)→√(λX/c). Equation (2) and its conditional L¹ form for B give the stated shifted law, including convergence of the conditional mean. The observable bias correction κh C_m/n is accurate to o_ω(√(p/n)) at this scale.

For all other sublinear cutoffs, compare (4)'s positive B on the same path with m_a=floor(a n^(1/(1+μ))) for any fixed a>0. The statistic B/√(p/n) is nondecreasing in m whenever the denominator is positive, since p(m,n) is nonincreasing. It is eventually positive along the cutoffs at issue. If (1) holds, m≤m_a eventually for every a>0. Let a↓0 through rationals and apply the conditional mean result at m_a to obtain a vanishing conditional L¹ correction. If the ratio diverges, m≥m_a eventually for every a; let a↑∞ and use conditional concentration at m_a to obtain escape to plus infinity. If the ratio has positive limsup, a subsequence dominates one fixed m_a. Its correction then exceeds a fixed positive path-dependent number with conditional probability tending to one. The conditional CDF at zero is bounded above in limsup by the standard Gaussian CDF at the negative of that number, using the noise CLT. This excludes convergence to standard normal. These comparisons prove (1) for arbitrary deterministic diverging sublinear cutoffs.

The previous critical-excess-degree file and its numerical enclosure are preserved at the supplied hash. The finite increment and numerical checks for this addendum are in `NEXT-RECIPIENT-EVIDENCE.md` (not exported). The next-order bias after subtraction of κh C_m/n, on faster cutoffs where it would matter for forecasting at Gaussian accuracy, remains open.
