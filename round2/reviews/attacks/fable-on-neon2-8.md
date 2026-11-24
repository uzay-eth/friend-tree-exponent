VERDICT: HOLDS

# Fable on Neon-2 claim 8 (all-label mean degree-share order E Z_1 = E Z_2 > E Z_3 > ...; per-tree two-birth gap at least delta*/(n(n+1)(n+2)), delta* = 2493671/62500000)

**Scope of the verdict.**
- Attacker: the Fable lead, Claude Opus 5.5. Pin: `834bb96` of the claimant's branch. The `mean-order` directory is unchanged from that pin to the claimant's present head (`git diff --stat` is empty).
- Re-derived by me: the analytic skeleton (section 1).
- Not re-derived by me: the six finite-state tables (`v, w, U, phi2, phi1, phi1_x`; 13,299 rows) and their hand-offs above the cutoffs. These were audited, replayed and attacked with planted errors by a GPT-6.2 Sol Ultrafast child of this lead (report `round2/checks/neon2-claims8-9-audit-gpt-family.md`, "Claim 8 (`834bb96`): HOLDS"). The claimant's own cold audit is also GPT-family. Astra's verdict is in `astra-on-neon2-8.md`.

## 1. Steps re-derived by the lead

**Yule embedding.**
- At size m the holding time is Exp(m) whatever the shape. So the jump chain, and with it every Z_v, is independent of the size amplitude S = lim e^{-t} N(t), which has law Gamma(N, 1) from size N.
- W_v = lim e^{-t} D_v(t) = S Z_v. Hence M(T, v) := E_T W_v = N E_T Z_v, and e^{-t} M(T_t, v) is a closed martingale.

**The sandwich (L_v − 1)_+ ≤ M ≤ D_v.**
- Upper bound: v gains at rate r_v = Σ_{x~v} 1/D_x ≤ D_v. So e^{-t} D_v(t) is a supermartingale, which gives M ≤ D_v.
- Lower bound, for g = (L − 1)_+. The leaf stock L rises by one at rate at least L, since each leaf neighbour sends at rate 1. It falls by one at rate L/D_v ≤ 1, since only v can promote its leaves.
  - For L ≥ 2: the drift is at least L − L/D_v ≥ L − 1 = g.
  - For L = 1: the drift is at least 0 = g.
  - For L = 0: g stays at 0.
  So e^{-t} g(T_t) is a submartingale. Since g ≤ D, uniform integrability from the Yule moments gives g ≤ M.

**Signed identity and endpoint split.**
- The identity is the same as in Astra claim 26, which I re-derived there: twins cancel, nonadjacent pairs cancel, and an edge {x, y} carries the coefficient s_y − s_x.
- Expanding (s_y − s_x)(M_{x|y} − M_{y|x}) and dropping s_y M_{x|y} + s_x M_{y|x} ≥ 0 leaves −s_x M_{x|y} − s_y M_{y|x}. This is valid for every tree.
- Weight sum: Σ_x r_x = n, since each caller distributes total mass one, and D_x + 1 ≤ n. So Σ_x r_x/(D_x + 1) ≥ 1.

**Consistency with my own certified mean bands** (Fable claim 3, REPORT.md; grouped exact upper bounds from `exact.log`).
- The claimed founder gap m_2 − m_3 ≥ 409872349/3000000000 = 0.136624 lies below my grouped upper bound 0.252499170. It is also below my own lower bound 0.1499925, so the founder bound is weaker than mine but consistent.
- The claimed later gaps delta*/((u−1)u(u+1)) are below my grouped upper bounds .052689658, .025582202, .015544936, .010622993 and .007833101 for u = 3, ..., 7 (child's comparison, checked by me against the same table).
- Comparison with Astra claim 27 at p = 1: Astra's constant 9/5000 multiplies Gamma(1+p)Gamma(n)/Gamma(n+p+2) = 1/(n(n+1)(n+2)), the same scale as delta*. Both are lower bounds, so there is no contradiction; Neon-2's constant is larger.

## 2. Checked by the GPT-family child only

- The finite tables and their boundary hand-offs (degree caps 40, 14 and 8, exit faces, and the parent-tracked variants).
- Optional stopping with finite-exit and zero no-exit payoffs.
- The replay: `ALL_EXACT_CHECKS_PASSED`, with 10 out of 10 bundled controls rejected.
- A planted one-tick change in the parent-exit cell (9, 2, 1): rejected.
- Exhaustive screens of all 47 free trees with n ≤ 8, using two births plus one or two more and the child's own moment bounds: no upper bound below the claimed floor. The child notes these upper bounds are loose, so this screen catches gross errors only.

Not covered: stochastic dominance; moments with p < 1 (see claim 9 and Astra claim 27).
