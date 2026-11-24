# Independent review: standalone Claim 7

**Recommendation: accept.** Both mass-conservation proofs work, and the label-TV and recipient-degree conclusions follow. I found no mathematical or systematic verification failure in the standalone claim. The proof is correct specifically because the localized process includes the crossing step, freezes the denominator, and charges *all* late promotions to dominate the late-label tail. Reviewed writer blob `0a759e72e238a35b43586c35e436385970e0deab`; writer file untouched.

## Exact inputs

- Both proofs use the labelled founder-edge random friend tree, the conditional oriented-edge weight `1/(n D_source)`, monotonicity of each born label's degree, and unit jumps of `C` starting with `C_2=0,C_3=1`. For `n>=3` every leaf has a core neighbour and the core induces a tree. No imported local limit is required: the bounded nonnegative supermartingale theorem applied from each deterministic birth time supplies the simultaneous `Z_v`.
- **Mean route:** only `E C_n <= K n^alpha`, with deterministic finite `K` and `alpha<1`, is imported. The upper half of Theorem 1 in `proofs/bracket-motifs.md` (`proofs/bracket-motifs.md#theorem-1-finite-size-bounds-for-the-mean-core-size`) matches the model, range `n>=3`, prefactor `14606753031/2989767809`, and exponent `1133/2000` exactly. Its lower bound, the existence of a mean exponent, and the motif verifier itself are outside this review. The polynomial estimate follows, for instance, from `sum_{j=3}^{n-1}1/j <= log(n/2)`; dyadic Markov–Borel–Cantelli and monotonicity extend the resulting core-fraction limit to all `n`.
- **Pathwise route:** only `C_n=O(n^beta)` a.s. for a deterministic `0<beta<1` is imported. [Addario-Berry et al., arXiv v1, Proposition 6.22](https://arxiv.org/html/2403.20185v1#S6.Thmthm22) states `n^(-lambda) X_n^{>=2}->0` a.s. for `1-(3-2sqrt(2))/2 < lambda < 1`, where `X_n^{>=2}=C_n`. This is the proposition-level `sqrt(2)-1/2` threshold, not the stronger `0.9` appearing elsewhere in the source. The proposition states its coupling conclusion with details omitted; this review checks its use as an input and does not recertify that proof. The alternate argument needs no integrability of the random constant in the pathwise bound. Neither route uses a hub-positivity theorem, leaf-neighbour convergence, or a recipient moment or empirical-law theorem.

## Promotion and stopped drift

Write `d_v(r)=D_v(r)/r`. For each born label, the exact conditional drift to size `r` is `(p_{r-1}(v)-d_v(r-1))/r<=0`. Conditioning on `sigma_v=r` uses the post-step sigma field `F_r`; the initial cost on that event is **`2/r`**, not `1/r`, `1/(r-1)`, or a pre-step weight with the promotion outcome already known. For each possible `r<=n`, pull the `F_r` indicator out of the conditional expectation and iterate the fixed-label supermartingale. Summing finite terms gives the mean tail; all `v>m` must satisfy `sigma_v>m`, since a newborn is a leaf. The converse is false: an old leaf can become core after `m`. In particular the first possibly contributing late-label promotion is actually `m+2`, but the `m+1` term remains a valid overcharge. Because `m>=3`, the bound `E[Delta C_j|F_j]<=C_j/j` never encounters its false founder-edge extension. For the final factor in (7), the sum is bounded by `sum_{j=m}^infty j^(alpha-2) <= m^(alpha-2)+m^(alpha-1)/(1-alpha)`; the factor of two is necessary already at the entry time.

For localization let `A_{r-1}` be the already-counted labels at size `r-1`; this random finite set is `F_{r-1}`-measurable. With `I_r=1_{tau>=r}`, the exact increment at `r>m` is

```
X_r-X_{r-1} = I_r sum_{v in A_{r-1}} [d_v(r)-d_v(r-1)]
               + 2 (B_r-B_{r-1})/r.
```

Here `I_r=1_{tau>r-1}` is predictable, and the counted-degree part stays constant once stopped. Consequently `X_n - 2 sum_{r=m+1}^n (B_r-B_{r-1})/r` is the supermartingale; `X_n` itself need not be one. All quantities are integrable at every finite size (in particular `X_n<=2`). This proves (8) without a stopping theorem at an unbounded random time and without conditioning the degree drift on future survival. If `tau<=m`, the charge and `X` both vanish. At a crossing `C_tau<=K(tau-1)^beta+1`, including precisely the unit overshoot. In the stated `K>=1,beta>0` domain `tau=3` cannot occur, but the writer's treatment of it is harmless. The Abel sum in (112)–(116) has the two correct boundary denominators `n` and `m+1`; dropping the negative boundary and integrating `x^(beta-2)` from `m` gives the stated constant.

On `{tau=infty}`, compare `X_n` to the late-*label* tail by inequality, not equality. In Fatou use the nonnegative random variables `1_{tau=infty} X_n` (or the indicated tails); `E[1_{tau=infty} X_n]<=E X_n`. There is no assertion that the drift stays negative under the conditional law given `{tau=infty}`. The finite-cutoff residual `1-sum_{v<=m}Z_v` is in `[0,1]` on this event. Taking `m->infty` at fixed `K` forces the residual to vanish there, and the increasing union over integer `K` has full probability. Mean route: the same finite-cutoff subtraction, bounded convergence for the finite block and then monotone convergence over labels prove conservation. Pointwise convergence by itself would not suffice in either route.

## Leaf versus core mass, and the degree law

At `n>=3`, core degree sum is `n+C_n-2` and leaf degree sum is `n-C_n`; the latter divided by `n` tends to **one**, despite each current leaf having degree fraction `1/n`. This is not recipient mass. In fact `ell_n=P(W_n is a leaf|T_n)=(1/n) sum_core L_u/D_u<=C_n/n`. A useful cross-check stronger than (9): `p_n(v)<=q_n(v)` on the core and `q_n(v)=0` on leaves, so the exact label L1 difference is `(C_n-2)/n+2 ell_n`, implying `(3C_n-2)/n`. The factor `1/2` converts L1 to TV only for the two probability vectors in the claim; `q_n` need not have mass one. In a six-vertex double star, `ell_n=2/9`, while uniform vertex sampling has leaf probability `2/3` and the L1 difference is `4/9>C_n/n=1/3`.

Conditional recipient-degree law is the pushforward of `p_n` under the time-dependent map `v -> d_v(n)`. Once `sum Z=1`, coordinate convergence implies label L1 convergence by truncation. For any continuous `f`, replacing `p_n` by `Z` costs at most `||f||_infty ||p_n-Z||_1`; then dominated convergence against the fixed summable weights `Z_v` moves the atom locations. The pathwise argument works for all `f` on the same event and establishes weak/Prokhorov convergence, not degree-law TV or almost-sure convergence of the sampled sequence. For example the uniform-vertex normalized-degree law instead tends to `delta_0`; the recipients in the claim have no limiting atom at zero. Conditional and unconditional mean corollaries use boundedness and need no further input.

## Fresh adversarial controls

`standalone_7_sol.py` independently enumerates labelled recipient histories through size seven using exact fractions and no writer imports. It conditions on each *post-promotion prefix*, checks pre-step drift prefix by prefix for two threshold stops, an early history stop, a recipient-dependent stop and no stop, and verifies the finite mean tail and recipient degree identities. The controls expose concrete failures of incorrect variants:

- On the size-four path `3--1--2--4`, label 2 has just been promoted and has recipient probability `3/8`. At the next step, stop on any recipient other than 2. If one omits the crossing step using this outcome, the conditional normalized old-term drift is **`+3/80`**. Conditioning on survival instead makes the old term `3/5 > 2/4`. The claim uses neither modification.
- On the same path, promoting label 3 at size five crosses `C_r>sqrt(r)` with `C_5=3`. Both tracked promotions after cutoff `m=3` have labels `<=3`, and the late-label tail is zero. The stopped sum is `4/5`, not `4/7` at horizon seven. Omitting the overshoot or treating old labels as synonymous with early promotions loses this branch.
- A newly promoted label at size five contributes `2/5` on that event. Charging `1/5`, or applying the `1/4` pre-promotion weight after selecting the promotion outcome, fails. The double-star check separates vertex counts, degree mass and recipient probability and falsifies the underfactored `C_n/n` bound for the label L1 error.

Run from the repository root:

```text
$ /usr/bin/time -f 'wall=%e rss=%M exit=%x' timeout 8s taskset -c 13 nice -n 19 python3 -B round2/checks/standalone_7_sol.py
PASS: 584 conditional promotion bounds, 3060 prefix drifts, 20 horizon budgets; live crossing sizes [5, 6, 7]
PASS: exact label/degree/TV identities and mean tail at sizes 3..7
MUTANTS: anticipative strict stop drift +3/80; moving denominator 4/7 != 4/5; entry 1/5 < 2/5; postselected survival gain 3/5 > 2/4; late promotion may be old-label
DOUBLE STAR: P(leaf recipient)=2/9, uniform leaf mass=2/3, label L1=4/9 > C/n=1/3
wall=3.18 rss=13952 exit=0
```

The finite checks falsify shortcuts; the asymptotic results rest on the proof and the stated input for either route. An optional clarity edit is to display the compensated supermartingale explicitly. No change is needed for correctness.
