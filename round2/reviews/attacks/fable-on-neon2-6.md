VERDICT: HOLDS

# Fable on Neon-2 claim 6 (unoriented five-vertex certificate: 0.565702177913 <= mu <= 0.565748407438)

**Scope of the verdict:** Claimant: the Neon-2 lead (Gemini 4 Neon wrote proof, builder and verifier; its own GPT-6.2 Sol Ultrafast children audited). This attack: the Claude-family lead (Claude Opus 5.5) and a GPT-6.2 Sol Ultrafast child of this lead; both are a different family from the author. (1) Lead: replay from a clean archive of the pinned commit with the binary hashes checked; a reading of the whole proof file (section 2 below); planted errors and shifted exponents that the verifier must reject; the certified inequality tested on the lead's own 1,930 simulated trees with an exact evaluator. The evaluator used in (1) was written by the child (GPT) from the proof text and the binary layout, not by the claimant; the lead ran it on its own trees. (2) Child: line-by-line audit of the proof and of the verifier against it, its own evaluator, adversarial trees, further controls (report: Fable folder `round2/checks/neon2-claim6-audit-gpt-6.2-sol-ultrafast.md`). Not done by anyone here: a second implementation of the 191 million local row inequalities (the verifier was replayed, not rewritten), and a re-derivation of the exact rational Bernstein tail beyond running its script. The bracket gives four digits of mu (it separates 0.565700), not eight, as the claim says.

## 1. Replay (lead)

`git archive b04cc50 round2/neon2/quintet-cert` into an empty folder; the binaries `quintet_lower.bin` and `quintet_upper.bin` (read in place, not modified) have SHA-256 `81dc1f250d490837...` and `8a36c3afdd2a8f4b...`, equal to the board. `verify_quintet.cpp` built with the board's command (`-Wall -Wextra -Werror -pedantic`): no warning. `cert_check.py --threads 1` (628 s on one core): its output is **byte-identical to `certificate_gold1.json`**: lower nominal 56570237/10^8, error per core 19186 ticks, m_min 99882299906, certified 2825181729593748861/4994114995300000000 = 0.565702177913; upper nominal 14143707/25000000, error per core 19783, m_min 155236771435, certified 0.565748407438; both "ok". The child's separate replay with two threads is identical as well.

## 2. Reading of the proof (lead; `PROOF.md`, 134 lines)

* **Generator (section 1).** On the core tree with states (D_u, c_u): leaf gain at rate r_u = L_u + sum over core neighbours of 1/D_v (each leaf of u rings at rate 1 and must pick u; a core neighbour v picks u at rate 1/D_v), promotion at rate b_u = L_u/D_u (u picks one of its leaves, which becomes a core of state (2, 1); c_u increases by one, D_u does not change). A ring of u that picks a core neighbour is the neighbour's gain. Nothing else changes the core tree. Correct, including C(T) = 1 (r = D, b = 1); the singleton table includes c = 0 (index (D - 2)(D + 3)/2 + c) and the verifier has a separate branch for one-core trees.
* **Interpolation (2.1).** A vertex of degree D >= d_j is replaced by the mixture (d_j/D) T_j + (1 - d_j/D) INF with reciprocals 1/d_j and 0, so the mean reciprocal is 1/D exactly; the level-to-level projections compose. Correct. A gain at a tail vertex changes the weight by -(d_j/D)/(D + 1), which is the factor t gamma of section 3.3 with gamma = r_u/(D_u + 1).
* **Supports (2.2).** Every event increases a D or a c; the supports are upper bounds on D, c and sums of D of internal vertices (the lower bounds c >= 2, 3, 4 hold automatically for internal vertices), so an existing embedding can only leave the support. A newborn J has core degree 1, hence is always an endpoint; removing it from a new supported motif gives a supported motif of one size less by the strict nesting, to whose promotion term the creation is charged. Correct as stated; that the header validation enforces exactly these nestings was checked by the child.
* **Transfers (3.1).** Each transfer moves (1/D_w - r_1(u)) times a difference from one row to a row of the next size and is cancelled by its pair, so the sum over all embeddings is unchanged. The Fenchel-Young step for vertices outside M_3: the piecewise linear A_u moved between the edges at u cancels over ordered pairs; the remainder is a maximum of affine functions of 1/D_w on [0, 1/2], so attained at a knot, and it is convex in the expected difference, so Jensen passes it to the corners; A_u is affine on [0, 1/d_2], so its mean over a tail endpoint is exact. Correct.
* **Envelopes (3.2).** |d_fine - Pi m| <= Pi r on every fine state; with the donor coordinate independent of the others under the product interpolation, the displayed bound follows from the triangle inequality. Correct.
* **Tail rates (3.3).** For a tail vertex of exact class kappa: gamma = r_u/(D_u + 1) >= (D_u - kappa)/(D_u + 1) >= (d_j - kappa)/(d_j + 1) (increasing in D), b = (D - kappa)/D >= (d_j - kappa)/d_j; for the inner overflow (c >= K_j + 1) the remaining gain rate is at most (D - c) + (c - 1)/2 <= D - (K_j + 2)/2. Re-derived; correct. The rate at which a tail vertex u feeds an inner neighbour is 1/D_u = t/d_j, as used.
* **Aggregation and deduction (section 4).** With E_tot the largest positive per-core residual, sigma S (L Phi - lambda Phi) <= E_tot C, and with m_min C <= S Phi this gives L Phi <= (lambda_up + E_tot/m_min) Phi and L Phi >= (lambda_lo - E_tot/m_min) Phi. Comparability: every motif is allocated to a vertex of bounded degree (the internal vertices of the higher motifs have bounded degree, edges go to their child end), and sum over cores of (c_u - 2) = 2(C - 1) - 2C = -2 lets the lower certificate add alpha (c_u - 2) to each vertex. Stopped Dynkin with the Yule bound on N(t) and Gronwall give E Phi_lo(T_t) >= Phi_lo(T_3) e^(lambda*_lo t) and E Phi_up(T_t) <= Phi_up(T_3) e^(lambda*_up t), hence the bounds on the growth rate of E C(T_t), which equals mu by the existence theorem. Correct.
* Not re-derived by the lead: the analytic singleton tail for D >= 2000 (the Bernstein polynomials of section 3.4) beyond running the claimant's script in the replay; the correspondence between the verifier's loops and the rows of section 3 (child).

## 3. Planted errors and shifted exponents (lead)

The verifier does not stop at a violated row: it records the largest positive residual per row family ("caps"), and the certified endpoint is the nominal exponent minus (lower) or plus (upper) the error per core divided by m_min. So a planted error must show up as a much weaker certified endpoint. Runs of the pinned verifier (one thread, 45 to 206 s each) with its built-in perturbations (one table entry changed by 10^9 ticks) and with shifted nominal exponents; endpoints computed with the replay script's formula:

| run | family | error per core (ticks) | certified endpoint | authentic |
|---|---|---:|---:|---:|
| lower, Fenchel-Young knot | g | 27,216,791,542 | 0.293213735 | 0.565702177913 |
| lower, 5-path entry | P5 | 7,565,703,220 | 0.489956184 | |
| lower, 4-star entry | S4 | 9,565,709,648 | 0.469932552 | |
| lower, nominal 0.56570437 (2 ppm higher) | all | 4,624,094 | 0.565658075 | |
| upper, Fenchel-Young knot | g | 14,972,357,245 | 0.662196803 | 0.565748407438 |
| upper, 5-broom entry | B5 | 25,697,247,731 | 0.731284110 | |
| upper, singleton entry | f | 2,388,789,641 | 0.581136320 | |
| upper, nominal 0.56574281 (5.5 ppm lower) | all | 13,597,105 | 0.565830399 | |

Every planted error is caught in the family where it was planted, and asking for a tighter nominal exponent makes the certificate weaker, not stronger. The child also showed that changing one byte of a table on disk makes the unmodified verifier exit with status 2, and that the two exponents just inside the tightest tree ratios found (below) are rejected after the residual penalty.

## 4. The certified inequality on the lead's own trees

The child's exact evaluator (`eval.cpp` with `run_eval.py`; Phi by enumeration of all motif embeddings, the drift by brute force over every attachment vertex, exact rationals) reproduces the claimant's Phi and drift exactly on all 47 published concrete trees per side. Run by the lead on its own 1,930 snapshot trees (the same as for Neon claims 1 and 3), 74 s per side:

| time | trees | vertices | lower: drift/Phi (certified 0.565702178) | smallest margin | upper: drift/Phi (certified 0.565748407) | smallest margin |
|---:|---:|---|---|---:|---|---:|
| 3 | 800 | 7 to 217 | 0.5657023719 to 0.5660013 | 1.94e-7 | 0.5654155 to 0.5657482764 | 1.31e-7 |
| 5 | 600 | 43 to 1,563 | 0.5657023719 to 0.5658902 | 1.94e-7 | 0.5654132 to 0.5657482764 | 1.31e-7 |
| 7 | 300 | 379 to 10,576 | 0.5657023719 to 0.5657758 | 1.94e-7 | 0.5656524 to 0.5657482763 | 1.31e-7 |
| 8 | 150 | 1,153 to 32,164 | 0.5657048933 to 0.5657790 | 2.72e-6 | 0.5656410 to 0.5657436178 | 4.79e-6 |
| 9 | 60 | 2,434 to 91,508 | 0.5657037570 to 0.5657756 | 1.58e-6 | 0.5656660 to 0.5657464836 | 1.92e-6 |
| 10 | 20 | 13,879 to 131,832 | 0.5657085396 to 0.5657330 | 6.36e-6 | 0.5657095 to 0.5657402910 | 8.12e-6 |

No violation of the drift inequality or of the comparability (Phi/C between 0.78 and 1.65) on any tree. The smallest margins occur at stars (a single core with all other vertices leaves), where the drift ratio equals the nominal exponent to within 2e-9: the certificate is exact on stars, and the certified endpoints sit 1.9e-7 and 1.3e-7 outside the nominal ones because of the residual penalty. On larger ordinary trees the ratios close in on each other (0.565709 to 0.565740 at time 10), around the numerical value of mu.

## 5. Not verified, and remarks

* Not verified by anyone here: a second implementation of the 191 million local rows per side (the verifier was replayed, not rewritten); the exact Bernstein tail beyond running the claimant's script; the existence of the exponent, which the claim takes from the earlier proof.
* Child's findings (none affects the bracket): the proof file omits the step from the one-edge seed (conditioning on the first event fills it); two sentences about the tail are misleading as written but the implemented bounds are right.
* Comparison: this bracket lies strictly inside Neon claim 3 ([0.565680, 0.565768]) and Astra claim 21 (current [0.565680279, 0.565768284]) and is about half their width. My numerical value 0.56572366 lies 2.15e-5 above its lower end and 2.47e-5 below its upper end. The rationals of smallest denominator in it are still 99/175 = 0.5657143 and 142/251 = 0.5657371, so the trivial bound "if mu is rational its denominator is at least 175" does not improve (my numerical value excludes 99/175 by 9.4e-6, the bracket does not).
* The bracket certifies the fourth decimal (mu > 0.5657) and gives four digits; eight digits are far away.

## 6. gold2 (commit b435b4c, [0.565713486333, 0.565736986679]; section added 16:14 UTC)

VERDICT (gold2): HOLDS

Scope as for gold1 above (attacker: Claude Opus 5.5; the exact evaluator is a GPT child's, run by the lead on its own trees; a GPT-6.2 Sol Ultrafast child of this lead audited separately: report `round2/checks/neon2-claim6-gold2-and-astra25-audit-gpt-6.2-sol-ultrafast.md`).
* Hashes first: `quintet-gold2/quintet_lower.bin` 41565d49c10809eb..., `quintet_upper.bin` f1dcd7ad3e6df82a..., equal to the board (6,344,902,168 bytes each).
* Clean archive of b435b4c: `verify_quintet.cpp` is byte-identical to the b04cc50 source (compared with `cmp`), so the verifier built and checked for gold1 was reused. `cert_check.py --threads 2` (237 s): output **byte-identical to `certificate_gold2.json`**; lower certified 35083855339718511/62017003637500000 = 0.565713486333 (error per core 19217 ticks), upper 127993535256830043/226242120050000000 = 0.565736986679 (20062), both "ok".
* Planted errors (pinned verifier, lower side): one Fenchel-Young knot entry changed by 10^9 ticks gives error per core 18,650,603,143 and a certified endpoint 0.377755116; one 5-path entry gives 7,565,715,239 and 0.489467300 (authentic 0.565713486333). The child ran five more on both sides, all weakening their endpoints.
* Own trees (the same 1,930 snapshots, exact drift and comparability, 30 s per side): no violation; smallest margins 1.96e-7 (lower) and 1.27e-7 (upper), again at stars and two-core trees, where the drift ratio sits at the nominal exponent; at time 10 the ratios lie in 0.5657159 to 0.5657284 (lower) and 0.5657189 to 0.5657333 (upper), around the numerical value of mu.
* My numerical value 0.56572366 lies 1.02e-5 above the lower and 1.33e-5 below the upper end of gold2. The rational of smallest denominator in gold2 is still 99/175 = 0.5657143 (inside by 7.99e-7); 142/251 is excluded.

