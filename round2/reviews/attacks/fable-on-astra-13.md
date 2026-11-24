VERDICT: HOLDS

# Fable on Astra claim 13 (positive renewal using current parent degrees)

**Scope of the verdict:** The Claude-family lead re-derived the bookkeeping of the construction (which physical change is paid by which emission line) and tested its first-moment identity by its own Monte Carlo of the particle system with the actual reciprocal 1/D, written from the claimant's statement alone: the identity holds at the 0.1% to 0.7% level for three scores, four horizons and two starting trees, and three planted errors in the construction are detected (the weakest of them at 4 to 9 standard errors). That is a test of the construction at short times. The weight inequality with the constant 0.53, the compactness in the initial mark, the identification of the critical root and the transfer from the permanent-mark renewal were audited and recomputed only by a GPT-6.2 Sol Ultrafast child, the same family as the claimant: that part is a check inside the GPT family. The claim gives no digits of mu, as it says itself.

Attacker: the Fable lead (Claude Fable 5.2; family different from the claimant's GPT). Program and log: Fable folder `round2/rival-attacks/code/astra13-current_mark_identity.py`, `astra13-run_all.sh`, `astra13-current_mark_identity.log`.

## 1. The construction in my words, and the bookkeeping

Every non-centre core vertex v is evaluated as f_a(S_v): S_v its fringe (with a silent slot for the parent), a = q(degree of the parent) read at the current time. The particle system: an active particle (a; S) lets S evolve by the physical dynamics inside S, dies at rate a, and, keeping itself, emits for every core child Y of its top (degree D, p = 1/D, d = q(D), d' = q(D+1), r = rate at which the top gains leaves from inside S) an active (d; Y) at rate d - p, an active (d; Y plus one leaf) at rate p, a dormant (d', d; Y) at rate r; for every grandchild Z below Y a dormant (q(D_Y + 1), q(D_Y); Z) at rate p; and a newborn active (d; J) at rate L/D. A dormant (l, u; S) dies at rate l and becomes the active (l; S) at rate u - l.

I checked the accounting term by term. In the physical drift of the sum of f_a(S_v): (i) the events inside S_v are the particle's own evolution; (ii) the one event S_v does not see, a leaf sent by the parent (rate p of the parent), changes S_v into S_v plus a leaf, which is the line "Y plus one leaf", and it raises the degree of v, hence the labels of the children of v, which is the grandchild line, at the sender's rate p and not at the rate 1/D_Y of v itself; (iii) a gain of v from inside S_v (rate r) changes the labels of its children from d to d', the root line; (iv) a promotion at v creates the evaluation f_d(J) and changes no degree. The artificial death at rate a of the child's particle is returned by the term d f_d(Y) hidden in (d - p) f_d(Y) + p f_d(Y+) = d f_d(Y) + p [f_d(Y+) - f_d(Y)]. A dormant stands for the difference f_l - f_u; this needs a score that does not depend on the mark, as the claimant says. With q(D) = 1/D the first line has rate zero. I found no term missing or counted twice.

## 2. Monte Carlo of the first-moment identity (own program, q(D) = 1/D)

Identity tested (the claimant's (1) in `current-mark/physical/README.md`, with the seed term): for a score f of the shape alone,

    E sum over non-centre cores v of f(S_v(t)) = sum over initial non-centre cores of E[score of the particle system from (q(D_parent); S_v(0))]
                                                  + int_0^t E[centre's emission lines at time u, applied to the expected score at time t] du.

Left side: plain simulation of the tree. Right side: my simulation of the particle system (u uniform, one emission line of the centre chosen in proportion to its rate, recursion over emitted particles). Scores: f1 = 1 (number of non-centre cores), f2 = leaves of the top, f3 = core children of the top.

| start | t | samples | f1: physical / particle system | difference in standard errors (f1, f2, f3) |
|---|---:|---:|---|---|
| three-vertex path | 1 | 400,000 | 1.1178 +- 0.0019 / 1.1170 +- 0.0026 | -0.25, -0.38, -0.25 |
| three-vertex path | 2 | 1,000,000 | 2.8642 +- 0.0023 / 2.8695 +- 0.0074 | +0.68, +0.52, -0.11 |
| three-vertex path | 3 | 300,000 | 5.8921 +- 0.0077 / 5.8914 +- 0.0412 | -0.02, -0.18, -0.22 |
| chain of four cores (degrees 5, 2, 3, 4) | 1.5 | 300,000 | 7.7155 +- 0.0051 / 7.7204 +- 0.0184 | +0.26, +0.55, +1.63 |
| the same chain, second run | 1.5 | 1,200,000 | 7.7203 +- 0.0026 / 7.7294 +- 0.0093 | +0.95, +1.79, -1.37 |

Planted errors (the same program with one line of the construction changed):

| planted error | start, t, samples | difference in standard errors (f1, f2, f3) |
|---|---|---|
| no dormant at a gain of the emitter's top (root line dropped) | path, 2, 400,000 | -37.7, -33.6, -31.9 |
| the same | chain, 1.5, 300,000 | -54.9, -50.7, -54.6 |
| no dormant for grandchildren (grandchild line dropped) | chain, 1.5, 300,000 | -3.4, -7.4, -4.4 |
| grandchild dormants at the child's own send rate 1/D_Y instead of the emitter's p | chain, 1.5, 300,000 | +2.0, +2.5, +3.3 |
| the same | chain, 1.5, 1,200,000 | +4.0, +4.0, +8.9 |

From the three-vertex path the two grandchild plants are not visible at these horizons (differences below 2 standard errors at t = 2, 400,000 samples): grandchildren appear late there, which is why the chain start was added. The root line and both grandchild plants are detected from the chain; the wrong-rate plant needed the larger sample.

The first version of my program disagreed by 16 standard errors; the cause was my own mistake (a particle that neither died nor had an event before the horizon was not counted), not the construction. After the correction the table above was produced without further changes.

## 3. Inside the GPT family: audit and recomputation by a Sol Ultrafast child

Report: Fable folder `round2/checks/astra13-audit-gpt-family.md` (first line "CLAIM 13: HOLDS"). The child

* audited `current-mark/physical/README.md`, `current-mark/physical/gap.md`, `current-mark/physical/truncation.md`, `current-mark/weight/README.md` and `UNCAPPED.md` line by line (table in the report) and supplied the missing arguments it needed; no blocking gap;
* rebuilt the weight certificate with its own program and exact rationals: d = 73/200, rho = 53/100, seam M = 400, 638,392 edge corners and 159,598 local end points, virgin normalisation 999609761/10^9 < 1, h_min = 86271979/10^9; also the smaller witness (rho = 0.553, M = 160) and the 21-label coarsening (rho = 0.56495);
* verified on all shapes of at most 9 vertices (485 shapes, 7,760 states), in exact arithmetic, the transfer identity between the physical generator and the five emission lines, the drift inequality, and the raw and scaled dormant rows; at sizes 6 to 8 also the resolvent identities (I + F = (I - K)^-1 and others) at discount 8; two planted errors (receiver's send rate in the grandchild line; wrong dormant diagonal) were rejected;
* reproduced the Taylor differences 0, 0, 0, 0, 1/96, -113/3456, 14227/124416 between the current-mark and the permanent-mark newborn responses: the two kernels are different functions of the discount with the same critical root, as claimed.

Remarks of the child, none of them affecting the verdict: the first replay command printed in `UNCAPPED.md` runs without flags and then replays the smaller witness (rho = 0.553), not the main one (the main one is the command with `--M 400 --d .365 --rho .53`); the sentence "decreases because dM >= 1" is true for integer degrees by a one-line computation that the text does not give.

## 4. Not verified

By me: the inequality R_c W <= 0.53 W (an exact certificate, recomputed only inside the GPT family), the compactness of the kernel in the initial mark, and the proof that the kernel has spectral radius 1 exactly at mu (it rests on claims 2 and 5 and on the renewal remainder proved for claim 6 in `fringe/almost-sure.md`). My Monte Carlo covers horizons up to 3 and cannot see the long-time behaviour. As the claim says, the kernel entries still sum over unbounded fringes: nothing here gives digits of mu.
