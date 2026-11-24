VERDICT: HOLDS

# Fable on Astra claim 2 (renewal characterisation r(K_mu) = 1; E[C_n] of exact order n^mu)

**Scope of the verdict:** The Claude-family lead rederived the mean identity, checked a finite rational weight grid, and simulated the kernel and three mutations at finite horizons. A DeepSeek child checked the identity exactly on small trees. The unbounded variation-of-constants, compactness, positivity, and exact-order lower-bound arguments were read, with the existing growth exponent and bracket taken as inputs. No part is credited only to GPT-family children.

Attacker: Claude Fable 5.2 (a different family from the GPT claimant). A Gemini child left partial work before a quota stop; the exact addendum is by DeepSeek. The claim in Astra folder `round2/fringe/renewal.md` (`bfcbbe5`) is that the compact positive kernel has a unique root of r(K_s)=1 for s>51/100, at μ, and that E[C_n] is of order n^μ.

## What was attacked and how

1. **The mean identity (sections 1–2 of Astra folder `round2/fringe/renewal.md`), re-derived by hand.** For f summed over the fringes
   S_v of all non-initial cores, the only input of the true dynamics that the intrinsic dynamics L0 of S_v omits
   is a leaf sent by the parent of v, rate p = 1/D_parent, effect S_v -> S_v^+. Every other event (including a
   leaf sent from v to its parent, and anything deeper in the tree) is internal to each fringe that contains
   it. Writing p[f(Y+) - f(Y)] = [(a_Y - p) f(Y) + p f(Y+)] - a_Y f(Y) and moving the last term to the row of the
   child gives exactly generator (1) plus the birth term; the terms a_Y f(Y) cancel pairwise for every
   non-initial child, and for children of the centre they are supplied by the immigration (2). The rates
   a_Y - p are non-negative because the promoter's degree never decreases. I found no error.
2. **Inequality (4) and the copy-weight identity, exact arithmetic, own code.** (1 - c/D)[h(c+1) - h(c) +
   1/(1 + dD)] <= q h(c) checked with Python fractions for all 0 <= c < 60, max(2,c) <= D < 400 (smallest slack
   2.3e-4 at (c,D) = (1,3)); outside that range the left side is at most 1/(1 + d c) or 1/61, below q h = 0.252.
   sum_Y a_Y W(Y) = (a_S + d) W(S) - h(c_root) also follows by factoring each first edge, as in the proof. So R W <= (51/100) W.
3. **Monte Carlo estimates of the kernel K_s and its spectral radius.** Simulator written from
   the statement alone (Fable folder `round2/src/residual_sim.cpp`; it shares only the random number
   and Yule-increment routines with the lead's friend-tree simulator): a particle is a planted marked fringe;
   it evolves by the ordinary dynamics with an inert parent slot, is killed at rate a_S, emits for each core
   child Y of its root a copy Y+ at rate 1/D_root and a copy Y at rate a_Y - 1/D_root; rewards = direct births,
   collected at each ring of a particle's root as L/D with mark 1/D. One family with an immortal original gives
   all rows a by the weights e^{-a t} (own rewards) and e^{-a tau} (rewards below a copy emitted at tau).
   Columns: marks 1/2..1/40 exact, five lumps beyond; rows also at the lump end points; eigenfunction
   interpolated linearly in the mark inside a lump. At the central discount, the estimated discounted reward
   falls by a factor of about 0.70 between the last two unit bins, corresponding to a finite-window reward
   growth slope of about 0.21 after undoing the discount. The analytic weighted bound is e^{0.51t}; the
   Monte Carlo evaluation instead uses a fitted geometric tail.

   | horizon T | batches × families | s = 0.5557 | s = 0.5657237 | s = 0.5757 |
   |---|---|---|---|---|
   | 14 | 16 × 5000 | 1.03766(200) | 1.00369(192) | 0.97208(185) |
   | 17 | 16 × 2500 | 1.03588(361) | 1.00214(350) | 0.97074(340) |

   Parentheses give delete-one-batch jackknife standard errors in the last digits, without tail or lumping
   bias. The T=17 row uses all 16 archived batches, as the command below specifies. At the central discount
   the truncated estimates are 0.99533(190) and 0.99933(350); the exact time-truncated kernel is a lower
   bound, but these Monte Carlo estimates are not certified bounds. Putting the whole lump mass on either
   endpoint moves the T=14 tail-corrected estimate between 1.0019 and 1.0052. The estimates are consistent
   with a crossing near the reference estimate 0.56572366 and with the certified interval [0.5650,0.5665].
   At s=5, the archived T=14 and T=17 evaluations give K(0,1/2) estimates of 0.07685 and 0.07684;
   the analytic first-sojourn term alone is 1/13 = 0.07692 (rounded).
4. **Planted errors.** At s=0.5657237, T=12, and 32,000 families per variant, the estimated radii with
   geometric tails are 0.94892(319) when copies never get the extra leaf, 1.98194(1484) when copies are
   not killed, and 0.91857(424) when newborn marks are 1/(D+1) instead of 1/D.

## What was not recomputed

The justification of the variation-of-constants formula (3) for the unbounded system, the compactness and
strong positivity of K_s, and the count of "good" vertices in section 5 (the lower bound in E[C_n] of exact
order n^mu) were read line by line and no error was found, but they were not re-derived independently.
The finite-time offset measurements cited in the attack on Claim 3 are consistent with a leading exponential, but do not establish a sharp asymptotic or a constant-order remainder.

## Commands (from Fable folder `round2/`; use its NumPy environment)

    mkdir -p bin outputs/residual
    g++ -O2 -std=c++17 -o bin/residual_sim src/residual_sim.cpp
    bin/residual_sim 41 16 5000 14 outputs/residual/T14.txt ; .venv/bin/python src/residual_eval.py outputs/residual/T14.txt
    bin/residual_sim 42 16 2500 17 outputs/residual/T17.txt ; .venv/bin/python src/residual_eval.py outputs/residual/T17.txt
    for v in 1 2 3; do
        RES_VARIANT=$`v bin/residual_sim 51 8 4000 12 outputs/residual/variant`$v.txt 0.5657237
        .venv/bin/python src/residual_eval.py outputs/residual/variant$v.txt
    done

Output lines (s = 0.5657237):

    horizon T=14, batches=16, families per batch=5000
    s=0.5657237: r(K_s) truncated at T = 0.99533 +- 0.00190;  with geometric tail (ratio per unit time 0.6992) = 1.00369 +- 0.00192
    horizon T=17, batches=16, families per batch=2500
    s=0.5657237: r(K_s) truncated at T = 0.99933 +- 0.00350;  with geometric tail (ratio per unit time 0.6975) = 1.00214 +- 0.00350
    variant 1: with geometric tail (ratio per unit time 0.6971) = 0.94892 +- 0.00319
    variant 2: with geometric tail (ratio per unit time 0.8706) = 1.98194 +- 0.01484
    variant 3: with geometric tail (ratio per unit time 0.6857) = 0.91857 +- 0.00424

## Addendum: exact check of the mean identity by a child of a third family

After the Gemini quota stop, the identity of sections 1-2 (physical generator of sum_v f(S_v) equals
sum_v (R + B) f(S_v) plus the centre terms) was checked in exact rational arithmetic by a child running
DeepSeek v4.1 (opencode; neither GPT nor Gemini nor Claude), with its own code written from a statement of the
construction: all 1039 marked trees reachable from the 3-path with at most 9 vertices (1, 2, 6, 18, 61, 208, 743
trees of sizes 3..9), five test functions including a hash-valued one: 5195 checks, 0 mismatches; and the
coefficient of every fringe value f(S) separately (the difference is linear in f): 0 of 1039 trees mismatch.
Planted errors restricted to R and B (centre terms unchanged): (a_Y - p) replaced by a_Y: 4010 of 5195 checks fail; f(Y+) replaced by f(Y): 2653 fail; newborn mark
1/(D+1): 2064 fail. Report and code: Fable folder `round2/checks/exact-identity/` (`REPORT.md`, `identity_check.py`). The command above regenerates the raw Monte Carlo inputs, which are working outputs rather than tracked files.
