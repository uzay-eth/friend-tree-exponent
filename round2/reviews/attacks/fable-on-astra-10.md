VERDICT: HOLDS

# Fable on Astra claim 10 (spatial limit; one-ended local limit; core-degree and branch-size laws)

**Scope of the verdict:** The Claude-family lead recomputed finite marginal counts and crossing ratios from simulations. The right-hand side of the crossing calculation uses the rate from the proof. The analytic spatial, rank-one and one-ended limits were read only by a GPT-family child, who also made the finite exact rate/flux checks and used Claims 2–7 as inputs.

Attacker: the Fable lead (Claude Fable 5.2; family different from the claimant's GPT) for the numerical tests of
sections 1 to 3, with own simulator and own counting code (Fable folder `round2/attacks/astra10-11/`:
`marginals.cpp`, `eval_marginals.py`, output `eval_16_20.log`). Section 4 is a proof audit by a GPT-6.2 Sol
Ultrafast child (same family as the claimant; used on the owner's order) and is labelled as a check inside the GPT
family.

## 1. Number of core neighbours: geometric tail with ratio 1/(1+mu)

3000 runs of my continuous-time simulator to t = 20 (mean 108,388 core vertices), laws time-differenced between
t = 16 and t = 20, non-initial core vertices, c = number of core neighbours (parent included):

    k            1       2       3       4       5       6       8       10      12             13             14
    P(c=k+1)/P(c=k)
                 0.6689  0.3127  0.4625  0.5379  0.5773  0.6004  0.6238  0.6336  0.6390(23)     0.6333(29)     0.6429(36)
    P(c>=k)(1+mu)^(k-1)
                 1       0.8218  0.5077  0.4135  0.3712  0.3486  0.3276  0.3197  0.3168(3)      0.3160(4)      0.3146(5)

At the reference estimate 0.56572366, the predicted limiting ratio is 0.63868. The displayed ratios rise from k=2 to k=12 and then fluctuate near that value; they are not monotone across the table. The scaled tail is near 0.315 at the largest levels shown. P(c = 1) = 0.475100(20), P(c = 2) = 0.317795(28), and P(c = 3) = 0.099379(19).

## 2. Size of the descendant branch: P(N >= n) ~ A_N n^-mu with 1.1618 <= A_N <= 1.5999

N = number of vertices (core and leaf) of the branch below a non-initial core vertex, the vertex included:

    n               2       4       8       16      64      256     1024       2^14       2^18       2^20
    n^mu P(N >= n)  1.4801  1.3613  1.3257  1.3117  1.3027  1.3012  1.3000(4)  1.2992(8)  1.2990(19) 1.2997(28)

From n=64 to 2^20 the displayed estimates range from 1.2990 to 1.3027; the table gives pointwise standard errors where shown, not an error on a fitted plateau. The interval in the heading gives the rounded theoretical bounds evaluated at the reference estimate 0.56572366.

## 3. The crossing identity mu P(N >= n+1) = E[(n + 1/D_parent - 1/D_v) 1{N = n}]

    n            2          3          5          10         20         100        1000
    right/left   0.9999(1)  1.0000(2)  1.0000(3)  0.9991(5)  1.0000(8)  0.998(3)   1.008(17)

The displayed ratios are consistent with one at the reported precision (largest deviation about two standard errors). The simulation inserts the analytic rate N + 1/D_parent − 1/D_v on the right-hand side. The finite exact rate check is separate (section 4).

## 4. Inside the GPT family: proof audit by a Sol Ultrafast child

Report: Fable folder `round2/children/astra10-11/REPORT.md` and tracked Fable folder `round2/checks/astra10-11-audit-gpt-family.md`.
Line-by-line audit of Astra folder `round2/fringe/spatial-limit.md` and `round2/fringe/other-marginals.md` (`eb8e738`), taking the cited earlier results
(claims 2 to 7) as inputs: every step is marked "follows" or "follows with a supplied argument"; no gap found.
Points worth recording:

* The rank-one limit is proved for each fixed bounded continuous score and uniformly over scores of a fixed
  radius; it is NOT convergence in operator norm on all bounded continuous scores (counterexample: the indicator
  that the whole rooted core has at most M vertices). The claim as posted ("on bounded local scores in the sense of
  the proof file") is within that scope.
* The change of measure for the proposal clock is applied at bounded stopping times and extended by monotone
  convergence; the constant B is positive, and the child supplies a uniform lower bound that the proof does not need.
* The rate of N -> N + 1 is exactly N + p_p - p_v with the root's parent slot not counted in N (the formula used in section 3); the Gamma(3/2) and Gamma(5/2) sandwiches and the constants of (4) check out.
* Exact enumeration through 12 vertices checked finite physical flux identities; two planted errors were detected.

## 5. Not verified

The earlier theorems (Claims 2–7) were not rederived here. The trajectories at t=20 have ordinary size of order 10^9 in expectation; the raw table does not record ordinary size. Standard errors describe the finite samples, not distance to an infinite-size limit.

To evaluate the prepared raw table from Fable folder `round2/`, run `.venv/bin/python attacks/astra10-11/eval_marginals.py attacks/astra10-11/run_16_20.txt`. The raw input is a working output; the tracked readout is Fable folder `round2/rival-attacks/code/astra10-11-eval_16_20.log`. The simulator source takes `SEED RUNS OUT T1 [T2 ...]` with 3000 runs and checkpoints 16 and 20 here; the seed is not recorded in the prepared input or run log.
