VERDICT: HOLDS

# Fable on Astra claim 26 (birth order of the hubs for every label: p_1 = p_2 > p_3 > p_4 > ...; conditional two-birth gap at least 371/[12500 n(n+1)] on every old tree with at least two cores)

**Scope of the verdict.** This is a Claude-family check by the Fable lead, running as Claude Opus 5.5.
- **Covered:**
  - release b450a31: RELEASE.md, its 12 archive commands, the edge-dual ledger with its low-end closure, and the founder pair;
  - in the synthesis 212dc9f: the seed and star constants and the telescoped bound (2).
- **Not covered:**
  - the positive-threshold statements of the claim-26 board entry: the window (n+2)x ≤ 2^-131400, and the intervals [0.49, 1) and [0.48, 1);
  - the claim-24 additions ce6d4cc and fa3ac4d. A GPT-6.2 Sol Ultrafast child of this lead audited these additions and found no error (`round2/checks/astra26-audit-gpt-family.md`). Astra is also GPT-family, so that is a same-family check only.
- **Hub:** "hub" means a positive limiting degree share (threshold zero), as in the claim.

## 1. Pins, hashes, replay (lead)

- **Archive and hashes:** `git archive b450a31` was unpacked into an empty folder. All 18 git blob IDs in the RELEASE.md inventory match `git hash-object` of the unpacked files. None of these files changed between b450a31 and Astra's HEAD b2d83e7 (`git diff --stat` is empty).
- **Replay:** the 12 archive commands were run one by one on one core. All exit 0, taking 0.9 to 24.1 s each. The lines that carry the result:
  - low-end: `PASS endpoint degree 2: minimum residual 0.000005478682` and `degree 3: 0.000020000000`;
  - edge-dual: `PASS all degree-two chain endpoint types and simultaneous low endpoint budgets`;
  - second review: `PASS 3 mixed original-tree ledgers with direct and suppressed ports; floor 371/12500`.

## 2. Reading of the proof (lead)

I read the following:
- RELEASE.md;
- the founder pair (`early/README.md`);
- the signed formula (`global-transport/README.md`);
- the rescue and return race (`hostile/README.md`);
- the vertex credit and paired debit (`discharge/README.md`);
- the directed excess, census identity, chain suppression and low-end closure (`edge-dual/README.md`, `interfaces/low-end.md`);
- the classifier code (`hostile/check.py`: `return_lower`, `rescue`, `serial_lower`; `discharge/check.py`: the rounding of G, U, L).

**Two-birth identity, re-derived by hand.**
- The first birth lands at u with probability r_u/n.
- Given that, the rate of the second birth at v ≠ u changes only when v ~ u, and then by −1/[D_u(D_u+1)]. The newborn's own rate is 1/(D_u+1).
- How the terms combine:
  - Two leaves at the same target are exchangeable and contribute zero.
  - Ordered pairs of nonadjacent targets cancel.
  - On an edge, the two orders give net mass s_v − s_u, where s_u = r_u/[D_u(D_u+1)].
  - The serial completion (the second birth attaches to the first) has mass w_u = r_u/(D_u+1).
- Hence n(n+1)·gap = Σ_u w_u I(u) + Σ_edges (s_v − s_u) H(u,v), as stated.

**Census identity.**
- The old core tree has K−1 edges, and Σ_u (c_u − 2) = 2(K−1) − 2K = −2.
- The correction term in J(u→v) is valid because s_v ≥ (1/d)/[e(e+1)]: v receives at least 1/D_u from u.
- Then debit(u,v) ≤ max(J(u→v), J(v→u)) ≤ b + (J(u→v) − b)_+ + (J(v→u) − b)_+.
- Summing gives N ≥ Σ_u R(u) + b + 2a.

**Chain suppression.**
- On a bare degree-two chain, interior sites have s = 1/6. The sites next to an endpoint have s = (1/2 + x)/6; a one-site chain has (x + y)/6.
- Consecutive interior sites have zero paired debit, and each site adds credit 2sG(3) > 0. So lengths one to three suffice, as claimed.

**Classifiers: transition rates checked.**
- `rescue(d,l)`: tracked births at rate l move (d,l) to (d+1,l+1); tracked promotions at rate l/d move (d,l) to (d,l−1).
- Return race R_d(k): success at rate 1/k, against competing growth at most k − 1 + 1/d.
- Cutoff payoff 1 − d^-l: a linear birth-death stock with death rate at most 1/d per leaf dies out with probability at most d^-l.
- `serial_lower` = B(2,1) + the old-neighbour bonus − U(2,0); and L(d) = R_1(d)·B(2,1).
- U enters only as an upper bound and L only as a lower bound.

**Not re-derived line by line (replayed only):**
- the 6,105-row first-send supersolution table U (`preplacement-transfer`);
- the analytic high-degree tails at degree 113 and above, and the joint high exposure .012 (`mixed-stock`, `edge-dual` transfers).

## 3. Constants, recomputed in rational arithmetic (lead)

**b + 2a.** 17/1000 + 2·634/100000 = 2968/100000 = **371/12500**.

**Founder pair.**
- Leaf of the three-vertex tree:
  - Parameters: one permanent slot (B = 1), immigration 1 from the tag's own clock, cap m = 4.
  - The chain gives q_0 = 2/3, q_1 = 3/4, q_2 = 6/7, q_3 = 56/61, with 1 − f(4) = 3/5.
  - So Pr(the leaf is never promoted) ≥ (3/5)(3/4)(6/7)(56/61) = **108/305**.
- Centre: its two leaf families (birth rate 1, death rate at most 1/2 because D ≥ 2) both die out with probability at most 1/4. So H_c ≥ 3/4.
- Result: p_2 − p_3 = (H_c − H_l)/2 ≥ (3/4 − 197/305)/2 = **127/2440** (≈ 0.052049).
- Consistency with my own certified bands (Fable claim 3):
  - p_2 − p_3 ∈ [.246544, .329961], so the founder bound is weaker but consistent;
  - p_3 ≤ .489775 < 197/305;
  - H_c = 2p_1 − p_3 ≤ 1.0742, which puts no constraint on the bound 3/4.

**Seed and stars (212dc9f).** My recomputation of the basic local bound, using Astra's rounded tables, gives:
- seed: 0.189409485;
- star of centre degree 2: 0.201819087;
- star of centre degree 3: 0.225772497.

These match the stated 0.18940948, 0.20181908 and 0.22577249. For a star of degree d ≥ 4:
- the adverse cap is at most 1/d < L(2), so every old-edge charge is zero;
- the leaf serial terms sum to G(2)/2 > 371/12500;
- G(p) ≥ 0.16694 > 0 for every p, so dropping the centre's own serial term is valid. Astra's table decreases towards 0.166942351; its minimum for p < 400 is 0.166943476.

**Telescoped bound (2).** Σ_{n=i−1}^{j−2} c/[n(n+1)] = c(1/(i−1) − 1/(j−1)), and the conditioning on F_m needs m ≤ i − 1. Checked.

## 4. Exhaustive tests on small trees (lead's programs)

**(a) Astra's ledger on every small tree.**
- Program: the cold reviewer's function `check_tree` (`round2/checks/preplacement-all-sol.py`). It asserts that every retained residual is at least 0 and that the floor 371/12500 holds, after the simultaneous chain bills, direct excesses and high-to-low vouchers.
- What I ran: it on every free tree with at least two cores and 4 ≤ n ≤ 19, i.e. 522,940 trees (204,986 of them with n ≤ 18). The release ran it on 13 fixtures.
- Result: no assertion fails.
- Over these trees, the smallest basic local lower bound
  LB(T) = Σ_u w_u G(D_u+1) − Σ_edges |s_u − s_v| (U − L)_+
  is **0.274042754**, at the 4-path. That is 9.2 times the floor. It agrees with the discharge's minimum root charge 27404273/10^8 up to that check's outward rounding on a 10^-8 grid.
- Program and logs: `round2/rival-attacks/code/astra26-ledger_vs_exact.py`, and `astra26-ledger_vs_exact-ledger18.log` and `-ledger19.log` in the same folder. That folder also holds:
  - my enumeration `astra26-claim26_exact.py` and its logs;
  - the comparison logs (`-exact11x2b`, `-exact8x3`, `-stars13x2`);
  - `astra26-replay.log` (exit codes, times and all 18 blob IDs);
  - `astra26-plants.log`;
  - the child's three programs (`astra26-child-*.py`).

**(b) Astra's lower bound against my own upper bound.**
- For each free tree with at least two cores I computed, in exact rational arithmetic, the conditional law after the two births and EXTRA further births. The states were merged on the shape with both tags marked.
- I then applied my own classifier bounds from Fable claim 3 (birth-death escape, first-promotion supersolution, bare-core bound, edge theorem). This gives a certified interval [lo, hi] for N_T = n(n+1)·gap.
- A contradiction would be either LB(T) > hi (Astra's lower bound above a certified upper bound) or hi < 371/12500.

| Trees | EXTRA | Contradictions | lo > 0 on | lo ≥ 371/12500 on |
|---|---:|---:|---:|---:|
| 425 trees, n ≤ 11 | 2 | **none** | 73 | 40 |
| 40 trees, n ≤ 8 | 3 | **none** | 34 | 28 |

- The smallest margin hi − LB is 1.033, at the 4-path: LB 0.274, my interval [0.081, 1.307].
- The trees my bounds leave undecided are mostly spiders with legs of length two. Astra's early notes treated these separately as hard cases.

**(c) Seed and stars, n ≤ 13 (two extra births).**
- My lower bounds are positive on all of them and at least 0.1196 > 371/12500. Examples: seed [0.1196, 0.8209]; star of degree 2 [0.1280, 0.9891]; star of degree 12 [0.1955, 1.3693].
- Astra's star and seed bounds lie inside my intervals.

**Power of these tests.** They can find gross errors only.
- My upper bounds sit at least 1.0 above Astra's lower bound, and N_T is of order 1. So an error in the constant of the size of 371/12500 would not be visible this way.
- They see only trees of up to 19 vertices (11 for the exact comparison). The delicate parts of the ledger (degrees of 10 and above, long chains, unbounded degree) rest on Astra's profile checks with concavity and convexity reductions. I replayed those checks but did not prove them again.

## 5. Planted errors

All plants were made in copies of the archive; the originals were not touched.

| Plant | Result |
|---|---|
| Terminal allowance a = .00634 → .00640 in `edge-dual/check.py` | **Rejected.** edge-dual fails (assertion, line 109); low-end fails (degree-two endpoint residual −0.000054521318). |
| Sign of the old-edge correction in `upper_out`: s − 1/[d e(e+1)] → s + 1/[d e(e+1)] | **Rejected.** Both fail; low-end degree-3 residual is −0.015594997842. |
| Terminal chain bill .021 → .0105 in `low-end.py` | **Rejected.** The chain-length-1 low-low residual is 0.021194 > b. |
| Terminal chain bill raised to .0211, .025 or .03 | **Passes.** Not an error: the degree-two terminal minimum 5.478682e-6 is reached without a chain port, so a larger bill is still affordable. It is listed to show what is not caught. |
| My own enumeration: the two labels exchanged | **Detected.** 11 violations among the 47 trees with at least two cores and n ≤ 8 (one extra birth); no tree has a positive lower bound. |
| My own enumeration: wrong attachment rule r(w) + 1/2 | **Not detected.** |

## 6. Conclusion

I found no error in the proof of p_1 = p_2 > p_3 > p_4 > ... at hubhood, nor in the stated constants 371/12500 and 127/2440, nor in the seed and star floors of 212dc9f.

What this check does not establish:
- the positive-threshold extensions;
- the claim-24 additions (GPT-family audit only);
- a line-by-line re-derivation of the U table and of the analytic high-degree tails (these were replayed, not reproved).
