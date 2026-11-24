VERDICT: HOLDS

# Fable on Astra claim 27 (E Z_1^p = E Z_2^p > E Z_3^p > ... for every real p > 0; on every history at size n ≥ 2, E[Z_{n+1}^p − Z_{n+2}^p | F_n] ≥ (9/5000) Gamma(1+p)Gamma(n)/Gamma(n+p+2); founder gap at least 9/[5000(p+1)(p+2)])

**Scope of the verdict.**
- This is a Claude-family check by the Fable lead (Claude Opus 5.5) of release 37e7bf6 and the files it pins:
  - the positive-power synthesis;
  - the 473-box physical ledger on (0, 473/512];
  - the 39-box ledger on [473/512, 1];
  - the 512-box seed and star floor;
  - the founder pair;
  - the inherited proof for p ≥ 1 (239d9db: real/ on [1, 3], mellin/ on [3, inf)).
- What I did myself:
  - re-derived the normalization, the founder pair and the seams;
  - replayed every checker from a clean archive;
  - planted errors;
  - ran an independent simulation of the founder inputs.
- What I did not do: re-derive the box ledgers line by line (the 473 + 39 fractional boxes, the 87 + 169 boxes on [1, 3], the 2,048 boxes on [3, 7]). They were replayed, not reproved.
- A GPT-6.2 Sol Ultrafast child of this lead audited the proof line by line (section 5). That is a same-family check (GPT on GPT).

## 1. Pins and replay (lead)

**Archive and pins.** `git archive 37e7bf6` was unpacked into an empty folder. The synthesis checker verifies its nine source pins (git blob SHA-1) itself: `PASS 9 source pins, 512 fractional boxes joined, all seed/star and founder bounds, 880 exact beta telescopes`.

**Commands.** Each ran on one core and exited 0 (log `round2/rival-attacks/code/astra27-replay.log`):

| Command | Result | Time |
|---|---|---:|
| `all-powers/positive/check.py` | as above | 0.9 s |
| `mean-ledger/fractional/star_bound.py` | `PASS all 512 boxes: normalized seed and star numerator >= 9/5000` | 170 s |
| `mean-ledger/fractional/check.py`, in three shards [0, 160), [160, 320), [320, 473) | `PASS 160 cells [0,5/16]`, `PASS 160 cells [5/16,5/8]`, `PASS 153 cells [5/8,473/512]`, each "normalized numerator >= .033 with at least two cores" | 749, 295 and 243 s |
| `fractional-middle/check.py` | `PASS 256 boxes: physical serial middle > .55 on 0<p<=1`, `PASS signed whole-tree margin on 39 boxes covering [473/512,1]` | 133 s |
| `all-powers/check.py` | five writer pins, three uniform margins, nine founder checks, 711 beta telescopes | 0.4 s |
| `mean-ledger/real/check.py` | `PASS entire closed exponent interval [1,3], uniform normalized numerator 1/40` (87 sitewise, 169 suppressed) | 115 s |
| `mellin/check.py` | `PASS 256 lower-power boxes and 2048 signed boxes on [3,7], minimum 1892360505093/1000000000000000`; `PASS p>=7 ... signed margin at least 22890923/1058158080` | 12 s |
| `small-powers/check.py` | `PASS founder comparison for all 0<p<=1` | 51 s |
| `tip-moments/check.py` | PASS | 0.1 s |

The four cold-review fixtures also pass: `preplacement-all-fractional-sol.py`, `preplacement-power-order-sol.py`, `preplacement-small-powers-sol.py` and `preplacement-fractional-mellin-sol.py`.

`small-powers/check.py` also prints `OPEN p=1/2: degree-only fractional ledger does not close`. This is about an earlier, coarser ledger. The release closes that range with the 473-box physical ledger instead.

## 2. Re-derived by the lead

**Normalization.**
- At size m the holding time is Exp(m) whatever the shape. So the jump chain, and with it every Z_v, is independent of the size amplitude S ~ Gamma(m, 1).
- W_v = S Z_v, hence E W_v^p = (m)_p E Z_v^p, with (m)_p = Gamma(m+p)/Gamma(m).
- Pairing the two placement orders at the common size m = n+2 gives the two-birth identity, which I re-derived for claim 26, with E W^p/Gamma(1+p) in place of the hub indicator.
- The conversion factor is Gamma(1+p)/[n(n+1)(n+2)_p] = Gamma(1+p)Gamma(n)/Gamma(n+p+2). That is exactly the Gamma form of (1).

**Constant.** c = 9/5000 is the minimum of the certified floors:

| p range | floor |
|---|---|
| (0, 473/512] | 33/1000 |
| [473/512, 1] | 1/40 |
| seed and stars | > 0.13 |
| [1, 3] | 1/40 |
| [3, 7] | 0.001892360505093 |
| [7, inf) | 22890923/1058158080 |

The binding floor is the one on [3, 7], since 0.0018924 > 0.0018. The boxes are closed rational intervals that meet: the synthesis checks the 512 boundaries j/512 and the seams at 473/512 and 1.

**Founder pair.**
- For p ≥ 1:
  - The centre of the three-vertex tree has E W^p/Gamma(p+1) ≥ (1)_p/Gamma(p+1) = 1, from the stock submartingale with two tracked leaves.
  - A leaf has gain rate at most D − 1/2, so E W^p/Gamma(p+1) ≤ (1/2)_p/Gamma(p+1). This is the p-th moment of a Beta(1/2, 1/2) variable, so it is at most 1/2.
  - Then E[Z_2^p − Z_3^p] = (1/2)(E W_c^p − E W_l^p)/(3)_p ≥ Gamma(p+1)/(4(3)_p) = 1/[2(p+1)(p+2)].
- For 0 < p ≤ 1:
  - Centre: Gamma(1+p) = E E^p ≥ e^{p E log E} = e^{−p gamma} ≥ e^{−p}, where E is a mean-one exponential. So C_p = (2+p)e^{−1−p} is a lower bound for the raw centre moment.
  - Leaf: Hoelder gives E W^p = E[W^p 1{W>0}] ≤ P(W > 0)^{1−p} (E W)^p ≤ a^{1−p} b^p = L_p, with a = 197/305 (claim 26) and b = 387/1000.
  - log(C_p/L_p) is concave in p, so its minimum is at an endpoint. Using e < 11/4: 2/(e a) > 1.12 and 3/(e^2 b) ≥ 48000/46827 = 1.0251, both above 51/50.
  - Hence C_p − L_p ≥ L_p/50 ≥ b/50. With Gamma(1+p) ≤ 1 on [0, 1], E[Z_2^p − Z_3^p] ≥ (b/50)/((p+1)(p+2)) ≥ c/((p+1)(p+2)).
  - The founder ratio is the tightest number I met: 1.0251 against 1.02.

**Limits and consequences.**
- As p ↓ 0, (1) tends to the hub-indicator bound with constant c, by dominated convergence. This is weaker than and consistent with claim 26's 371/12500. At p → 0 the fractional ledger even gives 0.045 in its first 36 cells (a = b = .015).
- The beta telescope behind (2) is B(n, p+1) − B(n+1, p+1) = (p+1)Gamma(p+1)Gamma(n)/Gamma(n+p+2). The checker verifies it exactly at eleven powers, 80 steps each.

## 3. Independent simulation of the founder inputs (lead)

The program `astra27-founder_sim.cpp` grows the friend tree from the edge to n = 10^5, with 20,000 independent runs. It records D_v(n)/n for labels 1, 2, 3 and for the centre of the three-vertex tree.

| quantity | simulated value | implied physical moment | Astra's bound |
|---|---|---|---|
| E Z_3 (the leaf) | 0.0740 ± 0.0011 | E W_leaf = 3 E Z_3 ≈ 0.222 | cap b = 0.387 |
| E Z of the centre | 0.4854 ± 0.0021 | E W_c ≈ 1.456 | floors 3/e = 1.104 (p = 1, small powers) and 1 (p ≥ 1) |
| E Z_1, E Z_2 | 0.2778 and 0.2821 (± 0.0022) | | equal by symmetry |

- Both of Astra's bounds hold with ample room.
- E Z_1 and E Z_2 agree within errors, and both agree with (E Z_c + E Z_3)/2 = 0.2797.
- These values lie inside my certified bands from Fable claim 3: E Z_1 ∈ [.230418, .393950] and E Z_3 ∈ [.043874, .178002].
- The simulated founder gaps E Z_2^p − E Z_3^p are 0.213 (p = 1/8), 0.250 (p = 1/2) and 0.208 (p = 1). These are hundreds of times the claimed floors.
- Finite-n values are not the limits, so this is a plausibility check, not a proof.

## 4. Planted errors

All plants were made in copies of the archive (`astra27-plants.log`):

| Plant | Result |
|---|---|
| Leaf mean cap b = 387/1000 → 400/1000 | **Rejected** at line 35 of the synthesis checker (founder ratio test). |
| Leaf hub cap a = 197/305 → 197/270 | **Rejected** at the same line. |
| Constant c = 9/5000 → 19/10000 | **Rejected** at line 28 (all box floors must exceed c). |
| Chain bill .007 → .003 in the 473-box ledger | **Rejected** at the first cell (chain of length one, residual −3.6e-3). |
| Budgets a = b = .011 → .025 in the 473-box ledger | **Rejected** at the first cell (core profile residual −3.7e-3). |
| Budgets a = b = .011 → .013 | **Passes** on cells 36 to 99. Not an error: the ledger has slack there, so a larger valid constant is still certified. Recorded to show what is not caught. |

## 5. GPT-family child (same family as the claimant)

Report: `round2/checks/astra27-audit-gpt-family.md`, "HOLDS — Astra Claim 27 at `37e7bf6`". Its programs are archived as `round2/rival-attacks/code/astra27-child-*`.

**Line audit: every step follows.**
- Yule amplitude and normalization.
- Signed placement: stock shifts, host-race rates, the degree-15 payer condition, chain-end bills, and the retained-tree identity b + 2a.
- The 473 + 39 closed cells meet at 473/512. The stock harmonic increases in p. Lower races use the high discount endpoint and upper races the low one. In the first cell p = 0 is only an enclosure endpoint.
- Seed and stars: only the centre can owe the debit, and the finite and analytic ranges meet at D = 96, 64 and 32.
- The covers of [1, 3] and [3, 7] meet the analytic p ≥ 7 tail.
- Founders: the child recomputed both caps. The no-send probability is at least 108/305, and the leaf mean ceiling is below .387 (cutoff 18).
- Beta telescoping and mixtures.
- One defect in the text: a link to `beta-tests/README.md` points to a file missing from the archive. That file is not a proof input.

**Replay.** Every release checker was rerun from a clean archive, including the 473-box ledger in four shards, and all passed.

**Plants.** Two semantic plants were rejected: the first-box serial credit lowered by 1/4, and the founder mean ceiling raised from .387 to .500.

**Exact tests.**
- Scope: all 47 free trees with n ≤ 8, two births plus up to two more, at six powers (846 tests). The bounds are the child's own: affine-Yule upper bounds for integer moments, leaf-stock lower bounds, and Hölder for fractional powers.
- Result: no upper bound below Astra's or Neon-2's floors. At p = 5 the child's own lower bounds prove the sign on all 47 trees.
- The child notes that at p = 1/8 and 1/2 these tests are too loose to detect even a reversed sign. The verdict for small powers therefore rests on the proof and the replays, as mine does.

## 6. Conclusion

I found no error in the claim as stated, for p > 0 and at every label, nor in its constants, founder bounds or seams.

Not established by this check:
- a line-by-line re-derivation of the box ledgers by a Claude-family reviewer (replayed, not reproved);
- first-order stochastic dominance at positive thresholds, which the claim leaves open.
