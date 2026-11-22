VERDICT: CORRECT

# Cold audit of `attempts/strogatz-astra/bounds/` (BOUNDS.md, hub.cpp, reproduce.sh)

Auditor: Claude Fable 5.2 (claude-code harness), child thread of the Astra lead, 2026-11-22.

**Claim audited.** For the random friend tree started from one edge and every n ≥ 3,

    (15493968692/17611447440) · Π_{m=3}^{n−1}(1 + 0.560/m)  ≤  E[C_n]  ≤  (15349885606/810439425) · Π_{m=3}^{n−1}(1 + 0.570/m),

so every limit point of log E[C_n]/log n lies in [0.560, 0.570].

**Result.** I found no error. The two tables are valid certificates for the stated inequalities.
The evidence is of three independent kinds, all run by me:

1. `reproduce.sh` reproduces `baseline.log` exactly.
2. My own exact-rational verifier re-checks the drift inequality at **every** finite state of both
   tables (3,550,712 and 5,960,583 states, 0 failures) from the dumped integer tables. It uses its own
   implementation of the event list and label transitions and no rounding.
3. My own tree-level test grows actual trees, re-derives every label from the written rules, checks the
   labels against the true tree at every step, checks that every true rate lies in its claimed interval,
   and checks the per-vertex drift inequality with the **true** rates: 3.58 million vertex checks,
   0 violations. Deliberately wrong rates are detected by both tools (10 mutants, below).

Files audited (sha256):

    6cd878002e04b693319dc4571a632ed0e4e80cf281dbd73c62a343ea0ec0ea37  hub.cpp
    20426f48cd2c8dd7b93b072ea4b9aab8b00829561bf84b9a8b29cd7d1a602a31  reproduce.sh
    8d32e28a88bf369cc4f460b0f3624d92174a2f199b223d023d760db3d504e447  BOUNDS.md

**What I did not read:** the board, any other check or opinion of this certificate, `transfer.cpp`,
`TRANSFER.md`, `SENSITIVITY.md`. I had read `attempts/x-strogatz/REPORT.md` (the earlier 0.543/0.586
work) for my other assignment; its event list is the same model, so I am not blind to the model, only
to this certificate.

**Limits of this audit** (read these before relying on it):

* My exact verifier was written after I had read `hub.cpp`. It is separate code in another language,
  keyed by label multisets instead of indices, but it is not a blind re-derivation.
* The tree-level test is a sample: about 18,000 distinct finite states per run out of millions, trees
  with n ≤ 600, core degree of a checked vertex ≤ 11. The exhaustive statement comes from item 2.
* The tables were dumped by a patched copy of `hub.cpp` (18 added lines after the `PASS` line). The
  patched build prints the same `PASS` lines as the original.
* Same machine and compiler as the author (g++ 13.3.0, x86-64).

---

## 1. Reproduction

    $ cp bounds/hub.cpp bounds/reproduce.sh <scratch>/ && cd <scratch> && time ./reproduce.sh
    ... PASS lower 56/100 product prefactor 15493968692/17611447440
    ... PASS upper 57/100 product prefactor 15349885606/810439425
    real 0m58.7s
    $ diff <output> bounds/baseline.log        ->  identical

The same binary at stronger exponents refuses, as it should:

    hub lower 113 200 128 12 64 18 8 6 5 4 4 3 3 3  -> newborn ratios [0.99091670819, 0.990916708206]  screen failed
    hub lower  57 100 ...                           -> newborn ratios [0.981228303599, ...]            screen failed
    hub upper 113 200 168 12 90 20 9 6 5 4 4 3 3 3  -> newborn ratios [1.00981695881, 1.00981695884]   screen failed
    hub upper  56 100 ...                           -> newborn ratios [1.0201452365, ...]              screen failed

From these ratios the scheme's own thresholds are about 0.5604 (lower) and 0.5698 (upper). The claimed
0.560 and 0.570 are close to what this state space can give; there is little room to spare.

## 2. Hand audit of the argument

I checked each item against the model (rate of a leaf gain at u is L+ρ, rate of a conversion is L/D,
in units of n × one-step probability).

**Event classification.** Partition the moves (V, W) by W. The state of a core vertex u records its own
(L, c), a label per core neighbour and, at c = 1, possibly a witness two steps away. It can change only
when W is: u itself (gain); a leaf of u (conversion); a core neighbour w (w gains); a leaf of w
(w converts); the witness a or a leaf of a. Nothing else touches a recorded quantity, because degree
and core degree of a vertex change only through its own gain and its own conversion. The list in
BOUNDS.md is complete.

**Rates.** All correct:
* neighbour (d, k), k exact: gain in [d−k+1/D, d−k+1/D+(k−1)/2], conversion (d−k)/d exactly;
* k capped: lower endpoint 1/D, upper endpoint the same as above (the maximum of d−j+(j−1)/2 over
  j ≥ k is at j = k);
* witnessed (2,2) neighbour: gain 1/D + 1/e, or 1/D + [0, 1/6] when the witness is pooled; witness
  gain in [e−j+1/2, e−j+1/2+(j−1)/2], witness conversion (e−j)/e;
* central credit: 1/d per exact neighbour, [0, 1/b] per mark b.

**History consistency of labels.** The labels are functions of the history, not of the tree alone
(a mark is a degree bound from the moment of pooling; only the greatest mark survives). Each rule only
discards information: resetting a mark to h_c, replacing a mark by the greatest threshold of the next
layer, capping k at K_{c+1} ≤ K_c, turning an exact degree ≥ h_{c+1} into a mark. Monotone cutoffs
(h_c nonincreasing, enforced; K = 20, 7, 5) make every projected label valid.

**Witness.** The two origins are exhaustive. A (2,2) vertex w had degree two for its whole core life, so
it was (2,1) just before its single conversion. A neighbour u of w with c_u = 1 is then either w's old
sole core neighbour (witness: the newborn (2,1)) or the newborn itself (witness: w's old sole
neighbour, read from w's own layer-1 label, exact below 6 because h_1 ≥ 6 and K_1 = 20). BOUNDS.md
asserts this without the one-line reason. In `hub.cpp` the index arithmetic matches: `label(2,2)` in
layer 1 is index 1 = witness (2,1), reached only from the conversion of a (2,1) sole neighbour;
`bridge(d,k)` = 1+(d−2)(d+1)/2+k−1 enumerates (e, j), 2 ≤ e ≤ 5, 1 ≤ j ≤ e, as indices 1..14 and the
pooled witness as 15 = `offset[3]−1`.

**Tail D ≥ H.** With q = L/D: Gq = (L+ρ)(1−q)/(D+1) − q/D ∈ [−q/D, 1−q] because 0 ≤ L+ρ ≤ D.
* Lower, f = Bq: drift ≥ q(s − B/D) ≥ δBq when B(δ+1/H) ≤ s.
* Upper, f = A(1+λq): drift ≤ Aλ(1−q) + qs ≤ λA(1+λq) when s ≤ Aλ(1+λ).
Both conditions hold for the actual pooled newborn weight s (exact check in §3): lower
17611447440·(0.56+1/128) = 9999999999.5 ≤ 10007920767; upper 11174432898·0.57·1.57 = 10000000000.4 ≥
9996896633. A vertex never returns below H, and every child of a tail vertex is the newborn with mark H.

**From drift to the product.** Φ_n = Σ_u f(x_u) satisfies E[Φ_{n+1} | history] ≥ (1+δ/n)Φ_n (≤ with λ).
Lower: f ≥ 0 and f ≤ sup f give E[C_n] ≥ E[Φ_n]/sup f. Upper: f ≥ inf f > 0 gives E[C_n] ≤ E[Φ_n]/inf f.
At n = 3 the centre (L, c) = (2, 0) is the only core vertex. The constants in the claim are exactly
root/sup and root/inf (§3).

**Integer rounding in `hub.cpp`.** `divide` rounds a signed quotient correctly in both directions
(truncation plus a correction by the sign of the remainder). Every term rate × increment is rounded up
for the upper table and down for the lower, with nonnegative integer rate numerators; the right-hand side is rounded the other
way; the boundary value at D = H is rounded in the conservative direction. So a nonnegative integer
residual implies the exact inequality. Products go through `__int128`, conversions back are range
checked, the build runs under UBSan. The size bounds quoted in BOUNDS.md are right
(total rate ≤ D+1+c(H+1)+H+1 ≤ 2.7·10^5, increments ≤ 3.2·10^11, product < 9.2·10^18).
The floating-point solver only proposes values; the final pass (`bound(..., actual)`) reads the completed
integer table and the actual newborn weights.

**Closure.** Every transition target is looked up with `index.at`, which throws on a missing state, and
the start state and all newborn states are enumerated, so every reachable state is in the table.

## 3. Independent exact re-verification of every state (`verify_full.py`)

Input: the integer tables dumped after `PASS`. For each row (c, multiset of labels) and each L with
2 ≤ D < H the script builds the events from the rules of BOUNDS.md, takes each uncertain rate at its
worst endpoint, and compares the drift with (exponent)·f in `fractions.Fraction`. At D = H−1 it uses the
exact tail function, not a rounded one.

    lower table H=128, exponent checked 14/25: 3550712 finite states re-verified exactly, failures 0, minimum slack 0.000 scale units (exact 0)
    tail coefficient 17611447440, pooled newborn weight s=10007920767, tail condition OK
    finite weights in [0, 17611447402], global inf 0, global sup 17611447440.0, root weight 15493968692
    prefactor: 15493968692/17611447440
    c<=2 subset minimum (finite part): 3448168319
    RESULT: PASS

    upper table H=168, exponent checked 57/100: 5960583 finite states re-verified exactly, failures 0, minimum slack 0.080 scale units (exact 2/25)
    tail coefficient 11174432898, pooled newborn weight s=9996896633, tail condition OK
    finite weights in [810439425, 17543680998], global inf 810439425, global sup 17543859649.9, root weight 15349885606
    prefactor: 15349885606/810439425
    c<=2 subset minimum (finite part): 4052646232
    RESULT: PASS

State counts, weight ranges, root weights, prefactors and the c ≤ 2 subset minimum 3448168319/10^10 all
agree with BOUNDS.md and `baseline.log`. Negative controls: the lower table checked at 0.561 fails at
3,549,839 states; the upper table checked at 0.569 fails at all 5,960,583.

## 4. Tree-level test with true rates (`tree_check.py`)

The script grows trees from the three-vertex path. Half of the moves follow the real dynamics; the rest
are forced conversions, feeding of low-degree vertices, chain extensions and feeding of the largest
vertex, to reach unusual shapes (any move sequence is a legitimate history). After **every** move it
checks all labels against the true tree: exact (d, k) equal the true degree and core degree, a capped k
is a lower bound, a mark is a legal threshold below the true degree, at most one large mark, a witness
is the true far endpoint, a (2,2) sole neighbour always has a witness, and the state exists in the
table. On checked trees it also verifies that every true rate lies in its claimed interval, and that
for every core vertex

    Σ_{true events} (true rate) · (change of f at u + [u is the creator] f(newborn))  ≥ δ f(u)   (≤ λ f(u) for upper).

The grouped-event bookkeeping was itself compared with a brute-force sum over all (V, W) moves on
8,174 small trees (exact equality).

| table | sampling | runs | trees checked | core-vertex checks | of which D ≥ H | violations |
|---|---|---|---|---|---|---|
| lower | mixed | 3 | 40,667 | 1,412,310 | 0 | 0 |
| lower | natural | 3 | 33,056 | 336,597 | 5,253 | 0 |
| upper | mixed | 3 | 42,620 | 1,476,203 | 0 | 0 |
| upper | natural | 3 | 33,950 | 351,090 | 2,811 | 0 |

The extreme value of (local drift)/f over all checked vertices is 0.560000 for the lower table and
0.570000 for the upper one: tight states occur in real trees. Hence the test is sensitive. The lower
table tested against 0.5601 gives 1,331 violations in 40 trajectories, the upper table against 0.5699
gives 105. Two deliberately wrong label rules in my own script (witness not updated when it gains;
neighbour's k not updated when it converts) are caught at n = 10 and n = 6.

## 5. Wrong-rate mutants of the verifier

For each mutant I changed one rate in a copy of `hub.cpp` (both in the solver and in the integer
checker), let the mutant certify its own table at the submission's parameters, and then tested that
table under the true rules with both of my tools.

| mutant | what is wrong | mutant's own verdict | exact verifier, true rules | tree test, true rates |
|---|---|---|---|---|
| M1 lower | neighbour conversion (d−k)/(d+1) | PASS | PASS (0 failures) | no violation |
| M1 upper | same | PASS | FAIL, 66,943 states | 8,139 violations |
| M1p lower | neighbour conversion (d−k)/(d−1) | PASS | FAIL, 69,188 | 14,130 |
| M2 upper | exact neighbour's credit 1/(d+1) | PASS | FAIL, 146,124 | 82,475 |
| M2p lower | exact neighbour's credit 1/(d−1) | PASS | FAIL, 128,131 | 111,448 |
| M3 lower | gain interval collapsed to upper endpoint | PASS, table byte-identical | PASS | no violation |
| M3 upper | same | PASS | FAIL, 118,674 | 56,113 |
| M5 upper | central conversion L/(D+1) | PASS | FAIL, 74,681 | 84,486 |
| M5p lower | central conversion L/(D−1) | PASS | FAIL, 22,575 | 114,310 |
| M6 lower | mark b gives credit exactly 1/b | PASS | FAIL, 140,090 | 6,791 |
| M6 upper | same | PASS, table byte-identical | PASS | no violation |
| M7 lower | witness's credit from the middle vertex 1 | PASS | FAIL, 5 | 2,731 |
| M7 upper | same | PASS | PASS (0 failures) | no violation |
| M8 lower | witnessed (2,2) gain without 1/e | PASS | FAIL, 52 | 19,832 |

(M1p, M2p, M5p, M8 upper and M2, M5 lower fail their own screen, so there is no table to test.)

Every mutant table that is invalid under the true rules is rejected by both tools, and the two tools
agree on all 14 tables. Four mutant tables are valid: two are byte-identical to the original and two
are wrong in the safe direction. They show which endpoints are active:
* lower table: the lower endpoint of an exact neighbour's gain interval is never used (a growing
  neighbour never raises the central weight), and slower neighbour conversions only hurt;
* upper table: a mark's credit is always taken at 1/b, and faster witness growth only hurts.

## 6. Plausibility of the final numbers

    exact E[C_n], 3 ≤ n ≤ 11 (enumeration, attempts/x-strogatz/exact_small.py): all inside the stated bounds
    n=10: 1.7762 ≤ 2.4892 ≤ 38.694
    simulated means from the earlier report (heuristic): n=1e8: 14956.1 ≤ 24678.3 ≤ 382783.1

## 7. Remarks for the author (none changes the verdict)

1. **Wording.** The last paragraph of BOUNDS.md writes E[Φ_{n+1} | T_n]. Φ depends on the labels, which
   depend on the history. The inequality is conditional on the history; the iteration is unaffected.
2. **Missing sentence.** Say why the two witness origins are exhaustive (§2 above).
3. **No negative control in `reproduce.sh`.** Adding the 113/200 runs (4 lines, about 1 minute) would
   show that the checker can fail.
4. **`minimum finite slack 0`** in the log is the directed-rounding residual. The exact minimum slack is 0
   for the lower table and 2/25 grid unit for the upper table.
5. **Margins.** The tail conditions hold with relative room 8·10⁻⁴ (lower) and 3·10⁻⁴ (upper), and the
   screen ratios are 1.00079 and 0.99969. Nothing is wrong, but any change of parameters needs a rerun.

## 8. How to rerun this audit

All audit code is in `attempts/strogatz-astra/branch-work/audit-scratch/` (scratch, not part of the
submission): `make_sources.py` (derives the dumping copy and the mutants from `bounds/hub.cpp`; every
textual replacement must match exactly once), `verify_full.py`, `tree_check.py`,
`tree_check_final.log` (the twelve runs of §4), and `run_audit.sh`, which chains the commands of
sections 1 and 3–5 (about 1 GB of temporary disk, removed at the end).

All numbers in this file come from running those same commands one by one in this session.
I started `run_audit.sh` once as a packaged rerun and stopped it during step 1: `reproduce.sh` pins the
first permitted CPU at nice 15, that CPU was busy with other projects' jobs, and the lead asked for no
long jobs. So the packaged script has **not** completed a full run. The scratch directory is 56 KB;
the table dumps (43 MB and 72 MB) and mutant binaries were deleted and are regenerated by the script.
