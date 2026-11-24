RESULT: RECOMPUTED AND EQUAL

Checked inside the GPT family. All programs here use independently written calculations; no claimant executable or solver logic was used. Witnesses were opened read-only. `git -C ASTRA_REPO diff --name-only e08ac70 HEAD --` for the specified documents, manifest and format header produced no output.

## 1. Quartet-donor table

Built `quartet.cpp` using `TMPDIR="$PWD" g++ -std=c++20 -O3 -Wall -Wextra quartet.cpp -o quartet`, then ran `taskset -c 13 nice -n 5 timeout 2400 ./quartet ASTRA_REPO/round2/outputs/mc-upper.bin` and the identical command with `mc-lower.bin`.

```text
sign 1 mutation none
arity 5 rows 17034500 positive sources 14460608 violations 0
arity 4 rows 6200922 positive sources 0 violations 0
arity 3 rows 1832220 positive sources 0 violations 0
arity 2 rows 310254 positive sources 0 violations 0
arity 1 rows 5325 positive sources 0 violations 0
newborn 7709975 total violations 0
sign -1 mutation none
arity 5 rows 17034500 positive sources 14500159 violations 0
arity 4 rows 6200922 positive sources 0 violations 0
arity 3 rows 1832220 positive sources 0 violations 0
arity 2 rows 310254 positive sources 0 violations 0
arity 1 rows 5325 positive sources 0 violations 0
newborn 7402100 total violations 0
```

The claimed values are 7,709,975 upper and 7,402,100 lower. Both match exactly. Construction repeatedly tightens a linear upper bound on the current score using the violated minimizing policy. The integer check then evaluates the original residual and tests maximality of each uncapped score by testing one tick more.

## 2. Edge and higher caps

Built with `TMPDIR="$PWD" g++ -std=c++20 -O3 -Wall -Wextra edge.cpp -lgmpxx -lgmp -o edge` and the same command for `higher.cpp -o higher`. Ran each binary for both witnesses via `taskset -c 13 nice -n 5 timeout 2400 ./PROGRAM ASTRA_REPO/round2/outputs/mc-SIDE.bin`.

```text
sign 1 scanned 6483402 impossible zero 183322 exact maximum ceiling 520 1 (2,2) 2299 (68,23) grouped maximum 524 1 (2,2) 2299 (68,23)
sign -1 scanned 6483402 impossible zero 183322 exact maximum ceiling 10646 36 (9,2) 229 (21,21) grouped maximum 10651 36 (9,2) 229 (21,21)
sign 1 family V rows 28329840 peak ceiling 167
sign 1 family C rows 10398780 peak ceiling 125
sign 1 family P rows 9643480 peak ceiling 35
sign 1 family F rows 5467736 peak ceiling 31
sign 1 family Y rows 11332200 peak ceiling 64
sign 1 family T rows 1591659 peak ceiling 87
deepest 337 per child 290 higher 627
sign -1 family V rows 28329840 peak ceiling 5647
sign -1 family C rows 10398780 peak ceiling 904
sign -1 family P rows 9643480 peak ceiling 705
sign -1 family F rows 5467736 peak ceiling 985
sign -1 family Y rows 11332200 peak ceiling 1723
sign -1 family T rows 1591659 peak ceiling 50
deepest 7440 per child 1641 higher 9081
```

The published edge caps 524 and 10,651 match the independently computed upward-rounded maxima when the two parts of each incoming singleton transfer are rounded separately. Exact rational ceilings are four and five ticks smaller. The six higher maxima and the aggregate 627/9,081 match. Those sums use a directed integer interval to find the exact ceiling, with GMP fallback for an ambiguous interval.

`TMPDIR="$PWD" g++ -std=c++20 -O3 -Wall -Wextra envelopes.cpp -lgmpxx -lgmp -o envelopes`; `taskset -c 13 nice -n 5 timeout 2400 ./envelopes ASTRA_REPO/round2/outputs/mc-SIDE.bin` for both signs also reproduced the envelope census and passed the exact projection inequalities. Real output:

```text
sign 1 negative-radius plant 0
negative stored radii 0
one-dimensional 2249983 failures 0
two-dimensional 7418600 total failures 0
sign -1 negative-radius plant 0
negative stored radii 0
one-dimensional 2249983 failures 0
two-dimensional 7418600 total failures 0
```

## 3. Full extension, including all arities and sources

Recomputed the singleton unused row, the edge unused row and edge donor reserve over physical reciprocals, the tightened intrinsic triple and quartet rows, and the uncovered triple-donor reserve on four-sets. The known-state physical rows use the underlying smaller motif weights directly. Each enlarged-set reward is projected over all its affected coordinates and floored once. Arithmetic on the row is directed integer: each term, including the negative discount, is floored in 1/4096 tick units. A floating candidate is corrected against the integer residual until it is feasible and one tick higher would fail or exceed the cap. All retained rows were checked.

Built `TMPDIR="$PWD" g++ -std=c++20 -O3 -Wall -Wextra extension.cpp -lgmpxx -lgmp -o extension`. Ran `/usr/bin/time -f 'wall %e peak KiB %M' taskset -c 13 nice -n 5 timeout 2400 ./extension ASTRA_REPO/round2/outputs/mc-SIDE.bin 4` for both signs:

```text
sign 1 highest 4 mutation none
arity 4 rows 10109560 positive sources 3734840 violations 0
arity 3 rows 21172384 positive sources 14473328 violations 0
arity 2 rows 510900 positive sources 245509 violations 0
arity 1 rows 55259 positive sources 55178 violations 0
newborn 17458365 total violations 0
wall 61.51 peak KiB 750064
sign -1 highest 4 mutation none
arity 4 rows 10109560 positive sources 3996010 violations 0
arity 3 rows 21172384 positive sources 14247516 violations 0
arity 2 rows 510900 positive sources 395272 violations 0
arity 1 rows 55259 positive sources 55178 violations 0
newborn 15289028 total violations 0
wall 64.35 peak KiB 750064
```

Claimed newborns: 17,458,365 upper and 15,289,028 lower, both equal. All eight row and positive-source counts match `../astra9-gpt/replay.log` lines 30 and 34. Additional capped-arity controls using `taskset -c 13 nice -n 5 timeout 2400 ./extension ASTRA_REPO/round2/outputs/mc-SIDE.bin R` with `R=1,2,3`: the command printed the same per-layer row/source counts wherever present, zero violations in each layer, and the following newborns. These are cumulative for *different truncated problems*, not an additive decomposition of the full credit. No claimant totals for these truncated runs were available for comparison.

| Maximum arity | Upper newborn | Lower newborn |
|---:|---:|---:|
| 1 | 271397 | 43573 |
| 2 | 281857 | 48478 |
| 3 | 11689413 | 10267854 |
| 4 | 17458365 | 15289028 |

## 4. Planted errors

The following used `taskset -c 13 nice -n 5 timeout 2400 ./PROGRAM ASTRA_REPO/round2/outputs/mc-upper.bin ...` with the arguments in the left column. For rejected occupation plants, failure counts below are rows with negative original residual evaluated on the perturbed table. For endpoint term rounding in upper-cap scans, the independent exact ceiling detects local underestimates; a global maximum alone is not the test.

| Program and arguments | Plant | Observed output |
|---|---|---|
| `quartet shift` | Quartet payload index +1 | First failing arity 5 shape 0 index 314947; 7,016,728 violations; newborn 7,861,803 |
| `quartet up` | Round the Bellman score up instead of down | First failure at arity 5 shape 0 index 319157; 15,847,197 violations; newborn 7,710,011 |
| `quartet donor` | Use `[0,1]` for outside gain and source maximization | First failure at arity 5 shape 0 index 319157; 9,129,226 violations; newborn 11,082,265 |
| `quartet duplicate` | Pay each selected incidence twice | First failure at arity 5 shape 0 index 319157; 14,460,608 violations; newborn 15,419,968 |
| `quartet overlap` | Book donor 4 a second time on donor 5's overlapping support | `shared quartet (0,1,2,3), contact 0, donors 4 and 5; legitimate 4639414 + 6185885; donor 4 entered again on donor-5 support accepted 0` |
| `edge shift` | Read old edge at next child label | Grouped cap 131,137,557,328,691 vs 524 |
| `edge down` | Round each upper term downward | `wrong-direction underestimates 6482700` |
| `edge donor` | Add the extra positive-drift allowance from extending FY gain beyond 1/2 | Grouped cap 560,736,142,409 vs 524 |
| `edge duplicate` | Book incoming FY contributions twice | Grouped cap 291,745,089,359 vs 524 |
| `higher shift` | Read old weight at next terminal label | Higher charge 2,050,283,472,495 vs 627 |
| `higher down` | Round each upper term downward | `wrong-direction underestimates 66763695` |
| `higher donor` | Maximize uncovered donor input on `[0,1]` | Higher charge 78,327,461,850 vs 627 |
| `higher duplicate` | Book each incoming cut contribution twice | Higher charge 260,908,832,706 vs 627 |
| `higher 3` | Add 100,000 ticks to first fork residual | Higher charge 100,492 vs 627 |
| `envelopes plant` | Negative projected radius | `one-dimensional 2249983 failures 1`; `two-dimensional 7418600 total failures 1` |
| `extension 4 shift` | Lookup enlargement with pre-promotion core degree | First failure arity 3 shape 0 row 10408841, residual −122996541 sub-ticks; 246,762 violations; newborn 2,456,904 |
| `extension 4 up` | Round all Bellman terms upward | First failure arity 4 shape 0 row 1805809, residual −4; 951 violations; newborn **unchanged** at 17,458,365 |
| `extension 4 duplicate` | Pay each edge/triple donor incidence twice | First failure arity 4 shape 0 row 1833777, residual −1283354894; 10,659,444 violations; newborn 30,740,706 |
| `extension 4 donor-source` | Use `[0,1]` for source gap of uncovered triple donor | First failure arity 4 shape 0 row 1833777, residual −1283354894; 1,713,981 violations; newborn 21,077,319 |
| `extension 4 donor` | Use `[0,1]` only in the lifetime minimization | **Zero row violations**, newborn 14,014,755; control exited 1 because the expected failure did not occur. This merely takes a minimum over a larger interval, so the wrong model gives a weaker still-valid bound. The newborn comparison detects the change. |

## 5. Proof-text reproducibility gaps

- The proof cites the format header rather than spelling out offsets and index order. Even the header does not identify the five centering and three FY sections in order; the array declarations in `motif-cert/export.c` lines 39–44 supply this byte-layout detail. No exporter or checker calculations were used.
- `motif-cert/README.md` describes upward rounding of groups without a full partition of the edge summands into the groups actually rounded. Combining `(i-r) Δf` before rounding gave caps 523 and 10,650; splitting `i Δf` and `-r Δf` reproduced the published 524 and 10,651. Separating the old-edge and newborn parts of the promotion alone did not reproduce the published caps. The exact rational result is independent of these groupings.
- `inside.md` specifies feasibility and score cap, but not the exact integer selection rule needed to make the logged newborn uniquely determined. Taking the largest feasible integer matches both logs.
- `continuation.md` likewise specifies acceptance, not a unique choice of accepted score. Maximizing the accepted integer row reproduces both totals. The text says to drop sources requiring interpolated receivers; in particular the finite degree-knot label is treated as a tail, not as an exact source-eligible state even when the physical mass on that label is one. The same condition drops the singleton source there. The quoted positive-source counts agree with this interpretation. The far-corner donor floor can be read as the sign-dependent minimum on the stated interval (using either endpoint).
- The extension replay available here gives row and positive-source counts per arity, and only one newborn per sign. It does not give independent newborn credits by arity.
- No remaining discrepancy was found in the requested numerical ledger. I did not regenerate witness weights, retest the previously assigned singleton caps/donor credits or lifetime ceiling, or rederive the final fractions.

`sha256sum ASTRA_REPO/round2/outputs/mc-{upper,lower}.bin`:

```text
cae73eb2d4241ace1e566e257c052e4b908a31b9346b624cf63bb38cec4f8f91  ASTRA_REPO/round2/outputs/mc-upper.bin
5030698d2ef127d1f57684b8f570021fa726e19e91199ecee0319b3b3d994c04  ASTRA_REPO/round2/outputs/mc-lower.bin
```
