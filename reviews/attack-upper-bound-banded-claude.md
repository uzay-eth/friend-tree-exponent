CLAIM 4 VERDICT: HOLDS

Referee: Claude Opus 5.5 (claude-code harness). Claim attacked: Astra's banded upper certificate
`E[C_n] <= (1539738114735/3908512529) prod_{j=3}^{n-1} (1 + (5657/9983)/j)` (BOUND.md, check.cpp, recipe.txt.gz).
I found no false step and no unjustified step. I could not break it. Every kind of event that can change a core
vertex's recorded state has a term in the checker's row. Each term's rate interval contains the true rate in every
tree, and each destination is the state the labelling rules produce. Both signed-transfer ends and both arrow ends
agree, and the tail and the rounding are in the safe direction. A term-by-term program audit of 1,532,992 finite
rows and 324 tail rows, from 2,336 simulated histories plus about 50 hand-built corner histories, found 0
mismatches. A mutation test shows the audit does detect corrupted terms.

All my files are in `checks/attack-astra4-opus/`: `terms_table.cpp`, `labelled.py`, `match.py`,
`run.py`, `corners.py`, `mutation.py`, `battery.sh` and `logs/`. Nothing outside that folder and this report was written.

---------------------------------------------------------------------------------------------------------------

## 1. Logic audit (BOUND.md, TRANSFER.md, BOUNDS.md against check.cpp)

Notation: the focal core vertex is u, with degree D, core degree c and L = D - c leaves. A core neighbour is w with
(d, k). The witness is x with (e, j). "Rate" means n times the one-step probability. The true rates are:
* gain of a vertex y: `sum_{v~y} 1/deg v`;
* conversion of a core vertex y: `(#leaves of y)/deg y`.

| event changing u's record | checker term (check.cpp) | rate interval vs truth | destination | verdict |
|---|---|---|---|---|
| u gains (L+1) | `LEAF` {target, L, L, 1} plus one `CG` term per label group, rate `[cnt/d, cnt/d]` (exact) or `[0, cnt/b]` (pooled) | true rate is `L + sum 1/d_w`. Exact neighbours give exactly 1/d. A pooled neighbour has true degree >= its mark b, so it gives 1/d_w in [0, 1/b]. Marks are lower bounds because degrees never decrease, resets go to h_c <= the old mark, and projections take the greatest mark <= the old bound. | same id at D+1; at a band end, `coarser` = the projection to the next band's cutoffs (witness kept at c=1, since h_1 = 620 in every band); at D+1 = H, `A + up(B(L+1), H)` | OK |
| u converts a leaf | `CC` {promoted + reward(newborn), L, L, D} | exactly L/D | `promoted` = projection to layer c+1, new (2,1) inserted, witness dropped, at most one high pooled mark kept (`index()`). Kept core degrees are capped by K_{c+1} <= K_c, so a saturated label never becomes "exact". | OK |
| birth (newborn charged to its creator) | inside `CC`: trial reward of `newborn(D,l,s)` | the rate is exact, so using reward >= actual newborn weight is safe whatever the sign; `run()` checks `rewards[i] >= weight(i, D=2)` for every newborn label (minimum reserve 4,340,525) | sole neighbour (D, c+1): pooled if D >= 620 (greatest of 620/1024/2048/2400); witness iff D <= 14 and 2 <= c+1 <= 4, taken from `groups.front()` (minimum label = minimum degree, then capped core degree) | OK |
| birth with witness | same | the creator had 1 to 3 core neighbours and D <= 14, so it is in band 0 with cutoffs 620/240/45 >= 6. A witness of degree e < 6 is exact, including (5,5) with K_3 = 5, since k <= d. | e >= 6 → pooled witness | OK |
| exact neighbour w gains | `NG`: target W(dest) - g(D,c;d,k) | without a witness: `[d-k+1/D, d-k+1/D+(k-1)/2]`. The k-1 other core neighbours of w have degree >= 2. If saturated (k = K_c), the lower end is 1/D and the upper end is max over j >= k of `d-j+(j-1)/2` = the value at j = k. With a witness: `d-k+1/D+1/e+[0,(k-2)/2]`, with `1/e in [0,1/6]` for a pooled witness. | `label(d+1,k,e,j)`: becomes pooled with mark h_c at d+1 = h_c; the witness is kept while still needed and dropped at d = 14 → 15 | OK |
| exact unsaturated neighbour w converts | `NC` (only if k < d and k < K_c) | exactly (d-k)/d | `label(d,k+1,2,1)`: the new (2,1) is the witness when one is needed (k: 1→2, 2→3, 3→4); at 4→5 the witness is dropped | OK |
| saturated or pooled neighbour converts; pooled neighbour gains | no term | the recorded state does not change | (none) | OK |
| exact witness x gains / converts | `FG` `[e-j+1/d, e-j+1/d+(j-1)/2]`; `FC` exactly (e-j)/e | x has e-j leaves, w gives 1/d, and its j-1 other core neighbours give at most 1/2 each. x differs from u and w because the graph is a tree. | e+1 (5→6 gives the pooled witness) / j+1 | OK |
| pooled witness (e >= 6) events | no term | no change | (none) | OK |
| events at distance >= 2 that are not witness events | no term | u's record depends only on deg u, the (deg, core degree) of its core neighbours, the stored marks, and the witness's (deg, core degree) | (none) | OK |
| signed transfer g(D,c;d,k) | added in u's central-gain target (`gain_out`), subtracted in the receiver's `NG` (`gain_in`) | The constructor rejects any key unless both ends are exact and unsaturated in each other's view (`d<u.cap && k<u.kcap && D<w.cap && c<w.kcap`). So both rows see the same true key, and both multiply it by the gain rate of the same vertex. Witness variants share the key. No key has D >= H or d >= H. | (none) | cancels exactly |
| arrow a(D,c), threshold T = 8 (18 entries, D <= 10, c <= 2) | receiver: `CG` target minus a when the label's degree/mark >= T; sender: `fee = up(sum cnt*receiver_max, D)` when D >= T | Exact neighbours: the subtraction is exactly a/d_w. Pooled with mark b >= T: the bound is `max(0, U-a)/b`, which is >= `(U-a)/d_w`. Pooled with b < T: the subtraction is omitted, which is safe since a >= 0. `receiver_max` is a max over all (d,k) consistent with the label (pooled: d >= b, any k; saturated: k >= K). | (none) | sender adds >= what receivers subtract |
| root, c = 0 | empty multiset; conversion inserts (2,1) | (as above) | (as above) | OK |
| c >= 12 (h_c = 2) | every neighbour is pooled with mark 2, except one high mark | 1/d_w <= 1/2 | newborn (2,1) becomes mark 2 | OK |
| tail D >= H = 2400 | analytic `f = A + Bq` | `Gq <= (L+rho)(1-q)/(D+1) <= 1-q`, conversions only lower q, transfers vanish, arrow fee <= a_*(1-q). The newborn is (2,1) with pooled mark H, whose trial reward is s. The checks `B+a_* <= lambda A` and `s <= lambda(A+B)` are enforced in integers in the constructor. | D = 2399 gain → `A + ceil(B(L+1)/H)` (safe); birth at 2399 keeps the creator at 2399 | OK |

Remaining points I checked:
* **Rounding.** Every term is `ceil(max(r_-(v-z), r_+(v-z)))`. `up()` is a ceiling for negative numerators too, because C++ truncates toward zero. The right-hand side is `floor(p z / q)`.
* **Completeness of the enumeration.** `label()` throws on a witness-requiring type without a witness, and `index()`/`indices.at` throw on states that are not enumerated. A true state that is missing would therefore abort my lookups. None did.
* **Construction order.** `get()` refuses successors that have not been computed yet.
* **Overflow.** Overflow is trapped by UBSan with `-fno-sanitize-recover`, and my build also ran under UBSan.
* **Final constants.**
  * The certificate is `Phi_3 = f(root) = 1539738114735` and `m = 3908512529`. m is the minimum over all finite rows and A, and A + Bq >= A > m. So `E C_n <= E Phi_n / m`.
  * Tail slack in integer units: `lambda A - (B+a_*) = 0.62` and `lambda(A+B) - s = 0.48`, both >= 0. Computed with A = 1380739398522, B = 383976882193, a_* = 398437500000, s = 10^12.

Small textual remarks, none of which is a gap:
* BOUND.md and TRANSFER.md describe the witness rule in slightly different words ("minimum recorded label" vs "sorted by degree and then label"). These are the same rule, and the checker uses the first group label.
* The remark that "the drift is nonincreasing in z" only matters for the search for z. The proof uses only the checked inequality at the final z.

## 2. Program audit, term by term

**Tagged checker (`terms_table.cpp`).**
* It is a copy of `astra_table.cpp` (= check.cpp plus table storage) with a `build()` that rebuilds each row's terms with tags: kind, label, count, target, lo/hi/den, destination id, base weight, transfer, arrow subtraction, reward and newborn label.
* `solve()` still checks the original terms. For every row it also calls `build()` and requires equality term by term.
* It adds the queries `TERMS`, `ID`, `REWARD`, `NBLABEL` and `STATE`, served on a UNIX socket.

```
$ g++ -std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all terms_table.cpp -o terms_table
$ ATTACK4_SOCK=.../server.sock nice -n 10 ./terms_table geometry.txt recipe.txt 5657 9983     (server.log)
verified rows 226545617, signed transfers 63309, minimum integer drift slack 3
scale 1000000000000, minimum newborn reserve 4340525, weight range [3908512529, 1764716280715], centre 1539738114735
PASS upper 5657/9983
ATTACK4 tagged rows compared 226545617
```
This is identical to `verified.log`, so the query reports exactly the terms that were checked. (`recipe.txt` = `gzip -dc recipe.txt.gz`, sha256 of the .gz = 584bb93c…0f60a as stated.)

**Labelled histories (`labelled.py`).**
* Written from the three .md texts, independently of brute.py and check.cpp.
* It keeps the true tree, the stored pooled marks and the witness identities.
* It asserts on every descriptor: marks <= true degree; marks are valid marks of the layer; at most one mark above h_c; exact neighbours are below the cutoff; the witness is adjacent to the near vertex, distinct from u and core; no stale witness.

**Matcher (`match.py`).** At each checked time it enumerates every event of the tree: a gain at each core vertex, and a conversion of a leaf of each core vertex. Each event is applied to a copy-on-write copy. The matcher records every old core vertex whose **recorded** state changes. Recorded means core degrees capped exactly as the checker's labels are.
* *Full enumeration* recomputes all cores after every event, and every leaf is converted separately. It was used at every check with n <= 60, which is 109,135 checks.
* *Local mode* is used otherwise. It recomputes only the radius-2 ball of the event vertex, with one leaf per creator multiplied by L.
* Full enumeration would flag any change outside the predicted neighbourhood as "not local". None was flagged.

For every finite row the matcher checks:
* the label multiset equals the checker state;
* the destination (band, layer, id) and its weight equal the checker term's; at the tail boundary, `base >= A+B(L+1)/H` with a gap below 1;
* the sum of TRUE rates of the events mapped to each term lies in `[lo,hi]/den`; the central gain is split into L and per-group `sum 1/d`;
* every state-changing true event has a term, and every term has its events;
* the transfer used equals g at the TRUE endpoint types, on both ends;
* the arrow subtraction equals `a(D,c) 1[d>=T]`, every sender charge is >= a at the true receiver type, and `fee*D >= true outgoing`;
* the newborn label equals the checker's, and the true newborn weight is <= the reward;
* the amended row with true rates, true transfers, true arrows and the true newborn satisfies `A_true <= checker integer bound <= floor(lambda z)`.

Tail rows are checked as `A_true <= lambda(A+Bq)`. Globally, the sum of the amended rows equals the exact drift `n E[Phi'-Phi | history]`, which is <= `lambda Phi`.

Commands (`battery.sh`, server on 1 core, client on 1 core). Real summary lines:
```
python3 run.py natural 1500 60 11 1 60 170   -> DONE natural histories=1500 seconds=160 queries=320529 violations=0   (87000 checks, all full enumeration, 407785 rows)
python3 run.py policy 700 120 12 2 60 180    -> DONE policy histories=700 seconds=156 queries=393664 violations=0    (42000 checks, 888521 rows)
python3 run.py hubtail 6 0 13 1 0 200 0.02   -> DONE hubtail histories=6 seconds=4 queries=5125 violations=0          (9315 finite + 56 tail rows)
python3 run.py hubtail 30 0 21 1 0 150 0.05  -> DONE hubtail histories=30 seconds=62 queries=20437 violations=0        (107267 finite + 261 tail rows)
python3 run.py natural 40 700 14 50 0 90     -> DONE natural histories=40 seconds=4 queries=56920 violations=0        (12168 rows)
python3 run.py policy 60 400 15 20 0 90      -> DONE policy histories=60 seconds=17 queries=91420 violations=0        (97672 rows)
python3 corners.py 1                         -> DONE corners seconds=3 queries=16150 violations=0                     (10264 finite + 7 tail rows)
```
Totals over all logs: 2,336 random or policy histories plus about 50 constructed ones; 133,864 checked times; **1,532,992 finite rows + 324 tail rows; 0 violations.**

Coverage counters (summed):
* bands 0–6: 1455232 / 52463 / 23188 / 1080 / 531 / 255 / 243 rows;
* c = 0…13+ all present, including 911 rows with c >= 12;
* matched terms: CC 1157833, NG 2421212, NC 1612223, FG 379048, FC 226352;
* rows with an exact witness 379048, with a pooled witness 134187;
* witness 5→6 events 36763; near neighbour 14→15 drops the witness 4291; NC installs a witness 76701;
* births with an exact witness 1219634, with a pooled witness 514818;
* neighbour crosses its cutoff 7256; band crossings 8446; gains entering the tail 58; rows at D = 2399: 58;
* rows with a saturated core degree 54192; rows with a pooled mark above the cutoff 60893;
* rows receiving an arrow from an exact sender 285496, from a pooled sender 66843; rows sending arrows 134857;
* rows with an incoming signed transfer 793411; NG events carrying a transfer 1769751.

Tightest values seen:
* `(bound - A_true)/z = 0` at rows such as (2,2,{(2,1),(2,1)}). These rows have only exact rates, so the interval bound is tight, never exceeded.
* `(lambda z - A_true)/z = 4.89e-11` at (2,2,{(2,2),(2,2)}).
* Global `(lambda Phi - drift)/Phi >= 1.75e-6`.

**Sensitivity (`mutation.py`).** I corrupted the served terms on 15 small histories:
* NG upper rate minus one unit: **999 violations**;
* FG destination shifted: **223**;
* incoming transfer plus one: **239**.

So the matcher would have caught such defects.

## 3. Constructed corner histories (`corners.py`, every step checked)

* A c=1 vertex whose sole neighbour has 19/20/21/26 core neighbours (saturation at K_1 = 20), first with all other neighbours of degree 2 (upper rate end), then with all other neighbours grown to 26+ (lower rate end).
* Two adjacent pendant cores, then random gains and conversions of the witness, the near vertex and the centre. This covers witness 2…5→6, near vertex 14→15, and k: 1→2→…→5.
* Neighbour cutoff crossings at c = 2/3/4/5 (240/45/13/8), a second neighbour crossing (two pooled marks, one reset), then central conversions that project the marks.
* Band crossings at c = 2 with neighbours of degree 100–160 becoming pooled (cutoff 85 at D = 80, then 46 at D = 200).
* Conversion projection of a c=1 neighbour of degree 300 (exact under 620) to c = 2 (cutoff 240 → mark 256), then that neighbour growing past 620/1024 while pooled.
* A centre with c up to 17 (cutoff 2) and one large neighbour.
* A vertex at D = 2399 with 3 core neighbours: birth at 2399, gain into the tail, then tail conversions and gains.

Result: `DONE corners ... violations=0`. The tightest row was `(lambda z - A_true)/z = 1.05e-10`.

## 4. Verdict

**HOLDS.** The certificate as checked proves the stated bound, on the following basis:
* every step of BOUND.md is implemented by check.cpp as written;
* every true transition of every core vertex is covered by a term with a valid rate interval and the correct destination;
* both cancellations (signed transfers and arrows) are exact or one-sided in the safe direction;
* the analytic tail and the newborn rewards are verified in exact integers;
* 1.53 million rows from real labelled histories matched term by term with no exception.

Limits: the program audit only covers states reached by the histories I ran, which are mostly D < 1000; the states at D >= 1000 come from about 36 hub histories. The argument for all other states is the logic audit in section 1. I did not re-derive the recipe values. They only enter through the checked inequalities, which I re-ran.
