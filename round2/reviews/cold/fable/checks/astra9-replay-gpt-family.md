RESULT: REPRODUCED AND RECOMPUTED

Checked inside the GPT family. The full replay finished; the independent recomputation covers the full singleton base family and both baseline singleton-donor lifetime tables.

## 1. Second-party replay

From this directory, `mkdir -p copy && git -C ASTRA_REPO archive e08ac70 round2 | tar -x -C copy`. From `copy` I launched the following privately in the background, with stdin from `/dev/null`:

```sh
/usr/bin/time -p taskset -c 13 nice -n 5 python3 -B round2/occupation-tail/replay.py --snapshots ASTRA_REPO/round2/outputs > ../replay.log 2>&1
```

`taskset -c 13 nice -n 5 sha256sum ASTRA_REPO/round2/outputs/mc-{upper,lower}.bin` and `stat -c '%s %n' ASTRA_REPO/round2/outputs/mc-{upper,lower}.bin` gave:

```text
cae73eb2d4241ace1e566e257c052e4b908a31b9346b624cf63bb38cec4f8f91  ASTRA_REPO/round2/outputs/mc-upper.bin
5030698d2ef127d1f57684b8f570021fa726e19e91199ecee0319b3b3d994c04  ASTRA_REPO/round2/outputs/mc-lower.bin
611658376 ASTRA_REPO/round2/outputs/mc-upper.bin
611658376 ASTRA_REPO/round2/outputs/mc-lower.bin
```

Both hashes and sizes match `certificate.json`. The full replay finished in **1175.44 s wall time (19 min 35.44 s)**. `taskset -c 13 nice -n 5 python3 -B compare_replay.py` checked the log against the claim and exited successfully:

```text
caps (printed upper,lower): (155608, 115075): match
first three credits (printed lower then upper): (7402100, 11253928, 301604, 7709975, 8148713, 3226994): match
ceilings (printed lower,upper): (2654889165321, 2655887184018): match
additional credits (printed lower,upper): (15289028, 17458365): match
final fractions: (('lower', '565662790716854568223/1000008252270910000000', '0.5656581227'), ('upper', '1698046074520368376961/3001152517940340000000', '0.5657979941')): match
completion: real 1175.44
```

Final lines from `replay.log`:

```text
lower mask 63 at 1133/2000 rows 55259 510900 21172384 10109560 positive sources 55178 395272 14247516 3996010 newborn 15289028/1099511627776
lower fraction 565662790716854568223/1000008252270910000000 outward 0.5656581227
mask 63 at 1133/2000 rows 324 6480 77778 477848 positive sources 304 4633 21232 76492 newborn 1040281/1099511627776
physical states 1212 one-ring transitions 11814 tail gain identities 52 positive sources 1205 caught flat,newborn,shift
upper mask 63 at 1133/2000 rows 55259 510900 21172384 10109560 positive sources 55178 245509 14473328 3734840 newborn 17458365/1099511627776
upper fraction 1698046074520368376961/3001152517940340000000 outward 0.5657979941
mask 63 at 1133/2000 rows 324 6480 77778 477848 positive sources 304 1280 14792 68337 newborn 1217604/1099511627776
physical states 1212 one-ring transitions 11814 tail gain identities 52 positive sources 1205 caught flat,newborn,shift
Outward physical interval [0.5656581227, 0.5657979941].
real 1175.44
user 318.34
sys 11.00
```

## 2. Independent singleton residual scan

`singleton.cpp` reads the layout directly, without including claimant code. In ticks, I used the three exact rational groups `gain_rate * σΔf`, `(L/D)*σ(f_promoted-f+f_J+g(post,J))`, and `-σλf`. The gain rate is `L+r_p+k*r_c` below 72, and selects its maximizing endpoint from `L,L+(k+1)/2` above that. I evaluated the edge using physical mass `72/D` and evaluated the analytic gain successor at the 1399 seam. The analytic calculation in `analytic.py` also evaluates the polynomial in `proof.md` using exact Fraction tensor Bernstein conversion at each endpoint and on all 16 panels.

Reproducibility detail: `proof.md` names the two singleton centers but does not map the first two of the five center blocks to `r_p,r_c`; neither the format header nor `snapshot.hpp` names the blocks individually. I used the order `r_par,r_ch,x_par,x_ch,x_cc` listed in the separate prose `round2/checks/motif-base-hy4.md:60` and confirmed the resulting caps. Likewise, the README describes per-group outward rounding but doesn't define singleton group boundaries explicitly. The three groups above match the published finite ledger. Using the degree of each individual analytic endpoint polynomial (ignoring identically zero leading coefficients) matches the analytic ledger; artificially elevating every polynomial to bidegree `(3,4)` still gives a valid upper-witness bound, −2542, but does not reproduce −2532.

```sh
taskset -c 13 nice -n 5 g++ -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all singleton.cpp -o singleton
taskset -c 13 nice -n 5 python3 -B analytic.py
for side in upper lower; do taskset -c 13 nice -n 5 ./singleton ASTRA_REPO/round2/outputs/mc-$side.bin; done
```

Output (coefficient dump omitted):

```text
upper: 384 Bernstein coefficients; max=-1623084896333/640000000; directed bound=-79531159920317/31404816000; ceiling=-2532
lower: 640 Bernstein coefficients; max=951645367447/10000000; directed bound=951645367447/10000000; ceiling=95165
mode 1 mutation none header DF,d,K 1400,72,24
centered 2555 exact 22643913778562053/14680064000000 ceil 1543 at (D,k) 70,68 grouped 1544 at 70,68
finite-tail 976744 exact 8464455513907388557/54801628000000 ceil 154457 at (D,k) 1399,0 grouped 154457 at 1399,0
mode 2 mutation none header DF,d,K 1400,72,24
centered 2555 exact 10976223337392334867/117440512000000 ceil 93462 at (D,k) 70,0 grouped 93464 at 70,0
finite-tail 976744 exact 341881489505947459/3585842000000 ceil 95343 at (D,k) 1339,0 grouped 95343 at 1339,0
```

All six directed ledger numbers match. The smaller single-ceiling centered maxima are explained by rounding the three groups separately. A preliminary draft used unsafe rational cross-products in mutation maxima; the reported build uses Euclidean ratio comparison and UBSan, and the reruns below supersede that draft.

## 3. Independent occupation inequality (3)

`occupation.cpp` reconstructs the two *baseline singleton-donor* tests, including their size-two source tables because the size-one rows have zero direct source. For each state I evaluated each endpoint choice for each gain and minimized the resulting exact positive-denominator Bellman quotient, flooring it and capping at `2^40`. I separately substituted the result in (3), choosing the minimizing endpoint from the sign of the actual successor minus current value. The donor pool has minimum degree 9, all outside slots are bounded by one half, and promotion incorporates both the old successor and the pair with the *post*-promotion parent. All pooled-target or other omitted zero rows satisfy the inequality by nonnegativity. For the younger-target test the only pair enlargement from a singleton has the newborn as target, i.e., as the child. No claimant checker sources were used for this calculation.

```sh
taskset -c 13 nice -n 5 g++ -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all occupation.cpp -o occupation
for side in upper lower; do taskset -c 13 nice -n 5 ./occupation ASTRA_REPO/round2/outputs/mc-$side.bin; done
```

```text
mode 1 mutation none
older pair positive-capable rows 2327600 violations 0
older singleton rows 66023 violations 0 newborn 8148713
younger pair positive-capable rows 2065415 violations 0
younger singleton rows 66023 violations 0 newborn 3226994
mode 2 mutation none
older pair positive-capable rows 2327600 violations 0
older singleton rows 66023 violations 0 newborn 11253928
younger pair positive-capable rows 2065415 violations 0
younger singleton rows 66023 violations 0 newborn 301604
```

These are exact matches to all four baseline singleton-donor newborn values. My reported counts exclude the singleton pool and the identically zero pooled-target pair rows. The replay counts those as well: 1 additional singleton row, and 100/85 additional pair rows for older/younger target (99/84 valid exact donor states plus one pooled donor state). The quartet-donor test also has a singleton starting layer, but its newborn value depends on size-five quartet sources; that different table was not independently recomputed. The continuation test has a separate singleton layer as well; its newborn value depends on the larger arities and was not independently recomputed.

## 4. Planted errors

Commands run individually in loops for both `side=upper,lower`: `taskset -c 13 nice -n 5 ./singleton ASTRA_REPO/round2/outputs/mc-$side.bin "$m"` for `sign donor shift`, and the same with `./occupation` for `sign donor discount shift`. For occupation, row substitutions always check the unmutated source, rates, discount and index convention. The first upper occupation loop hit the shell timeout after all upper outputs; the lower commands were rerun separately. Relevant real stdout values:

| Residual scan mutation | Upper centered / finite-tail grouped maxima | Lower centered / finite-tail grouped maxima |
|---|---:|---:|
| Reversed sign | 3971 / 664692952930 | −10006 / 23610459225784 |
| Donor bound 1 | 1544 / 664693049629 | 93464 / 95343 |
| Shift `f` table one slot | 1217853960515 / 333228524375 | 8026779705796 / 1287009429712824 |

The lower tail's donor-bound mutation is inert for the maximized residual (signed gain jumps are nonpositive); the upper tail clearly rejects it. The occupation donor mutation changes both the source maximization and the exterior gain interval. In the lower case its overwide rate interval is conservative and leads to smaller newborn scores, so the row check still passes but matching the ledger rejects the result.

| Occupation mutation | Upper older / younger newborn | Upper violations, pair; singleton, older / younger | Lower older / younger newborn | Lower violations, pair; singleton, older / younger |
|---|---:|---|---:|---|
| Reversed source sign | 10957962 / 294559 | 109595 / 68616; 0 / 0 | 8364770 / 3304806 | 1994629 / 1806888; 0 / 0 |
| Donor bound 1 | 29559144 / 8558727 | 2191900 / 1944120; 0 / 0 | 9814819 / 69394 | 0 / 0; 0 / 0 |
| `113/200` in place of `1133/2000` | 8210972 / 3259437 | 2068354 / 1822729; 10390 / 59683 | 11335374 / 303615 | 1662922 / 1441957; 16839 / 8 |
| Pair reward indexed at pre-promotion `c` | 4369813 / 1362170 | 0 / 0; 2565 / 23 | 8062806 / 1975 | 0 / 0; 2164 / 0 |

For example, the upper wrong-discount older-receiving singleton at `(D,c)=(919,17)` had residual `-521784/3676000`; the upper pre-promotion-index older-receiving singleton at `(997,24)` had residual `-3934331202/3988000`. The lower wrong-index younger-receiving rows still hold with the resulting conservative proposal, but the newborn `1975` differs from `301604`.

## 5. Credit and boundary reread

- **Singleton:** the residual is booked once per observed vertex, and the continuation reserve uses `cap_f − R_max`. The baseline donors use `R_max − R_actual`, so the same event may be involved without paying its same summand twice. Older/younger tables choose distinct *directed* incidences even if both use the same pair of vertices. Verified the baseline source and lifetime row numerically as above.
- **Edge and higher motifs:** the donor losses subtract the actual input from the *same* maximum used in the unused-row credit; after both payments the row is at most `max(0,R_max)`. Edge FY bookkeeping requires both removal and addition and consistent physical-vs-continuous exterior domains; the prose expressly requires that. Triple donors can only use uncovered slots, so are disjoint from centered transfers to quartets. The quartet intrinsic reserve still maximizes its uncovered gains, while baseline quartet donors pay the gap to the actual gains. The proof's donor-leaf deletion argument identifies an incidence by the retained set and the deleted leaf, even when different deletions of one chronological set yield rewards. Verified these algebraic interfaces in the linked prose, not the actual edge/high-row source census, transfer-support gates or one-to-one identification in their implementation. No stationary LP contribution appears in the two claimed proofs; adding one would require its own disjointness analysis.
- **Founder:** any motif/transfer affected by the founder is within a fixed core radius, including the recipient of a lost transfer. All those row remainders must be bounded in *absolute* value (not treated as nonnegative credit). The interpolation increment bound and analytic singleton formula bound `D|Δf|`; all other affected jump contributions are bounded times polynomials in nearby core degrees. This plus `physical-test.md` gives negligible normalized expected error for either sign. The published text gives a bounding recipe, not an explicit index-by-index clipping list or constants for `E`/`E_1`; I did not check those against the implementation or reprove the mean asymptotic it cites. The lower comparison shift belongs only to value comparison: `sum(k−1)=−1` on the full core tree, hence a constant correction with vanishing normalized contribution, not a new drift credit.
- **Birth and omitted boundaries:** all promotions (including founder promotions) produce exactly one newly counted noninitial singleton. A nonsingleton containing the founder is not in the physical observation sum; every included new nonsingleton has a unique old set obtained by deleting the newborn. Pooled old coordinates never revert because `D,c` are nondecreasing. Dropping a pooled member's possible positive enlargement is safe; dropping its own old-score change is safe only because its pooled label stays unchanged at its own transition. An omitted zero-score row can only acquire nonnegative score. The strict ordering of old transitions and componentwise nesting are plausible from the stated cutoffs; I did not enumerate all chronological shapes and successor cases outside the baseline singleton-donor grids. The initial and terminal discounted endpoint treatments are conditional on the global moment growth estimate; they do not require a sign for the omitted founder terms.
- **Continuation boundary convention:** `continuation.md:27` says receiving motifs with an interpolated position are discarded, but is not fully explicit that this also applies to the singleton source. A centered/exact tuple is needed here: one cannot simply evaluate `cap_f − R_max` at a finite tail knot and linearly extend its value without justification. I checked the replay implementation's eligibility predicate *after* the independent calculations: `occupation-tail/solver.cpp:167` explicitly requires `inner(x[0])`, so both tail corners receive zero for that source. For triple donors its near-corner exact value at `1/d'` and far-corner lower bound at 0 imply an affine lower bound over the whole reciprocal interval, since the real raw gap is affine before flooring. I did not independently check every continuation source eligibility/floor, every projection or the row ledger/counts. The text alone does not print the exact discrete floor sequence or the complete tables needed to compare row by row.

Not independently redone: `y`, final endpoint arithmetic and the claimant's physical-tree oracle (lead already did these); the quartet-donor lifetime table, additional affine lifetime table, other base motif families, and full source eligibility census. No candidate mathematical discrepancy identified in the independent singleton work.
