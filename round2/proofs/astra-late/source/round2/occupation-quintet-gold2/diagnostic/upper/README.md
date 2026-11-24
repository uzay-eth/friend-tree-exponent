> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# Original GOLD2 upper: complementary mother-high H3

The new stationary floor is **`2.462450387956266e−7` normalized potential units per core**. Added to the frozen upper D24 uniform-prior ledger, it gives **`μ < 0.565736054406458`**, rounding the rational bound outward. This is an inward step of about `1.9905e−7` relative to the previous upper bound, or **2.33%** of its distance to `241/426`. It does not exclude `241/426`.

**Joint upper H3 update:** Adding the two corrected [upper-frontier floors](../../post-fy/upper-frontier/README.md) gives **`μ < 0.565736027670973`**. The two add-ons cover mother `D=19…79,c=9…12` or mother `D=80…120,c≤10` with middle `3≤D≤5`. They use middle `c≤4` and child `D≤7,c≤3`. The first avoids the frozen upper H3 arity cap of 8, and the second avoids the degree-two middle mask in this receipt. `compose.py` checks both successful floor logs and their exact sum of `13266167036667732636019911/364805903750000000000` ticks/core. Together the three fresh floors give `1356759587941864580840943426321/4417767772160000000000000` ticks/core. The remaining gap to `241/426` is `8.328140456272e−6`, requiring another `1.03025269735839e−5` normalized numerator units/core at the same denominator.

## Physical target and measured evidence

On the original GOLD2 upper witness the broad priority H3 negative charged pool was 3.992 / 3.861 ppm pooled / held out, with 3.590 / 3.384 in the mother-high orientation. H3 was larger than any other individual higher-row family. The compact D24-paid H3 already covers mother degrees 19–79, so the new cut uses a nonfounder mother `h` at **degree 80–480, arity ≤12**, her child `m` with **degree and arity 2**, and `m`'s child `y` with **degree ≤12, arity ≤3**. The single middle state has zero ordinary leaf stock at payment; degrees and arities are actual physical counts. The exact same mask measures **0.5065 / 0.4582 ppm** pooled / held out on the preserved trees. This measurement motivates the cut; it supplies no stationary input. The smaller 80–420, young degree ≤8 and arity ≤2 tested cut has floor `1.85468e−7` and is nested, not additive.

The reward is the positive part of the **negative original signed charged H3 row**, including the two edge transfers and the original covered-centering corrections. It can be added to the existing upper H3 price because the mother-degree supports are disjoint. It is also distinct from singleton and edge row payments; none of the three lower credits was recalculated or carried over. At each newly paid occurrence, `P+[-P]_+=max(P,0)` remains below the retained nonnegative H3 family cap, so D24 can keep its existing local cap and host savings, even at overlapping middle marks. Do not add a box gap on the same occurrence.

## Exact proof and composition

`price.py` SHA-pins the existing chronological H3 generator and runs it only on upper GOLD2, snapshot `neon2-gold2-GE7v6RUF`, source `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`, upper payload SHA-256 `f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`. The existing reader freshly verifies this digest. A nonfounder birth starts the one-vertex prefix. Each designated actual promotion continues the current prefix at its new arity and starts the longer one; paths leave the tracked set only at monotone degree/arity exits. The gain hazard and charged affine row share the same exterior rate. Integer Bellman residuals verify the two lower branches of the clipping function, on all rate boxes, with the original Gold2 upper bound as discount and lower bound for the stationary birth flux. No sampled rate or occupancy appears. The exact floor is

```
8803350896087646439582830411/35750368572477400958998937600000000  units/core
8803350896087646439582830411/32514770803097600000000              ticks/core
```

`compose.py` reads the frozen **upper uniform-prior** D24 trial and adds the exact new ticks to its exact uniform margin before dividing by the certified D24 upper mean potential. This reuses a margin valid on the original full prior; it does not extrapolate the narrower fixed-trial contradiction. The remaining distance to `241/426` is `8.35488e−6`; with that denominator the additional sufficient numerator is about `1.03356e−5/core`. For scale, the entire upper sampled higher-row negative pool was 11.87 / 12.05 ppm, inclusive of already-priced H3, so a single compact floor is not evidence of reaching the fraction. Other sampled higher-row families, notably S3 and the B5/S4 pair, are disjoint logical next candidates if separately certified against their original retained caps.

## Evidence

Run from Astra root (`R=round2/occupation-quintet-gold2/diagnostic/upper`):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 1 nice -n 19 timeout 135s \
  /usr/bin/time -f 'wall=%e user=%U rss=%M exit=%x' \
  round2/.work/venv/bin/python -Bu "$R/price.py" > "$R/mother80-480.log" 2>&1
python3 -IB "$R/compose.py"
bash "$R/sample.sh" > "$R/sample.tsv"
```

Exact replay exited 0 (`wall=56.40 user=46.49 rss=426948`); 42 literal labelled rows, 166,848 terminal states, 141,152 potentially paying rows. Composition exited 0 and printed `0.565736054406458`. The preserved-tree replay also exited 0 (`wall=2.42 user=1.15 rss=243236`), with no new simulation. Its TSV lists path, core count, nested-mask row count/reward sum, and main-mask row count/reward sum.
