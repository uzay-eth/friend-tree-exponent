# Independent discrete cross-check

This folder tests the packaged simulators against code written separately from the discrete growth rule. Claude Opus 5.5 wrote it without reusing any code from [`sim/`](../sim/), [`fixed-size/`](../fixed-size/), [`route-exact/`](../route-exact/) or the motif-bracket verifier. Every comparison below agrees within Monte Carlo error. Every negative control fails by tens to thousands of standard errors.

| Check | What is compared | Result |
|---|---|---|
| Exact enumeration | Exact $`\mathbb E[C_n]`$ from all reachable unlabelled trees vs the packaged exact table | identical rationals and tree counts for $`n=3,\ldots,21`$ |
| Direct chain at fixed size | $`10^8`$ direct trajectories vs exact $`\mathbb E[C_n]`$, $`\mathbb E[S_n]`$ and $`\mathrm{Var}(C_n)`$ | max $`\lvert z\rvert`$ 1.36 over $`n\le 21`$ |
| Fast CT generator vs direct chain | $`10^6`$ fresh fast-generator runs vs the direct chain stopped at an exact Yule size draw, $`t=2,4,\ldots,12`$ | max $`\lvert z\rvert`$ 0.94 for core and $`S`$ |
| Fixed-size generator vs direct chain | Packaged discrete generator, 9,600 fresh seeds, and the archived 2,400 trajectories vs the direct chain at sizes 100 to $`10^6`$ | max $`\lvert z\rvert`$ 0.72 (fresh) and 1.62 (archived) over nine sizes |
| Negative controls | Attach to $`V`$ itself; degree-biased neighbour; deterministic size $`3e^t`$ | rejected, $`\lvert z\rvert`$ from 39 to 7,445 |

## Programs

[`naive_chain.c`](naive_chain.c) stores every vertex, its parent and its children. Each step draws a uniform vertex $`V`$, then a uniform entry of the actual neighbour list of $`V`$. The new leaf goes to that neighbour. Each run uses its own SplitMix64 stream, seeded from the run index, so results do not depend on scheduling. It records $`C`$, the number of vertices of degree at least two, and

```math
S=\sum_{u\,:\,D_u\ge2}\frac{L_u}{D_u},
```

where $`L_u`$ counts the leaf neighbours of $`u`$ and $`D_u`$ is its degree. Its `yule` mode draws the size $`N=G_1+G_2+G_3`$ from three independent geometric variables on $`\{1,2,\ldots\}`$ with success probability $`e^{-t}`$. It then grows one discrete tree to size $`N`$. This is the law of the continuous-time tree at time $`t`$, because the population is a Yule process independent of the jump chain and the path on three vertices is deterministic. The `fixedN` mode replaces $`N`$ by $`\mathrm{round}(3e^t)`$ as a negative control. Model 1 attaches the leaf to $`V`$ itself, which gives the random recursive tree. Model 2 chooses the neighbour of $`V`$ with probability proportional to its degree.

[`exact_small.py`](exact_small.py) propagates the exact distribution of unlabelled trees from the single edge. It stores each tree once under a canonical centre-rooted encoding and uses rational arithmetic throughout. A leaf attached to $`w`$ has probability $`\frac1n\sum_{v\sim w}1/D_v`$. The script compares the number of reachable trees and $`\mathbb E[C_n]`$ with columns 2 and 3 of the [packaged exact table](../../verify/motif-bracket/exact_rationals_n21.txt).

[`compare_fixed.py`](compare_fixed.py), [`compare_ct.py`](compare_ct.py) and [`compare_sizes.py`](compare_sizes.py) pool the outputs and print $`z`$-scores. A $`z`$-score is the difference of means divided by its standard error. That error is the sample standard deviation over the square root of the run count, combined in quadrature for two independent samples. The variance check uses the normal-theory error and is approximate. Sizes and times observed along the same trajectories are correlated, so the rows are not independent tests. [`agg_ct.awk`](agg_ct.awk) reduces fast-generator output to per-checkpoint sums.

## Results

**Exact enumeration.** Every tree count and every rational $`\mathbb E[C_n]`$ for $`n=3,\ldots,21`$ equals the packaged table. At $`n=21`$ there are 2,144,505 unlabelled trees. The run took 987 s with a peak of 813 MB under Python 3.12. [Output](results/exact_small_n21.txt)

**Direct chain against exact values.** Four batches of 25,000,000 trajectories each, all grown to $`n=21`$, give the [pooled table](results/compare_fixed.txt). The largest $`\lvert z\rvert`$ over the means of $`C_n`$ and $`S_n`$ is 1.36. The sum of the 18 squared core $`z`$-scores is 4.26. Sample variances of $`C_n`$ also agree with the exact variances.

**Fast CT generator against the direct chain.** Four fresh batches of 250,000 runs of the packaged [`friend_ct`](../sim/friend_ct.cpp) ran to $`t=12`$ with fixed run counts. Each took 115 to 122 s and used about 4 MiB. Their mean population differs from $`3e^t`$ by at most 0.59 standard errors at $`t=1,2,4,8,12`$. A separate run of $`10^7`$ trajectories to $`t=1`$ gives mean population 8.15605 against $`3e=8.15485`$, $`z=1.02`$, and sample variance 14.0053 against $`3e(e-1)=14.0123`$ ([output](results/ct_t1_population.txt)). The [table](results/compare_ct.txt) compares their core and $`S`$ means with the direct chain at matched Yule draws.

| $`t`$ | fast $`\mathbb E[C]`$ | direct $`\mathbb E[C]`$ | $`z_C`$ | fast $`\mathbb E[S]`$ | direct $`\mathbb E[S]`$ | $`z_S`$ |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 3.86431 ± 0.00231 | 3.86235 ± 0.00073 | 0.81 | 2.26235 | 2.26133 | 0.94 |
| 4 | 12.25959 ± 0.00763 | 12.26382 ± 0.00540 | −0.45 | 7.02040 | 7.02251 | −0.43 |
| 6 | 38.35563 ± 0.02409 | 38.34484 ± 0.03814 | 0.24 | 21.79025 | 21.78494 | 0.21 |
| 8 | 119.25929 ± 0.07520 | 119.19020 ± 0.23777 | 0.28 | 67.56021 | 67.52564 | 0.25 |
| 10 | 370.05850 ± 0.23374 | 369.57440 ± 1.64592 | 0.29 | 209.43730 | 209.18179 | 0.27 |
| 12 | 1147.52115 ± 0.72542 | 1157.61650 ± 11.55850 | −0.87 | 649.26656 | 654.95748 | −0.87 |

The direct runs number 10,000,000, 2,000,000, 400,000, 100,000, 20,000 and 4,000 at these times. At $`t=12`$ the mean tree has about 488,000 vertices.

**Fixed-size generator against the direct chain.** Four direct batches of 3,000 trajectories each were grown to $`10^6`$ vertices and observed at the nine published sizes. Each took 274 to 286 s and peaked at about 77 MiB. The packaged [`fixed_size`](../fixed-size/fixed_size.cpp) generator was also rerun on 9,600 fresh seeds, 20261130000 to 20261139599, which do not overlap the archived seeds. Its means agree with the direct chain at every size, with $`\lvert z\rvert\le0.72`$ ([table](results/compare_sizes_fresh.txt)). The archived 2,400 trajectories lie about one standard error below the direct means at all nine sizes, with $`z`$ from $`-1.62`$ to $`-0.81`$ ([table](results/compare_sizes_archived.txt)). The same trajectories contribute to every size, so this is one fluctuation of about one standard error carried along the trajectories. It does not indicate a bias in the generator. A common relative offset leaves log–log slopes unchanged. A separate 1,000,000-trajectory direct run at sizes 100, 316 and 1,000 gives 9.5984 ± 0.0048, 18.5656 ± 0.0093 and 35.7914 ± 0.0181 ([output](results/sizes_small.txt)).

**Negative controls.** The same comparisons reject every wrong model. Against the exact values at $`n=21`$, attaching to $`V`$ gives $`z=7026`$ for the core mean and the degree-biased neighbour gives $`z=-1010`$ ([table](results/compare_fixed_negatives.txt)). Against the fast generator, the deterministic size $`3e^t`$ gives $`z=-56.6`$ at $`t=4`$ and $`-39.0`$ at $`t=6`$. The degree-biased chain at a Yule size gives $`z=899`$ at $`t=4`$ ([table](results/compare_ct_negatives.txt)). The degree-biased model agrees at $`n=4`$, where both neighbours of the centre of the path have degree one, and departs from $`n=5`$ on.

## Reproduce

Run from `simulations/cross-check/`. These results used GCC 13.3 and Python 3.12. Rebuilt binaries reproduced the saved `yule_t6.txt` and `neg_fixed_n21_model1.txt` byte for byte.

```bash
mkdir -p build
gcc -O2 -o build/naive_chain naive_chain.c -lm
g++ -O2 -o build/friend_ct ../sim/friend_ct.cpp
g++ -O2 -std=c++17 -o build/fixed_size ../fixed-size/fixed_size.cpp
E=../../verify/motif-bracket/exact_rationals_n21.txt

# exact enumeration
python3 exact_small.py 21 $E

# direct chain against exact values
for i in 1 2 3 4; do build/naive_chain fixed 2026112310$i 25000000 21 > build/fixed_n21_$i.txt; done
python3 compare_fixed.py $E build/fixed_n21_?.txt

# fast CT generator against the direct chain
for i in 1 2 3 4; do build/friend_ct 2026112360$i 250000 12 1 | awk -f agg_ct.awk | sort -k2n > build/ct_t12_$i.txt; done
runs=(0 0 10000000 0 2000000 0 400000 0 100000 0 20000 0 4000)
for t in 2 4 6 8 10 12; do build/naive_chain yule 202611235$(printf %02d $t) ${runs[$t]} $t > build/yule_t$t.txt; done
python3 compare_ct.py "build/ct_t12_*.txt" build/yule_t{2,4,6,8,10,12}.txt
build/friend_ct 777001 10000000 1 1 | awk '!/^#/ { R++; s += $3; ss += $3 * $3 }
    END { m = s / R; v = 3 * exp(1) * (exp(1) - 1); print R, m, 3 * exp(1), (m - 3 * exp(1)) / sqrt(v / R), "sampvar", ss / R - m * m, v }'

# fixed-size generator against the direct chain
for i in 1 2 3 4; do build/naive_chain sizes 2026112391$i 3000 100 316 1000 3162 10000 31623 100000 316228 1000000 > build/sizes_$i.txt; done
build/naive_chain sizes 20261123950 1000000 100 316 1000 > build/sizes_small.txt
for i in 0 1 2 3; do build/fixed_size $((20261130000 + 2400 * i)) 2400 > build/fresh_$i.tsv; done
(head -n 1 build/fresh_0.tsv; for i in 0 1 2 3; do tail -n +2 build/fresh_$i.tsv; done) | gzip -n > build/fixed_size_fresh_9600.tsv.gz
python3 compare_sizes.py build/fixed_size_fresh_9600.tsv.gz build/sizes_?.txt
python3 compare_sizes.py ../fixed-size/core_counts.tsv.gz build/sizes_?.txt

# negative controls
build/naive_chain fixed 20261123200 2000000 21 1 > build/neg_fixed_n21_model1.txt
build/naive_chain fixed 20261123300 2000000 21 2 > build/neg_fixed_n21_model2.txt
python3 compare_fixed.py $E build/neg_fixed_n21_model1.txt
python3 compare_fixed.py $E build/neg_fixed_n21_model2.txt
build/naive_chain fixedN 20261123704 2000000 4 > build/neg_fixedN_t4.txt
build/naive_chain fixedN 20261123706 400000 6 > build/neg_fixedN_t6.txt
build/naive_chain yule 20261123804 2000000 4 2 > build/neg_yule_t4_model2.txt
python3 compare_ct.py "build/ct_t12_*.txt" build/neg_fixedN_t4.txt build/neg_fixedN_t6.txt build/neg_yule_t4_model2.txt
```

Each saved file in [`results/`](results/) is the output of the matching command, with `results/` in place of `build/`.
