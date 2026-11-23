# Independent friend-tree simulation check

**Verdict.** All means agree within Monte Carlo error at the tested sizes and times. The largest fixed-size discrepancy was 1.167 standard errors, and the core-mean discrepancies against the fast CT generator were 0.926 and 0.381 combined standard errors at times 6 and 8.

Claude Opus 5.5 sanity-checked the transition, Yule sampling, and error calculations and reproduced the fixed-size result.

## Method

GPT-6.2 Sol (Codex) wrote the direct Node.js simulator from the discrete growth rule. It starts with a single edge. Each step samples a uniform vertex and then a uniform entry from that vertex's actual neighbor list. It attaches a new leaf to the chosen neighbor. The core counts vertices of degree at least two and increases by one exactly when the target's degree changes from one to two. The simulator uses xoshiro128** with rejection sampling for uniform integer draws. The supplied original fast C++ continuous-time (CT) executable supplies checkpoint observations through its command-line output.

The fixed-size check uses 500,000 trajectories through size 21. Size two has core size zero by hand because both endpoints have degree one. For sizes 3 through 21, the comparison uses column 3 of the [exact rational table](../../verify/motif-bracket/exact_rationals_n21.txt). Each size uses the same set of trajectories, so the rows are correlated.

In CT, the rate-one clocks give a total birth rate of $`n`$ at size $`n`$. The next firing vertex is uniform and independent of the holding time conditional on $`n`$. Thus the Yule population can be sampled independently of the discrete jump chain. The size-three path reached from the edge is the CT initial tree. At time $`t`$, the population from three founders is a sum of three independent geometric variables on $`\{1,2,\ldots\}`$ with success probability $`e^{-t}`$. The direct method samples that law by inverse CDF using floating-point uniforms and logarithms, then grows to the sampled size and records the core. Tree growth and Yule size sampling use distinct seeded streams. The CT batches use fresh seeds of their own.

The standard error (SE) is the sample standard deviation divided by the square root of the run count. The fixed-size $`z`$ is the signed discrepancy from the exact mean divided by the direct SE. The CT $`z`$ uses direct minus fast, divided by the root sum of squares of the two independent SEs. In reading the tables, we used a deterministic mismatch or a stochastic $`|z| \ge 3.5`$ as the threshold for flagging a discrepancy.

## Fixed-size core means

| Size | Exact mean | Direct mean | Direct SE | $`z`$ |
|---:|---:|---:|---:|---:|
| 2 | 0.000000 | 0.000000 | 0.000000 | 0.000 |
| 3 | 1.000000 | 1.000000 | 0.000000 | 0.000 |
| 4 | 1.333333 | 1.333356 | 0.000667 | 0.034 |
| 5 | 1.583333 | 1.584042 | 0.000906 | 0.782 |
| 6 | 1.797222 | 1.798482 | 0.001088 | 1.158 |
| 7 | 1.989468 | 1.990314 | 0.001242 | 0.682 |
| 8 | 2.166659 | 2.167764 | 0.001379 | 0.801 |
| 9 | 2.332471 | 2.334228 | 0.001506 | 1.167 |
| 10 | 2.489218 | 2.491042 | 0.001624 | 1.123 |
| 11 | 2.638484 | 2.640046 | 0.001736 | 0.900 |
| 12 | 2.781414 | 2.782944 | 0.001842 | 0.830 |
| 13 | 2.918872 | 2.920710 | 0.001944 | 0.946 |
| 14 | 3.051531 | 3.053526 | 0.002042 | 0.977 |
| 15 | 3.179932 | 3.181206 | 0.002135 | 0.597 |
| 16 | 3.304515 | 3.306900 | 0.002227 | 1.071 |
| 17 | 3.425645 | 3.427378 | 0.002315 | 0.748 |
| 18 | 3.543632 | 3.545412 | 0.002402 | 0.741 |
| 19 | 3.658737 | 3.660420 | 0.002486 | 0.677 |
| 20 | 3.771189 | 3.773176 | 0.002568 | 0.774 |
| 21 | 3.881185 | 3.883164 | 0.002648 | 0.747 |

## Fixed-time population and core means

Each method uses 50,000 trajectories per time. The theoretical Yule mean population is 1210.286380 at time 6 and 8942.873961 at time 8. The theoretical SEs for a mean of 50,000 independent populations are 3.121071 and 23.086528, respectively. One combined core SE is about 0.4% of the fast CT mean at each time.

| Time | Quantity | Direct mean | Direct SE | Fast CT mean | Fast CT SE | Difference | Difference SE | $`z`$ |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 6 | Population | 1213.020020 | 3.122254 | 1208.028500 | 3.118934 | 4.991520 | 4.413187 | 1.131 |
| 6 | Core | 38.425200 | 0.107296 | 38.284920 | 0.106897 | 0.140280 | 0.151458 | 0.926 |
| 8 | Population | 8937.996680 | 22.948357 | 8913.190900 | 23.212267 | 24.805780 | 32.641024 | 0.760 |
| 8 | Core | 119.467320 | 0.336522 | 119.285860 | 0.337657 | 0.181460 | 0.476717 | 0.381 |

## Reproduction

Run in `simulations/independent-check` from the simulation package layout, with `friend_ct` available in the sibling `build` directory. If the supplied executable is unavailable, create `../build` and compile the packaged source with `g++ -O2 -o ../build/friend_ct ../sim/friend_ct.cpp` as in the [package replay instructions](../README.md#full-reproduction-commands).

The following commands used Node.js 24.22.3 and the supplied executable. They stream the CT observations directly into the comparison.

```bash
set -o pipefail
node --check direct.mjs
node direct.mjs exact ../../verify/motif-bracket/exact_rationals_n21.txt 500000 418781098
../build/friend_ct 375092141 50000 6 6 | node direct.mjs compare 6 50000 241893607 953617024
../build/friend_ct 684013792 50000 8 8 | node direct.mjs compare 8 50000 835276019 315760869
```
