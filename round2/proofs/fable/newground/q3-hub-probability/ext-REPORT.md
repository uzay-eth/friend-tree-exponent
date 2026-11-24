# Exact unlabelled certificate for random friend tree hub ordering

**Author:** Fable's GPT-6.2 Sol Ultrafast orbit-certificate author. **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** Astra (GPT) audited exact receipt agreement with the separate Claude implementation and the single-engine grid step; no outside-family review covers the terminal grid row. Uses the original edge-theorem classifiers and their strengthened variants; see [F3](../../README.md#f3).

`dp.cpp` computes with GMP integer limbs and rationals. `full.log`, `narrow.log`, `late.log`, `lastpair.log`, and `merged23.log` contain exact rational endpoints; `SUMMARY.md` gives the best bounds found and the marginal enclosures. Human-readable decimals are not used in any decision. Disregard preliminary `extension.log` and `local.log`, which predate the robust core-lower recurrence. The baseline is the lead's default classifier **including the edge theorem**, with `K=40` and the identical `B` series tail. `validate.py` instruments only the display of the lead's code in memory and compares both endpoints of each reported interval as Python `Fraction`s: all 30 endpoints agree at each of sizes 9 and 10.

## State space and arithmetic

Represent each unrooted tree by AHU parenthesis encoding at its centre, with a virtual root above the two halves in the bicentred case. Integer ordering of child encodings suffices as a canonical total order. The rooted shapes above a fixed unrooted tree correspond exactly to its vertex orbits under automorphisms. Orbits are found by folding equal encoded sibling branches onto the first equal sibling, position by position. This is exact by induction on subtree size; automorphisms can permute equal sibling branches and act independently inside them. Automorphisms preserve the centre or swap the two centres, which the virtual root allows.

At size `n`, set `J_n=lcm(1,...,n-1)` and `D_n=prod_{j=2}^{n-1} j J_j`. The attachment integer weight is `sum_{x~w} J_n/deg(x)`. For each old orbit, propagate its marked mass with this weight at each attachment vertex, then map the representative of the old mark to the new orbit. No orbit multiplicity for the mark is needed: taking the representative instead of the actual mark is an automorphism of the old tree and permutes attachment sites with unchanged rates. The unrooted mass is the sum over mark orbits of any one label. Use this mass to initialize the newborn mark, assigning it to its orbit in the new tree. Symmetry identifies the first two label laws. The code checks rooted and unrooted census counts, total integer rate and per-label unit mass. The difference classifier applies sign before it multiplies bounds, and collects integer coefficient buckets indexed by exact terminal rationals.

## Proofs of strengthened classifiers

Use the equivalent rate-one per-vertex clock embedding. At a core vertex with `k` leaves and `c` core neighbours, births of leaves occur at rate `β≥k`, and promotion of each of its leaves occurs at rate `1/(k+c)`, increasing `c` by one and leaving total degree invariant. Core degree is otherwise constant; an attached leaf never changes core degree except at promotion. Stop for the comparisons below on the specified boundaries. The population is Yule so no explosion occurs. Optional stopping for the nonnegative promotion supersolutions bounds the probability of success by any finite time, and hence ever. For lower-bound recursions below, the stopping time is a.s. finite: in a positive finite `(k,c)` region the changing events have rate at least one, `β≤k+c` in that region, and a uniformly positive chance of a sufficiently long consecutive string of births reaches the upper boundary.

1. **Leaf lower bound.** Suppress the parent clock until its first exponential(1) ring, independent of the suppressed evolution. The marked leaf remains a leaf; the parent's degree started at `d` and has expectation at time `t` at most `d e^t`, since its increase rate is at most its degree. Jensen and independence give probability at least `1/(2d)` that this first ring selects the marked leaf. On promotion the marked vertex has one core neighbour and one leaf, so `B(1,1)/(2d)` is a lower bound. This is the first-ring argument of the lead's d1 child, used here pointwise for terminal leaves.
2. **Core lower bound.** Set `A(0,c)=0` and at total degree 24 use the baseline rational `B(k,c)`, except at `(24,0)` where the boundary is conservatively zero. In order of decreasing total degree and then increasing `k`, let `a=A(k+1,c)` be the birth endpoint and `b=A(k-1,c+1)` the loss endpoint; set `A(k,c)=max(B(k,c), min(a, [(k+c)a+b]/(k+c+1)))`. Validity follows by induction: the conditional probability of birth at the next changing event lies in `[(k+c)/(k+c+1),1]` by the pointwise hazard ratio, so the minimum payoff is at one of those endpoints, without assuming monotonicity. Keeping `c+1` after a loss is the improvement over freezing `c` permanently.
3. **Leaf upper bound.** At degree 24 use baseline `F_c(k)`. In decreasing `c` and then decreasing `k` set `P(k,c) = min(F_c(k), [k(k+c)P(k+1,c) + (k-1)P(k-1,c+1)+1]/[k(k+c)+k])`; the loss term is zero when `k=1`. The stop is promotion or reaching degree 24. The code verifies the supersolution inequality at rate `β=k` directly after taking the minimum, and checks `P(k,c)≥P(k+1,c)` including the boundary. Thus increasing the birth rate to its true value only decreases the generator. On promotion the payoff is 1.
4. **Upper bound for a core with no leaves.** For any one currently adjacent core neighbour, let `k` count its unmarked leaves and let `c≥1` be its starting core degree. While it has not selected the marked core, its leaf population has births at rate at least `k`, losses at rate `k/D`, and success at rate `1/D`, with `D≥k+c`. Stop at `k=0` with bound 1, `k=41` with bound `1/40`, or on success with payoff 1. Solve `(k(k+c)+k+1)H_c(k)=k(k+c)H_c(k+1)+kH_c(k-1)+1`. The checked monotonicity makes extra births decrease the generator; the combined success/loss increment is nonnegative, so larger `D` decreases it as well. For the upper boundary, `f(k)=1/(k-1)` for `k≥2` and `f(0)=f(1)=1` is a supersolution: for `k≥3` birth drift is at most `-1/(k-1)`, and the combined success/loss term is `(k²-3k+4)/[(k-1)(k-2)D]≤1/(k-1)` since `D≥k+1`; `k≤2` is trivial or has nonpositive drift. Union bound over neighbours, taking the minimum of `H` and the old bound for each neighbour, then clipping at 1. Finally use the edge theorem with these strengthened uppers to raise the lower bound for the opposite endpoints.

## Scaling and the asymptotic gap

The number of rooted trees grows asymptotically like a constant times `2.9558^m m^(-3/2)`; the cost includes polynomial factors for the labels, canonicalizations, and exact arithmetic. The observed first-passage `m-u` is in `SUMMARY.md`: the improved bound already succeeds at the earliest possible size `m=u+1` throughout the certified range, whereas the baseline first succeeds at `m=u+2` for `3≤u≤18` and fails for `u=19` at `m=21`. These data do not show growing look-ahead for the improved classifier, but neither prove it stays constant nor address the exponential cost even at constant look-ahead. A positive asymptotic without a computable onset proves eventual positivity, but cannot be joined to a finite prefix unless the first asymptotically covered comparison is no later than the first uncertified comparison. Either a quantitative remainder with an overlap, or a different structural proof for the gap, is needed.

## Commands and outputs

Run in this directory. Never compile with `NDEBUG` (checks guard exact arithmetic and several of them also execute the arithmetic).

```sh
taskset -c 12 nice -n 5 g++ -std=c++20 -O3 -Wall -Wextra dp.cpp -lgmpxx -lgmp -o dp
taskset -c 12 nice -n 5 ./dp 10 9 1 > quick.log
taskset -c 12 nice -n 5 bash -c 'python3 validate.py 9 quick.log && python3 validate.py 10 quick.log' > validation.log 2>&1
taskset -c 12 nice -n 5 bash -c 'ulimit -c 0; ulimit -v 11718750 && /usr/bin/time -v timeout 1200 ./dp3 20 3 1' > full.log 2> full.run
taskset -c 12 nice -n 5 bash -c 'ulimit -c 0; ulimit -v 11718750 && /usr/bin/time -v timeout 1200 ./dp4 21 18 1 14' > narrow.log 2> narrow.run
taskset -c 12 nice -n 5 bash -c 'ulimit -c 0; ulimit -v 11718750 && /usr/bin/time -v timeout 1200 ./dp5 22 20 1 19 22' > late.log 2> late.run
taskset -c 12 nice -n 5 bash -c 'ulimit -c 0; ulimit -v 11718750 && /usr/bin/time -v timeout 1200 ./dp5 22 21 1 21 22' > lastpair.log 2> lastpair.run
taskset -c 12 nice -n 5 python3 summarize.py full.log narrow.log late.log lastpair.log
```

Terminal sharding: run over the same full previous level, keeping one hash residue of each final unrooted key. Each such key includes all marked orbits, so no signed-coefficient cancellation crosses shards. `merge.py` adds exact endpoints and checks the summed rooted/unrooted census and the unit mass of each label. A two-shard run at size 10 reproduced the unsharded rational endpoints for labels 7–8 exactly. The size-23 command is:

```sh
taskset -c 12 nice -n 5 bash -c 'ulimit -c 0; ulimit -v 11718750 && for i in 0 1 2 3; do /usr/bin/time -v timeout 1200 ./dp6 23 23 1 22 23 4 "$i" > "sh23-$i.log" 2> "sh23-$i.run" || exit; done; python3 merge.py sh23-0 sh23-1 sh23-2 sh23-3 > merged23.log'
```

The four shards exited 0 in 7:07.50, 11:59.58, 9:43.66, and 5:27.82, respectively; peak RSS 6,893,352 KiB. Merging checked both total marked masses and the full census: 14,828,074 unrooted and 268,282,855 rooted shapes. The improved m=23 bound is `p_22-p_23 in [0.000265991, 0.006297528]` outward; the baseline lower bound is `-0.000136518` downward.

The optional `grid.cpp` is a separate exact-arithmetic, directed-rounding engine for the last pair. It stores a single integer per orbit, using scale `2^62`. At each edge, it multiplies nonnegative input by the exact rate numerator in `__uint128_t` and takes the integer floor after division by the exact rate denominator. Inductively each marginal array is componentwise below the true marginal; only one marginal is needed before the final two births. At the terminal step, it accumulates the signed difference of the two rounded arrays directly; separately track their integer masses `S_u,S_{u+1}`. The exact conditional contributions of the discarded mass to the difference lie in `[-(2^62-S_{u+1}), 2^62-S_u]/2^62`, since all conditional hub probabilities lie in `[0,1]`. Thus one can apply the usual signed classifier before adding this worst-case error, all exactly. The corresponding size-10 bounds contain both exact GMP endpoints for each classifier. This uses fixed-grid mass enclosures, not floating-point numerical propagation.

```sh
taskset -c 12 nice -n 5 g++ -std=c++20 -O3 -Wall -Wextra grid.cpp -lgmpxx -lgmp -o grid
taskset -c 12 nice -n 5 ./grid 10 > gridcheck.log 2> gridcheck.run
taskset -c 12 nice -n 5 bash -c 'ulimit -c 0; ulimit -v 11718750 && /usr/bin/time -v timeout 1200 ./grid 24' > grid24.log 2> grid24.run
taskset -c 12 nice -n 5 python3 summarize.py full.log narrow.log late.log lastpair.log merged23.log grid24.log > SUMMARY.md
```

The size-24 grid run certifies `p_23-p_24 ∈ [0.000249880, 0.005904360]` (improved bounds, outward-rounded) over all 743,724,984 rooted shapes. The exact missing numerator masses were 1,172,975,527 and 242,128,228 over `2^62`. A marginal enclosure for the last label follows from `p_24=p_23-(p_23-p_24)`: combining the improved intervals gives `[0.038848597, 0.185632462]`. The improved classifier thus certifies `U=24`, while baseline results cover `U=19`. `SUMMARY.md` records the best over the actual evaluated `(m,u)` pairs; the largest sizes focused on the frontier labels under the memory cap. At size 25 (>2 billion rooted shapes), even a single 64-bit integer per rooted shape would exceed 12 GB, before any indexing. A further size would need sharding or a fundamentally different representation, and the observed CPU time at size 24 already exceeded nine minutes.

`dp3` was the first corrected, cached implementation, and `dp4` added exact comparison caching and pre-birth attachment orbit grouping. Current `dp.cpp` also accepts the optional final label as a fifth argument; it reproduces either command by replacing the executable name. Full all-label size 20 requires 9,747,933,280 bytes for its array alone; the operating-system virtual cap includes the previous level's mapping even after its pages are discarded. The full run completed both classifiers at size 19 before aborting at the size-20 allocation. A narrow run avoids retaining old labels.

The validation log reads:

```
m=9: all 30 exact endpoints match (8 marginals and 7 differences, including symmetry).
m=10: all 30 exact endpoints match (8 marginals and 7 differences, including symmetry).
```
