# Pre-placement adjacent gaps and a larger exact early prefix

**Certified:** `p1=p2>p3>…>p11`, for the original edge seed. The conditional gap before the two placements is strictly positive on **every tree of size 2 through 7**. The general conditional claim remains open: the size-eight scan proves positive gaps on 19 of its 23 shapes and neither sign on the other four. A straddling interval is not evidence of a reversal.

## Exact split and what a transfer would need

Write `I(u)` for the ultimate-hub probability difference of middle and tip after the serial placement at `u`. Write `H(u,v)` for first minus second after the two old-target placements at `u,v`. Conditional on the pre-placement tree, multiply the probability gap by `n(n+1)`. The numerator is the sum of

```
S(u) = w(u) [ I(u) − (1/D_u) Σ_{v~u} H(u,v) ],
w(u) = r_u/(D_u+1).
```

Indeed `H(v,u)=−H(u,v)`, while the adjacent target correction is `r_u/[D_u(D_u+1)]`. Nonedges and repeats cancel exactly. For a core-grouped score, assign every old leaf to its parent, then combine the two appearances of each within-group edge *before* multiplying by an interval. In particular a leaf `a` with core parent `o` contributes, jointly with the parent's share on that edge, `(r_o δ_o−1/(2D_o)) H(a,o)`. Each between-group edge is kept with its two different weights; on summing groups those shares reduce to the already paired unordered-edge coefficient. The certificate for the whole tree uses that pairing, not a sum of independently rounded group lower endpoints.

A sufficient pointwise transfer would compare `I(u)` with the neighbour average of `H(u,v)`. The code tests the separate vertex and core-grouped hypotheses as well as the whole-tree hypothesis. They are strictly stronger tests: even for small trees where the whole gap is certified positive, the vertex and group intervals often straddle zero. A fixed-tag reversal after specified placements does not sign `S`, and a negative `S` alone would not give a whole-tree counterexample. The present bounds yield neither a general transfer inequality nor a negative conditional gap. The first difficult size-seven shape is the tripod with three arms of length two; five jumps leave it unresolved, but six certify `[.000385990574,.036818752400]`.

## All-time tails and exact recursion

Reuse the [physical terminal bounds](../README.md) for cores and for no-send probabilities. There is one additional lower tail: if the focal tag is still a leaf and its parent has current degree `D`, its promotion rate before the next degree gain of the parent is `1/D`. The competing parent gain rate is at most `D`; all other transitions leave that degree and the leaf status in this race unchanged. By predictable-rate comparison, the probability of promotion before a parent gain is at least `1/(1+D²)`. At promotion the new core has one leaf and degree two, so multiply this by the conservative tracked-core value `V(2,1)` at cutoff 12. This does not require the exterior to be frozen. For the leaf upper tail use the original parent's permanent slots and keep the tag's own immigration. For a bare core upper tail use a union bound on the first incoming arrow, with zero immigration for each neighbour comparison. These bounds are on *ultimate* success, and both upper and lower bounds are subtracted with their proper signs.

For an exact finite pullback canonicalize the tree by rooting at tag one and marking tag two. Sorting rooted child encodings compresses all equivalent labelled histories; rerooting at tag two identifies exchanged queries and flips interval endpoints. If the two encodings are the same, symmetry proves gap zero. Each physical transition is integrated over its holding time to the exact `r_i/n` jump law. At each step take the interval expectation, round both endpoints outwards to integer ticks of size `10^-12`, and intersect with the direct terminal interval. The recursion uses a common denominator from the lcm of sending degrees. The outer tree-shape law and signed placements are computed as rational numbers, not sampled. The pre-placement shape law is independent of old labels, since both tags have yet to arrive. Every free tree has a positive weight by pruning one leaf at a time back to the initial edge; shape generation therefore also covers every reachable conditional configuration up to isomorphism.

## Results and replay

Combining the earlier result with five more births per placement certifies the following additional gaps (outward-rounded):

| Gap | Lower | Upper |
|---|---:|---:|
| `p5−p6` | .005021411567 | .055570905903 |
| `p6−p7` | .002613919817 | .042503762893 |
| `p7−p8` | .001461596901 | .033941586453 |
| `p8−p9` | .000818172405 | .027991301188 |
| `p9−p10` | .000438713999 | .023641075524 |
| `p10−p11` | .000205667071 | .020338619058 |

For the whole-tree scan with five extra births, positive counts by sizes `3,4,5,6,7,8` were `1/1, 2/2, 3/3, 6/6, 10/11, 19/23`. The separate six-birth tripod completes size seven. Size two is the deterministic edge used for the earlier `p3−p4` certificate.

```
timeout 55s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/early/further/check.py --horizon 5 --n 4 5 6 7 8
timeout 95s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/early/further/check.py --horizon 5 --n 9
timeout 85s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/early/further/check.py --control --scan --horizon 5 --n 3 4 5 6 7 8
timeout 45s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/early/further/check.py --fixture --n 4
```

The control cross-checks exact shape weights against labelled parent-vector enumeration through size six and compares leaf endpoint bounds through parent degree 64. The size-ten five-birth prefix attempt did not finish within 95 seconds; no `p11−p12` claim is made. A route beyond this scan needs a sharper signed comparison across the serial and parallel completions, rather than treating each post-placement reversal as a pre-placement counterexample.
