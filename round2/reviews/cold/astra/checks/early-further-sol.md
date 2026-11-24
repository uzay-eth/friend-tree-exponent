# Cold audit: extension through label 11

**Verdict: holds.** The extension (`open/age-order/early/further/README.md`) certifies the six new gaps through `p10−p11` and positivity conditional on every prebirth tree of size 2–7. The size-eight whole-tree result is 19 positive of 23, with four unresolved. I checked the signed recurrence, independently tested it on small labelled cases, and reproduced the larger gaps on pinned source bytes.

## Multiplicity and cancellation

- **Unmarked seeds.** A free-shape representative is a state of the lumped Markov chain, not one of several equiprobable labellings. For each concrete target in that representative, add the full `shape mass × target rate / n` to the key of its completion. Repeated orbit images are summed. There is no division by the automorphism group, stabilizer, or number of rooted types. For the four-vertex star the center contributes `3/4` while each of three leaves contributes `1/12` to the same completion: that orbit has mass `1/4`. Every tree has an increasing construction after pruning to the edge and reversing the prunings, and every specified target in that construction has strictly positive rate. The shape enumeration thus covers every conditional prebirth state up to old-label symmetry.
- **Tagged states.** Rooting at tag one and marking tag two gives a complete invariant of the ordered two-tagged tree: the sorted recursive child data reconstructs the tree with both marks. Changing which tag supplies the root flips the difference; equality of both codes gives an actual tag-reversing automorphism. Decoding produces one vertex per occurrence of a child encoding, even if several children share the same encoding. The recurrence loops over every such target and keeps its full mass. In the independent fixture the canonical key is instead center-rooted with both tags colored, and the comparison recurrence retains labelled states while averaging successors.
- **No cancellation across different tagged graphs.** For distinct old targets `i,j`, the opposite order gives an isomorphism of the full final trees *with the two tag identities exchanged*, holding the old vertices fixed. This permits cancellation of the common product `r_i r_j`; it does not permit cancellation with a third completion just because that completion has the same free shape. The residual for an old edge is `(r_j δ_j−r_i δ_i)/[n(n+1)]` on the `i,j` placement; for a nonedge it is zero. Repeats are exact twins; a chain contributes `r_i/[(D_i+1)n(n+1)]` with no reversed order. The fixture builds the entire two-jump law by directed-arrow enumeration and compares its signed coefficient map after canonicalizing **ordered** tags and accounting for reversal. In a concrete fork completion of `T3`, the center-to-newborn chain and center-to-old-leaf placement give the same free graph but disjoint tagged reversal pairs: one contrast compares a middle and tip, the other compares two tips. The test fails if those contrasts are merged by their free shape alone. When the writer leaves different terms separate, multiplying each by its full signed interval is conservative even if additional tagged isomorphisms could tighten it.
- **Rounding.** The `lcm` makes each `r_i` an integer multiple of `1/lcm`; division by size gives the exact jump law. Taking weighted sums of integer ticks and using floor for the lower endpoint and ceiling for the upper remains outward even on negative values. Reversal swaps and negates the endpoints, also outward. The direct interval and pullback bound the same ultimate signed difference, so intersecting is valid. Grouped diagnostic rows are not used in the whole-tree lower endpoint. In a leaf–parent group the formula uses the direction `H(leaf,parent)` and the leaf weight `1/(2D_parent)`; between groups the two oriented weights are retained, then combined once for the whole tree.

## All-state terminal checks

After the two births, every graph has size at least four, so a degree-one tag always has a core parent. Until promotion or the first parent gain, its promotion clock has predictable rate `1/D`, exactly, and the parent-gain hazard is at most `D`; disjoint transitions can be split into independent proposal channels. Couple a rate-`D` dominating gain clock and require the exact promotion clock to beat it. This proves `1/(1+D²)` without freezing other rates. At promotion the tag starts at degree two with one leaf; multiplying by the audited `V(2,1)` is valid by the strong Markov property. The degree-one upper tail is the no-send bound with one immigrant per unit time from the still-leaf receiver, and `B=D−(L−1)` counts the permanent slots outside the tracked stock. A bare core cannot be a hub without another gain, whose first source is among its existing core neighbours; each upper charge uses zero immigration and `B=D−L`. For a core with leaves use the audited tag lower bound and the trivial upper bound one. If stock is already at the no-send cap, the code simply uses the all-time barrier at that stock; if degree is at or above the core cap it uses the family survival bound directly. These are true ultimate-event probability bounds after any finite pullback. The only domain caveat is that the new absolute leaf routine is not a two-vertex-edge terminal, and it is never required there.

## Exact replay

From the repository root, `timeout 20s taskset -c 13 nice -n 10 python3 -Bu round2/checks/early-further-sol.py` passed:

```text
200 free-shape weights through size 10; labelled through size 8; star orbit 1/4,3/4
141 physical endpoint rows and 38 exact whole-gap controls through size 7
fork negative control: same free shape, distinct tagged contrasts
```

The fixture (`early-further-sol.py`) pins SHA-256 of the extension script, extension note and imported terminal checker. The imported checker also matches the bytes in `408aa28`. The two extension files were untracked when reviewed, and I found no separately authenticated receipt of their larger outputs, so I used a single combined replay for the missing larger means and the size-eight count. The short writer command `timeout 25s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/early/further/check.py --horizon 5 --n 4 5 6` took `wall=2.92 rss=58332 exit=0`. The larger command was:

```sh
/usr/bin/time -f 'wall=%e rss=%M exit=%x' timeout 95s taskset -c 13 nice -n 10 python3 -uB -c 'from runpy import run_path; from fractions import Fraction as Q; w=run_path("round2/open/age-order/early/further/check.py"); expected={7:"[0.000818172405, 0.027991301188]",8:"[0.000438713999, 0.023641075524]",9:"[0.000205667071, 0.020338619058]"}; exec("for n in expected:\n    lo=hi=Q(0); positive=negative=count=0\n    for t,m in w[\"classes\"](n).items():\n        a,b=w[\"compare\"](w[\"decode\"](t)[0],5)[0]; lo+=m*a; hi+=m*b; positive+=a>0; negative+=b<0; count+=1\n    out=w[\"fmt\"]((lo,hi)); assert lo>0 and out==expected[n]; print(n,out,positive,count,negative,flush=True); assert n!=8 or (positive,count,negative)==(19,23,0); w[\"bound\"].cache_clear()")'
```

It returned `wall=88.90 rss=901024 exit=0`. Both passes gave the following exact-arithmetic intervals, displayed rounded outward:

| Gap | Interval |
|---|---|
| `p5−p6` | `[.005021411567,.055570905903]` |
| `p6−p7` | `[.002613919817,.042503762893]` |
| `p7−p8` | `[.001461596901,.033941586453]` |
| `p8−p9` | `[.000818172405,.027991301188]` |
| `p9−p10` | `[.000438713999,.023641075524]` |
| `p10−p11` | `[.000205667071,.020338619058]` |

The writer's `--scan --n 3 4 5 6 --horizon 5` returned `1/1, 2/2, 3/3, 6/6` positive (`wall=4.14 rss=59272 exit=0`); `--scan --n 7 --horizon 5` returned `10/11`, the sole unresolved interval `[-.000231327500,.037405850164]` on the three length-two-arm tripod (`wall=8.69 rss=130896 exit=0`). `--fixture --n 2 --horizon 5` gave the tripod at depth six `[.000385990574,.036818752400]` and the size-two tree at depth five `[.027991318747,.122321658633]` (`wall=6.90 rss=105020 exit=0`). The combined larger pass counted `19/23` size-eight shapes positive and none negative. No writer files were edited.
