LOCAL CRITERION: CORRECT

Hostile referee: GPT-family (GPT-6.2 Sol Ultrafast), reviewing the Claude-family lead's §4.2. The finite implication and the C++ PASS decision are correct. Two non-certificate claims need correction: the zero-mass ratio display and the sufficiency of absolute uniform convergence in the final stationary extrapolation.

## Line-by-line audit of `COUPLING.md`

I use **core** to mean a vertex of degree at least 2 in the base tree, and **parent** of a leaf to mean its unique neighbour, regardless of birth labels. Sections 1–3 provide the alternating-injection order and the extension from a coupling at size `u+1` to stochastic domination of the limits.

| Lines | Verdict and derivation |
|---|---|
| 181–187, rates | An original vertex `w` receives one term `1` from each new leaf directly attached to it; in the first step this affects only `x`. The changed old degree is `D_x+1` at `x`; each old neighbour loses `1/D_x−1/(D_x+1)=eps_x`. No double change at `x` since there are no loops. The new vertex has sole neighbour `x`, so its rate is `1/(D_x+1)`. Multiplying the four possible second rates (`r(x)+1`, `1/(D_x+1)`, `r(y)−eps_x`, `r(y)`) by `r(x)` gives the four weights, divided by `n(n+1)` from the two choices of a uniform vertex. Rates sum to tree size by summing the contributions from every chooser. A base leaf of a core `p` has `D=1`, rate `1/D_p`, and `rho=1/(2D_p)`. |
| 188–192, shape and piece 1 | For distinct base attachment sites `S_T(x;y)` is rooted at the only new leaf on `x`. If `y=x`, there are two sibling leaves at `x`; their transposition is an automorphism. Their common full outcome mass is used once. |
| 193–194, piece 2 | `(x,y)` viewed from the older root and `(y,x)` viewed from the younger root have identical unlabeled edges and root site after exchanging the names of the two new labels. For nonadjacent distinct sites neither second rate changes, so the weight product is symmetric. The reverse rooted shape is treated too (the argument over ordered pairs). |
| 195–199, piece 3 | The same exchange works for adjacent sites. Subtraction gives `r(y) eps_y − r(x) eps_x = rho(y)−rho(x)`. There is one diagonal minimum for each orientation. For `rho(z)>rho(x)`, the excess lies at `S_T(x;z)` on the older side, and at `S_T(z;x)` on the younger. Equality leaves no residual, whether a core–core or core–leaf edge. Shape coincidences across outcomes do not double-count: the residuals and diagonals are accounted for as masses of individual outcomes, then added by shape. |
| 200–202, piece 4 | When the second site is the new vertex, the final tree gains the path `x–a–b`. The older root `a` is degree 2 with leaf child `b`. The younger root `b` is a leaf with degree-2 parent `a` attached to `x`. The shared outcome weight is `U_x=r(x)/(D_x+1)=D_x rho(x)`. These events are disjoint from pieces 1–3. |
| 203–209, low endpoint a leaf | At high endpoint `p`, `S_T(x;p)` replaces the old leaf branch `p–x` by `p–x–(root)` and adds a new leaf at `p`. `K_T(p)` retains `p–x` and adds `p–a–(root)`. Map `p` and the remainder of `T\{x}` identically, old `x` to `a`, new leaf at `p` to old `x`, and root to root. Summing the edge excess over old leaves of `p` costs `L_p(rho(p)−1/(2D_p))^+ ≤ L_p rho(p) ≤ D_p rho(p)=U_p`. This spends only younger `K_T(p)` mass and never touches the older universal shape. In the younger residual rooted on `p`, the other new leaf extends `x`. Exactly one old leaf at `p` stops being a leaf, so there are `L_p−1` other leaf neighbours. |
| 210–212, high endpoint a leaf | Older root at core `x`: the parent has its `D_x` original neighbours plus the root, so there are `D_x` others. The other new leaf is attached to `z`, a neighbour of `x`, and is not a neighbour of `x` itself; the promotion of `z` does not remove `z` from the count. Younger root at old leaf `z`: its parent `z` has one other neighbour `x`, which is a core, so there are zero other leaf neighbours. The sum over leaves is `L_x(1/(2D_x)−rho(x))^+`. |
| 213–214, both cores | The older root at `x` again gives `D_x` other neighbours. The younger root at `z` gives `L_z` other leaf neighbours: `x` receives the other new leaf, but was already non-leaf and stays non-leaf. `z` gains only its own root directly. If both new leaves are directly adjacent to the root's parent, that is the separate case `x=y`, already coupled diagonally. |
| 216, order test | Match the younger leaf root to the older leaf root and their unique children to each other. At the odd-depth parent pair inject the `k` older children (other neighbours of its parent) into `k` younger leaf children. These pairs are at even depth; the younger member has no children, so the next injection is empty. The order follows irrespective of the unmatched branches. |
| 218–220, formula | The demand sum counts every positive core–core difference once, from its low endpoint, and every low-core/high-leaf edge once. Supply counts every positive core–core difference once, at the high endpoint's leaf-count level, plus precisely the low-leaf/high-core residuals at `L_p−1`. A low-core/high-leaf residual is not counted as supply, but is retained for the universal step. For `L_p=0` the low-leaf supply is zero and no negative level is needed. More generally, for `n≥3` no edge joins two leaves. Founders are just vertices of `T`; neither their degrees nor their contribution to any of these cases depend on labels. |
| 221–222, mixing and Hall | Weight each conditional bucket by `P(T)` and divide by the common `n(n+1)`. At fixed `n` there are finitely many shapes and levels. Process demand from high levels down. Before a demand at level `k` is processed, every previous unit at level `≥k` has used supply from the tail `j≥k`. The tail inequality ensures there is enough left for all demand at level `k`. This argument allows fractional splitting and use of supply from different base trees; the order comparison is between final rooted shapes and does not require common base trees. Conversely, the tail inequalities are necessary for routing along these level arcs, as every demand from the tail needs supply from the same tail. |
| 223–224, universal balance | After pieces 1–3's diagonal matches and the `K_T(p)` reserve, equal masses have been removed from both total-one marginals. On the older side only level demand and universal shapes remain. Routing the entire demand removes the same mass from each side. The younger residual is the unused reserve, low-core/high-leaf residuals, and unused level supply: all leaf-rooted. Its mass equals the universal mass exactly. A universal root matches its own leaf child to the younger root's only child and the odd injection is empty. Any mass-matching allocation therefore works; §3 extends the resulting coupling. |
| 225, `K≤1` | Correct: every demand has `D_x≥2`, so its tail is unchanged below 2, while the supply tail can only increase. |
| Boundary `n=3` / `u=4` | `T_3` is a path, with `rho(center)=1/3` and `rho(ends)=1/4`. The two residuals cost `2(1/3−1/4)=1/6` in local units; `U_center=2/3`, so reserve suffices. There is `1/6` of supply at level 1 and nothing at levels `≥2`; no demand. In actual probability units the reserve is `1/72` and the older universal mass is `7/72`. `u=1` is the initial edge (symmetric rooting); `u=2` uses the size-3 path and `u=3` uses the direct size-4 coupling at the start of §4. §4.2 starts at `u=4` and does not invoke a nonexistent core at size 2. |
| 227–234, exact claims | Confirmed exact formulas, tree law, and PASS for the archived `u=4,...,21` run. For `u=4` there is no defined ratio with positive denominator. For `u=5` the maximum defined ratio is 0 at `K=2`, but the program prints a sentinel `K=0`; for `u=6,...,21` the maximum is reported at 2. File names and placeholders in the draft need updating. |
| 236–239, simulations | Read `attacks/sd/stationary_hall.py` but did not replay the snapshots or independently verify numeric simulation claims. The script pools local masses over trees within a file, clips levels above `2^20` into the last bucket, and excludes zero-supply levels from its ratio scan. The representative printed ratios cannot by themselves imply a bound at every level, and `dem>0, sup=0` would be missed by that scan. |
| 240–242, stationary implication | Absolute uniform convergence of the normalized means, plus strict inequalities at each fixed level (even with a fixed *relative* limiting margin), does not generally prove simultaneous finite-size inequalities. If limiting gaps tend to zero as `K` grows, the absolute errors can still exceed those gaps at moving levels. The finite criterion is unaffected. See replacement text below. |

## Independent exact computation and comparison

`referee.py` uses `Fraction`. An unrooted code is the minimum of the sorted-child parenthesis codes over **all** choices of root. It starts from the single edge, enumerates every attachment site in a representative at every step, and adds the probability to its resulting free shape. Tree isomorphisms preserve all site rates, justifying this aggregation even if the original labelled trees have different numbers of histories. At each size it checks total probability one and compares the free-tree count and the displayed values against the archived C++ log. All tail values through `K=n` are in `exact.json`; both tails are zero above `n`. The comparison below uses the lead's display convention for `u=4,5` only.

| Pair | `n` | Free trees | Python = log `E dem_2` | Python = log `E sup_2` | Python = log largest ratio | Display `K` |
|---|---:|---:|---:|---:|---:|---:|
| (4,5) | 3 | 1 | 0.000000000 | 0.000000000 | 0.000000 | 0 |
| (5,6) | 4 | 2 | 0.000000000 | 0.166666667 | 0.000000 | 0 |
| (6,7) | 5 | 3 | 0.031250000 | 0.150000000 | 0.208333 | 2 |
| (7,8) | 6 | 6 | 0.054490741 | 0.185949074 | 0.293041 | 2 |
| (8,9) | 7 | 11 | 0.077596826 | 0.211563533 | 0.366778 | 2 |
| (9,10) | 8 | 23 | 0.095623021 | 0.235327033 | 0.406341 | 2 |
| (10,11) | 9 | 47 | 0.112433668 | 0.257921609 | 0.435922 | 2 |
| (11,12) | 10 | 106 | 0.127598587 | 0.278730118 | 0.457785 | 2 |
| (12,13) | 11 | 235 | 0.142335163 | 0.299195565 | 0.475726 | 2 |

## Explicit coupling and planted error

Built the transport using the proof's diagonal, swapped, adjacent-minimum, reserve, level, and universal allocations, **not** a generic max flow. The exact laws come from evaluating the actual rates afresh after each first attachment; a separate traversal of all labelled histories from the initial edge independently checks both marginals. The program checks the predicted adjacent difference and counts levels in the **final rooted trees**, and evaluates `≥` on every positive plan entry by memoized alternating recursion with bipartite augmenting paths. Full sparse plans and full marginal laws are in `exact.json` (fraction strings, root codes as above).

| `u` | Older / younger shapes | Positive plan entries | Labelled histories for marginal cross-check | Plan mass |
|---:|---:|---:|---:|---:|
| 1 | 1 / 1 | 1 | direct | 1 |
| 2 | 2 / 1 | 2 | direct | 1 |
| 3 | 3 / 2 | 3 | direct | 1 |
| 4 | 6 / 4 | 8 | 24 | 1 |
| 5 | 13 / 9 | 15 | 120 | 1 |
| 6 | 29 / 20 | 42 | 720 | 1 |
| 7 | 68 / 48 | 103 | 5,040 | 1 |
| 8 | 163 / 115 | 256 | 40,320 | 1 |

**Planted error:** drop the reserve match, retaining the original demand formula and treating the corresponding younger shapes as available to the universal step. At `u=4` the scalar inequalities at `K≥2` still pass, but the explicit construction only transports `71/72`; both exact marginals are short by `1/72`. The older deficit code is `(((()())))`, the younger deficit under the program's iteration order is `(((())()))`. This is failure of the modified construction, not failure of stochastic domination or a claim that no other transport exists.

## `local_scalar.cpp` code audit

| Lines | Check |
|---|---|
| 23–30, 54–95 | The tree key is the recursively sorted child code at a centre (or sorted pair of centre branches under a virtual root). Simultaneous leaf peeling leaves either one final-layer vertex or two adjacent ones. The virtual-root case has two more bits and reconstructs the edge when decoded. The DP sums equal keys; the fixed limbs are sufficient over the checked range (next total denominator bounded below the explicit signed-width guard). |
| 96–115, 122–131 | With `J=lcm(1,...,n−1)`, each real degree divides `J`. The `wt` numerator is `J*r(w)` and propagation uses scaling `Den[n+1]=Den[n]*n*J`. It iterates **all** sites, so counts orbit multiplicities even when several vertices lead to the same free shape. The initial key 110100 is a virtual root above the two founders; the reconstructed starting tree is an edge. The census and total checks are additional to the exact arithmetic. |
| 133–148, 161–162 | Write `L1=lcm(1,...,n)`, `L2=lcm(1,...,n+1)`. Each degree divides `L1`, and adjacent integers `D,D+1` are coprime and individually divide `L2`, so `D(D+1)` divides `L2`. Line 144 produces `L1*L2*rho`. Line 147 produces `L1*L2/(2D)`. Line 148 applies the sign and uses levels `D` and `L_f−1` respectively. Skipping `L_f=0` prevents negative indexing; equality adds zero. Dividing the final weighted numerators by `Den[n]*L1*L2` gives the expectations of the local terms, before the further conditional probability factor `1/[n(n+1)]`. |
| 149–154 | The check `z<x` excludes only the reverse orientation of each core–core adjacency (no loops). Either strict sign puts the magnitude in the smaller-`rho` vertex's degree bin and the larger-`rho` vertex's old leaf-count bin. Integer bins are nonnegative; each bin is bounded by `(n−1)L1 L2` since, expanding the leaf multiplicities, each of the `n−1` edges contributes at most 1 before scaling. This is well inside signed 128 bits. The two 64-bit GMP multiply-adds reassemble the value on this Linux LP64 build. |
| 155–165 | All nonzero bins are below `n+1`. The descending loop accumulates the tails and directly checks `cd>cs` as integers. A positive numerator with zero supply sets `ok=false` and reports an infinite ratio. The reporting defect is that zero ratios with positive denominators are not registered, and all-undefined ratios have the same sentinel display. This does not alter PASS. |
| 120–129, input envelope | Outside the checked range, `UHI` is unchecked: `FR[have+1]` has no entry beyond base size 25. Requesting `u≥27` can index it out of bounds if enough limbs and memory permit reaching that point. Six limbs suffice for sizes through 23; seven suffice for sizes 24 and 25, but not for every higher size (the formal size-26 denominator already needs 471 bits). The large-size capacity is unrelated to the archived `u≤21` run. |

## Exact replacement text

1. At the first use of “core” in §4.2, add: “A core is a vertex of degree at least 2 in the base tree `T`. The parent of a leaf means its unique neighbour in the tree being discussed, regardless of birth labels.” Clarify the ambiguous line 212 parenthetical as “its parent `z` has zero other leaf neighbours.”
2. For lines 227–234, as far as the specified log is concerned: “`attacks/sd/local_scalar.cpp` computes the two expectations exactly from the unrooted law. In `attacks/sd/local_scalar_4_21.log` every pair `(4,5)` through `(21,22)` passes. At `u=4` both tails are zero for every `K≥2` and there is no defined ratio. At `u=5` the maximum defined ratio is zero at `K=2`. For `u=6,...,21` the maximum is attained at `K=2`.” Retain the existing numeric table with its placeholders removed (or extend it and cite the corresponding additional log separately).
3. For lines 239–242: “The sampled ratios below 1 suggest a stationary version of the criterion. To deduce eventual domination one must control all finite-size levels simultaneously. For example, it would suffice to prove a uniform relative gap at stationary levels with positive supply and prove that the normalized finite-size demand and supply converge there with errors uniformly `o(s_K)`, where `s_K` is the stationary supply tail. Levels with `s_K=0` require a separate finite-size inequality. Neither the gap nor the required uniform tail control is established here.” For why absolute convergence is insufficient, let `s_K=2^{-K}`, `d_K=s_K/2`, and at stage `m` use `s_{m,K}=s_K`, `d_{m,K}=max(d_K,2s_m)` for `K≤m` and `d_{m,K}=d_K` otherwise. These are nonincreasing tails with uniform absolute errors at most `2s_m→0`, but `d_{m,m}=2s_m>s_{m,m}` at every stage. The same obstruction works with finite-support tails by setting both stage tails to zero for `K>m`.
4. For the ratio display in the C++ tail loop, after the exact PASS check, replace line 159 with:

   ```cpp
   if (cs == 0) {
       if (cd > 0 && worst >= 0) { worst = -1; worstK = Kk; }
   } else if (worst >= 0) {
       mpq_class r(cd, cs); r.canonicalize();
       if (worstK == 0 || r > worst) { worst = r; worstK = Kk; }
   }
   ```

   Print “undefined (no level has positive supply)” without a `K` when `worstK==0`, “infinite at K = …” when `worst<0`, and the numeric maximum and level otherwise. To guard the separate input-envelope issue, immediately after the declaration of `FR` add:

   ```cpp
   if (uhi > static_cast<int>(std::size(FR))) {
       fprintf(stderr, "UHI too large: the free-tree census ends at base size %zu\n", std::size(FR)-1);
       return 2;
   }
   ```

   Replace the build comment about limbs with “`-DNLIMBS=7` is needed for base sizes 24–25; this program's census stops at base size 25.”

## Commands actually run

From `round2/children/refL`:

* `taskset -c 13 nice -n 5 round2/.venv/bin/python referee.py` — exit 0 on the final run; nine `MATCH` rows as above, all eight exact plans of mass one, and planted-error line `plan mass 71/72, exact marginal deficits 1/72 each`. An earlier diagnostic incorrectly asserted that the two same-edge residuals must be incomparable; that assertion failed and was replaced with the actual incomplete-plan marginal check. The reserve mutation makes no claim about whether a different routing could repair the plan.
* `TMPDIR=round2/children/refL taskset -c 13 nice -n 5 g++ -O2 -std=c++17 ../../attacks/sd/local_scalar.cpp -lgmpxx -lgmp -o local_scalar` — exit 0, no output.
* `taskset -c 13 nice -n 5 bash -o pipefail -c './local_scalar 4 21 | diff -u ../../attacks/sd/local_scalar_4_21.log -'` — exit 0, no diff (entire 18-line log reproduced byte for byte).
* A pinned Python integer-bound calculation printed denominator bit lengths `365, 398, 434, 471` for sizes `23, 24, 25, 26`, respectively. No large simulation or large DP beyond the specified log was rerun.
