> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; Claims 3 and 10; Gemini reviewed the base sequential law, with separate scopes for addenda. See the [release index](../../../../../README.md#claim-20) for exact scope and citations.

# Repair of the moving-cutoff marked forest comparison

Fable's Claim 20b objection (`fable-on-astra-20b.md (see the release review index)`) to the signed-union step of [the frozen diagonal note](../../forest-age/diagonal/README.md) at `6f26ab5` is correct. We replace that step, without modifying the original. The conclusions are:

* **Uniformly counted components:** for every deterministic `N=N(n)→∞` with `N≤(log n)^α`, `0<α<1`, the unnormalized marked count error is `o(d_N)` almost surely, where `d_N∼A N^-μ log N`. This proves the normalized full marked TV comparison to the moving **exact-lattice** stationary law. In particular it restores the board's comparison on its stated range.
* **Uniform surviving cores:** full marked TV to the moving exact-lattice stationary law holds on the smaller range `N≤(log log n)^α`, `0<α<1/2`; more generally it holds if a deterministic envelope `N_*≥N` satisfies `N_*² log(N_*+1)=o(log log n)`. On this range a replacement error is `O_path(n^-η/200)+O((log log n)^-4)` per original core. The frozen bound (3) and full surviving-core TV on the original log-power range are not proved here. Local convergence there still follows, in fact for any `N(n)→∞`.
* Fixed-shape fine limits are unaffected by the alphabet gap. See §5 for a capped one-coordinate proof of the transfer.

We use the mass norm `||λ||=sup_{|f|≤1}|λf|`, which is twice the usual TV distance for probability differences. All physical marks, including parent degree, entire ordinary top-rooted shape, and original core kinds or count `K`, are as in the frozen note. Optionally distinguish the sampled surviving core inside the component. No original degree stub other than the parent mark is required. A parent ratio `D_parent/N` or its logarithm is equivalent to the exact integer parent degree at fixed `N`.

## 1. Counterexample to the alphabet estimate

Set `h=floor(2 log₂ log n)`, `L=2^h`. Take a complete binary core skeleton of height `h`; at each bottom vertex attach independently one through four ordinary leaves. The top has an external, deleted parent of degree `N+1` (give it `N` other leaves). Every skeleton vertex is an original core of degree at most five. The top height is `h`, up to the choice of counting levels from zero or one; the frozen height cap eventually exceeds it. The ordinary size is at most `6L−1`, below the frozen size cap `U`; the parent is below its cap. Every fixture can be built chronologically starting with the parent–top edge: add a leaf at a target already present, using any one of its existing neighbours as a sender. Each such directed event has positive probability. If insisting on a particular founder convention at size three, the unique centre can be the high parent by adding one of its other leaves first.

The binary skeleton is distinguished by the original kinds. Its rooted automorphism group has exactly `2^(L−1)` elements (the recursive swap at each internal vertex). Therefore the `4^L` colourings yield at least `2^(L+1)` distinct admissible types, regardless of stabilizers. In fact their orbit counts satisfy `a_0=4`, `a_(h+1)=a_h(a_h+1)/2`, by unordered pairs of child types. The checker enumerates ordinary rooted representations and the independent colour orbits at heights 0–2. Since `L` lies between `(log n)^2/2` and `(log n)^2`, `T` grows faster than every power of `n`; the logarithm of the number of sign patterns is `T log 2`, not `log T`. The claimed union in the original proof is invalid; using individual coordinate bounds and summing them over this same alphabet does not fix it.

## 2. Cylinder tool and the corrected counted-component census

Use the frozen [downward cylinder estimate (2)](../../forest-age/diagonal/README.md#1-uniform-physical-cylinder-estimate-at-deterministic-sizes) with `η=2^-348`. For depth `k=o(log n)`, norm `n^{o(1)}`, and a family of at most `exp(n^{o(1)})` tests, its right side is `O_path(n^-η/200)` simultaneously. This uses the centered estimate, deterministic size, and normalization by the actual `C_n`; it does not require a convergence rate of `C_n/n^μ`. The relevant positive rate exponent `μ−11/20` is larger than `η/100`.

Fix `α<γ<1`. Admit a nonfounder component if its **ordinary** size is at most `J=floor((log n)^γ)` and its deleted parent's degree is at most `B=floor((log n)^b)`, with any fixed `b>8α` (increasing `b` is harmless). From a high parent read its low core children, explore each downward connected low piece, and stop a piece as soon as its ordinary size is known to exceed `J`. Include all real ordinary leaves of surviving cores in this size. A search reaching core depth `J` from a child has already found `J+1` vertices, so a downward depth `J+3` suffices, including full original degrees at the boundary. Reject the founder component. The alphabet has at most `B exp(O(J log J))=n^{o(1)}` types: labelled encodings bound rooted ordinary trees, and kinds add at most `2^J`; the count `K` is determined by the marked vertices. Optional exact original degrees of retained vertices add at most `N^J`, which is still harmless. For the uniform count, any signed sum has score norm at most `B`. The true sign union thus fits the cylinder theorem and the admitted marked error is `O_path(n^-η/200)`.

Under `π`, the unnormalized omitted count for volume `>J` is

```
≤ d_N ( stationary shape-TV error + 2N/J ) = o(d_N).
```

Indeed the existing [inert uniform-age approximation](../../forest-age/README.md) is TV on the whole ordinary shape, and the inert size at age `≤log N` is dominated by a two-founder Yule count, with mean at most `2N`. For the high-parent cap, [core moments](../../../../../../../reviews/cold/astra/checks/age-effective-sol.md#effective-local-convergence-versus-a-limit-without-a-modulus) give `π(c_o 1{D_o>B})≤C B^-1/8=o(d_N)`. Count unknown shapes by their tops at the high parent, not at all surviving cores; no multiplication by `K` enters either stationary estimate.

To get the finite remainder, first compare the total component counts. At each high parent count the low core children; set this score to zero if its core degree is greater than `c_* log n`. It reads downward depth at most two, has norm `O(log n)`, and its stationary cap error is polynomially small by the exponential core-degree moment. For deterministic `c_*` sufficiently large, eventually almost surely no physical vertex is affected: the own rate-one initiation count since birth dominates the number of new core neighbours, plus a fixed initial allowance. The Yule size epoch is at most `2 log n` eventually, and a union over the first `n` vertices of the Poisson upper tails on that deterministic time window is summable for sufficiently large `c_*`. Equivalently one may use the existing maximum core-degree bound. Add the possible founder component, at cost `1/C_n`. The total difference is `O_path(n^-η/200)` (choosing `c_*` sufficiently large for the stationary polynomial error as well). Subtract the admitted masses from the totals on both sides: the finite discarded count is its stationary counterpart plus `O_path(n^-η/200)`. The full mass-norm difference is then `o(d_N)`. In particular the finite component total divided by `C_n d_N` tends to one. Normalizing costs at most another constant factor because the masses are asymptotically equal. All test families are deterministic; the cylinder union can include all integer cutoffs below the envelope. The relative statement holds on every diverging schedule under that envelope, not with an error modulus at bounded `N` inferred from the stationary large-`N` approximation.

Combining with the existing [stationary marked coupling](../../forest-age/marks/README.md) gives the stated inert mixture and Pareto parent comparison in the original test metric: bounded tests uniformly Lipschitz in log relative degree and arbitrary in shape and `K`. The exact lattice law and a continuous Pareto law are not close in full TV on the parent coordinate. The shape and `K` projection does admit the TV comparison. The founder mark has vanishing conditional counted mass.

## 3. Why this is not a survivor argument

With subscript `N`, the unnormalized stationary weight of a component is multiplied by `K`; the total survivor mass is `π(2≤D_o≤N)→1`, not `d_N`. The cutoff `J` would need an estimate on `∫ K 1{|F|>J} dM_N`, not on `M_N(|F|>J)`. Neither ordinary-size escape from the stationary count approximation nor its unweighted TV error controls that weighted integral. For example a mass `1/K` at an omitted type tends to zero in count norm while retaining one unit of survivor mass. Comparing finite truncations of the **limiting** decorated tree does not by itself fix the moving whole-component statement: under a single `π` sample its components exhaust the infinite decorated tree as `N→∞`. For each fixed finite type set the stationary survivor mass in that set tends to zero. The exact marked law on finite types has no finite limiting probability on which to apply the fixed-alphabet pointwise-plus-mass argument. This does not refute TV between two *moving* laws.

Local convergence on any deterministic `N(n)→∞` follows by a different sandwich. For any fixed radius the uncut empirical marked-ball law tends to `π`. Restrict to balls of bounded maximum original degree first; for sufficiently large `n` such balls are retained. The omitted fraction can be made arbitrarily small under `π` and hence empirically. The surviving-core denominator tends to one by the same argument at radius zero. This applies to local rooted graph marks, not the full component or its distant cut-parent mark.

## 4. A survivor TV repair on a larger log-log range

The original radius proof uses equal lengths for old and recent segments. We separate them. Work under `π` with the same spatial-eigenmeasure backward sampling at time `−t` as in [the low-forest radius proof](../README.md#2-backward-window-a-quantitative-radius-bound-for-all-fixed-n). The focal old cell has count `Y_t` with `E[Y_t²|old]≤2e^(2t)`. Let `H≥1` be integer, `t=8N²+H`, `r=32N` and

```
ell = ceil(8 (t+r(N−1)+H)),       k = r+2 ell.
```

Use one low path of length `k` from a focal cell member. It has a single old-tree segment (possibly empty) and at most two recent segments. If the old segment has at least `r` edges, it contains an old low path of length `r` starting at the focal old anchor. That anchor is on the path whenever it actually exits the focal cell through old edges. Old degrees are nondecreasing. As in the original proof, a specified old path survives low with probability at most `exp(−rt/(4N))`, using the stopped Poisson domination and `t≥4N²`. The union is over at most `N^r` initially low paths. Its contribution to stationary mass, with the cell weighting paid by Cauchy–Schwarz, is at most

```
sqrt(2) exp((1−μ)t + (r/2)log N − rt/(8N)) ≤ sqrt(2) exp(−H).
```

Here `t≥8N log N` makes the `log N` term at most `rt/(16N)`, and `r/(16N)=2`. Otherwise a recent segment has at least `ell` edges, hence some occupied cell reaches depth `ell` from its old anchor. This also holds for a path contained wholly in the focal cell: each of its two segments relative to the turning point is no longer than a root-to-descendant depth. A relevant old anchor is either the focal old root, included even if high, or is within `r` through old low vertices. There are at most `2N^r` candidates. The uniform rate-one depth bound is `t^ell/ell!`; no independence between the cells or with `Y_t` is used. By Cauchy–Schwarz the stationary contribution is at most

```
2 exp((1−μ)t) (N^r t^ell/ell!)^(1/2) ≤ 2 exp(−H).
```

Indeed `t^ell/ell!≤(3t/ell)^ell≤(3/8)^ell≤exp(−ell/2)`, with `ell≥8(t+r(N−1)+H)` and `log N≤N−1`. Thus `π(R_N≥k)≤4 exp(−H)`, with **`k=O(N²+H)`** instead of a cubic initial radius. This is a radius bound on the stationary survivor measure per original core, including the deletion convention. It does not assert a bound on the maximum component in the finite tree.

For a deterministic envelope `N_*≥N` satisfying `N_*² log(N_*+1)=o(log log n)`, set `H=ceil(10 log log log n)` for large `n`, take `k` at `N_*`, and admit nonfounder components of core height less than `2k` and parent degree at most `B=U^32`, with `U=(N_*+1)^(2k+4)`. A component excluded for height has core diameter at least `2k`, and every core in it has radius at least `k`. All admitted components have ordinary size at most `U`. Now `log U=O((N_*²+H)log(N_*+1))=o(log log n)`: the hypothesis also bounds `log(N_*+1)=O(log log log n)`. Thus `U=(log n)^{o(1)}`. The corrected alphabet bound for pointed components is `BU exp(O(U log(U+N_*)))=n^{o(1)}`, and the signed score has norm at most `BU`. Depth `2k+3=o(log n)` and the true sign union therefore fit (2). This step is a separate truncation for survivor weights, not a reuse of `J`.

Stationary discarded survivor mass costs `4 exp(−H)` for height and `U π(c_o 1{D_o>B})=O(U^-3)` for cap by transporting each surviving core to the deleted parent. For the latter estimate each of its `c_o` child ports accounts for at most `U` possible retained cores. These errors are `O((log log n)^-4)` because `H≥10 log log log n` and `U≥3^(2k+4)`. To bound the finite discarded survivor mass, subtract the admitted survivor census from the total empirical surviving-core score `1{2≤D≤N}`. This includes the entire possible founder component in the discarded mass, not merely one core. The mass norm difference is bounded by the admitted error, the total error, and twice the stationary discarded mass. Divide by the respective survivor masses which tend to one; this yields the stated full moving TV, optionally pointed and with exact parent mark. All chosen scales can again be used simultaneously for `N≤N_*`.

## 5. Fine fixed shapes

One fixed ordinary shape `R`, including the core singleton, can be counted by one downward parent score: stop the exploration above `|R|` and cap the parent's core arity at `c_* log n` as in §2. The same argument works for the fixed-shape score multiplied by a prescribed bounded parent test varying with `n,N`. Its norm is `O(log n)`, depth is fixed, and the one-count founder cost is negligible. Thus the stationary-to-finite error is `O_path(n^-η/200)=o(N^-μ)` on any fixed log-power range. Combining with the frozen [fixed-shape sequential estimates](../../forest-age/fixed-shapes/README.md) retains the pair coefficient `2A/3`, the other inert occupations for `|R|≥2`, and the separate positive singleton coefficient. This argument does not pass to all shapes at cutoff-dependent volume by summing fine fixed-shape errors.

## Exact replay

```
nice -n 10 python3 -IB round2/open/low-forest/marked-repair/check.py
PASS binary alphabet: exact orbit counts 4, 10, 55; symbolic lower bound and reachable realizations
PASS 7848 exact size-capped censuses; 5340 strict survivor-versus-component losses
PASS asymmetric path split: 18 integer exponent controls
```

The finite check compares parent discovery to direct connected components, verifies ordinary-size admission, the weighted and unweighted remainder identities, and tests integer exponent inequalities in §4. The stationary tails and growing-alphabet convergence are analytic.
