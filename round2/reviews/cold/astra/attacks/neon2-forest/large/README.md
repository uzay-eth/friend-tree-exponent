# Independent stress test: Engine A with an arbitrarily slow size/cutoff ratio

**Recommendation to the official reviewer: HOLDS**, for Engine A of `ft2-neon2@6950241589978a7064ade0f6a2ee97cbbdcc6a34`, conditional on the cited Astra bounds. This covers the count and log-volume law on every deterministic large-cutoff subsequence `N≥n^γ`, `N/n→0`, for fixed `γ>0`. With the separate small-cutoff engine this is the required large-cutoff half of the all-deterministic-sequence claim. I do not certify Engine B here. No fixed-shape or `o(hub count)` error follows. Below is an independent reconstruction with a shorter proof of the old-child estimate.

## Clock and degree gates

Use the actual Yule size clock from `T₃`, with `n e^(−τ_n)→Y>0`. For each fixed width `ω>0` and integer `r`, use the deterministic bracket `[t_-,t_+] = [log n−(r+1)ω, log n−rω]`; impose the event `e^(rω)<Y<e^((r+1)ω)` only at the pathwise conclusion, never in the laws of the clocks. Gamma size mass makes the excluded boundaries null. Put `L=log N`, `b=4 log t_-`, `s=t_-−L+b`, and `δ=t_-^(−1/4)`. On the pinned regime `L≥γ log n`, `b=o(L)` and `s−b=log(n/N)−(r+1)ω→∞`, with no rate requirement. In particular `s≤t_-` eventually, even for oscillating cutoffs.

Localize at `s` on maximum core degree at most `J t_-`, for each fixed integer `J`. For an existing degree `d≤t_-²`, normalized finite-time Yule upper domination after elapsed `h≤t_+−s` gives probability at most `2^d exp(−e^(b−ω)/2)` of exceeding `N` at `t_+`. Here `E exp(p Geom(p)/2)≤2` follows from `(1−p/2)e^(p/2)≤1`. For births after `s`, use the rate-one Yule upper comparison starting at one **at the birth stopping time** and sum by compensation, not by conditioning on all future birth times. These bad-event bounds times `O(n)` are summable.

For `d>t_-²`, the finite-duration tagged construction gives normalized degree between a sum of `d−c` independent normalized single-founder Yule sizes minus a normalized immigration population, and a sum of `d` such sizes. Here `c≤J t_-`; the normalized immigration population plus the normalized unit is distributed as a normalized single-founder Yule size, because its unnormalized size plus one has rate equal to that size. Thus it has the same exponential-moment envelope; independence from the lower sum is unnecessary. The centered log mgf for normalized single-founder size is `O(θ²)` uniformly in elapsed time near zero, by the common exponential-moment envelope and mean one. Since `c=o(δd)`, applying Chernoff separately to the sums and subtraction term bounds relative error `>δ/4` by `O(exp(−c₁δ²d))`, summable after multiplying by `O(n)` at both endpoints. Exhaust the localization integers. With `H(t)={D(t)>N}`, the resulting sandwich is

```text
H_-={D(s)>(1+δ)e^b} ⊆ H(t_-) ⊆ H(τ_n) ⊆ H(t_+) ⊆ H_+={D(s)>(1−δ)e^(b−ω)}.
```

Every selected vertex is already a core at `s`. By the **common-event** intermediate degree theorem, applied at random size `|T_s|∼Ye^s` with threshold rounding, `K_±=|H_±|` divided by `A X (n/N)^μ` tend to `(Y e^(-(r+1)ω))^μ` and `(Y e^(-rω))^μ`, respectively. The conditions for that input hold because `e^b→∞` and `e^b/|T_s|→0`. Thus `K_±→∞` at whatever rate the size/cutoff ratio allows.

## Old children: no retrospective Poisson law needed

Here is a direct alternative to pinned Step 3.2. On the common pathwise event, `C(T_u)=O(e^(μu))` and `c_max(u)=O(1+u)` for **all** `u≥0` (extend the eventual bounds over finite compact time by nonexplosion). Take `M=ceil(4b/μ)=o(L)`. When `s≥M`, all hubs promoted by `s−M` number at most `C(T_(s−M))`; their total old core degree is bounded pathwise by

```text
C(T_(s−M)) c_max(s) = O((1+s)e^(μ(s−M))) = o(L K_+),
```

since `(1+s)/L=O(1)` and `e^(μ(s−M))/K_+=O(e^(μ(b−M)))→0`. If `s<M`, there are no such nonfounder cores. The founder's own possible component always costs at most one; any separately charged founder core-degree term is `O(1+s)=o(LK_+)`. No conditional inference on the already observed old selected count is used.

For young cores, do not use a bound only polynomial in `log n` in an all-`n` union. Group `j≤t_-<j+1`; bin `s` upward to integer `m` with `s≤m<s+1` and use a duration cap `M_j=ceil(4(4log(j+1)+1)/μ)+2`. There are `O(j)` bins. A nonfounder promoted in `(m−M_j,m]` starts at core degree one; at its promotion stopping time the subsequent ring count is dominated by `Pois(M_j)`, and thus so is excess core degree at `m`. This bounds the summed expectation over **all** young candidate cores, by the Markov property and promotion compensation (their expected number is `O(e^(μm))`). The truncated tail above `1+10M_j` has expectation `O(e^(μm−10M_j))`, from the Poisson tilted first-moment bound. Hence a threshold `e^(μ(m−4log j))` has probability `O(j^(4μ−40/μ))`; the union over bins and integer `j` is summable. Returning to the selected young hubs adds at most `(1+10M_j)K_+` below this truncation and `O(K_+)` above, both `o(LK_+)`. If a duration extends below time zero, treat the founder separately. This proof covers very slow `s−b` and requires no regularity in `n↦N`.

## Marks, mismatches, and the random epoch

Given the physical history at deterministic `s`, attach independent inert branch drivers to the rate-one clock points of each selected core. A newly promoted leaf branch contains no selected core and branches from different promotions are disjoint. One can use the inert internal directed channels on the newborn subtree and a separate child-to-parent channel with rate `1/D_child` for the external change. The latter does not alter the inert mark; a parent-to-child send can disrupt it. Preassign fresh trajectories also for unused points, and use independence under the **unmasked** conditional law at `s`. Any joint construction for each `n` can be put over the same physical law by sampling its extra data from its conditional distribution given the physical trajectory, independently across the countable uses if desired. This does not require the actual physical processes at different `n` to be independent.

On the eventual core-degree cap up to `t_+`, stop at its first violation using a deterministic multiple of `t_-`. Given the initial history, all selected hub-to-core sends up to the stop have integrated intensity bounded by `O_J(K_+ L t_-^(−3))`: `D_v≥(1−δ)e^(b−ω)` and `c_v=O_J(t_-)`. The count is stochastically dominated by the corresponding Poisson count via thinning or the counting exponential supermartingale with the deterministic-in-that-history bound. Each such send costs **two** marks at worst: the unused mark at its own clock point, and the first disruption of an older mark. The oversized-mark count is independently marked Poisson as a marginal given the starting history, of mean `O(K_+ L t_-^(−4))` by the ideal Yule size upper bound at `t_+`. Neither count needs to be Poisson conditional on future localization, survival, final size, or the hub sandwich, and no independence between these two counts is needed.

The old-child bound and the upper hub sandwich imply that **no** post-`s` promoted branch is internally cut at `τ_n`, including at its root: its root was leaf or absent at `s` and every descendant was absent. In fact this applies through `t_+` by monotonicity. So a good arrival is one exact ordinary component when its parent is deleted; even a bad arrival can produce at most one component in that branch. Pre-`s` promotions whose parent is deleted cost at most `Σ_(H_+) c_v(s)`, plus the possible founder component. For the upper log-volume count use qualified early marks at their **younger** age at `t_-` and count **all** arrivals in `(t_-,t_+]`; for the lower bound take only arrivals by `t_-` from `H_-`, qualified at their **older** age at `t_+`. Subtract the combined mismatch budget in the lower bound and add it in the upper. This takes all decisions on membership/cutoff at the actual random size, and applies monotonicity to the *unpruned, coupled* volume only.

For each of these mark counts the conditional mean is the length integral of the appropriate inert event probability times the number of parents. For totals the length is `L−b+O(ω)`. For the CDF at any fixed `0<x<1` the integral is `xL+o(L)`: under the upper and lower Yule count comparisons, at age `u=zL`, `log V_u/L→z` in probability for each fixed `z>0`. Change variables, use bounded convergence off `z=x`, and note `b=o(L)`. The late-arrival upper term has mean `K_+ω=o(K_+L)`. Uniform additive Poisson deviation at accuracy `ρ K_± L` has probability at most `2 exp(−c_ρ K_±L)` for each fixed `ρ>0`; the small-mean error counts satisfy the same type of estimate. For every `ρ,γ` fix deterministic `K₀` large enough that `c_ρ γ K₀>3`. Apply the conditional bound **on the starting-history event** `K_±≥K₀`, then use `K_±→∞` pathwise to remove that event eventually. This is summable over all `n` even when `n/N` is slower than every iterated logarithm; no tail in `K_±` alone is used. Take countable intersections of bracket widths, tolerances and interior rational `x`; squeeze widths down to zero. Monotonicity of the CDF in `x` gives uniform convergence on `[0,1]` from the continuous uniform limit, including endpoints by interior bounds. No monotonicity in the cutoff or observation size is assumed. The precise quantifier proved is for each deterministic sequence, on an event of probability one. This does not silently exchange the two quantifiers.

## Small local qualifications

- The pinned statement that `s` has a bin among `{0,…,j+2}` is fine eventually since `L≥γ log n≫b`. The single notation `ε` for bracket and error tolerance may be split into `ω` and `ρ` as above; either the fixed-bracket little-oh follows by separate tolerances or a fixed `O(ε)` error suffices for the final squeeze.
- “No proper descendant” undersells the useful fact: the newborn branch **root also** is absent from the deleted set. Conversely a branch created by an *old* promotion can already contain a deleted core, and it is not entitled to a full-volume comparison. Using the top-core count and the old-child budget handles these cases.
- Finite-shape or false conditioning controls below deliberately fail. They are failures of shortcuts, not counterexamples to the pinned count or CDF statements.

## Exact replay

```text
$ taskset -c 8 nice -n 10 timeout 8s python3 -I -B round2/attacks/neon2-forest/large/check.py
PASS: 41113 exact nested-deletion, fresh-branch and old-top rows
CONTROLS: a late birth is needed in the upper CDF; one send can spoil two marks
CONTROL: retrospective no-return selection changes an ideal mark from 1/2 to 9/17
```

No rival code is imported. The enumeration checks the top-core decomposition including singleton survivors, under the full old-core deletion gate. In the rational selection fixture the two first-arrow choices have identical physical total size after the first birth, so specifying the second source does not itself change their odds. No-return selection then weights these choices differently because their return rates are `1/3` and `1/4`. The almost-sure assertions are analytic.
