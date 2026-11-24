> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../../README.md#claim-18) for exact scope and citations.

# Early adjacent hub probabilities

For the original edge seed, with `p_i=Pr(Z_i>0)`, this proves

```
p1 = p2 > p3 > p4 > p5.
```

These are comparisons of ultimate hub probabilities. The finite passes below only expose a bounded number of physical births; their terminal bounds account for all future time. They do not, together with eventual adjacent ordering, fill all the intervening labels.

## A no-send bound, without an age cutoff

Fix a directed edge `x→u` with `x` a core. Track `j` of the leaf neighbours of `x` other than `u`. Let `B` be the number of permanent incident slots outside those initially tracked (including `u`); then `D_x≥j+B` throughout the comparison, also after tracked leaves are promoted, and `B≥1`. Before the first send along `x→u`, track the new leaves produced by the living tracked leaves. They give gain rate `j`. If `u` is a leaf, it remains one until that first send and its own clock adds tracked immigrants at rate one. If `u` is a core, use zero immigration. Track loss is at most `j/max(2,j+B)` and the specified send rate is at most `1/max(2,j+B)`.

First take the probability barrier `f(j)=min(1,2/(j+1))`, with value one on the send event. For `j≥1`, its drift is bounded above by

```
j [f(j+1)−f(j)] + j/(j+1) [f(j−1)−f(j)] + [1−f(j)]/(j+1).
```

This is nonpositive: it is `−1/3` at `j=1` and `−(j−2)/[(j+1)(j+2)]` at `j≥2`. At `j=0` the barrier is one and its drift is nonpositive. Bounded stopped Dynkin and monotone convergence bound the probability of **ever** getting the specified send by `f(j)`. This works with arbitrary changes to the exterior.

Sharpen by stopping the stock at `m`, with freezer payoff `1−f(m)`, or at the specified send, with payoff zero. For the surrogate process at stock `i`, use gain `g=i+immigration`, loss `b=i/max(2,i+B)`, and send `k=1/max(2,i+B)`. The discounted-free probability ratio of advancing from `i` to `i+1` before sending is `q_i=g/(g+b+k−b q_(i−1))`, with the downward term zero at `i=0`. The freezer value at `j<m` is `(1−f(m)) ∏_(i=j)^(m−1) q_i`. Each `q_i∈[0,1]`, so this value is increasing in stock. Against the true rates its drift is nonnegative until stopping: actual losses and sends are no faster, and gains are at least those used. The comparison process hits `m` or is killed a.s.; alternatively the physical stopped process leaves its finite stock interval or sends a.s. by the positive rates of a string of stock gains whenever gains are available, or by first hitting stock zero, where the zero-immigration payoff is zero. More directly, while the physical path remains inside a finite stock interval with positive immigration the chance of a run of gain events over a unit time has a uniform positive lower bound (gain bounded below and all monitored loss/send hazards bounded above). With zero immigration this also works before absorption at stock zero. Optional stopping gives the asserted lower bound on no send. Use the better of it and `1−f(j)`. A leaf therefore has ultimate hub probability at most one minus this bound. For a bare core, union-bound the next gain over the first-send probabilities from its present neighbours, all cores. **Use zero immigration for this latter calculation**: the receiving core does not have a leaf's deterministic return send.

## A lower bound for a core and a direct proof for 2/3

Track living leaf families at a focal core. A conservative tracked state `(d,l)` has physical degree at least `d` and `l` living tracked leaves. Their initiations create a new tracked leaf and raise `d` by one at rate `l`. Use independent `1/d` per-leaf removal proposals and thin to match the actual promotion hazard `1/D≤1/d`, discarding the track even on rejection. Ignore other degree gains. At `l=0` pay zero. At a chosen `d` cutoff, surviving families with individual rates `1,1/d` give probability at least `1−d^(−l)` of unbounded degree, hence positive mass by the existing [tagged-family argument](../../../../../../astra/source/round2/fringe/hub-children.md). The backward equation below the cutoff is

```
V(d,l) = [d V(d+1,l+1)+V(d,l−1)]/(d+1),   l>0.
```

The comparison process has no other transitions and reaches one boundary in finitely many jumps. Its descendant tags are only a subset of physical leaves. No autonomy of the actual core stock is assumed.

From the edge seed the center of `T3` is equally likely to be label 1 or 2, and label 3 is always a leaf. Let `H_c,H_l` be the center and leaf probabilities in that deterministic shape. Rerooting and exchangeability of the two leaves give `p1=p2=(H_c+H_l)/2`, `p3=H_l`. The center has two leaves and `H_c≥3/4` from two tracked families at death rate `1/2`. For a selected leaf, retain its initiating clock; the parent starts with one other leaf, and `B=1`. With cap `m=4` the no-send certificate is `1−H_l≥108/305`. Consequently `p2−p3≥127/2440>0` without any enumeration of later tree shapes.

## Signed finite-clock pullback

For a given labelled state and two tags, form an interval for their **signed difference** by subtracting the upper tail of the second from the lower tail of the first, and conversely for the upper endpoint. An exact tag-swap automorphism sets the interval to zero. Pull the interval back through actual rate-one initiation clocks, integrated over their holding times: the next target law is `r_w/n`, where `r_w=Σ_(x~w)1/D_x`. Recursing for `h` further births and intersecting with the direct all-time terminal interval is valid by the Markov property. The computation uses exact rational arithmetic throughout; printed decimals are rounded outwards. A count of further births is a finite stopping horizon, not a deterministic bound on elapsed time.

For `p3−p4` start from the deterministic edge and for `p4−p5` from the deterministic `T3` shape (the old labels do not affect the future rates or these two new tags). In the two-old-target sector, pair reverse orders before applying any interval; nonedges and repeats cancel exactly. An edge `{u,v}` contributes `[r_vδ_v−r_uδ_u] Δ(u,v)` with `δ=1/[D(D+1)]`; the serial target contributes `r_u/(D_u+1)` times the middle-minus-tip interval. A negative coefficient swaps interval endpoints. Divide by `n(n+1)`. This preserves the signed contributions and does not discard the adverse possibility that one of the tags becomes a hub only beyond the replay horizon.

```
timeout 25s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/early/check.py
```

The default uses four extra births per placement, a core cutoff of 12 and a first-send stock cutoff of 32. Certified intervals:

| Gap | Lower | Upper |
|---|---:|---:|
| `p2−p3` | .186507074525 | .392551907917 |
| `p3−p4` | .022645238637 | .134418576039 |
| `p4−p5` | .005316116377 | .087696910622 |
