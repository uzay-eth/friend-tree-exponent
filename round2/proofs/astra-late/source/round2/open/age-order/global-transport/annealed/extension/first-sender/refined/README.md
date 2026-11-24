> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../../../README.md#claim-26) for exact scope and citations.

# Red-family conditioning closes the half-threshold gap

For **every old tree with older birth label `b≥20`**, the two next birth labels have strict tail order at every **`1/2≤x<1`**, conditional on the tree before either birth. Together with [Fable's paired receipts for `b=2,…,19`](../../../../../../../attacks/fable5/flow/README.md) and their [strictness argument](../../../../../../../attacks/fable5/README.md), this gives all-label order on that entire interval. Founders 1 and 2 are tied, and at `x=1` all tails are zero. The single `20→21` receipt is unnecessary. This closes all nine previously missing half-threshold pairs; it does not give those pairs all-threshold transport receipts.

## 1. Expose the first lost family, rather than charging its expected size

Use the physical serial middle at common size `b+1`, its one red leaf and blue mark, and the [first blue sender index `F`](../README.md). Given colour endpoint `θ`, put `u=1−θ`; `F` is independent of the entire colour word. Given `F=f`, immediately before that send there are `R=1+S` red vertices, all good leaves of the mark, where `S` is the sum of `f+1` geometric waits with survival `θ^k` (geometric counts start at zero). Conditioned also on the tree at that point, the probability of actual red loss is exactly `R/D`, where `D` is the mark's degree. On a red hit the lost leaf is uniform among the `R` leaves. For a uniform lower bound, on a miss choose a uniform dummy red leaf and charge its entire family too. In either case label the `R` extant red vertices as distinct family founders; red newborns inherit the family of their red sender, irrespective of the physical attachment.

Conditional on the exposed past and the chosen founder, this red family word is exactly a Pólya word: each red step picks a uniform red sender, and adds one member to that sender's family. The chosen family's limiting fraction `W` of all reds is `Beta(1,R−1)` if `R≥2`, and `W=1` if `R=1`. Its first two moments for all `R≥1` are `1/R` and `2/[R(R+1)]`; for the *initial actual* lost-family share these are multiplied by `R/D`. We use the full law of the conservatively charged family, which avoids needing an estimate of `D` or a second-moment inequality. Given `W=w<1`, each future red sender lies outside the charged family with probability `1−w`, independently across red steps; within that class the sender remains uniform. The future colour word and blue drivers retain their laws. These disintegrations follow by multiplying finite colour and family-word probabilities; no future attachment target is exposed to choose the founder.

## 2. Conditional damage to the surviving red families

For `w<1`, the density of `W=w` given `S=s≥1` is `s(1−w)^(s−1)`. Multiplying by the negative-binomial weight of `s` gives, **conditional on `w`**, that `S−1` is the sum of `f+2` geometric counts with survival parameter `z=θ(1−w)`: use `s binom(s+f,f)=(f+1)binom(s+f,f+1)`. The initial uncharged red population is `S`. Before blue choice `j>f`, each of the `j−f` new red strings adds an independent geometric number of uncharged reds with survival parameter `z/(1−θw)≥z`. Indeed deleting charged-red outcomes from the iid word leaves uncharged-red versus blue probabilities proportional to `z` and `u`. Thus the uncharged red population before choice `j` dominates `1+` the sum of `j+2` geometrics of parameter `z`. Its mean reciprocal is at most `(1−z)/[z(j+1)]`, by integrating the generating function on `[0,1]`.

At the forced send all uncharged reds remain good. At subsequent blue choices only the mark can spoil one such leaf, with sender probability `1/(b+j)`; at an uncharged-red step the bad fraction cannot increase in conditional expectation (good leaves reproduce good leaves). Multiply the reciprocal loss estimate by the uncharged asymptotic share `z`, and use bounded convergence for that share and Fatou for the deficit, as in the parent proof. This gives

```
E[(z−Z)_+ | θ,F=f,W=w] ≤ (u+θw) Σ_(j>f) 1/[(b+j)(j+1)]
                         ≤ (u+θw)/b · log(1+b/(f+1/2)).
```

The last estimate replaces `j+1` by `j` and uses the midpoint integral of the convex function `1/[v(b+v)]` over each `[j−1/2,j+1/2]`. The blue hazard is unchanged even conditional on the red family schedule, because neither the colour word nor family assignment uses subsequent blue sender uniforms.

## 3. Integrate the exact family-size distribution

For `θ>x` set `t=b(θ−x)/u`, `a=bθW/u`, and `C_b(f)=log(1+b/(f+1/2))`. The success floor conditional on `θ,f,W` is zero for `a≥t`, and otherwise is

```
G(t,a,C) = (1 − C(1+a/b)/(t−a))_+ .
```

This is Markov **after** exposing the initial family's full random share, and applies only to subsequent losses. Given `θ,f`, for `0≤a<bθ/u`, the survival of the scaled family size is exactly `(1+a/b)^(-(f+1))`: average `(1−w)^S` with the geometric PGF. There is an atom at the upper cap `bθ/u` when `R=1`, but that cap exceeds `t`. Thus in computing the payoff we can replace `a` by the uncapped variable with this survival, independent of `t` conditional on `f`. Integrating over `θ>x` contributes `(1−x)^b`, and conditional on that event `t` has the parent proof's survival `(1+t/b)^(-b)`.

For a uniform bound over `b≥B`, use the same quantile `q∼Uniform(0,1)` for `F=floor((b−1)(1−q)/q)`, and set `g=(B−1)(1−q)/q−1`. Drop `g≤0`. Since `(g_b+1/2)/b` and `g_b/b` increase with `b` (where `g_b=(b−1)(1−q)/q−1`), use `C≤log(1+B/(g+1/2))` and `1+a/b≤1+a/B`. The uncapped `a` has CDF bounded below by

```
1 − exp(−(g/B) a/(1+a/B)),
```

using `log(1+y)≥y/(1+y)` and `f+1≥g_b`. This lower CDF bounds the payoff from below because it decreases in `a`. Replace `t` by an independent mean-one exponential, which decreases the increasing payoff. This leaves just uniform/exponential rectangles and a one-dimensional family CDF. In each `q` bin use its right endpoint (both changes lower the payoff); in each `t` bin use its left endpoint. Round the lower CDF down at the family grid endpoints, assigning its increments to the right endpoints and the residual to infinity. The exponential upper lookup for the CDF uses `e^-h≤1−h+h²/2` at mesh `1/64`, upward fixed-point products and a downward index; all payoff and `t` bin masses are rounded down. The checker certifies a serial lower coefficient `>41/200` at `B=20`, multiplying `(1−x)^b`.

Use the parent proof's unchanged paired beta debit with `θ_*≤(3/4)((1−x)/(bx))^(1/6)`. Its rational corner test at `B=20,x=1/2` uses the outward sixth root `61/100`; monotonicity of the cost extends the result to all larger `b,x<1`. The additional `B=21` replay uses a coarser grid. Neither computation enumerates physical trees.

```
nice -n 10 python3 -IB round2/open/age-order/global-transport/annealed/extension/first-sender/refined/check.py
PASS exact red-family conditional weights and cap atom
PASS b>=20, x>=1/2: refined coefficient > 41/200
PASS b>=21, x>=1/2: refined coefficient > 507/2500
```
