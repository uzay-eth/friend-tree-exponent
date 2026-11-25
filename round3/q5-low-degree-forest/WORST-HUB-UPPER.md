> **Public checkpoint export note:** Edited export of `forest-astra/WORST-HUB-UPPER.md` (frozen original SHA-256 `9b2b7b0a02f37757857ce80fe7d565685ef0653ed9dc90850fc3b8dff9c82cd4`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

RESULT: proof candidate — at fixed supercritical polynomial exponent, the worst single-parent log-volume KS distance converges almost surely to Neon's two-bin value $`D^*(B)`$. Old children that can still contribute admit a stronger exponential census than raw old arity.

# Sharp upper bound for the worst counted-volume parent

Astra lead. Let $`N(n)=n^{\alpha+o(1)}`$ be deterministic integer cutoffs, with $`\mu/(1+\mu)<\alpha<1`$ fixed, and set $`B=\mu(1-\alpha)/\alpha\in(0,1)`$. In the physical size-$`n`$ tree delete vertices of original degree greater than $`N`$. In each deleted parent's row count every component with an original surviving core whose top has that parent, including singleton cores. Omit the founder component. Write $`F_{u,n}(x)`$ for the row fraction with ordinary volume $`V\le N^x`$, with the zero function for an empty row. Then

```math
\max_{u\in H(\tau_n)}\ \sup_{0\le x\le1}|F_{u,n}(x)-x|
\ \longrightarrow\ D^*(B)\quad\text{almost surely}, \qquad (1)
```

where

```math
D^*(B)=\sup\{|p-x|:\ 0<x,p<1,\ 1-e^{-D_{\rm KL}({\rm Ber}(p)\|{\rm Ber}(x))}< B\}.
```

This is a statement about empirical rows of the observed tree for each fixed deterministic schedule. The new assertion is the upper bound. The matching lower bound is Neon's D1: two recent/early Poisson counts realize a tilt of cost strictly below $`B`$ at many predictable low-arity candidates. Its threshold gaps and private mark tests transfer those counts to actual rows. D2 concerns $`\log N/\log n\to1`$ and supplies geometric-block bookkeeping; (1) uses no endpoint vanishing theorem. The two-parent extremal constant is outside (1).

**Sources.** I verified Neon's review before opening the frozen bytes. Its SHA-256 is 8d876c8be0842903e4288143a480ecda3198fd3f3c46b6d0a4211abf96794ea8 (28,390 bytes). The reviewed Astra §5 input matches c36055fa6649ae5d779c72215ba8e8e1c6765ab36e02bc7f882ab571d7c5fcb8. The public baseline is c250a1f2fdd53ced91b7b485e94ac8a16106b9c7. The model inputs used below are the Yule clock, $`\mathbb E C(T_v)\le C e^{\mu v}`$, the a.s. maximal core-arity bound $`O(t)`$, and the tagged-degree and Poisson–inert branch couplings in Claim 5a, §§3.1 and 3.3. `REVIEW-ASTRA-WORST-PARENT-KS.md` (SHA-256 `8d876c8b...`, matching lower-bound input, not exported in this checkpoint) · [Baseline couplings](../../round2/proofs/neon2/q5/README.md) · `REVIEW-NEON-HUB-DEPENDENCE.md` (SHA-256 `c36055fa...`, not exported in this checkpoint)

## 1. One old forest for a whole block

For $`e^m\le n<e^{m+1}`$ write $`L=\alpha m`$, $`L_-=\min\log N(n)`$ and $`L_+=\max\log N(n)`$ over that block. Both are $`L+o(L)`$. Localize with a deterministic integer $`C`$ so that eventually every $`\tau_n`$ in the block is in $`[t_-,t_+]=[m-C,m+1+C]`$; this follows from the a.s. positive Yule size limit. Use the deterministic old time

```math
s=t_+-L_-+4\log m,\quad h_\pm=t_\pm-s,\quad
D_f=2e^{L_+-h_-},\quad w=\log D_f+\sqrt m .
```

Thus $`h_\pm/L\to1`$, $`D_f\ge m^4`$, $`w=o(L)`$, $`s=(1-\alpha)m+o(m)`$, and $`s>w`$ eventually. The margin $`w`$ handles subpower variation of $`N`$ within a block. Define the old candidate set $`\mathcal H_m`$ of cores with $`D(s)\ge m^4/2`$.

Eventually every deleted parent in the block belongs to $`\mathcal H_m`$, and every old vertex with $`D(s)>D_f`$ is deleted at **all** observations in the block. For the first assertion use upper Yule degree domination to $`t_+`$, since $`e^{L_--h_+}=m^4`$; include vertices born after $`s`$ using the expected birth compensator up to $`t_+`$. For the second use the tagged-stock lower comparison at $`t_-`$ with $`e^{L_+-h_-}=D_f/2`$. In detail, on the old-measurable gate $`\max c(s)\le K_* (m+O_C(1))`$, an old vertex with $`D>D_f`$ has $`D-c=(1-o(1))D`$ tagged leaves. Their normalized $`h_-`$-time Yule sum has exponentially small lower tail at $`3D/4`$. The removed families are bounded by the normalized rate-one immigrating Yule comparison of Claim 5a §3.1, whose tail above $`D/8`$ is $`O(e^{-cD})`$. The resulting degree exceeds $`e^{L_+}`$. All marginal failures union over $`O(e^{t_+})`$ vertices in expectation against $`e^{-c m^4}`$ (or the weaker summable bounds of §3.1). Exhaust $`K_*,C`$. No future deletion set is used to condition a Poisson count.

For each $`u`$ already a core at $`s`$, define the **observable**

```math
R_u=\#\{y:y\text{ is a core neighbor of }u\text{ at }s,\quad D_y(s)\le D_f\}.
```

For every actual deleted $`u`$, the number of row entries with a top already a core at $`s`$ is at most $`R_u`$. A deleted old child contributes no component with that child as the surviving top. Any component hanging off that deleted child has its own immediate parent. Including the core neighbor toward the founder in $`R_u`$ is harmless.

## 2. Exponential census of potentially contributing old children

For $`\theta=3/4`$, with $`A_m=\sum_{u\in K(T_s)}e^{\theta R_u}`$, we claim

```math
\mathbb E A_m\le C_\theta\exp\{\mu s+(e^\theta-1)w\},\qquad
\limsup_{m\to\infty} L^{-1}\log A_m\le B\quad\text{a.s.} \qquad (2)
```

Constants may absorb a bounded initial time or seed. The sum is over **all** old cores; no bias from conditioning on becoming a hub is involved.

At a promotion, the parent's new core child starts with one leaf. Use a lower stock chain $`J`$ started at one, absorbed at zero, with upward rate $`j`$ and downward rate $`j/(j+1)`$. At equal positive stock the true child receives stock at rate at least its stock, and loses it at rate $`j/D\le j/(j+1)`$ since its parent is a core. An order coupling uses the child's leaf-source birth measure and its own marked rate-one initiating clock (put leaf choices first). Other incoming arrows can add to the true stock. Starting above the lower chain cannot reverse the order at a single unit jump. For **one specified parent** its proposal times have rate one, and one may preassign to every proposal independent private drivers for a prospective child's lower chain, even if no promotion occurs. Use these fresh drivers for a newly promoted leaf from that stopping time. Its own clock and the aggregate leaf-source births at its site can be realized by independent Poisson drivers. Distinct such roots are disjoint, and the parent's proposal clock is independent of these private drivers. Incoming physical arrows can depend on earlier child histories. The order comparison holds pathwise despite that dependence. This construction is needed only marginally for each parent when taking the expected sum in (2).

Here is a bound with enough numerical slack. From one, the probability to reach $`k`$ before zero is $`1/\sum_{j=1}^k1/j!`$, by the finite birth–death difference equation. The nonnegative supermartingale $`e^{-v}J_v`$ converges. Started at $`k`$, its limit has mean at least $`k-1`$ and second moment at most $`k(k+1)`$: integrate the scaled death drift, bounded in absolute value by $`e^{-v}`$, and use Yule$`(k)`$ domination for uniform integrability. Conditional on reaching $`k`$, the probability of positive limit is at least $`(k-1)^2/[k(k+1)]`$. Let $`k\to\infty`$. The positive-limit probability from one is at least $`1/(e-1)>1/2`$. Consequently

```math
q_m=\Pr_1(J_w\le D_f)\quad\text{satisfies}\quad
\limsup q_m\le1-1/(e-1)<1/2 ,
```

because $`w\to\infty`$ and $`D_f e^{-w}=e^{-\sqrt m}`$. If an accepted child proposed by $`s-w`$ has $`J_w>D_f`$, its **degree** is greater than $`D_f`$ already at $`s`$ by monotonicity, even if its stock fell later. This is the reason to test the old degree before freezing the fresh kernel.

Fix a parent at its own core promotion time $`v\le s`$ (handle the initial core separately). It begins with at most one core neighbor. Every subsequent new core neighbor requires one of its own proposals. Up to $`s-w`$, count each proposal whose private chain fails the test, including dummy proposals; in the last $`w`$ units count all proposals. This dominates $`R_u-1`$. Conditional on the post-promotion filtration, independent marking of the unthinned proposal clock gives

```math
\mathbb E[e^{\theta R_u}\mid\mathcal F_v]\le
e^\theta\exp\left\{(e^\theta-1)
[q_m(s-w-v)_+ + \min(w,s-v)]\right\}.
```

This does not claim that outcomes of physically accepted children are independent. The future tests for the early marks finish by $`s`$. For large $`m`$, $`q_m\le1/2`$, and $`\tfrac12(e^{3/4}-1)<9/16<\mu`$ since $`e^{3/4}<17/8`$. Integrating the bound against promotion times with expected intensity at most $`Ce^{\mu v}`$ proves the first part of (2). Explicitly the integrand after extracting $`e^{\mu s+(e^\theta-1)w}`$ decays as $`e^{-[\mu-q_m(e^\theta-1)](s-v)}`$. Conditional bounds at promotion times followed by compensation justify this integration, without requiring the future statistic to be predictable. Markov and Borel–Cantelli over integer $`m`$ give the second part since $`\mu s/L\to B`$, $`w/L\to0`$.

One can expose a still earlier forest, before the old counts have been realized. For deterministic $`0\le v_0\le s`$, put $`a_\theta=e^\theta-1`$, $`\lambda_m=q_ma_\theta`$, and $`C_0=C(T_{v_0})`$. The same argument conditional on $`\mathcal F_{v_0}`$ yields

```math
\mathbb E[A_m\mid\mathcal F_{v_0}]
\le C_\theta e^{a_\theta w}
\left[e^{\lambda_m(s-v_0)}
\sum_{u\in K(T_{v_0})}e^{\theta c_u(v_0)}
+ C_0 e^{\mu(s-v_0)}\right]. \qquad (2a)
```

For the first term charge all neighbors already core at $`v_0`$ and thin only later proposals. For the second use the uniform all-seed core mean, conditional at $`v_0`$, to bound the expected subsequent promotion intensity by $`C C_0e^{\mu(v-v_0)}`$. The integral converges with a uniform factor because $`\lambda_m<\mu`$ with slack. The original process starts at $`T_3`$ in this time coordinate. Taking conditional expectations of (3) or (4) below and substituting (2a) gives the joint tail while both the intervening old population and the fresh kernel are still random. This stronger conditioning does not require simultaneously constructing independent surrogate marks for every existing parent.

For comparison, the identical argument counting all old neighbors gives (2) with $`c_u(s)`$ replacing $`R_u`$ only for the sufficient range $`\theta<\log(1+\mu)`$. That raw exponential estimate also admits the conditional joint bound below with its smaller $`\theta`$. The improvement in (2) is applied specifically to physical row contributors.

## 3. Conditional joint tail and transfer to the observed parent row

Given $`\mathcal F_s`$, run the proposal/inert coupling for the old candidates. Use only proposals on $`(s,t_-]`$ in the two-bin ideal count $`(X_u,Y_u)`$; $`X_u`$ counts the favorable marks, $`Y_u`$ the others. For an upward deviation at a fixed threshold $`x`$ use the favorable event $`V_{\rm ideal}(t_--\sigma)\le e^{xL_+}`$. For a downward deviation use $`V_{\rm ideal}(t_+-\sigma)> e^{xL_-}`$ instead. These monotone endpoint tests overcount the respective favorable physical marks everywhere in the block. Given the old forest each row's two bins are independent Poissons with means $`L(z+o(1)),L(1-z+o(1))`$, where $`z=x`$ and $`z=1-x`$, respectively. The errors are deterministic for each fixed grid and block. Indeed at ages below $`(x-\epsilon)L`$ Markov and $`\mathbb E V_a\le2e^a`$ make the probability of overshooting exponentially small, while above $`(x+\epsilon)L`$ the lower $`1+\mathrm{Yule}(1)`$ comparison bounds undershooting by $`e^{-c_\epsilon L}`$. The intervening width is at most $`2\epsilon L+o(L)`$. Use $`h_\pm/L\to1`$ and a fixed endpoint shift, including the clipped age intervals at $`x=0,1`$.

For every $`a\ge0`$, $`b\in\mathbb R`$, $`t\ge0`$ fixed and $`p\in[0,1]`$, we obtain the joint upper tail

```math
\Pr\{\exists u\in\mathcal H_m:R_u\ge a L,\ (1-p)X_u-pY_u\ge bL\mid\mathcal F_s\}
\le A_m\exp\{L[-\theta a-tb+z(e^{(1-p)t}-1)+(1-z)(e^{-pt}-1)+o(1)]\}. \qquad (3)
```

More usefully, for $`t(1-p)\le\theta`$ the old cost can be kept in the **same row**:

```math
\Pr\{\exists u\in\mathcal H_m:(1-p)(X_u+R_u)-pY_u\ge-\epsilon L
\mid\mathcal F_s\}
\le A_m\exp\{L[t\epsilon-1+z e^{(1-p)t}+(1-z)e^{-pt}+o(1)]\}. \qquad (4)
```

Both follow by exponential Markov and summing over the predictable set; in (4) use $`\sum_u e^{t(1-p)R_u}\le A_m`$. At time $`s`$ the old counts are known. The old randomness in (3)–(4) is expressed by the almost-sure census (2). No conditioning on the selected size epoch, or on survival of a mark, takes place.

For the physical comparison, stop a hub's bad sends (to existing cores) at the common arity exit $`c_{\max}>2K_*(m+1+C)`$. Before that exit their predictable rate at a candidate is $`O(m^{-3})`$ since $`D_u\ge m^4/2`$. Over the whole horizon the count $`S_u`$ is conditionally dominated by $`\mathrm{Pois}(O(m^{-2}))`$. The late proposals $`B_u`$ on $`(t_-,t_+]`$ have constant mean $`1+2C`$. On $`\{A_m\le e^{(B+\delta)L}\}`$, a union bound over all candidates costs at most that factor: conditional Poisson tails above any fixed $`\epsilon L`$ decay as $`\exp[-c_\epsilon L\log L]`$. Borel–Cantelli and the eventual localization imply $`\max_u(2S_u+B_u+1)=o(L)`$ for these charges. Comparisons apply to actual hubs at the observation where the degree sandwiches hold.

For a deleted $`u`$, every accepted post-$`s`$ branch is free of deleted descendants at observation: all deleted vertices were already old cores at $`s`$. A coupled early branch with no separation contributes one entry with the same parent, with favorable indicator bounded by its ideal endpoint test. A failed proposal contributes at most one phantom; a bad send can also first separate one previously accepted branch. For the score favorable count minus $`p`$ times total count, omitting a phantom costs at most $`p`$, and changing a matched binary indicator at most one. Late births cost at most $`B_u`$, and the initial founder direction can cost one exceptional branch. Pre-existing tops contribute at most $`(1-p)R_u`$. Thus a physical favorable fraction at least $`p`$ implies the event in (4) up to an $`o(L)`$ tolerance. For a nonempty row the downward event uses the complementary fraction. For an empty row, the downward event under the zero-CDF convention still implies the necessary **unnormalized** score inequality $`0\ge0`$. This avoids dividing by a possibly small denominator.

## 4. Why the old-count restriction does not change the KS optimum

For a deviation $`d\in(0,1)`$ take $`0\le z\le1-d`$, $`p=z+d`$, and define

```math
I_\theta(z,p)=1-\inf_{t\ge0:\,(1-p)t\le\theta}
[z e^{(1-p)t}+(1-z)e^{-pt}],\qquad
J(d)=\min_{0\le z\le1-d}D_{\rm KL}({\rm Ber}(z+d)\|{\rm Ber}(z)).
```

Then for $`\theta=3/4`$,

```math
\inf_{0\le z\le1-d} I_\theta(z,z+d)=1-e^{-J(d)}. \qquad (5)
```

Without the restriction on $`t`$, the bracket's minimum is $`e^{-D_{\rm KL}}`$, at $`t_*=\log[p(1-z)/(z(1-p))]`$ in the interior. Its maximum over $`z`$ of these minima is $`H_*=e^{-J(d)}`$. Clipping can only increase a minimum of the bracket. It suffices to show it never increases one past $`H_*`$. Where the restriction binds, put $`y=1-p`$. The clipped bracket is

```math
H_\theta=e^\theta[1-d-y+(d+y)e^{-\theta/y}].
```

Its second derivative in $`y`$ is $`e^\theta e^{-\theta/y}\theta[(d+y)\theta-2dy]/y^4`$, nonnegative since $`d+y\le1`$ and $`2dy/(d+y)\le1/2`$. On each interval of clipping its maximum is at the boundary. At an interior switch it agrees with the unrestricted value. At $`p=1`$ the restriction disappears and the infimum is $`z`$. At $`z=0`$ the clipped value is $`e^{-\theta d/(1-d)}`$. Using the unrestricted trial point $`z=(1-d)/2`$, put $`v=(1+d)/(1-d)`$ to get

```math
\frac{1-d}{d}J(d)\le (1-d)\log v=\frac{2\log v}{1+v}
\le 2/e < 3/4.
```

The $`z=0`$ value is therefore also below $`H_*`$. This proves (5), with the endpoint infima interpreted by limits. The fixed-threshold rate can genuinely be changed by clipping. Equation (5) is about optimizing over the KS threshold.

The raw old-arity bound alone would leave a real gap in this argument. At $`B=3/8`$, $`d=1/2`$, $`z=0`$, using $`\theta=\log(1+\mu)`$ gives rate $`1-e^{-\theta}=\mu/(1+\mu)<3/8`$. Yet $`J(1/2)\ge 1/2`$ by $`D_{\rm KL}({\rm Ber}(z+d)\|{\rm Ber}(z))\ge2d^2`$ (integrate the second derivative in $`p`$), and $`1-e^{-1/2}>3/8`$. This is a failure of that bound, not a physical counterexample.

For every $`d>D^*(B)`$ with $`d<1`$ the right side of (5) strictly exceeds $`B`$. Indeed reflection covers both deviation signs. If the minimum were at most $`-\log(1-B)`$, a strictly smaller positive deviation still exceeding $`D^*(B)`$ at the same $`z`$ would have strictly smaller divergence. Here $`z=0`$ cannot minimize, and moving $`p=1`$ inward handles the other endpoint. Compactness gives a minimizer. Also $`D^*(B)<1`$ for $`B<1`$, since deviation tending to one forces the Bernoulli divergence to infinity.

Fix such a rational $`d`$ and a finite rational CDF grid. For each direction and grid point with $`z+d\le1`$, choose a **finite** tilt in (4) for which the cost exceeds $`B`$ by a positive margin. This is possible even when the infimum is at infinity. Take the smaller margin over the finite choices, then choose the old-census and transfer tolerances small enough. Intersect the conditional bound with the old-measurable census gate first. On it failure probabilities are $`O(e^{-cL})`$, summable in integer blocks. Then use the eventual physical sandwiches and remove localization. Grid monotonicity adds at most one grid spacing, uniformly over rows and sizes in a block. Countably refining the grid and sending $`d\downarrow D^*(B)`$ proves the upper bound of (1). Together with frozen D1 this proves the claimed limit. Even for $`\alpha`$ very close to one the summability is on blocks, not on every integer $`n`$. D2 has $`s=o(L)`$ and can bound every raw pre-existing arity by $`o(L)`$ there; (2) replaces precisely that step at a fixed exponent below one.

## Verification and handoff

I authored [a finite checker](worst_hub_upper_check.py) for the birth–death hitting recurrence, exact exponential slack, the raw-arity negative control, and a numerical stress test of the clipped variational calculation. It uses no Neon replay. The conditional tail, curvature argument and almost-sure transfer above are the asymptotic proof.

~~~sh
python3 -B round3/q5-low-degree-forest/worst_hub_upper_check.py
~~~

Exit 0. Actual stdout

~~~text
exact stock hitting-equation checks: 780
exact concentration slack and raw-count negative control: passed
numerical envelope grid: 202851 cases, 43097 clipped, max minimum-cost difference 2.22e-16
numerical D* on the finite grid at B=0.1: 0.22817
numerical D* on the finite grid at B=0.375: 0.47165
numerical D* on the finite grid at B=0.9: 0.90004
~~~

**Review status (historical pre-review line).** Neon independently audited the earlier lower bound and authored D1/D2. The eligible-child thinned census, conditional tilt and sharp upper proof in this addendum were written as Astra's new work pending independent review (*subsequently cross-family audited `HOLDS` in [`REVIEW-WORST-HUB-UPPER.md`](REVIEW-WORST-HUB-UPPER.md); still awaiting human refereeing*). Sampling from fixed observed rows gives no process-level independence statement.
