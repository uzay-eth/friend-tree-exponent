# A uniform upper time block from one edge

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 1: Claude and Gemini; inputs include finite seam and uniform one-port comparisons. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


Working proof. The exponent is the continuous-time mean exponent of `proofs/existence.md`. All logarithms are natural. The interval starts just after an existing core promotes one of its leaves. The separate transition from the founding edge to the first core is handled in Section 8. The new core is called $`v`$, its parent is $`p`$, and its initial branch consists of $`v`$ and one leaf. Both across-edge endpoints are core forever. The number $`C_B(t)`$ counts cores on the $`v`$ side, including $`v`$. Put $`F(t)=\mathbb E_{T_2} C(t)`$, starting from the single edge.

## 1. Statement and parameters

Uniformly over the entire finite exterior present at the birth of $`v`$,

```math
 \mathbb E C_B(t)\le A(t)e^{\Lambda t} F(t+1).                                      \qquad\text{(1)}
```

The following constants are explicit. Choose an integer $`d\ge3`$ (the Hölder parameter), an even integer $`J\ge2`$, integers $`B,H\ge2`$, $`K\ge1`$, $`M>K`$, a number $`0<\kappa<1`$, and an integer $`b\ge2`$ such that $`s=1-1/b>\kappa`$. Define

```math
\begin{split}
 \gamma&=1/B,\quad \delta=(B-1)/H,\quad z=(1+1/(2\gamma))^d,\\ c_d&=d(d-1)(3/2)^{d-2}/8,\qquad V=d\gamma+c_d/K+(z-1)K/M,\\ D_0&=(1-\kappa)^{-b}/(s-\kappa),\\ \Lambda&=1/J+(1-\kappa)/2+V/d+2/d+\gamma+\delta .
\end{split}
```

Set $`n_*=3+(B-1)H`$,

```math
 \pi_*=(n_*-1)e^{-2}(1-e^{-1})^{n_*-2}\prod_{j=2}^{n_*-1}\frac1{j(j-1)},
```

and

```math
 A(t)=\pi_*^{-1}\left[1+D_0 e^{1/2}2^{1/d}z^{(M+K)/d}(2+t)^{J/2}(t+\gamma^{-1})\right]. \qquad\text{(2)}
```

For integer $`j\ge4`$, take $`d=B=j,\ J=2j,\ H=j^2,\ \kappa=1-1/j,\ b=2j,\ K=\max(1,\lceil jc_d\rceil),\ M=K+1+\lceil j(z-1)K\rceil`$. This gives $`\Lambda\le 6/j+2/j^2`$. At fixed $`j`$, $`\log A(t)=O_j(\log(2+t))`$. Section 8 converts (1) into an upper exponent bound for any fixed time block.

## 2. Joint density across the cut

Keep the edge $`pv`$ as a degree slot on each side. Each oriented core edge $`x\to y`$ adds one leaf at $`y`$ at rate $`1/D_x`$. Each leaf of $`x`$ gives a new leaf at $`x`$ at rate one. Budding at $`x`$ occurs at rate $`L_x/D_x`$. These disjoint transition channels may be driven by independent Poisson measures, thinned at their predictable rates. An outgoing core-to-core arrow changes only its receiver. The sum is exactly the friend-tree generator. Newly created cores can be named separately on each side by order of promotion.

For an interior trajectory $`x`$, write $`a_x(u)=1/D_v(x_{u-})\le1/2`$, $`A_x(u)=\int_0^u a_x(w)\,dw`$, and

```math
 r_x(u)=L_v(x_{u-})+\sum_{y\text{ core child of }v}1/D_y(x_{u-}).
```

For a prescribed, deterministic input $`a_x`$, let $`Q^{a_x}(dy)`$ be the law of the exterior alone. It has the original exterior transition rates, the retained slot at $`p`$, and independent Poisson immigration at $`p`$ of intensity $`a_x`$. Write $`b_y(u)=1/D_p(y_{u-})\le1/2`$. This function is nonincreasing. The arrows $`p\to v`$ have no effect on $`y`$.

Here is the joint path density, with no conditioning of the true exterior on an observed interior. Use counting measure on interior jump types and Lebesgue measure on their ordered times. Let $`I(x)`$ be the common factor consisting of the integrated survival for all intrinsic interior events, including the root's intrinsic gain rate $`r_x`$, times the jump rates for all events other than root leaf gains. If $`G(x)`$ is the list of root leaf gains, the joint path measure is

```math
 I(x)\,dx\ Q^{a_x}(dy)\ e^{-\int_0^t b_y(u)\,du}
           \prod_{u\in G(x)}[r_x(u)+b_y(u)].                                       \qquad\text{(3)}
```

Indeed, expand the usual time-ordered continuous-time jump density on both sides. The factors for exterior jumps and their survival rates, including input from $`v`$, give the normalized prescribed-input law $`Q^{a_x}`$. The remaining rates and survival factors give (3). One may first sum histories with bounded total jumps and then use nonnegative summation. Prescribed-input exteriors are nonexplosive since their total birth rate is bounded by their vertex count plus one. The feedback in the joint law survives in the dependence of $`Q^{a_x}`$ on $`x`$.

Use $`R_0`$ for the interior with zero incoming arrows and its parent slot retained. For deterministic $`0\le\beta\le1/2`$, $`R_{\beta,\tau}`$ is the same interior with immigration of rate $`\beta(u)`$ before $`\tau`$, one forced leaf gain at $`v`$ at $`\tau`$, and zero immigration afterwards. Superscripts on expectations below identify these laws. A forced gain adds no core.

## 3. An exterior estimate for one late arrow

For every deterministic $`0\le a\le1/2`$, every finite exterior and every $`u\ge0`$,

```math
 \mathbb E_{Q^a}[1/D_p(u-)]\le D_0 e^{-\kappa\int_0^u a(w)\,dw}.                    \qquad\text{(4)}
```

The first $`b`$ prescribed immigrants each give a new leaf at $`p`$. By the last of them $`D_p\ge b`$ permanently. Each subsequent immigrant starts a tracked family of leaves at $`p`$. Each tracked leaf rings at rate one and adds another leaf there. Its promotion rate is $`1/D_p\le1/b`$. Give each tracked leaf an independent death proposal of rate $`1/b`$, remove it from the tracked family on every proposal and thin the proposal for the true promotion with probability $`b/D_p`$. The tracked families are independent linear birth-death processes with rates $`1,1/b`$, contained in the actual leaves. Their births, death proposals and the prescribed immigration Poisson process use independent drivers. Untracked leaves and all other transition channels use their correct rates. This constructs the exterior marginal under $`Q^a`$.

A family survives forever with probability $`s=1-1/b`$. In operational time $`A(u)=\int_0^u a`$, the arrival $`W`$ of the first surviving family after the first $`b`$ immigrants is distributed as an independent sum $`\mathrm{Gamma}(b,1)+\mathrm{Exp}(s)`$, using rates. If necessary extend $`a`$ after the observation window. Its exponential moment is $`\mathbb E e^{\kappa W}=(1-\kappa)^{-b}s/(s-\kappa)`$.

For a family of age $`h`$, conditional on positivity, its size $`Z_h`$ is geometric on the positive integers with success parameter $`q=s e^{-sh}/(1-(1-s)e^{-sh})\le e^{-sh}`$. This follows by solving $`\partial_h f=(f-1)(f-(1-s)),\ f(0,z)=z`$ for the probability generating function. For a geometric variable, the expected reciprocal is $`-q\log(q)/(1-q)\le\sqrt q`$, including $`q=1`$ by continuity. The inequality follows from $`-\log q\le(1-q)/\sqrt q`$. Therefore $`\mathbb E[1/Z_h\mid\text{survival forever}]\le s^{-1}e^{-sh/2}`$, since survival implies positivity. Conditional on $`W\le A(u)`$, the selected family's law is its law given survival and its physical age is $`h\ge2(A(u)-W)`$. On $`W>A(u)`$ use $`1/D_p\le1`$. Both contributions together are at most $`s^{-1}e^{-\kappa A(u)}\mathbb E e^{\kappa W}`$, giving (4).

## 4. The last exterior arrival and finitely many profiles

In (3), mark a root gain as intrinsic or coming from $`p`$ by expanding its factor. If all are intrinsic, deleting $`e^{-\int b_y}\le1`$ and integrating $`Q^{a_x}`$ bounds the contribution by $`\mathbb E^{R_0} C_B(t)`$. Otherwise distinguish the last exterior gain at time $`\tau`$. Its factor is $`b_y(\tau)\,d\tau`$. Earlier gains have factors $`r_x+b_y`$, later ones have factors $`r_x`$. Drop the outside survival factor after $`\tau`$.

For each exterior path round $`b_y`$ up to multiples of $`1/J`$ at the left end of each unit interval, holding the result throughout that interval. Values lie in $`\{0,1/J,\ldots,1/2\}`$. The resulting deterministic candidate $`\beta`$ belongs to a list $`\mathcal P_t`$ of size at most $`(2+t)^{J/2}`$: each of the $`J/2`$ downward grid thresholds has a drop index or is never crossed. On the prefix, $`\beta\ge b_y`$ and $`\int_0^\tau(\beta-b_y)\le t/J+1/2`$, the last term coming from total variation at most $`1/2`$. Left limits at jump times give the same bound. Thus the prefix joint factors are at most $`e^{t/J+1/2}`$ times the corresponding factors for this profile. Sum over all profiles to remove their dependence on the exterior path. What remains inside the exterior integral is at most $`\int b_y(\tau)Q^{a_x}(dy)`$. Applying (4) gives

```math
 \mathbb E C_B(t)\le \mathbb E^{R_0} C_B(t)
 +D_0 e^{1/2+t/J}\sum_{\beta\in\mathcal P_t}\int_0^t
       \mathbb E^{R_{\beta,\tau}}[C_B(t)e^{-\kappa A_x(\tau)}]\,d\tau .             \qquad\text{(5)}
```

In this integration $`x`$ includes the forced gain at $`\tau`$. Its prefix survival and root jump factors are exactly those of $`R_{\beta,\tau}`$, and the intrinsic suffix starts at the post-gain state. Marking the last gain counts each marked history once.

## 5. A change of immigration rate at the root

For every such profile and every $`\tau\le t`$, with $`P=R_{\gamma,\tau}`$, let $`Z=dR_{\beta,\tau}/dP`$. Then

```math
 \mathbb E_P Z^d\le z^{M+K} e^{Vt}.                                                 \qquad\text{(6)}
```

Only the prefix root gain rate changes, from $`r+\gamma`$ to $`r+\beta`$. The likelihood is the product of $`u=(r+\beta)/(r+\gamma)`$ at these unmarked root gains times $`e^{-\int_0^\tau(\beta-\gamma)}`$. All other prefix rates, the forced gain and the suffix kernels agree at the same path. If $`L_v\ge K`$, factor the $`d`$-th powered likelihood on this predictable set as the likelihood for changing the rate to $`(r+\gamma)u^d`$, times the exponential integral of

```math
 (r+\gamma)(u^d-1)-d(\beta-\gamma)\le c_d/K .
```

Indeed $`r\ge K`$, $`|\beta-\gamma|\le1/2`$ and $`u\le3/2`$; Taylor's formula gives the displayed constant. This auxiliary probability law changes no rate when $`L_v<K`$, and preserves all promotion intensities. It is nonexplosive since its total rates are bounded by a constant depending on $`d`$ times the current vertex count plus one.

On $`L_v<K`$, the compensator term is at most $`d\gamma t`$, and each root gain contributes at most $`z`$ to the powered product. While $`D_v<M`$ there are at most $`M`$ root gains since each increases $`D_v`$. After $`D_v\ge M`$, each gain at $`L_v<K`$ increases $`\min(L_v,K)`$ by one. Decreases occur only on promotions with $`1\le L_v\le K`$, of rate at most $`K/M`$. Thus the low-leaf gain count is bounded by $`M+K+N`$, where $`N`$ counts these promotions after $`D_v\ge M`$ on the prefix. Under the auxiliary law its predictable intensity is at most $`K/M`$, giving $`\mathbb E z^N\le e^{(z-1)Kt/M}`$ by Poisson domination. This proves (6). No large-leaf gain count enters the estimate.

## 6. Hölder exponents

Let $`q=d/(d-1)\in(1,2)`$, $`\theta=(2-q)/q=1-2/d\in(0,1)`$, $`f=C_B(t)e^{-A_x(\tau)}`$, and $`h(\tau)=\mathbb E^{R_{\gamma,\tau}}f`$. Under $`P=R_{\gamma,\tau}`$, $`C_B`$ starts at one and its birth rate is at most $`C_B`$. The forced leaf gain at $`\tau`$ does not change this count. A one-ancestor Yule process thus gives $`\mathbb E_P f^2\le2 e^{2t}`$. Since $`A_x(\tau)\le t/2`$, Hölder with exponents $`d,q`$, then interpolation between moments one and two, gives

```math
\begin{split}
 \mathbb E^{R_{\beta,\tau}}[C_B(t)e^{-\kappa A_x(\tau)}]
 &\le e^{(1-\kappa)t/2}(\mathbb E_P Z^d)^{1/d}(\mathbb E_P f^q)^{1/q},\\ (\mathbb E_P f^q)^{1/q}
 &\le (\mathbb E_P f)^{(2-q)/q}(\mathbb E_P f^2)^{(q-1)/q}
 \le h(\tau)^\theta 2^{1/d}e^{2t/d}.
\end{split}                                                                          \qquad\text{(7)}
```

For the middle inequality write $`f^q=f^{2-q}(f^2)^{q-1}`$ and use exponents $`1/(2-q),1/(q-1)`$. This incurs only $`2t/d`$ in the logarithm.

## 7. A physical finite seed and the linear bound

Take $`T_*`$ with $`v`$ joined to its leaf and to $`p`$, $`p`$ joined to $`B-1`$ additional core vertices, and each of these vertices joined to $`H-1`$ leaves. Write $`S(t)=\mathbb E_{T_*}C(t)\ge1`$. Use an independent rate-$`\gamma`$ proposal process for $`p\to v`$, thinned to $`1/D_p`$. Then

```math
 \mathbb E^{R_0} C_B(t)\le e^{\gamma t} S(t),\qquad
 \int_0^t h(\tau)\,d\tau\le \gamma^{-1}e^{(\delta+\gamma)t}S(t).                     \qquad\text{(8)}
```

For the first inequality exclude every proposal up to $`t`$. For the second, require that none of the $`B-1`$ exterior neighbors of $`p`$ makes a proposal to $`p`$ on $`[0,t]`$. Each such edge has an independent rate-$`1/H`$ upper proposal, since its source degree is permanently at least $`H`$. This event has probability $`e^{-\delta t}`$. Integrate over the last $`p\to v`$ proposal at $`\tau\le t`$; its density is $`\gamma e^{-\gamma(t-\tau)}`$, with an independent unconditioned Poisson prefix of rate $`\gamma`$. Require no $`v\to p`$ arrow before $`\tau`$. Given the imposed interior path from that prefix, the independent reverse-channel Poisson measure meets this requirement with probability $`e^{-A_x(\tau)}`$. Up through the last proposal these events keep $`D_p=B`$: initially all its neighbors are cores, no arrival can add a leaf there, and its own outgoing arrows do not change its degree. Every prefix proposal and the last one are accepted. After $`\tau`$ there are no forward proposals, so the suffix needs no restrictions. The resulting branch is exactly the path in $`R_{\gamma,\tau}`$. Integrating its core count proves (8).

The standard process can realize $`T_*`$ at time one with probability at least $`\pi_*`$. Its Yule size at one has mass $`(n_*-1)e^{-2}(1-e^{-1})^{n_*-2}`$ at $`n_*`$. A fixed parent-before-child construction from one of its edges succeeds at every size $`j`$ with probability at least $`1/[j(j-1)]`$: choose an existing neighbor of the attachment target, then the target. The jump chain is independent of its holding times. By the Markov property $`S(t)\le\pi_*^{-1}F(t+1)`$.

Finally $`h^\theta\le1+h`$ and $`S(t)\ge1`$. Hence (8) gives

```math
 \int_0^t h(\tau)^\theta\,d\tau
 \le t+\gamma^{-1}e^{(\delta+\gamma)t}S(t)
 \le (t+\gamma^{-1})e^{(\delta+\gamma)t}S(t).                                        \qquad\text{(9)}
```

Substituting (6)–(9) into (5), using the profile count and $`S(t)\le\pi_*^{-1}F(t+1)`$, gives (1) with (2). In particular the resulting right-hand side is linear in $`F(t+1)`$.

## 8. Uniform upper and exponent bounds

For any initial finite tree $`T`$ with $`C(T)\ge1`$, each core born in $`(0,t]`$ lies in exactly one branch whose root was born directly from an initial core. Different such branches are disjoint, by core ancestry. Conditional at the birth of each root, the branch starts with $`v`$ and its one leaf, so (1) applies to its true exterior. The rate of direct promotions of each initial core is at most one. The tower property and compensation for those promotion times give

```math
 \frac{\mathbb E_T C(t)}{C(T)}\le 1+\int_0^t A(w)e^{\Lambda w}F(w+1)\,dw
 \le 1+tA(t)e^{\Lambda t}F(t+1)=:U(t).                                               \qquad\text{(10)}
```

Iteration at deterministic times yields $`\mu\le t^{-1}\log U(t)`$ for every $`t>0`$: start from $`T_3`$ and use the Markov property. The single-edge process delays $`T_3`$ by an independent exponential time of rate two; its nondecreasing mean has the same exponential rate, or just use $`F(u)\le\mathbb E_{T_3}C(u)`$ for this upper bound. For $`t\ge1`$, $`F(t+1)\le e(F(t)+1)`$: from any nonempty core the promotion rate is at most its size, and the still-single-edge state has probability $`e^{-2t}`$ and promotion intensity two. Integrating $`F'(u)\le F(u)+2e^{-2u}`$ proves the inequality. Also $`F(t)\ge1-e^{-2t}`$.

Consequently the upper endpoint has an explicit error above $`t^{-1}\log F(t)`$ bounded by $`\Lambda+O_j(\log(2+t)/t)`$, with (2) and (10) giving the constants without big-O. For lower endpoints one can use `proofs/existence.md`, Proposition 6.2 and Section 7, with their explicit parameter and threshold choices. Each $`F(t)`$ is computable by truncating the independent Yule size sum with the negative-binomial first-moment tail and enumerating the finite tree jump chain. Taking $`j`$ large and then $`t`$ large gives an effective two-sided approximation, without a local motif ansatz. These constants are intended as a computability bound.

## 9. An explicit error as a function of block time

Let $`t\ge1`$, $`f(t)=t^{-1}\log F(t)`$ and $`c=1-e^{-2}`$. Use the parameter family indexed by $`j`$ in Section 1 and denote (2) by $`A_j(t)`$. Equations (2) and (10) imply

```math
 \mu-f(t)\le E_+(t):=
 \min_{4\le j\le\max(4,\lfloor t\rfloor)}
 \left\{\frac6j+\frac2{j^2}
    +\frac1t\log\!\left[c^{-1}+e(1+c^{-1})t A_j(t)\right]\right\}.                 \qquad\text{(11)}
```

Indeed $`F(t+1)\le e(1+c^{-1})F(t)`$, and $`1\le F(t)/c`$. The integer range eventually contains every fixed $`j`$, so $`E_+(t)\to0`$.

The other direction follows quantitatively from Proposition 6.2 and Section 7 of the [existence proof](../../proofs/existence.md). To make the combined scale explicit, for $`j\ge4`$ set $`M_0=j`$, $`\epsilon_0=e^{-j}`$, $`k=\lceil 6j e^j\rceil`$, $`h=k+1`$, $`\eta_j=(7\log 2+9)/j`$, and let

```math
 p_j=\frac{e^{-2}(1-e^{-1})}{2}\,e^{-1}\Pr(\mathop{\mathrm{Pois}}\nolimits (1)\ge j)\,
        e^{-3}\frac{j}{j+2}
        \left(\frac{e^{-(k+1)}}{k}\right)^{k-2}.
```

At $`t\ge 2jk`$ the bad-history threshold in its Lemma 5.1 is met: with $`x=t/(jk)\ge2`$, $`1+\log(6x)\le 2x`$. Its additive bad-history error $`\sqrt{k(k+1)(2k+1)}e^{-t}`$ is at most $`1`$. The instantaneous promotion rate is at least one in every tree, including the initial edge. Therefore $`F(t)\ge t\ge 2(k+1)`$ and the bracket in Proposition 6.2 is at least $`F(t)/2`$. Iteration yields

```math
 \mu\ge\frac{\log F(t)-\eta_jt-\log(8k/p_j)}{t+h}.
```

Since $`F(t)\le2e^t`$, $`f(t)\le2`$ for $`t\ge1`$. This gives the computable lower-side error

```math
 f(t)-\mu\le E_-(t):=
 \min\!\left(\{2\}\cup
 \left\{\eta_j+\frac{2h+\log(8k/p_j)}t:
   4\le j\le\lfloor t\rfloor,\ 2jk\le t\right\}\right).                              \qquad\text{(12)}
```

Again $`E_-(t)\to0`$. Equations (11) and (12) give the interval $`[f(t)-E_-(t),\, f(t)+E_+(t)]`$. For a direct finite calculation one may use integer times $`t`$. At fixed $`t`$ write $`q=e^{-t}`$, $`a_n=\mathbb E C_n`$, and truncate $`F(t)=\sum_{n\ge2} (n-1)q^2(1-q)^{n-2}a_n`$ at $`n=N`$. The finite coefficients $`a_n`$ follow exactly by iterating the rational transition probabilities over labeled trees. The omitted part is bounded by the first moment of the negative-binomial size tail,
```math
 (1-q)^{N-1}\big[N(N+1)q+2(N+1)(1-q)+2(1-q)^2/q\big].
```
This tends to zero with a known bound at each fixed $`t`$. Strictly computable choices such as $`k=\lceil 6je^j\rceil`$ may equivalently use $`k=6j\,3^j`$ in (12) to avoid computing a transcendental ceiling; this also meets the same threshold.

One asymptotic simplification is
```math
 \left|t^{-1}\log F(t)-\mu\right|
      =O\!\left(\frac{\log\log t}{\log t}\right).
```
Indeed $`\log z\le j\log j`$, $`K\le 1+j^3(3/2)^j`$, $`M+K\le 5jzK`$ for $`j\ge4`$, and $`\log(1/\pi_*)=O(j^3\log j)`$. In (2) these give $`\log A_j(t)\le\exp\{j\log j+O(j)\}+O(j\log(2+t)+\log(t+j))`$ with absolute implicit constants. In (12), $`\Pr(\mathop{\mathrm{Pois}}\nolimits (1)\ge j)\ge e^{-1}/j!`$ gives $`2h+\log(8k/p_j)=\exp(O(j))`$. For $`j=\lfloor\log t/(4\log\log t)\rfloor`$ and sufficiently large $`t`$, both overhead terms divided by $`t`$ are $`o(1/j)`$ and $`2jk\le t`$. The explicit formulas (11)–(12) apply at smaller times as well.
