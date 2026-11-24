> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../../README.md#claim-18) for exact scope and citations.

# Ultimate hub probability of a deterministic late birth

Let `v` be the **birth label**, `Z_v` its limiting fraction of edges, and `S_v` the time from birth to first promotion, using the rate-one continuous embedding (`S_v=∞` if it never promotes). Using the repository's [core mean and L² limit](../../../../../../astra/source/round2/fringe/almost-sure.md) and [uniform-core local limit](../../../../../../astra/source/round2/fringe/spatial-limit.md), the pointwise asymptotic is

```math
 \Pr(Z_v>0)\sim c\,v^{\mu-1},\qquad
 c=\mathbb E[X] \,\pi(F),\quad X=\lim C_n/n^\mu,
 \qquad \mu/2\le \pi(F)\le7.
```

Here `π` is the stationary shape law in the sense `π Q_t=e^{μt}π`, and `F` is the birth score below. This gives a positive, finite specification of the leading constant in terms of a local limit and the ordinary physical dynamics. No independence of sibling futures or of local mass and the Yule size limit is needed. No numerical value for `c` is claimed.

## 1. The pointwise birth flux

Write `K_n` for the cores at size `n≥3`. Given the whole prebirth tree, each of its `n` vertices is equally likely to initiate the next birth and chooses a uniform neighbour. For a core `o` define the score `F(G,o)` by adding:

* **Core initiation:** choose a uniform neighbour (core or leaf) of `o`, do the physical birth there, and measure the probability that the newly added degree-one vertex is ultimately a hub. A leaf target is promoted **by this birth**, but it is the *newly added vertex*, not this target, whose future is counted.
* **Leaf initiations:** multiply by `L_o` the ultimate-hub probability for a marked new leaf attached to `o`. Each of its `L_o` leaves chooses `o` deterministically.

The probability after placement is the ordinary Markov law from the modified marked graph. Use `F_T` if the marked leaf must also promote by age `T`. Use `B_T` instead when counting promotion by age `T` regardless of ultimate hub status, and `B=B_∞`. Leaf identities are interchangeable until chosen, so these are well-defined scores on rooted marked core trees. This is the exact conditional formula, for any one of these endpoint/delay tests `E`:

```math
 n\Pr(E_{n+1}\mid T_n)=\sum_{o\in K_n} F_E(G_n,o).
```

In particular the number of ultimate hubs among labels at most `n` has mean `Σ_(v≤n) p_v`; its increment is the expected score above divided by `n`. In contrast `C_n` counts promotions **by size** `n`, including leaves born at earlier labels, and `#\{Z_v\ge x\}` is a mass-threshold count over all labels. Neither can simply be differenced to get `p_v`.

## 2. The tail estimate: the tagged leaf must stay in the parent process

After placement the tagged leaf has a core parent `x` of degree at least two. Suppress only `x→tag`, retaining the tag itself and its rate-one initiations (which add other leaves at `x`). In this suppressed process let `D(t)` be the parent degree. The missing channel can be restored with an independent exponential threshold at accumulated hazard `A(t)=∫_0^t ds/D(s)`, until its first ring. Conditional on the suppressed path the density of a **finite** promotion delay is `e^(−A(t))/D(t)` and there is a separate atom `e^(−A(∞))` at infinity. This representation is for the first promotion only; its continuation must use the new postpromotion state and fresh future rates.

Let `k` be the other leaf stock at `x` just after placement. Track a subset of those leaves and all new leaves born from the tag as immigrants. Each tracked leaf makes a new tracked leaf at rate one. Use independent rate-`1/2` proposals on each tracked leaf to cover its physical removal, whose rate is `1/D≤1/2`; discard the track on every proposal, even a rejected one. The suppressed degree dominates `1+J_t`, where `J` is linear birth–death with immigration, starting at `k`, with per-member birth/death rates `1,1/2` and immigration rate one. This construction splits the directed channels and keeps their physical rates. In particular the tracked initial families and the immigrant families in `J` can use independent drivers. Write `I_t` for the immigrant part.

The following uniform bound is useful both at small and large parent degrees:

```math
 (k+1)\,\mathbb E[(1+J_t)^{-1}]\le K(1+t)e^{-t/2}.
```

Indeed a family `Y_t` from one has `E Y_t=e^{t/2}` and `E Y_t²=4e^t−3e^{t/2}` (apply the generator to `j,j²`). Paley–Zygmund gives probability at least `1/16` of size at least `e^{t/2}/2`. If `M` of the `k` initial families are this large, then on `M>0` use `2e^(−t/2)/M`, with `1/M≤2/(M+1)` and `E[1/(M+1)]≤16/(k+1)` by the binomial integral. On `M=0` use the independent `I_t`; `(k+1)Pr(M=0)` is bounded. For clarity, the single-family pgf `g_t(z)` solves `∂_t(1−g)=(1−g)/2−(1−g)^2`, with `g_0=z`. Poisson immigration gives `E z^I_t = exp(−∫_0^t(1−g_s(z))ds)=1/[1+2(e^{t/2}−1)(1−z)]`, so integration over `z∈[0,1]` gives `E(1+I_t)^−1=log(1+b)/b`, `b=2(e^{t/2}−1)` (value one at zero). This proves the estimate.

In the core-initiation score each postbirth parent has some `k≥0`; in the leaf-initiation score `k=L_o`, exactly the potentially enormous multiplier. Thus in **either** score, both the finite-delay density and its hub-restricted part are bounded, after summing the two initiation sectors per core, by `K'(1+t)e^(−t/2)`. This also bounds the density per core for the remaining promotion delays of all leaves present at an observation time: for the `L` exchangeable existing leaves use `k=L−1`. Consequently

```math
 0\le F-F_T\le B-B_T\le K''(1+T)e^{-T/2}
```

uniformly in the rooted environment, excluding the atom at infinity. The atom is positive for every fixed finite-degree placement, even under the infinite-core law: `E A(∞)<∞` by the bound, so `A(∞)<∞` a.s. The bare fact that the typical newborn goes to a large hub is consistent with the result; the many births there are offset by small probabilities per birth.

There is also a shorter quantitative bound for the total finite mass. For `f(j)=1/(j+2)` and the immigration–birth–death generator `L_J`,

```math
 -(j+1)L_J f(j)=\frac{j^2+j+2}{2(j+2)(j+3)}\ge\frac16
 \quad(j=0,1,2,\dots),
```

where the inequality reduces to `j(j−1)≥0`. Dynkin's formula and nonnegativity give `E A(∞)≤6/(k+2)`. Consequently `B≤1+6 L_o/(L_o+2)≤7`, even including the core-initiation sector (at most one). More strongly, for each finite `n≥3`, the resulting pointwise estimate is `p_(n+1)≤7 E C_n/n`.

## 3. Local continuity of the ultimate-hub mark

Use the [old-cell cutoff construction](../../../../../../astra/source/round2/operator.md) for any fixed finite time after placement, retaining the marked vertex and its promotion history. After a finite modification, the same pathwise cutoff controls those histories; keeping track of the tag does not create any new route between old cells. Thus bounded tests of a fixed-radius output history are continuous at each fixed input. The factors `L_o`, `1/D_o`, and the finite list of choices in the birth score are locally fixed. Here is why reading the ultimate mark does not require an assumption about continuity at infinite time.

In finite or locally finite initial states, define the hub event by unbounded degree. This is equivalent in the physical finite tree to `Z>0` by the [tagged-family argument](../../../../../../astra/source/round2/fringe/hub-children.md). At the first hit of `d≥2` from below, the new leaf starts a birth–death subfamily with rates `1,1/d`; its survival has probability `1−1/d`. Survival implies that the stock tends to infinity: the supercritical linear birth–death process on survival tends to infinity, for example by the embedded biased walk and nonexplosion. The conditional probability of unbounded degree without leaf stock tending to infinity is therefore at most `1/d` for every `d`. If two permanent adjacent vertices both had bounded degree, both degrees would eventually be constant, yet a Poisson proposal channel of positive fixed rate between them would then eventually give a gain. Thus on a nonhub path, the marked vertex has eventually constant degree and each of its neighbours has diverging stock. These arguments apply in the infinite local construction too.

For `F_T`, observe at horizon `h≥T`. A lower finite-history test requires promotion by `T` and a hit of `d` by `h`, multiplied by `1−1/d`. An upper test requires promotion by `T` and, at `h`, reads the capped sum over neighbours of the [first-arrow charges](../README.md) `b(k)` using the other leaf stocks. Indeed a future hub must gain again and the next gain can only come along an edge already present. On a nonhub path that upper charge tends to zero; on hub paths it has liminf in expectation at least the hub probability by its conditional upper-bound property. More explicitly, its expected excess above the true probability is bounded by its expectation on nonhub paths, which tends to zero. Let `h→∞`, then `d→∞` for the lower test. These continuous, bounded tests per fixed placement squeeze the probability at each input, establishing continuity. `B_T` is already a finite-history test. Section 2 gives boundedness and uniform convergence in `T` for both scores, so `F` and `B` are also bounded continuous. There is no assertion of uniform truncation by postpromotion degree on the nonhub paths.

## 4. Asymptotic and delay laws

By the local limit, `C_n^−1 Σ_o F(G_n,o)→π(F)` almost surely, with the same statement for `F_T`, `B_T`, and `B`. The L² bounded core eigenmartingale and the comparison with `C_n` in the [core-limit proof](../../../../../../astra/source/round2/fringe/almost-sure.md) give uniform integrability of `C_n/n^μ`. Multiplying and taking expectations in the exact fixed-label formula proves `v^(1−μ)p_v→E[X]π(F)` (`n=v−1`). At promotion a fresh leaf initiates a birth–death family with rates `1,1/2`; survival implies hubhood, so `F≥B/2` by conditioning at promotion. The discounted flux identity below gives `π(B)≥μ` and hence positivity and the stated lower bound. As a finite-size check, `Σ_(v≤n) p_v≥ E C_n/2` by the same promotion bound, so the pointwise limit could not have a coefficient smaller than `μ E[X]/2`.

Let `B_o^{hub}(ds)` be the finite-delay hub score measure at a root (cumulative `F_T`), and `B_o^{prom}(ds)` the analogous measure without the hub requirement. For fixed `T`,

```math
 \Pr(S_v\le T\mid Z_v>0)\longrightarrow \pi(F_T)/\pi(F),
 \qquad \Pr(S_v<\infty)\sim\mathbb E[X]\pi(B)v^{\mu-1}.
```

The hub-conditioned limit is a probability law, since the missing tail is `O((1+T)e^(−T/2))`. The unconditioned atom instead tends to one, with the complement just displayed. For a size cutoff `⌊rv⌋`, `r>1`, the waiting time from size `v` to `⌊rv⌋` concentrates at `log r` (variance `O(1/v)`) and is independent of the whole discrete graph path, including hubhood and eventual promotion. Sandwich between ages `log r±ε` and use the density bound to pass `ε↓0`. The limiting conditional tail beyond that size multiple is at most `O((1+log r)/√r)`. The same argument works conditional on eventual promotion using `B_T/B`. This controls arbitrarily late *finite* delays after taking the late-label limit; it never counts `∞` as a finite delay.

In fact this tail bound holds uniformly for all sufficiently large `v` and `r≥e²`. Let `K_v` be the promotion size (`∞` if no promotion), and `U` the elapsed clock time to `⌊rv⌋`. Its mean is at least `log r−1/v`, and its variance at most `1/(v−1)`. Hence `Pr(U≥log r−1)≥1/2` uniformly for large `v`. By independence of `U` from the jump chain,

```math
 \tfrac12\Pr(\lfloor rv\rfloor<K_v<\infty)
 \le \Pr(\log r-1<S_v<\infty)
 \le K\frac{\mathbb E C_{v-1}}{v-1}\frac{1+\log r}{\sqrt r}.
```

The same inequalities hold with the additional condition `Z_v>0` on the events. Dividing by either positive asymptotic probability gives the corresponding uniform conditional bounds, even for growing size multiples.

## 5. Stationary birth-to-promotion balance

The stationary outgoing promotion rate is `π(L/D)=μ`. Express that same rate by time since true birth:

```math
 \int_{[0,\infty)} e^{-\mu s}\,\pi(B_o^{\mathrm{prom}}(ds))=\mu.
```

Indeed start the physical clock at `T_3`, and set `m_u(f)=E Σ_(o core at u) f(G_u,o)`. The eigenmeasure convergence gives `e^(−μu)m_u(f)→aπ(f)`, where `a=lim e^(−μu) E C(u)>0`. To count cores at time `t`, count every postinitial newborn that promotes by `t` once, using the conditional birth score at its birth time. The initial three labels contribute at most three, and the remaining expected count is `∫_0^t m_(t−s)(B_s) ds` by birth compensation and the Markov property at placement. Divide by `e^(μt)`: for each fixed `s`, the integrand tends to `a e^(−μs)π(B_s)`, dominated by a constant times `e^(−μs)` since `B_s≤7` and the normalized mean core count is bounded. Therefore `∫_0^∞ e^(−μs)π(B_s)ds=1`; exchanging the nonnegative integrals gives the claimed identity. The discount is essential: promotions seen at a late observation time come from earlier, smaller birth cohorts. The constant for a **deterministic birth** uses the undiscounted hub-restricted measure `π(F)`, not the promotion flux `μ`.

## Finite replay

`taskset -c 13 nice -n 10 timeout 10s python3 -uB check.py` matched the split flux to the direct target-location weights on 5,912 trees, and checked the reciprocal drift and the two family moment equations at 1,001 stocks. The asymptotic and the infinite-delay truncations use the analytic arguments above.
