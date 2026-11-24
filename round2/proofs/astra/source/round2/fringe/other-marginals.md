# Core degree and descendant size in the limiting fringe

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 10: HY4 Flash, with Claude finite checks and GPT proof review in Fable; inputs include Claims 2–7. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


Write $`c`$ for the number of non-leaf neighbours, $`D`$ for the ordinary degree, and $`N`$ for the number of ordinary vertices in the finite descendant fringe of a core, counting the core itself. Under the [limiting neighbourhood law](spatial-limit.md), with the unique end as parent direction, there are constants $`B>0`$ and $`A_N>0`$ such that
```math
 \Pr(c\ge k)\sim B(1+\mu)^{-(k-1)},\qquad
 \Pr(c=k)\sim B\frac{\mu}{1+\mu}(1+\mu)^{-(k-1)},              \qquad\text{(1)}
```
```math
 \Pr(N\ge n)\sim A_N n^{-\mu},\qquad
 \Pr(N=n)\sim\mu A_N n^{-1-\mu}.                               \qquad\text{(2)}
```
In particular the exponential decay rate of the degree in the tree *induced by the non-leaves* also identifies $`\mu`$. The ordinary degree has the power law obtained in [the degree argument](../degree-law.md). Each probability in (1)–(2) equals the almost-sure limiting proportion among physical cores for its fixed threshold or point.

## Core degree

Use the birth measure $`\beta=\pi B_{\rm birth}/\mu`$ of the spatial mean operator. The [stable-age formula](../degree-law.md) follows by compensation for every bounded nonnegative observable of the rooted environment. Under it the focal vertex is newborn at age zero with $`c=1,D=2`$, and is then observed at an independent age of exponential law with rate $`\mu`$. Give this focal vertex an independent rate-one Poisson clock of *promotion proposals*, accepting each with the current probability $`b=L/D`$. All the outgoing arrows to old cores can be driven by separate clocks as in the spatial construction; the proposal clock here contributes only promotions. The other physical transition rates keep their dependence on the actual current tree.

Let $`T_j`$ be the time of the $`j`$th accepted proposal and $`R_j`$ the number rejected up to $`T_j`$ when this time is finite. Put $`r=1+\mu`$. Change only this one Poisson clock to rate $`r`$ and keep the acceptance rule and all other conditional rates. On a finite interval $`[0,t]`$ the density of the original law relative to the changed law is $`\exp((r-1)t)r^{-P(t)}`$, where $`P(t)`$ counts its proposals. The two laws can be constructed using identical remaining independent drivers. Stop first at $`T_j\wedge t`$ and then increase $`t`$. The identity on the event $`T_j\le t`$ and $`r-1=\mu`$ give
```math
 r^j\Pr_\pi(c\ge j+1)
 =r^j\mathbb E_{\beta}[e^{-\mu T_j};T_j<\infty]
 =\mathbb E_{\beta}^{(r)}[r^{-R_j};T_j<\infty].                  \qquad\text{(3)}
```
These quantities decrease to $`B=\mathbb E_{\beta}^{(r)}r^{-R_\infty}`$, with $`r^{-\infty}=0`$. The proposal count tends to infinity almost surely, so finitely many acceptances implies infinitely many rejections.

Here $`B>0`$ uniformly over the initial newborn environment in the sense that every such environment gives positive probability of finitely many rejections. Indeed the focal vertex can acquire $`m>r`$ leaves in a fixed first time unit by silencing its promotion proposals and requiring rings from its tagged leaf. Start from any such state. Write $`\rho=\sum_{y\sim_{\rm core}o}1/D_y`$, so that its leaf drift is $`L+\rho-rb`$. Its core degree is at most its initial value plus $`P(t)`$, and $`0\le\rho\le c/2`$. Its ordinary degree rises at rate $`L+\rho\le D`$ and stays fixed at promotions. Thus it is dominated by a rate-one Yule population started at its current degree. The local martingale part of $`e^{-t}L(t)`$ has finite expected quadratic variation over $`[0,\infty)`$, since the unscaled total leaf jump rate is at most $`D+r`$. Its integrated remaining drift $`\int_0^\infty e^{-u}(\rho-rb)\,du`$ is absolutely integrable in mean. Consequently $`e^{-t}L(t)`$ has a nonnegative almost-sure and $`L^1`$ limit $`W_L`$ whose mean is at least $`m-r>0`$. On $`W_L>0`$, the eventual exponential lower bound for $`L`$ and the almost-sure linear bound on $`c`$ give $`\int_0^\infty r c(t)/D(t)\,dt<\infty`$. This is the rejection compensator. The number of rejections is finite almost surely on that event, by stopping when the compensator reaches successive integer levels. This proves strict positivity in (3). Taking successive differences gives the point estimate in (1).

## Ordinary size of the downward branch

In a one-ended initial tree the end remains unique over every bounded time interval. Each old vertex gains a finite attached cell and the old edges are the only edges between such cells. A focal newborn has a finite downward branch and one boundary edge to its promoter. Let $`p_v=1/D_v`$ and $`p_p=1/D_{\rm parent}`$. Every initiation inside the branch adds a vertex there except for an arrow from $`v`$ to its parent. The reverse arrow is its only external input. Thus $`N(0)=2`$ and the rate of the increase $`N\mapsto N+1`$ is precisely $`N+p_p-p_v`$. Both boundary degrees are at least two. In particular $`e^{-t}(N(t)+1/2)`$ is a nonnegative supermartingale and has an almost-sure limit $`W_N`$.

The increasing process is stochastically sandwiched on whole paths by pure-birth chains with respective rates $`n-1/2,n+1/2`$ at count $`n`$, both started at 2. Couple by the same Poisson marks and order the introduction times of the next particle, or use thinning on a common random-time measure. At a tie the lower rate is below the physical predictable rate and the upper rate is above it; unit jumps cannot cross distinct ordered integer states. For the chain of rate $`n+\epsilon`$, $`e^{-t}(N_\epsilon(t)+\epsilon)`$ is a nonnegative martingale bounded in $`L^2`$. Its limit has the gamma distribution of shape $`\alpha=2+\epsilon`$ and unit scale: the finite-time excess over 2 is negative binomial, with probability generating function $`[e^{-t}/(1-(1-e^{-t})z)]^{\alpha}`$ by its forward equations. Its rescaled Laplace transforms tend to $`(1+x)^{-\alpha}`$. Therefore $`W_N>0`$ almost surely, and
```math
 \frac{\Gamma(3/2+\mu)}{\Gamma(3/2)}\le A_N:=\mathbb E_{\beta} W_N^\mu\le
 \frac{\Gamma(5/2+\mu)}{\Gamma(5/2)}.                            \qquad\text{(4)}
```
The upper martingale maximal inequality also bounds $`\sup_t(e^{-t}N(t))^\mu`$ by an integrable variable. The hitting times $`H_n`$ of $`N=n`$ are finite and tend to infinity. At those times $`ne^{-H_n}\to W_N`$. The exponential age mixture and dominated convergence give $`n^\mu\Pr_\pi(N\ge n)\to A_N`$.

Finally the stable-age formula applied to the intensity of crossing from $`n`$ to $`n+1`$ gives, for $`n\ge2`$,
```math
 \mu\Pr_\pi(N\ge n+1)=
    \mathbb E_\pi[(n+p_p-p_v)\mathbf1_{\{N=n\}}].
```
Indeed the newborn starts at size two and crosses the given level exactly once, with density of its discounted crossing equal to its predictable jump intensity integrated up to the independent age. The factor on the indicated level lies between $`n-1/2`$ and $`n+1/2`$. This gives the point estimate in (2).
