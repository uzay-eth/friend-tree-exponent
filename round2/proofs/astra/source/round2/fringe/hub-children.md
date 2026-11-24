# The children of a fixed hub and the final degrees of its finite-degree children

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 11: HY4 Flash for ranked/intermediate limits; GPT for the fixed-hub extensions; inputs include Claims 3, 6 and 7. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


Write $`D,L,c`$ for the degree, the number of leaf neighbours and the number of core neighbours of a vertex. Use independent rate-one initiating clocks in continuous time. For a vertex present at time $`t`$, its local mass is $`W=\lim_{u\to\infty}e^{-u}D(t+u)`$. In the original tree this is $`e^tYZ_v`$, where $`Y=\lim e^{-s}|T_s|>0`$ and $`Z_v=\lim D_n(v)/n`$. A *silent planted newborn* starts with one leaf, has one exterior parent slot, and receives no arrows from that slot. Its mass in its own time scale is denoted $`W^0`$. Set $`q=\mathbb P(W^0>0)`$. All other edge directions have the physical rates, including incoming arrows from children. Addario-Berry et al. establish fixed-label degree convergence, a hub at one end of every edge, and positive probability to freeze at any specified current degree. [Random Friend Trees, Theorems 3.1, 3.2 and 3.4](https://arxiv.org/abs/2403.20185)

## 1. A vertex with unbounded degree has positive mass

The gain rate at a core is $`a=L+\sum_{w\sim_\mathrm{core}v}1/D_w`$, with a zero term for a silent slot. Thus $`L\le a\le D`$, and $`e^{-u}D(t+u)`$ is a nonnegative supermartingale. Whenever $`L=m`$, the gamma-minus-exponential coupling gives the conditional lower bound $`W\ge G_m-E`$ in a coupling, where $`G_m`$ is gamma of shape $`m`$ and unit scale and $`E`$ is exponential of mean one. Tagging leaves uses only their own birth clocks and the focal rate-one clock. The same argument applies at stopping times and with a silent parent. Consequently the conditional chance of zero mass at the first time $`L=m`$ is at most $`\Pr(G_m\le m/2)+\Pr(E\ge m/2)`$, which tends to zero. Unbounded $`L`$ implies positive mass almost surely. [Tagged-leaf comparison](hub-spectrum.md)

On first reaching degree $`d`$ from below, a newly acquired leaf is present. Its own rate-one initiations grow a tagged family of leaves at the focal core. A tagged member can leave only by being selected by the focal core, at hazard $`1/D\le 1/d`$. Give each tag an independent proposal clock of rate $`1/d`$ and accept a physical promotion with probability $`d/D`$. Discard the tag at every proposal, whether or not physical promotion occurred. A marked subpopulation then follows the linear birth-death process of individual rates $`1,1/d`$. Other directions preserve the physical generator by thinning their own channels. Survival of this family implies unbounded $`L`$. The conditional probability of bounded $`L`$ together with eventual infinite degree is therefore at most $`1/d`$ at the first hitting time of $`d`$. Letting $`d`$ increase proves the dichotomy. A finite-degree core eventually has no leaves: once its last incident edge is born, any permanently remaining leaf would ring again and cause another gain. In a countable initial core this argument holds simultaneously for all its vertices. In the silent newborn it gives $`q\ge 1/2`$.

## 2. A factorial bound on the final degree of a nonhub

For $`D\ge2`$ put
```math
 A_D=\frac{D(D+1)!}{12\,4^{D-2}},\qquad
 V(D,L)=A_D\left(\frac2{D+1}\right)^L .
```
At a gain the ratio of the new $`V`$ to the old one is
```math
 \frac{D+1}{2D}\left(\frac{D+1}{D+2}\right)^L\le\frac{D+1}{2D}<1.
```
At a promotion the ratio is $`(D+1)/2`$. Its rate is $`L/D`$. Substitution of $`a\ge L`$ shows that the drift is nonpositive at every physical state. Stopping on finite degrees and applying Fatou makes $`V`$ a nonnegative supermartingale. On the event of zero mass it eventually equals $`A_{D_\infty}`$ by Section 1. Since $`A_d`$ increases, the following bound holds for any current core and any exterior environment:
```math
 \Pr(W=0,\ D_\infty\ge d\mid \mathcal F_t)\le V(D,L)/A_d.       \qquad\text{(1)}
```
A newborn has $`V(2,1)=2/3`$. The initial centre has $`V(2,2)=4/9`$. In particular the expected number of cores promoted by size $`n`$ with zero mass and final degree at least $`d`$, including the initial centre, is at most $`2\mathbb E C_n/(3A_d)`$. Apply (1) at each of the stopping times of promotion up to size $`n`$. With $`\mathbb E C_n\sim c_2 n^\mu`$, Markov's inequality and Borel–Cantelli on dyadic sizes give almost surely
```math
 \limsup_{n\to\infty}
 \frac{\log\log n}{\log n}
 \max_{v:\,D_n(v)\ge2,\,Z_v=0}D_\infty(v)\le\mu .              \qquad\text{(2)}
```
The empty maximum is zero. Indeed $`\log A_d=d\log d+O(d)`$, and using $`d=(\mu+\epsilon)\log n/\log\log n`$ makes the exceptional probabilities summable on this mesh. The maximum is nondecreasing between mesh points.

The silent newborn has a matching logarithmic lower bound. Before the first promotion it is a star. At degree $`j`$ its gain and promotion rates are $`j-1`$ and $`(j-1)/j`$. Hence it reaches a star of degree $`d`$ without promotion with probability $`2/d`$. Thereafter temporarily suppress all arrows from new core children into the root. In this modified process the probability that all $`d-1`$ leaves are promoted before any one produces a root gain is $`(d+1)^{-(d-1)}`$. At $`i`$ leaves the two competing rates are $`i/d`$ and $`i`$, while arrows into existing children have no effect at the root. On this event the modified root stays forever at $`D=d,L=0`$. Its sends into each child, from that child's promotion, can be driven by an independent rate-$`1/d`$ Poisson clock. These and the child-internal drivers are independent of the competing clocks in the event. This construction uses the separate directed channels of the same generator.

For each new child embed the birth-death process $`B`$ of individual rates $`1,1/2`$ from its initial leaf as in Section 1. The processes can be chosen independently for these disjoint child branches. On survival, $`e^{-u/2}B_u`$ converges to a strictly positive limit, as follows from the elementary linear birth-death generating function, and $`\int_0^\infty du/B_u<\infty`$. Define
```math
\theta=\mathbb E\!\left[\mathbf1_{\{B\ \mathrm{survives}\}}
         \exp\left(-\int_0^\infty du/B_u\right)\right]>0, 
```
taking the integrand to be zero on extinction. In the modified process the rootward reinstatement channel at a child has conditional no-arrow probability $`\exp[-\int du/D_{\rm child}(u)]`$. It can be realised independently conditional on the modified paths. Its probability, with survival of the embedded family, is at least $`\theta`$. On all these no-arrow events the full and modified trees agree. Thus for the silent newborn and every $`d\ge2`$,
```math
 \frac2d\left(\frac{\theta}{d+1}\right)^{d-1}
 \le \Pr(W^0=0,D_\infty=d)
 \le \Pr(W^0=0,D_\infty\ge d)\le\frac{2}{3A_d}.                 \qquad\text{(3)}
```
In particular $`q<1`$ and the logarithm of either probability in (3), divided by $`d\log d`$, tends to $`-1`$.

## 3. Physical offspring of a hub

Fix a physical vertex $`v`$. Start the following construction at its promotion time $`\tau`$ (at time zero for the founder). Its future initiating clock is a rate-one Poisson process with times $`\tau+t_j`$. Independently at every clock point preassign a fresh silent planted newborn process and let $`W_j^0`$ be its limiting root mass. If a point promotes a leaf, use its preassigned branch for that newborn. Arrows from this branch back into $`v`$ have predictable rate $`1/D_{\rm child}`$ and leave the branch unchanged. Generate these arrows using a separate channel. The branch agrees internally with its assigned trajectory until $`v`$ first sends an additional leaf into that core child. On such a send resume the branch with its new state and fresh evolution. A clock point of $`v`$ choosing a core uses no new trajectory. The strong Markov property at promotions allows fresh drivers for the former leaf. This gives the physical process and a Poisson process with independent trajectory marks on a common probability space.

On $`\{Z_v>0\}`$ the outgoing rate into all existing core neighbours is $`c_v/D_v`$. Here $`c_v\le c_v(\tau)+\#\{j:t_j\le t-\tau\}=O(1+t)`$ and $`e^{-t}D_v(t)\to YZ_v>0`$. Its compensator is finite. Localising at bounded accumulated compensator shows that only finitely many such sends occur on this event. Every other clock point promotes a new core, and only finitely many of the resulting branches ever receive a send from $`v`$. Up to finitely many entries, the positive child masses $`Z_w`$ therefore coincide with
```math
 (e^{-\tau}/Y)e^{-t_j}W_j^0,\qquad W_j^0>0.                      \qquad\text{(4)}
```
All children of positive degree fraction eventually promote, so this comparison counts all hub children. The random overall scale is the same across the list. The independent marked Poisson process is constructed before restricting to the event that $`v`$ is a hub. Its conditional law given that event need not itself be Poisson.

By the Poisson mapping theorem the positive points $`e^{-t_j}W_j^0`$ have intensity $`\Pr(W^0\ge x)\,dx/x`$. The mean number above any $`x>0`$ is finite since $`\mathbb E W^0\le2`$; it is asymptotic to $`q\log(1/x)`$ as $`x\downarrow0`$. The strong law for Poisson counts, (4) and the common scale give, simultaneously at every hub,
```math
 \lim_{x\downarrow0}
 \frac{\#\{w:\ w\text{ a child of }v,\ Z_w\ge x\}}{\log(1/x)}
       = q.                                                        \qquad\text{(5)}
```
Also, the fraction of hub children among the first $`m`$ promoted children tends to $`q`$. Among those same $`m`$ children the largest eventual degree of a nonhub is asymptotic to $`\log m/\log\log m`$. For the last assertion assign the value zero to hub marks. The two bounds (3) give log tail $`-(1+o(1))d\log d`$. For independent marks, at $`m=2^j`$ the probability of exceeding $`(1+\epsilon)\log m/\log\log m`$ is at most $`m^{-\epsilon+o(1)}`$, while the probability of missing $`(1-\epsilon)\log m/\log\log m`$ is at most $`\exp[-m^{\epsilon+o(1)}]`$. Both are summable on that mesh; monotonicity and (4)'s finite changes to entire trajectories give the assertion.

The logarithmic count (5) describes the direct children of one fixed hub. The count over all vertices has the power law $`AXx^{-\mu}`$ of [the ranked spectrum](hub-spectrum.md). Finite exceptional children can matter at every vertex in a growing genealogy, so the one-hub construction supplies no equation for $`\mu`$ by itself.
