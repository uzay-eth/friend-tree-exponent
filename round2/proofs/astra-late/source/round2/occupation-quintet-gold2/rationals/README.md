> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# A rational friend-tree exponent would have denominator at least 426

For the original edge-seeded undirected random friend tree, let $`\mu`$ be the core-count growth exponent. The original Neon-2 Gold2 upper witness and three physical occupation bounds on its lower witness give
```math
 \frac{99}{175}<
 \frac{126905728962161948969119589499374909399}
      {224328301684105796116575000000000000000}
 \le \mu \le
 \frac{127993535256830043}{226242120050000000}
 <\frac{142}{251}.
```
The two outer fractions are Farey neighbours of order 425. Consequently, if $`\mu=a/b`$ in lowest terms, then $`b\ge426`$. The fraction $`241/426`$ is inside the displayed enclosure. Irrationality and transcendence remain open. [Original Gold2 production audit](../../../../../../reviews/cold/astra/attacks/neon2-quintet/gold2/README.md) · [Independent physical audit](../../../../../../reviews/cold/astra/checks/gold2-h3-sol.md) · [Independent arithmetic audit](../../../../../../reviews/cold/astra/checks/gold2-rational-sol.md)

## Three credits on the original lower witness

A core is a vertex of ordinary degree at least two. Write $`C(T)`$ for the number of cores, $`D`$ for ordinary degree and $`c`$ for the number of core neighbours, including the parent of a nonfounder. All three prices use the stationary law at a uniformly counted physical core. One potential unit is $`S=2^{40}`$ ticks. The original certificate and its signed comparison constants are pinned at `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`. [Production receipt](../../attacks/neon2-quintet/gold2/production.json) · [Physical spatial law](../../../../../astra/source/round2/fringe/spatial-limit.md)

| Original allocation | Certified ticks per core | Retained occurrences |
|---|---:|---|
| Signed three-core path (H3) | $`1084423030952180046412413909/1550425090937500000000`$ | The path runs from a nonfounder mother, $`19\le D\le79,\ c\le8`$, through her child, $`D\le5,\ c\le4`$, to that child's child, $`D\le18,\ c\le3`$. |
| Fenchel–Young term, small paying neighbour | $`39526146756920663734311/10^{17}`$ | Paying neighbour degree at most five, host degree 19–251 and arity at most 17. The absorbing donor class includes every degree from 12 upwards. |
| Fenchel–Young term, larger paying neighbour | $`1046834122639903380009/(2\cdot 10^{16})`$ | Paying neighbour degree 6–34, arity at most seven, donor degree at most five, host degree at most 150 and arity at most ten. At least one paying-edge endpoint has degree above 18. |

For H3 let $`P`$ be the outward signed physical row after the transfers between edge, path and higher families. The original verifier charges each occurrence to its middle with a nonnegative local cap $`K\ge P`$. Thus $`P+(-P)_+=\max(P,0)\le K`$. Each paid three-label set has a unique chronological middle and a unique mother. The two Fenchel–Young credits are on original ordered edge terms at exact G2 hosts outside M3. At the actual reciprocal of the distinct donor they recover part of the difference from the original knot maximum. The current degree of the paying neighbour separates the two sets of ordered terms. Each physical gain has the same receive rate in the signed reward and in the corresponding transition. These properties permit summation even when several occurrences share vertices. [H3 support and chronological registration](../post-fy/README.md) · [Small-payer computation](../edges/README.md) · [Larger-payer computation](../payer/README.md) · [Independent FY audit](../../../../../../reviews/cold/astra/checks/gold2-fy-sol.md)

The occupation proofs register a new core and its actual incident labels at promotions, evolve the saved degrees and stocks on the physical tree, and check the discounted generator inequalities at both endpoints of every ranged receive rate. Finite-step comparisons use the original Gold2 upper endpoint for the discount and its lower endpoint for the newborn rate. Summing the prefix, edge and completed-path continuation inequalities cancels the intermediate entrance terms. One source remains at each actual core birth. Its mean per core is $`\mu`$, since in continuous time each promotion increases $`C`$ by one. The exact fractions in the table already include this normalization. [Cold replay of all three physical prices](../../../../../../reviews/cold/astra/checks/gold2-h3-sol.md)

## Converting ticks to an exponent

Let $`\lambda=7071421/12500000`$, $`E=55+792+18370=19217`$ and $`m_{\max}=1808603193693`$, all from the original lower receipt. Denote its total potential in ticks by $`P(T)`$, the summed credited reward rate by $`R(T)`$ and the rate-one-per-vertex generator by $`L`$. The original signed comparison with the three payments restored reads
```math
 LP(T)-\lambda P(T)\ \ge\ R(T)-E\,C(T),\qquad
 0<m_{\min}C(T)\le P(T)\le m_{\max}C(T)
```
from the three-vertex tree onwards. The original unsigned potential allocation is used throughout. Each paid reward is a bounded per-core function of a physical neighbourhood, with bounded registration multiplicity at its host.

Put $`F(t)=\mathbb E P(T_t)`$ and $`K(t)=\mathbb E C(T_t)`$. The stable-age mean law gives $`e^{-\mu t}K(t)\to k>0`$ and a limiting per-core expected reward at least the sum $`r`$ of the three floors. Integrate the displayed drift after multiplication by $`e^{-\mu t}`$ and divide by time. The boundary term $`[e^{-\mu t}F(t)]_0^T/T`$ vanishes by comparability. Thus
```math
 (\mu-\lambda)\,T^{-1}\!\int_0^T e^{-\mu t}F(t)\,dt
 \ \ge\ k(r-E)+o(1).
```
Here $`r-E>0`$, so $`\mu>\lambda`$. The limsup of the positive average on the left is at most $`m_{\max} k`$, giving $`\mu\ge\lambda+(r-E)/m_{\max}`$. This argument uses a single comparison for the total potential. The size-two seed enters this process by its first birth, a finite delay that preserves the exponent. [Core-count mean theorem](../../../../../astra/source/round2/fringe/asymptotics.md) · [Signed accounting audit](../../../../../../reviews/cold/astra/checks/gold2-rational-sol.md)

Even rounding the three prices down to 699,000, 395,000 and 52,000 ticks leaves $`1,126,783`$ ticks/core after subtracting $`E`$. At $`99/175`$ the trial margin is then
```math
 1,126,783-\left(\frac{99}{175}-\lambda\right)m_{\max}
 =\frac{2737543234271}{87500000}>31,286.
```
Using the exact three fractions gives more than 32,325 ticks/core and the lower endpoint above. On the upper witness the original adjusted endpoint lies below $`142/251`$ by exactly $`3697635659207/56786772132550000000`$. This calculation uses three occupation floors and the two original certificates. The numerical value of a finite-window estimator and the choice of a stationary linear program play no role.

## Denominator exclusion

Since $`142\cdot175-99\cdot251=1`$, a fraction $`a/b`$ strictly between these two numbers has positive integer gaps $`k=175a-99b`$ and $`\ell=142b-251a`$. Substitution gives $`b=251k+175\ell\ge426`$. Equality forces $`k=\ell=1`$, giving $`241/426`$. The exact checker also enumerates the fractions in both closed certificate intervals. Gold2 alone admits only $`99/175`$ through order 425. After the three credits the least possible denominator is 426.

Run `python3 -IB round2/occupation-quintet-gold2/rationals/check.py` from the Astra root. It authenticates the original JSON receipt and all three recorded reward inputs, verifies their units, recomputes both sets of margins, and enumerates the fractions. The paths and commands for full physical replays are in the independent audit. [Checker](check.py) · [Physical replay commands](../../../../../../reviews/cold/astra/checks/gold2-rational-sol.md)
