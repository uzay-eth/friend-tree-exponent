# First differential sends along an existing path

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. Claim 14: Claude for transfer/scalar identities; GPT for full contraction and sign transfer; inputs include Claim 12 and the cap-48 continuation certificate. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


In the paired intrinsic friend trees, mark a particular existing outgoing
core edge at a vertex that has just acquired its first extra neighbour.
Let $`\tau`$ be the first differential send through that edge. With $`d`$ the
degree in the smaller copy at the initial time, the following bound is
uniform over the surrounding trees and their previous evolution.

```math
 (d-1)\,\mathbb E[e^{-s\tau};\tau<\infty]\le q_s.
```

The rational values $`q_{.39}=.986`$, $`q_{.5}=.710`$ and $`q_{.565}=.605`$ work.
Each value also works at larger discounts. The tagged recipient must be
uninfected, and the source may have one previously infected neighbour at
the initial time. All other old and future contacts continue to evolve.

## A two-count comparison

Use the common directed arrows at the smaller of the two actual reciprocal
rates, with a separate channel for the surplus. This gives the physical
marginals and preserves the nested neighbour sets proved in the
parity coupling. Write $`k`$ for the number of private
neighbours of the source, $`i`$ for the number of its infected common
neighbours, and $`B`$ for its degree on the smaller side. A silent parent
slot counts toward $`B`$. Before the tagged contact, $`B\ge d`$ and the
following conditional intensity bounds hold.

* Each private neighbour sends back at rate at most one. Each infected
  common neighbour contributes a difference at most one. Thus the rate
  of an increase of $`k`$ is at most $`k+i`$.
* A common neighbour can become infected only through a differential
  send from this source. The total rate of the contacts other than the
  tagged one that increase $`i`$ is at most $`k/(B+k)\le k/(d+k)`$.
  This includes neighbours attached after the initial time.
* The hazard of the tagged contact is $`k/[B(B+k)]\le k/[d(d+k)]`$.

Dominate the two counts by $`(K,I)`$ from $`(1,1)`$, with increments at rates
$`K+I`$ and $`K/(d+K)`$ respectively. Add a separate proposed tagged contact
at rate $`K/[d(d+K)]`$ and stop the comparison at its first occurrence.
All three bounds increase with the counts. Thinning the three dominating
clocks separately couples the systems until the comparison stops. The
physical tagged contact cannot precede that stop. Other physical events,
including real root rings, are included at their own marginal rates.
Starting the second count at one also covers the initial leaf insertion,
which has no infected common neighbour. Both pure-birth counts are
nonexplosive since their combined rate grows at most linearly.

For computation retain $`K\le300`$ and $`I\le24`$. At either upper exit the
discounted first-hit transform is at most $`1/(1+sd)`$ since the remaining
tagged hazard is at most $`1/d`$. At an interior point put
$`a=K+I`$, $`b=K/(d+K)`$ and $`h=b/d`$. The exact first-event numerator is
$`h+aP(K+1,I)+bP(K,I+1)`$ with denominator $`s+h+a+b`$. An integer ceiling
at scale $`10^{12}`$ on each row gives upper values, which may also be
intersected with the same uniform boundary. This checks $`2\le d\le255`$.

For larger degrees $`\mathbb E I_t\le1+t`$ and
$`\mathbb E K_t\le 3 e^t-t-2`$. The first-hit transform is bounded by the
discounted integrated hazard. After multiplying by $`d-1`$ its integrand
is at most $`e^{-st}\min(1,3e^t/d)`$. Integrating across $`\log(d/3)`$ gives
an upper bound $`(3/d)^s/[s(1-s)]`$ for $`0<s<1`$. The script checks each
large-degree comparison by taking integer powers of rationals.

## Cutting incoming arrows at depth

Let $`V_s b(S)`$ be the intrinsic discounted occupation of the root leaf
fraction. If one leaf is added at a core $`v`$ of depth $`h\ge1`$, let
$`v=v_h,v_{h-1},\ldots,v_1,o`$ be its path to the root. Then

```math
 |\Delta_v V_s b(S)|\le
 {1\over s}\prod_{j=1}^{h}{q_s\over D_{v_j}(S)-1}.
 \qquad\text{(1)}
```

Indeed the two root rewards agree until the first root infection.
Before that time the rooted core has the same root degree and the
same leaf neighbours on both sides. A current structural leaf attached
to the root could first become different only through its root's arrow.
The initial insertion is at an old core. The infected common sites form
a connected subtree, so first infections on this old path occur in order.
At each one apply the tagged bound conditionally on the whole past. The
current degree at that time bounds the initial degree from below.
Once infection reaches the root, the absolute remaining reward is at
most $`1/s`$.

Consider the intrinsic generator obtained by omitting all arrows from
real core children below level $`h`$ into their parents at level $`h`$.
It keeps track of the real child slots at the boundary, including
promoted children, and can ignore evolution below those slots. Let
$`\bar V_s`$ denote its discounted integral and let $`c_o`$ be the number
of immediate real core children. A state-dependent upper error is

```math
 |(V_s-\bar V_s)b(S)|
 \le {1\over s}\bar V_s\left[
 \sum_{|v|=h}\frac{c_v^{\rm child}}{2}
 \prod_{j=1}^{h}\frac{q_s}{D_{v_j}-1}\right](S).
 \qquad\text{(2)}
```

To see this, use discounted Dynkin for the bounded function $`V_s b`$
under the modified dynamics. Every omitted gain has rate at most
$`1/2`$ and adds a leaf at the receiver. Equation (1) bounds its
contribution. Finite stopping can be removed because the integrand
is at most $`q_s^h c_o/(2s)`$, by summing successively through each
parent's children and using $`c_v^{\rm child}\le D_v-1`$. Moreover
$`\mathbb E c_o(t)\le c_o(0)+t`$. The modified process is dominated in
population size by a Yule birth process. These bounds also give the
uniform-in-degrees estimate

```math
 |(V_s-\bar V_s)b(S)|\le\frac{q_s^h}{2s}
 \left(\frac{c_o(S)}s+\frac1{s^2}\right).
```

For a virgin $`J`$, the same sum after removal of $`q_s^h/(2s)`$ is at
most the number of level-$`h`$ cores. Every new core comes from a
promotion at a rate at most one per parent, so its expected value is
at most $`t^h/h!`$. One may integrate the smaller of $`t`$ and $`t^h/h!`$.
For instance, on $`[0,T]`$ use
$`e^{-st}\le e^{-sT}(T/t)^{sT}`$ and $`sT<h+1`$. The resulting error is
at most

```math
 \frac{q_s^h e^{-sT}}{2s}
 \left(\frac{T^{h+1}}{(h+1-sT)h!}
 +\frac{T}{s}+\frac{1}{s^2}\right).
```

[The integer replay](certificate.py) takes $`T=2h/5`$ and rounds the
exponential upward by $`(1+sT/256)^{-256}`$. These estimates cover the
direct intrinsic root integral. Pricing $`V_s C H_s`$ in the signed
newborn equation also requires control of the whole-child endpoint score.
