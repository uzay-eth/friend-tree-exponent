> Produced by GPT-6.2 Astra.

# Existence of the mean growth exponent

Start with a single edge. Write $`C_n`$ for the number of non-leaves after the tree reaches size $`n`$. We show that $`\lim_n \log \mathbb E C_n/\log n`$ exists. The comparison below does not identify its numerical value.

Give each vertex an independent rate-one Poisson clock. When its clock rings, it chooses a uniform neighbour and a fresh leaf attaches to that neighbour. Write $`F(t)=\mathbb E C(t)`$ for the mean core size from the edge. From a finite tree $`T`$ of size $`k`$ use $`C_T(t)`$ and $`N_T(t)`$ for the core and total population. The latter is a Yule process, since the rate of the next birth at size $`j`$ is $`j`$. In particular
`E[N_T(t)^2] = k(k+1) exp(2t) - k exp(t)`, as follows by applying the birth generator to $`j`$ and $`j^2`$.

## Finding places to start a new branch

The core vertices induce a tree whenever their number $`C`$ is positive. For a core vertex let $`c`$ count its core neighbours and $`L`$ its leaf neighbours. Call a vertex of type B one with $`c=2,L=0`$. Call a core vertex eligible if $`c\le2`$ and $`L\ge1`$, or if it is of type B with a neighbour of type B.

There are at least $`C/4`$ eligible vertices. For $`C=1`$ the centre is eligible. For $`C\ge2`$ let $`I`$ count the B vertices with no B neighbour. They form an independent set, and their $`2I`$ incident core edges are distinct. Thus $`I\le(C-1)/2`$. The vertices of core degree at least three are at most the number of core-degree-one vertices, by summing $`c-2`$ on the core tree. Each of the latter has a leaf and is eligible. Every ineligible vertex is either among the $`I`$ or has core degree at least three. The claimed fraction follows.

Fix integers $`M,k\ge2`$ and consider all rooted unlabelled trees $`(T,v)`$ on $`k`$ vertices. There exist $`h<\infty`$ and $`p>0`$, depending only on $`M,k`$, with the following property. For each eligible vertex $`u`$ at time $`t`$, one can arrange, with conditional probability at least $`p`$, that at $`t+h`$ an edge $`uv`$ separates a branch isomorphic to the prescribed $`(T,v)`$, and $`D_u\ge M`$. The branch starts at a leaf of $`u`$ promoted after $`t`$ and contains no vertices which were core at $`t`$.

Here are explicit finite clock requirements.

1. If $`u`$ has a leaf, keep $`u`$ silent for the first unit of time. Otherwise take a B neighbour $`w`$ of $`u`$, keep $`u`$ and the other neighbour of $`w`$ silent for that unit, and require the first ring of $`w`$ to occur and choose $`u`$. Before that ring $`w`$ has degree two. In either case $`u`$ now has a leaf, and still has at most two core neighbours.
2. For a second unit keep $`u`$ silent and require a chosen leaf of $`u`$ to ring at least $`M`$ times. That vertex stays a leaf because only $`u`$ can redirect to it. These rings give $`u`$ at least $`M`$ leaves.
3. In the next unit require exactly one ring of $`u`$. It chooses a leaf with conditional probability at least $`M/(M+2)`$. The chosen leaf becomes $`v`$ with one new leaf on its side of the edge. Keep its two-vertex branch silent for the rest of this unit (probability at least $`e^{-2}`$).
4. Keep $`u`$ silent for a fixed further interval, and grow the two-vertex branch to the specified rooted tree, with $`v`$ as the specified root. Start with any edge incident to that root in $`T`$, then adjoin the other vertices in parent-before-child order. For each addition prescribe exactly one ring of an internal neighbour of the desired attachment vertex, its choice, and silence of the other internal clocks. The chosen source has degree at most $`k+1`$. Specifying unit subintervals for the prescribed births, and silence of the newborn until its subinterval ends, gives a positive lower bound uniform in the exterior. The empty list of additions handles $`k=2`$.

These requirements use only finitely many unit-rate clocks and finitely many choices with probabilities bounded away from zero. Extra arrivals at $`u`$ from outside the branch can only help stages 1–3. A fresh choice in this construction can be made using any fixed ordering of the vertices. Branches successfully planted at distinct core vertices present at $`t`$ are disjoint, including their subsequent growth on the far side of the cut edge: each is attached at just its own old core vertex and contains no old core vertex. This is a statement about any one realised tree. No independence between different successful plantings is required.

## Comparing a planted branch with a free tree

Consider a successful planting. On the far side of $`uv`$, compare the next $`s`$ units of time with a freestanding copy of $`T`$. The clocks, birth choices and new vertex identities inside the branch can be coupled up to the first discrepancy. At $`v`$ a reference ring chooses an internal neighbour uniformly. Retain this choice with probability $`D_v/(D_v+1)`$, where $`D_v`$ is the current **internal** reference degree. Otherwise choose the exterior neighbour $`u`$. Use independent acceptance marks for the rings at $`v`$. All other internal choices agree.

For rings of $`u`$, put the edge to $`v`$ first in its ordering, using an independent uniform mark. If none of these marks belongs to $`[0,1/M]`$ during the comparison, $`u`$ never chooses $`v`$, since $`D_u\ge M`$ forever. This event has probability $`e^{-s/M}`$ and is independent of the reference history. On this event and on acceptance at every reference ring of $`v`$, the two internal trees agree throughout. Any reference core vertex is then an actual core vertex. Thus a reference history, including its clock sources, contributes its core count with at least the weight
```
exp(-s/M) product over rings of v of D_v/(D_v+1).
```
The exterior tree may evolve freely subject to the stated mark condition. No independence of the unconditioned sides of the original edge was used.

## A high-degree vertex among the starting vertices

Fix $`0<\epsilon<1`$ small. Choose a constant $`A`$ with $`A\log A-A+1>5`$ (for example $`A=8`$). Choose $`b=b(\epsilon)>\epsilon`$ so that $`b\log(b/\epsilon)-b+\epsilon>5`$ and $`b\to0`$ as $`\epsilon\downarrow0`$. For instance $`b=7/\log(1/\epsilon)`$ works for all sufficiently small $`\epsilon`$. Fix $`k\ge 6M/\epsilon`$. All these quantities are fixed before sending $`s`$ to infinity.

For a reference run from any $`k`$-vertex tree, call a history good if each of the $`k`$ original clock sources rings at most $`bs`$ times by $`\epsilon s`$ and at most $`As`$ times by $`s`$, and at least one original vertex has degree at least $`M`$ by $`\epsilon s`$. Its complement has probability $`O_{k,M,\epsilon}(e^{-4s})`$ for large $`s`$. The ring counts follow from the Poisson Chernoff bound $`\Pr(\operatorname{Pois}(as)\ge xs)\le e^{-s(x\log(x/a)-x+a)}`$ for $`x>a`$. For the degree condition, rank one fixed original neighbour first at each original source. Until an original degree reaches $`M`$, every original clock mark in $`[0,1/M]`$ adds one to an original degree. There can then have been at most $`kM`$ such marks if all original degrees stay below $`M`$. The marks up to $`\epsilon s`$ form a Poisson count of mean $`k\epsilon s/M\ge6s`$. Its lower tail at $`kM`$ is at most $`e^{-4s}`$ for large $`s`$. Both orderings in the comparison may be used: at an internal reference ring the choice among internal neighbours uses a separate mark from the acceptance mark for the extra edge.

By Cauchy–Schwarz and the Yule second moment, bad histories contribute $`O_k(e^{-s})`$ to $`\mathbb E C_T(s)`$, uniformly over the finitely many starting trees, with constants allowed to depend on $`k,M,\epsilon`$. Assign each good history to one original vertex $`v`$ that has reached $`M`$ by $`\epsilon s`$ (break ties by a fixed numbering). For histories assigned to $`v`$ the comparison weight is at least $`e^{-\eta s}`$, where
`eta = b log(2) + (A+1)/M`. Before $`\epsilon s`$ each acceptance costs a factor at worst $`1/2`$. Subsequently $`D_v/(D_v+1)\ge e^{-1/M}`$. Consequently, for each $`T`$ and $`s`$ some original $`v`$ accounts for at least $`(\mathbb E C_T(s)-O(e^{-s}))/k`$ in unweighted good core count. The choice of $`T,v`$ is made from the laws, before planting.

## Iteration and return to discrete size

Stop the tree from an edge when it reaches size $`k`$. If it has not reached $`k`$ by $`s`$, then $`C(s)\le k`$. Otherwise the strong Markov property and monotonicity give
`F(s) <= k + max_T E[C_T(s)]`,
where the maximum runs over the finite set of $`k`$-vertex trees. Choose a maximizing $`T`$ and a root $`v`$ as above. Plant that same rooted seed at each eligible old vertex in the preceding first-moment bound, sum over the disjoint branches and average over the tree at $`t`$. For all $`t\ge0`$ and sufficiently large $`s`$,
```
F(t+h+s) >= F(t) [p/(4k)] exp(-eta s) [F(s)-k-O(exp(-s))].
```
The comparison for each vertex is a marginal lower bound on the actual descendant count in its own branch; there is no requirement that every planting succeed.

Here $`F(s)\to\infty`$. Indeed, after the first birth, the total rate of promoting a leaf is $`\sum_u L_u/D_u\ge1`$: for a single core it equals one, and otherwise the core tree has at least two end vertices, each contributing at least $`1/2`$. Now take logarithms. Write $`B_s`$ for the entire multiplier of $`F(t)`$ on the right of this inequality. Iteration over blocks of length $`h+s`$, and monotonicity between block endpoints, give $`\liminf_{t\to\infty}t^{-1}\log F(t)\ge (h+s)^{-1}\log B_s`$ for fixed large $`s`$ (start at any time with $`F(t)>0`$). Taking a sequence of $`s`$ realising the limsup gives
`liminf log F(t)/t >= limsup log F(t)/t - eta`.
The exponents are finite by the Yule bound. Take $`M\to\infty`$ and $`\epsilon\downarrow0`$ to obtain the common limit $`\mu`$.

Finally let $`a_n=\mathbb E C_n`$. The holding time at size $`j`$ is exponential of rate $`j`$, independent of the uniform clock source and its neighbour choice, so $`F(t)=\mathbb E a_{N(t)}`$. From two vertices, $`N(t)`$ is the sum of two independent geometric variables on $`\{1,2,\ldots\}`$ with success probability $`e^{-t}`$. At $`t=\log n`$, $`\Pr(N(t)\ge n)`$ is bounded away from zero, so monotonicity gives $`a_n=O(F(\log n))`$. At $`t=(1-\alpha)\log n`$ for any $`0<\alpha<1`$, $`\mathbb E[N(t);N(t)>n]`$ decays faster than any power of $`n`$ by the geometric tail. Thus $`F((1-\alpha)\log n)\le a_n+o(1)`$. Since $`F((1-\alpha)\log n)\to\infty`$, letting $`\alpha\downarrow0`$ gives $`a_n=n^{\mu+o(1)}`$.
