> Produced by GPT-6.2 Astra.

# Local certificates for exponents 0.560 and 0.570

Let the random friend tree start from one edge. At each step choose a uniformly random vertex, choose one of its neighbours uniformly, then attach a new leaf to that neighbour. Write $C_n$ for the number of vertices of degree at least two. For every $n\ge3$,

```
(15493968692/17611447440) product_{m=3}^{n-1} (1 + 0.560/m) <= E[C_n]
E[C_n] <= (15349885606/810439425) product_{m=3}^{n-1} (1 + 0.570/m).
```

Consequently the liminf and limsup of log E[C_n]/log n lie in [0.560, 0.570]. This improves [0.543, 0.586]. The interval contains 17/30. The positive-weight check below provides an input for a separate sample-path argument.

## State and rates

Time units below are $n$ times one-step conditional probability. A core vertex $u$ has $L$ leaf neighbours, $c$ core neighbours, degree $D=L+c$, and credit `rho = sum_{core w ~ u} 1/d_w`. A leaf attaches to $u$ at rate $L+\rho$. One of its leaves becomes a new core vertex at rate $q=L/D$. At such a conversion, $u$ changes to $(L-1,c+1)$ and the new vertex has one leaf and one core neighbour, namely $u$, whose degree is still $D$ and whose core degree is now $c+1$. These rules follow by summing the choice probabilities over all neighbours of the attachment site. They apply from $n=3$.

Below a central degree cutoff $H$, the state of $u$ records the multiset of labels of its core neighbours. The cutoff $h_c$ for a neighbour's degree decreases with $c$. The command-line list gives $h_1,h_2,\ldots$, followed implicitly by 2 for all later entries. A vertex with $c=0$ has the empty multiset. Put $K_1=20,K_2=7,K_c=5$ for $c\ge3$.

* A neighbour of degree $d<h_c$ records $(d,k)$, with $k$ its core degree capped at $K_c$. Thus $k<K_c$ gives the core degree exactly, and $k=K_c$ gives a lower bound. The contribution to the central credit is exactly $1/d$.
* A pooled neighbour records a lower bound on its degree, using the greatest available threshold at most the known bound. The thresholds are $h_c$, powers of two strictly between $h_c$ and $H$ starting at 4, and $H$. Among pooled neighbours only one retains the greatest mark. All remaining pooled neighbours reset to $h_c$. A mark $b$ contributes an unknown credit in $[0,1/b]$. Degrees never decrease, so the marks remain valid when unobserved degrees grow.
* When a central conversion increases $c$, project the existing multiset with its known labels to the new cutoffs and insert the newborn type $(2,1)$. Exact degrees below the new cutoff retain their projected core degree. Known degree lower bounds suffice for the pooled labels since the cutoffs are nonincreasing. Keep only the largest pooled mark as above. At birth the sole neighbour is labelled from the creator's post-conversion type.

For a recorded neighbour $w$ of type $(d,k)$ with exact $k$, its gain rate lies between $d-k+1/D$ and $d-k+1/D+(k-1)/2$. Its $d-k$ leaves contribute one each, the central vertex contributes $1/D$, and its other core neighbours each contribute at most $1/2$. A gain changes its degree to $d+1$ and leaves $k$ fixed. A conversion increases $k$ by one with rate $(d-k)/d$. If $k=K_c$, the growth interval has lower endpoint $1/D$ and the same upper endpoint. This follows by maximizing $d-j+(j-1)/2$ over the possible true core degrees $j\ge k$. Conversions of saturated or pooled neighbours preserve their recorded state. The central vertex's own gain increases $L$ by one and preserves its multiset. Repeated copies of one type contribute with their multiplicity.

There is one refinement at $c=1$. If the sole neighbour $w$ has degree and core degree $(2,2)$, also record the other core neighbour $a$ of $w$. This is distinct from the central vertex because the graph is a tree. Record the exact degree $e$ and core degree $j$ of $a$ for $e<6$, or the mark $e\ge6$. The gain rate of $w$ is $1/D+1/e$, with $1/e\in[0,1/6]$ in the pooled case. A gain of $w$ discards the witness and gives type $(3,2)$. While the witness is exact, its gain rate lies in $[e-j+1/2,e-j+1/2+(j-1)/2]$, and its conversion rate is $(e-j)/e$. The first change increases $e$, the second increases $j$. The middle vertex $w$ contributes the known $1/2$.

A core vertex of type $(2,2)$ had degree two throughout its core lifetime and converted its single leaf while of type $(2,1)$. Thus this witness has two possible origins. When the sole neighbour was $(2,1)$ and converts its leaf, its other neighbour is the new $(2,1)$ vertex. When the central vertex is born to a creator formerly of type $(2,1)$, the witness is the creator's previous sole core neighbour. The creator's pre-conversion multiset supplies its label. Here $h_1\ge6$ and $K_1\ge6$, so the required witness degrees below 6 are always exact. The convention also covers the first conversion of the three-vertex star, which creates an edge between two pendant cores. Central conversions away from $c=1$ discard the witness.

For $D\ge H$ the weight ignores the multiset. A vertex can never re-enter the finite region, and every child it creates has the same newborn state with one pooled mark $H$. Denote that state's actual weight by $s$.

## Drift and the infinite region

For a weight $f$, the local drift is the sum over the above events of the rate times the change in the central vertex's weight, plus $q f(\text{newborn})$. In the sum over old core vertices each existing vertex accounts for its own change and each birth is charged once to its creator. Thus the sum of local drifts is exactly $n$ times the conditional expected change of `Phi = sum_u f(u)`. The tree may induce dependencies between rates. Enclosing each term using its own interval is valid by additivity.

For a gain of the central vertex the change in $q$ is $(1-q)/(D+1)$. A conversion changes it by $-1/D$. The generator of the existing vertex therefore obeys

```
-q/D <= Gq = (L+rho)(1-q)/(D+1) - q/D <= 1-q
```

since $0\le L+\rho\le D$. For the lower tail take $f=Bq$, where $B>0$ and $B(\delta+1/H)\le s$. Its drift including births is at least $q(s-B/D)\ge\delta Bq$. For the upper tail take $f=A(1+\lambda q)$, where $A>0$ and $s\le A\lambda(1+\lambda)$. Its drift is at most $A\lambda(1-q)+qs\le\lambda A(1+\lambda q)$. The boundary gain from $D=H-1$ uses this same analytic function at $D=H$. A birth at $D=H-1$ leaves the creator at $D=H-1$.

## Exact checks and reproduction

Run `./reproduce.sh` on Linux with a C++17 compiler supporting signed 128-bit integers and the undefined-behaviour sanitizer. The script pins the compiler and each run to one permitted CPU at nice level 15. It needs no data, Python packages, network access, or saved table. Full output is in `baseline.log`.

The solver proposes table values by backward substitution for fixed newborn rewards followed by a few updates of those rewards. The verifier quantizes to integer units of $10^{-10}$. It repairs each value using only backward successors, then checks every row again with the *actual* newborn weights extracted from the table. Each rational term in the drift is rounded upwards for an upper bound and downwards for a lower bound, including the boundary gain. The right-hand side is rounded in the opposite direction. A zero reported residual satisfies the non-strict inequality. The program checks the analytic tail coefficient against the actual pooled newborn weight separately. Floating-point threshold estimates do not enter the proof of a passing table.

The ordering is decreasing $D$, then decreasing $c$, then reverse lexicographic order on the sorted labels. A central gain advances $D$, a central conversion advances $c$, a neighbour transition advances its label, and a witness transition advances the corresponding refined bridge label. Independent of this construction order, the final pass reads the completed table. The multiset enumeration permits at most one above-cutoff pooled mark and includes every projected state.

The checked builds use the undefined-behaviour sanitizer with immediate termination on error. The calculator also bounds $H\le512$, rational denominator at most 10000, and exponent between $1/4$ and 1. Quantized table weights and the analytic tail coefficient are in $[0,16\cdot10^{10}]$; this range is enforced on proposals, tail coefficients, and repairs. Differences including a newborn or upper tail are at most twice that weight limit in magnitude. The total local transition rate is at most $D+1+c(H+1)+H+1$, allowing one extra witness. Directed term rounding adds fewer than $4H+10$ grid units in a row. In particular the signed residual and the repair step (at most four times its magnitude plus 5) fit in 64 bits. Products for rate fractions, the drift accumulator, tail tests and repair division use signed 128-bit intermediates and checked conversion after division. The largest unreduced displayed low-degree minimum numerator is bounded by $16\cdot10^{10}\cdot 2\cdot512\cdot10000<2^{63}$.

The finite sums or suprema below use all enumerated states, including any inconsistent local configurations. For the lower certificate $h_2=18>2$. The finite minimum on $c\le2$ and either $L\ge1$, or $c=2,L=0$ with a tracked $(2,2)$ neighbour, is exactly $3448168319/10^{10}$. In the lower tail on this subset $q\ge(H-2)/H$, so the global subset minimum is the same. The lower global supremum is $17611447440/10^{10}$ and the root weight at time 3 is $15493968692/10^{10}$. The upper weights are in $[810439425,17543859650]/10^{10}$, with root weight $15349885606/10^{10}$. The counts are 3,550,712 finite states for lower and 5,960,583 for upper.

Let $\mathcal F_n$ be the full growth history through size $n$, including the choices used to retain labels. The respective drift inequalities imply `E[Phi_(n+1) | F_n] >= (1+delta/n) Phi_n` and its upper counterpart. At time 3 the centre is the unique core vertex. Dividing the iterated expectation bounds by the lower global supremum or the upper global infimum gives the products stated above.
