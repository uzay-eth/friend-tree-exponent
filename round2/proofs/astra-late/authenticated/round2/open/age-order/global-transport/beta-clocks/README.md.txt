# Finite-size waiting laws and a conditional birth-scale cutoff

Write $b$ for the older of the next two birth labels, $m=b+1$ for their common completion size, and $Z_i$ for label $i$'s limiting degree share. For every old tree of size $b-1$ with $b\ge20$,
\[
  0<x<1,\qquad mx\ge23/10
  \quad\Longrightarrow\quad
  \Pr(Z_b>x\mid T_{b-1})>\Pr(Z_{b+1}>x\mid T_{b-1}).
\]
The bound holds at every finite size in its stated range. The critical local comparison is between the serial middle, born with one fresh leaf, and the two adverse leaf placements on old edges. We improve the middle's tail by accounting for the smaller chance that it loses its tracked leaf before the next gain. This uses the physical discrete-size process directly. [Poisson stock comparison](../../hub-ranking/poisson-tail/README.md) · [Signed placement identity](../../adjacent/README.md) · [Earlier cutoff of nine halves](../../hub-ranking/beta-tests/physical-tail/poisson/README.md)

## Exact law of a waiting stage

Write $B_{m,r}(x)=\Pr(\operatorname{Beta}(r,m-r)>x)$ for $0<r<m$, and use mass one for the random variable when $r=m$. A monitored event has probability $h/s$ at size $s$ until its first occurrence, $0<h\le m$. Let $N$ be the size *after* that event. If the continuation then pays a random variable with law $\operatorname{Beta}(r,N-r)$, its unconditional law is
\[
 U_h\,Y_{m,r},\qquad U_h\sim\operatorname{Beta}(h,1),\quad
 Y_{m,r}\sim\operatorname{Beta}(r,m-r)
\]
with the two variables on the right independent. Indeed the $p$-th moment of the proposed answer is $T_m=(r)_p h/[(m)_p(h+p)]$. It satisfies $mT_m=h(r)_p/(m+1)_p+(m-h)T_{m+1}$. The no-event probability tends to zero since the hazards have divergent sum. Backward iteration gives uniqueness of bounded solutions. Integer $p\ge1$ and equality of total mass determine the law on $[0,1]$.

Successive fixed races of intensities $h_1,\ldots,h_j$ therefore insert independent factors $U_{h_i}$. A branch of weight $a$ in an $h$-race has conditional probability $a/h$, independent of its exit size. One can first fix an embedded sequence of branch choices when its rates depend only on the already reached states. Throughout the calculation below its terminal shape is at most the initial size $m$. Equivalently the logarithm of the inverse product has the law of a sum of independent exponentials of rates $h_i$. These are distributional identities for the constructed races.

## Extra serial credit before the first tracked loss

Let $q=1-x$. From a serial completion choose the fresh leaf and track only its family. At stock $j$ before the first loss, retain gains initiated by the $j$ designated leaves, with rate weight $j$. Each creates a new leaf at the tag. The tag has degree at least $j+1$, since its original neighbour persists. Its send to a designated leaf has weight at most $d_j=j/(j+1)$. Supplement actual losses to that rate using unused births and discard one tracked leaf at a supplementary loss. The designated leaves stay real and distinct. The total virtual weight $h_j=j+d_j$ fits into the size $s$. Other births leave the virtual state fixed. At a loss, or after the chosen number of pure gains, use the existing Poisson stock payoff $H_{s,j}$ on the surviving subset.

The Poisson minorant is exactly harmonic with gain weight $j$ and loss weight one. Each size spent before the first loss gives a positive source
\[
 \frac{1-d_j}{s}\,[H_{s+1,j}(x)-H_{s+1,j-1}(x)].
\]
For a visit to class $j$, summing over the times spent there gives $(1-d_j)/h_j$ times this payoff difference at the next exit size. Reaching that class requires $j-1$ gains with probability $a_j=\prod_{i<j}i/h_i$. Bounded optional stopping at the finite sequence of almost surely finite exits gives the original $H_{m,1}=e^{-1}q^{m-1}$ plus all these sources. The construction also works for a current tag with one tracked leaf and another permanent neighbour.

For the first class alone there is a simple bound valid at every $m\ge3$ and $0<x<1$:
\[
 \frac{\Pr(Z_{\rm middle}>x)}{q^{m-1}}
 \ \ge\ e^{-1}\left(1+\frac{q}{2mx+5q}\right).
\]
The number of size steps to the first event has the law of a geometric variable with a random success probability $P\sim\operatorname{Beta}(h,m-h)$, for $h=3/2$. Its survival for $l$ steps is $(m-h)_l/(m)_l$. Thus $\mathbb E[q^{N-m}]=\mathbb E[A/(t+A)]$, where $A=mP$, $t=mx/q$, $\mathbb E A=h$, and $\mathbb E A^2\le h(h+1)$. Cauchy–Schwarz gives $\mathbb E[A/(t+A)]\ge h/(t+h+1)$. Multiplication by $(1-d_1)/h$ proves the formula.

For each source in the longer calculation, $H_{N,j}-H_{N,j-1}$ is $e^{-1}\sum_{k<j}(B_{N,r}-B_{N,r-1})/k!$, $r=j-k$, with $B_{N,0}=0$. Use the exact race law for both terms with the same factors. Given $v=(\prod_{i\le j}U_{h_i})^{-1}-1$, their difference at initial size, divided by $q^{m-1}$, is the positive binomial mass
\[
 {m-1\choose r-1}(x/q)^{r-1}(1+v)^{r-1}
       \big(1-(t/m)v\big)_+^{m-r}.                                      \tag{1}
\]
For $y=mx\in[l,u]$, $m\ge M=21$, use $l\le t\le u/(1-u/M)$. The binomial coefficient multiplied by $(x/q)^{r-1}$ is at least ${M-1\choose r-1}(l/M)^{r-1}$. For $0\le z<M$,
\[
 (1-z/m)^{m-r}\ge \min\{(1-z/M)^{M-r},e^{-z}\}.
\]
Indeed the sign of the derivative of its logarithm is the sign of
$[w+(1-w)\log(1-w)]/w^2-r/z$ at $w=z/m$. The first term is increasing in $w$ by its positive-coefficient series. As $m$ grows the derivative can change only from positive to negative. At $z\ge M$ zero is a valid lower bound. These two inequalities give lower rectangles for (1), including both threshold endpoints.

The mesh is $1+v=(129/128)^i$. For each $j$ the survival of the sum of the first $j$ exponentials at log-time $s$ is
\[
 \sum_{i=1}^j e^{-h_i s}\prod_{a\ne i,\;a\le j}\frac{h_a}{h_a-h_i}.
\]
This follows by partial fractions or convolution of independent distinct-rate exponential densities. The checker bounds each exponential independently before subtracting endpoints. Negative lower bounds on cell masses are discarded. It pays the increasing factor of (1) at the left end and the decreasing factor at the right end. All unfinished cells beyond the integration cutoff contribute zero to the lower bound.

## Finite-size adverse tails

A pending leaf with parent degree $k$ has the locked upper $B_{m,r}(x)$, $r=1/k$. With $D=m-r$, $U\sim\operatorname{Beta}(1,D)$ and $s=1-r$, the beta integral and $\Gamma(D+r)/\Gamma(D)\le D^r$ give
\[
 q^{1-m}B_{m,r}(x)\le \frac{1}{\Gamma(r)}
       \mathbb E\big[(Dx/q+D U)^{-s}\big].
\]
Here $Dx/q\ge y$ for $y\ge r$. If $A$ is uniform on $[0,1]$, $DU$ has law $D(1-A^{1/D})$. This expression increases with $D$ since $1-(1+t)e^{-t}\ge0$, where $t=-\log A/D$. For every size in this result $D\ge20$. The 128 equal-probability rectangles therefore bound the right hand side above with $D=20$ inside the expectation and $Dx/q$ replaced by the lower edge $l$ of the threshold box. The gamma denominator is rounded down using an odd Taylor integral on $[0,1]$ and decreasing right-end rectangles on $[1,9]$. For $k\ge32$ the simpler bound $l^{1/k}/[(k-1)l]$ follows from $\Gamma(1/k)\ge k-1$.

## Signed whole-tree comparison

The source prices enter the signed tree ledger at the common completion size. All prices use the same factor $q^{m-1}$. For an old core of degree $d$, $c$ core neighbours, $L=d-c$ leaves and reciprocal receive $R$ from its core neighbours, put $s=(L+R)/[d(d+1)]$, $t_0=1/(2d)$ and $g=M_{\rm serial}-U_2$. Its own serial placements and those of its leaves, after paying all old-leaf edges, yield at least
\[
 C=g(ds+Lt_0)-L[(t_0-s)_+U_2+(s-t_0)_+U_{d+1}].
\]
Stars are checked by $C>0$. Contract each maximal chain of bare degree-two cores. With base $B$, census constant $a$, ordinary chain fee $f$ and terminal fee $f_t$, the retained-vertex inequality at $k$ chain ports is
\[
 C+a(c-2)-B\ge k f_d+(c-k)(s U_{d+1}-B)_+,\qquad
 f_d=f_t\ (d=2),\quad f_d=f\ (d\ge3).
\]
The balance is concave in $R$, so $R=k/2,c/2$ suffice. A direct edge spends at most one base and its outward endpoint's positive excess. A chain spends its two fees and one base. A tree with $V$ retained vertices has $V-1$ contracted edges, so the vertex base covers all these base charges and leaves one base in reserve. At its boundaries the checker maximizes the two possible directed differences in $s$; its inner edges pay $U_3|s-s'|$. Each chain vertex buys $2sg$. Lengths one through three suffice because adding a further homogeneous bare vertex buys positive credit without changing a boundary or an interior cost. Summation leaves $2a+B>0$ by the tree identity for the retained degrees. [Original chain accounting](../annealed/extension/third-physical/README.md)

For $d\ge32$ the own cap is at most $z=17/(16\cdot32l)$. Here $(17/16)^{32}>9/2$. All its own outward charges fit into $dsz$, with each direct port below the base. When $L\le d/2$, after paying this own charge the available stock credit is at least $\min(0,[32(g-z)/33+(g-U_2)/2]/2)$; in addition there are at least $16$ core ports paying $a-f$ apiece. When $L>d/2$ the credit is at least $32(g-z)/66+g/4-U_2/66$, since a positive leaf-side difference is at most $1/[2d(d+1)]$. Both bounds are compared with $2a+B$. A high-degree star receives positive credit because $g>z$ and $s\ge t_0$. These estimates cover the arbitrarily high degrees of the original tree.

The closed boxes cover $[23/10,9/2]$. The earlier symbolic comparison applies throughout $y\ge9/2$. The new common cutoff is sufficient. The sign at uncovered thresholds remains open. All calculations use rational inequalities on the physical upper and lower prices. [Finite-size checker](check.py) · [Signed ledger](ledger.py)

Run `python3 -IB round2/open/age-order/global-transport/beta-clocks/check.py` from the Astra repository. The five closed boxes have endpoints $23/10,5/2,11/4,3,7/2,9/2$. Each uses eight positive stock sources. The verifier evaluates every integer chain-port count and both receive corners below degree 32, all ordinary/terminal boundary types, and the unbounded-degree inequalities. The fixed rational allowances used in the five checks are recorded in `BOXES`.
