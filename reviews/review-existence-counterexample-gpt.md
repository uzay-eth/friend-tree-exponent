VERDICT: COUNTEREXAMPLE to the unrestricted claim that the mean core count of a planted branch lies below the free mean at an internally non-leaf root.

# A planted branch with larger mean core count

The root $`v`$ of the free tree has two neighbours $`a,b`$. Each has nine leaf neighbours in addition to $`v`$. In the planted tree, join an exterior parent $`u`$ to $`v`$ and give $`u`$ three exterior leaves. The free tree has 21 vertices, and $`D_u=4`$ in the planted tree. Count vertices on the $`v`$ side of $`uv`$ using their degrees in the full tree. Count all non-leaves in the free process. Both counts start at three.

This configuration satisfies the planting hypotheses with $`M=2`$. Start at the centre $`u`$ of a three-vertex star, keep its two leaves through the first stage, add exactly two leaves at $`u`$ in the second, and promote one in the third. The promoted vertex is $`v`$, and the newborn is its first internal neighbour. Grow the remaining rooted tree by the prescribed internal births while the exterior stays silent. Every restriction fixes finitely many allowed choices or Poisson counts on bounded intervals and has positive probability. See [the planting lemma](../proofs/existence.md).

## The first two derivatives

Write $`K`$ for the count and $`S=\mathcal L K`$ for its instantaneous drift. An ordered move from $`x`$ to a neighbour $`y`$ has rate $`1/D_x`$ and adds a leaf at $`y`$. The root is already core and every inside leaf has an inside parent. Thus $`S`$ sums the fraction of leaf neighbours over the inside cores. At time zero both drifts are $`9/10+9/10=9/5`$.

| Change at the first birth | Effect on the difference of the drifts of $`S`$ |
|---|---|
| A new leaf at $`v`$ | The free rate is $`2/10`$ and the gain in $`S`$ is $`1/3`$. The planted rate is $`2/10+1/4`$ and the gain is $`1/4`$. |
| A new leaf at $`a`$ or $`b`$ | Each rate changes by $`1/3-1/2`$. The gain in $`S`$ is $`10/11-9/10=1/110`$ in both trees. |
| Promotion of an inside leaf | Its parent has degree ten in both trees. The promotion rate and the drift increment agree. |
| Attachment at an exterior vertex | The inside drift stays fixed. |

Subtracting the free curvature from the planted curvature gives

```math
(\mathcal L^2 K)_{\rm planted}-(\mathcal L^2 K)_{\rm free}
=\frac{2/10+1/4}{4}-\frac{2/10}{3}
 +2\left(\frac13-\frac12\right)\frac{1}{110}
=\frac{113}{2640}.
```

For the Taylor remainder, stop each process after its second birth. This gives a finite state space and preserves the first two generator derivatives at the initial state. If $`B_s`$ is the unrestricted number of births by time $`s`$, the error in the mean is at most $`\mathbb E[B_s\mathbf 1_{B_s\ge3}]`$. For $`b\ge3`$, $`b\le b(b-1)(b-2)/2`$. Starting from $`m`$ vertices, the Yule birth count is negative binomial and has third factorial moment $`m(m+1)(m+2)(e^s-1)^3`$. The error is therefore $`O(s^3)`$ in both trees. With the common count and slope cancelled,

```math
\mathbb E K_{\rm planted}(s)-\mathbb E K_{\rm free}(s)
=\frac{113}{5280}s^2+O(s^3)>0
```

for every sufficiently small $`s>0`$. Lemma 4.1(c) of the [mean-exponent proof](../proofs/existence.md) bounds the planted mean below by a weighted free expectation. Proposition 6.2 uses precisely that lower bound, which is compatible with this example.

## Exact reproduction

The [Python checker](../verify/tree-drift/existence_counterexample.py) evaluates the count after every ordered move and then applies the generator a second time using `Fraction`. It decides whether a newborn is inside from the receiver of the move. Running from the repository root gave

```text
$ python3 verify/tree-drift/existence_counterexample.py
free count, slope, curvature: 3 9/5 1583/1650
actual count, slope, curvature: 3 9/5 13229/13200
actual minus free curvature: 113/2640
```

As a negative control, I changed the inside-set update in memory to exclude newborns attached at vertex zero (the root). The curvature difference became $`-1/330`$ and the script's exact assertion failed. This checks sensitivity to how births at the boundary root are classified.
