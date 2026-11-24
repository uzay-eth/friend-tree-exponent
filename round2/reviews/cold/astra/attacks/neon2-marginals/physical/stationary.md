# Independent physical marginal constraints

Write `D=L+c` for ordinary degree, leaf stock and core degree; let `rho=sum_{u~core o}1/D_u`, `a=L+rho` and `b=L/D`. The physical generator of the existing focal core sends `(D,c,L)` to `(D+1,c,L+1)` at `a`, and to `(D,c+1,L-1)` at `b`. The second event also adds a newborn core `(2,1,1)` to the census. Accordingly `mu*pi(f)=pi(Af)+mu*f(2,1,1)` for bounded mark tests of integrable generator. Below, no witness reward is used.

## A leaf-stock cut linked to the arity and degree cuts

For integer `j>=0`, let `s_j=pi(L=j)`, `u_j=pi(b; L=j)` and `r_j=pi(rho; L=j)`. The cut `1[L<=j]` gives

```
mu*(pi(L<=j)-1[j>=1]) = u_(j+1)-j*s_j-r_j.
```

In particular `r_0+mu*s_0=u_1`, where the incoming rate is exactly `pi(1/D; L=1)`. In a lifted `(D,c)` census, every `u_j` is fixed by the degrees, while `r_j` is a projection of the same reciprocal edge flux used in the degree-crossing rows. Even if a table only exposes a few marginal intervals, the `j=0` row correlates the bare mass with the stocked mass and the bare-edge reciprocal flux, without estimating rewards, assuming independence or truncating a neighbour reciprocal to zero. At a trial `mu`, these are linear feasibility rows; a certified infeasibility could rule out that trial. Joint mark-level balances already imply these rows with complete aggregation, so their value is in strengthening coarse marginal models, especially overflow bins. Here `rho<=c/2` and `pi(c)=2` by sending one unit to the unique parent, so the cut has integrable crossing intensity even at unbounded arity.

Another exact consistency identity is `pi(rho)=pi(c/D)=1-mu`, by reversal of oriented edges and the birth-rate identity. It couples all of the `r_j`; no per-bin symmetry of conditional means follows. Use `mu*pi(c>=k+1)=pi(b; c=k)` for `k>=1` and `mu*pi(D>=k+1)=pi(a; D=k)` for `k>=2`. These cuts must share the same `(D,c)` and reciprocal edge census.

## An oriented degree-two edge cut

Write `J=(2,1)`, `B=(2,2)`, and `X=1[root=J, parent=B]` with the indicator rooted at the child. A physical `B→J` edge is created at rate `1/2` by every `J` promotion, at which the promoter becomes `B` and its child is newborn `J`. No external gain can create `B`, and a `J` parent (one core neighbour already accounted for by its own parent) cannot have an existing child. The exit rates for the edge are: child gain `1+1/2`, child promotion `1/2`, parent gain `1/2+1/D_grandparent`; only distinct arrows are counted. The stationary equation is therefore

```
pi(J)/2 = (mu+5/2)*pi(X) + pi(X/D_grandparent).
```

This is a bounded test of the physical spatial generator; the entrance occurs at the new root, not as an extra arrival into an existing `J` row. It gives `pi(J)/(2*(mu+3)) <= pi(X) <= pi(J)/(2*(mu+5/2))` without any unobserved moment enclosure. Also `pi(1/D_parent;J)` includes exactly `pi(X)/2` from this edge type, since a core of ordinary degree two with a core child necessarily has label `B`. For a feasibility model, keep the remaining parent contributions in bins `D_parent>=3` and the grandparent reciprocal in a three-site census. This oriented flow can constrain exponent candidates jointly with an independently bounded incidence/reciprocal column; plugging in just the published marginal boxes does not certify exclusion. A stationary model that already includes exact oriented edge balances gets this equation by projection.

## Downward size: a boundary flux with sign, not just a box

Let `N` be the ordinary size of the finite downward branch of the focal core, excluding its parent. A newborn has `N=2`; its exact unit-increase rate is `N+p_parent-p_root` with `p=1/D`: all internal initiations count except the root's outgoing arrow to its parent, and the latter's arrow into the root is the only exterior input. Thus at every `n>=2`,

```
mu*pi(N>=n+1)=n*pi(N=n)+pi(p_parent-p_root; N=n).
```

This cut is bounded in jump intensity. The support satisfies `N>=D+c-1`: the root, its `D-c` leaves and `c-1` child-core subtrees each of size at least two. Under the transport moving the root's `N=n` mark to the parent, `pi(p_parent;N=n)=pi(p_root * #{children w:N_w=n})`. Also `N=1+L+sum_{w child}N_w` pathwise. These can be used jointly without claiming independence between parent and child. For example the child count with size in any common threshold set `[n,infinity)` is at most `floor((N-1-L)/n)` and is zero at `N<=n`; at small `N` retain the exact integer partitions of `N-1-L`. For large `N` the count bound may be capped by `c-1` to keep it integrable. This tightens a model that only allows an unrestricted interval for the parent term and imposes no cost derived from the witness. A numerical improvement requires a joint census or a dual infeasibility proof; none is asserted here.
