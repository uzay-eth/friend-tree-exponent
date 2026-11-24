> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Stationary joint occupations of core patterns

A mother and a core child are endpoints of the same physical edge. Their degrees control different arrow rates along that edge. This calculation constrains both rates with one joint occurrence measure, and extends the same constraint to a path, fork or star and any of its exterior neighbours.

On the six original D18 higher families with **every motif degree at most three**, the stationary sum of lower-witness charge minus upper-witness charge has the following outward enclosures. Each entry is per uniformly sampled core. Both computations use the same generator balances, witness charges, discount interval and objective.

| Exterior neighbour information | Lower bound | Upper bound |
|---|---:|---:|
| Individual reciprocal in [0, 1/2] | −4.100739 × 10⁻⁶ | 4.136497 × 10⁻⁶ |
| Common parent and child occurrence measures through five cores | −2.892607 × 10⁻⁶ | 2.698686 × 10⁻⁶ |

Sharing the exterior neighbours narrows this interval by 32.1%. The unpriced degree classes and the signed potential term still enter the exponent test. These two rows give a comparison on identical support. The independent referee checked 418,720 generator rows and the projection on 1,117 physical ports. [Audit](../../../../../../reviews/cold/astra/checks/embedding-flow-sol.md)

## Variables and balances

Write $`\pi`$ for the oriented stationary uniform-core probability and use $`\eta=\pi/\mu`$. At a nonfounder core with original degree $`d`$ and $`k`$ core children, the number of leaves is $`l=d-k-1`$. Promotions occur at $`b=l/d`$. Gains occur at $`g=l+1/d_{\rm mother}+\sum_{y\text{ child}}1/d_y`$. The founder correction disappears in the stationary per-core limit. In a finite tree the mother of the selected apex can be the founder, which needs a separate endpoint correction.

For each finite labelled rooted pattern with specified $`(d,k)`$ at every role, $`m(P)`$ is its expected number of injective embeddings per unit of $`\eta`$, counted from the apex. Canonical ordering of child subtrees identifies equal expectations. A pattern with identical sibling subtrees still counts the choices for their distinct slots. To price an original unordered physical motif, divide the labelled count by the automorphism number of the coloured pattern. For each census, balance incoming gains, promotions and the creation of an embedding at a virgin tip against its outgoing rates and the discount. The one-vertex virgin has a birth source of one. If several labelled virgin tips can be the newborn, each legal deletion appears once in the incoming sum.

At a port with $`n`$ available neighbours, write $`z`$ for their aggregate reciprocal-weighted occupation. For one specified further neighbour $`y`$ whose degree is at most $`J`$, let $`m_y`$ be its extended labelled pattern mass. For the exterior mother $`n=1`$ and she becomes the apex of the extension. Mass transport counts the same incidence at either endpoint. The exact physical inequalities are
```math
 \sum_{d_y\le J}\frac{m_y}{d_y}\ \le z\ \le
 \sum_{d_y\le J}\frac{m_y}{d_y}+
 \frac{n\,m(P)-\sum_{d_y\le J}m_y}{J+1}.
```
Both row directions and all admissible deletions refer to the same extension columns. This permits interaction between different geometric patterns without assuming an independent donor. The full physical measure satisfies the bounds, including at the degree cutoff.

Each exact balance uses the common $`\mu`$. In one finite relaxation its term $`\mu m`$ is enclosed by splitting $`m`$ between the two prior endpoints. The total measure is in $`[1/\mu_{\rm hi},1/\mu_{\rm lo}]`$. The remaining scalar constraints use $`\pi k=1`$, $`\pi b=\mu`$ and the degree tail from the stable-age comparison. Rational dual repair gives certified extrema on the physical law. For each fixed labelled topology with $`r`$ roles of degrees at most $`D`$, the sum of all its decorated masses is at most $`(D-1)^{r-1}/\mu_{\rm lo}`$. Each port group has a further factor at most $`rD/2`$. Every reduced cost is checked and its group's largest shortfall is charged to this total. The signs are converted back from $`\eta`$ to $`\pi`$ with the same prior endpoints.

## Replay and scope

From the repository root, use the original two read-only D18 files and a Python environment with NumPy and SciPy.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -EBu round2/occupation-neon/embedding-flow/solve.py --lower "$D18_LOWER" --upper "$D18_UPPER" --cap 3 --share-through 0
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -EBu round2/occupation-neon/embedding-flow/solve.py --lower "$D18_LOWER" --upper "$D18_UPPER" --cap 3 --share-through 4
OPENBLAS_NUM_THREADS=1 python3 -EBu round2/checks/embedding-flow-sol.py "$D18_LOWER" "$D18_UPPER"
```

The historical `--unpaid` switch filters the two earlier complete-fringe size reserves. Later four-core path, bent-branch and live three-core donor payments can occur on its support. A new positive payment requires subtracting the relevant earlier charge or proving disjointness at the level of each physical incidence.

For convergence as the support grows, the [identification theorem](identification/README.md) imposes incidence-mass tails at every retained threshold and rational time-one cylinder inequalities on narrowing discount intervals. Its reconstruction orders child slots and then identifies the orientation by the unique end of the stationary tree. The finite calculations in the table use the displayed generator and reciprocal constraints. The theorem has a [separate cold audit](../../../../../../reviews/cold/astra/checks/embedding-identification-sol.md).
