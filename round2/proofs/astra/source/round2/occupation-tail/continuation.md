# Occupation with a reciprocal-degree tail

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 9: Claude and HY4 Flash on specified components; GPT for full replay; inputs include the Round 1 physical exponent and prior, signed raw witnesses and disjoint occupation credits. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


The proposed additional nonnegative occupation credits, in units of $`2^{-40}`$ per newborn, are 17,458,365 for the upper witness and 15,289,028 for the lower witness. In combination with the nominal trial points and the three credits in [the first calculation](../occupation/inside.md), they give the additional candidate interval $`[0.5656581227,0.5657979941]`$. The occupation and source calculations in this extension await a cold review.

## Additional physical sources

Keep the same raw signed witness and raw error $`e`$ as in the first calculation. At a singleton the per-vertex raw cap for that family minus the maximized singleton row is a nonnegative gain. At an oriented edge with known ordinary endpoint states use the positive part of the negative of its maximized signed row. A donor outside the edge can also contribute the gap between the maximal and actual noncentered or piecewise-linear reciprocal credit. This third contribution attaches to its three-vertex support, and uses the same maximization domain as the edge row. The integer implementation uses the set of physically possible reciprocals. All three assignments are described in the [edge and singleton construction](../motif-cert/reserve/README.md).

For each triple and quartet with known states allow the positive part of the negative of its tightened signed row, leaving its unknown external rates maximized. A donor adjacent to a triple whose four-set receives no transfer for that slot contributes an additional nonnegative gap, attached to its four-set. The proof of these two gains enumerates the transferred lower-order terms at the actual endpoint states. [Triple and quartet construction](../motif-high/occupation/source.md)

These sources use distinct summands in the raw physical generator. Within each family the row slack is calculated with its outside inputs still maximized, and the actual donor pays a separate difference. The baseline five-vertex credit uses external donors of quartets. The baseline edge credit uses external donors of a *singleton* in its uncentered range. Neither is one of the additional sources listed here. Denote the sum of the additional physical sources by $`F`$. Then $`F`$ can be added to the right of (1) of the first calculation. In particular this calculation does not add a further stationary relaxation of the same singleton or edge rows.

## The interpolation and the lifetime test

Again use all chronological connected sets of noninitial cores up to an indicated arity, with an integer score at each table row. Choose a degree knot $`d`$ and a total core-degree ceiling $`m\le d`$ for each position. With total core degree $`c\le m`$ use the true state for $`D<d`$. For $`D\ge d`$ use a two-point mixture, with mass $`t=d/D`$ on $`(d,c)`$ and $`1-t`$ on $`(\infty,c)`$. The table value at a physical set is its expectation under independent choices of these labels, a multiaffine interpolation. At $`c>m`$ set the value to zero. The new enlarged-set table, reached on promotion of an old member, uses weakly smaller degree and core-degree cutoffs at each old position.

For old coordinates, the following corner rate intervals give a *lower* bound on the physical generator whenever each rate is multiplied by its actual signed successor change and minimized on the stated interval.

1. At an exact receiver, start with $`D-c`$ plus the reciprocals of its exact internal neighbours. Each of its unobserved core slots can add any rate in $`[0,1/2]`$. An internal neighbour with a tail label at knot $`d'`$ can add any rate in $`[0,1/d']`$ at either corner. The unknown physical rate is inside the interval in every instance. For a fixed physical configuration, the same gain rate multiplies the mixture of the signed jumps, so taking the smaller endpoint contribution at each corner is valid even when the neighbour is itself being interpolated. Its exact promotion rate is $`1-c/D`$.

2\. At a tail receiver, write $`H_d,H_\infty`$ for the two score values with all other labels fixed. A physical gain with rate $`a=D-c+\sum i_y`$ contributes

```math
a\,\Delta H=t\,\frac{a}{D+1}(H_\infty-H_d).
```

The fraction $`a/(D+1)`$ lies between $`(d-c)/(d+1)`$ and one. Put this interval on the transition from the finite tail label to infinity, with zero gain at infinity. At either tail corner the promotion rate can be taken in $`[1-c/d,1]`$. Since $`D`$ stays fixed at a promotion, the pre- and post-promotion interpolation masses are the same. This also covers a promotion that passes the core-degree ceiling.
3. To evaluate an enlarged-set reward at a coarse knot $`d'\le d`$, project an old finite label of degree $`D`$ to the appropriate mixture, using mass $`d'/D`$ if $`D\ge d'`$. Send infinity to infinity. Projection composes because $`(d/D)(d'/d)=d'/D`$. The newborn has $`(D,c)=(2,1)`$. The code integrates the product of up to four two-point projections and rounds the nonnegative enlarged-set contribution down once. If a coarse total-core-degree cap is exceeded, the enlarged-set score is zero.

For a source on an exact tuple, the code uses the eligible integer row or incidence. It drops a source whenever an interpolated position in its receiving edge, triple or quartet would be required. For a donor outside an edge it can use the minimum on $`[0,1/d']`$ at both corner labels. For a donor outside a triple the physical credit $`J(g,p)=\max(g,0)/2-gp`$ is affine in the reciprocal $`p`$ of that donor, once the three receiving members have specified states. Use its floor at $`1/d'`$ for the near corner and at the far corner use a floor valid throughout $`[0,1/(2^{31}-2)]`$, including zero. Their linear interpolation lies below the real physical credit, by affinity. Every possible donor is identified by deleting a leaf of the quartet. The remaining triple inherits its chronological order. A four-set with two tail members cannot collect even one such credit using this prescription. The proof uses the real raw credits before flooring, which allows interpolation of endpoint floors. Thus the physical $`F`$ dominates the interpolation of the sources. A configuration outside a core-degree cutoff has zero old score and cannot have a negative old-particle drift there. The gains and promotions of each old set, and the unique enlargement obtained by adjoining its youngest member, now give the same inequality (3) and the same birth account (4) as in the first calculation. Boundedness of the expected normalized sums follows from the core-degree moment estimate there.

In each corner row, for every admissible rate interval $`[a,b]`$ the integer checker includes the smaller of $`a(V'-V)`$ and $`b(V'-V)`$. The promotion destination includes the rounded-down enlarged-set score. The checker rounds each term, including the negative discount, downward to units of $`1/4096`$ of a table tick and requires the resulting sum to be nonnegative. All successors within an old table are ordered strictly after the current index, except for missing successors of value zero. No interpolation approximation is accepted without these one-sided row inequalities.

## Values calculated

The degree knots in age order, at the four arities, are $`(1400)`$, $`(72,56)`$, $`(56,24,16)`$ and $`(24,10,7,6)`$. The corresponding core-degree ceilings are $`(40)`$, $`(16,10)`$, $`(14,9,7)`$ and $`(8,5,4,3)`$. The numerator and denominator of the discount are 1133 and 2000. All five source groups are enabled. The calculation visits respectively 55,259, 510,900, 21,172,384 and 10,109,560 corner rows in each direction.

The new scores add to $`v`$ in (4) of the first calculation. All previous $`e`$ and $`y`$ remain the same. The resulting fractions are
```math
 \frac{565662790716854568223}{1000008252270910000000},\qquad
 \frac{1698046074520368376961}{3001152517940340000000}.
```
The corresponding raw computations and direct subtraction against physical initiator rings are recorded separately. [Upper rows](../outputs/occupation-tail-affine-upper.txt) · [Lower rows](../outputs/occupation-tail-affine-lower.txt) · [Physical reconstruction](physical.py)


Run `python3 -B round2/occupation-tail/replay.py` from the repository root. The command rebuilds the frozen snapshots and all original occupation inputs before running the new tables. The physical check also evaluates actual one-leaf and promotion transitions on 1,212 configurations, including degrees above one billion. Each sign exercises 52 tail gain identities on nonconstant computed table entries and rejects flattening the interpolation, shifting its numerator, and inflating the newborn score. [Replay](replay.py)
