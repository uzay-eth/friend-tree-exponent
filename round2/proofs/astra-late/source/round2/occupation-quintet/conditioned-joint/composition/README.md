> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Additional edge and three-core payments on the original five-vertex witnesses

The original Neon-2 five-vertex witnesses at b04cc50 give the two certified stationary trials
```math
 \frac{5657046097649}{10^{13}}\le\mu\le\frac{56574562817}{10^{11}}
 \quad=\quad[0.565704609764900,\;0.565745628170000].
```
These fractions use the original gold1 files. Neon-2 has also issued larger gold2 witnesses with different potentials. The calculation here quantifies the unused drift in gold1. [Original proof](../../../../../../../reviews/cold/astra/attacks/neon2-quintet/README.md) · [Previous conditioned trials](../README.md) · [Cold review of the composition](../../../../../../../reviews/cold/astra/checks/quintet-composition-sol.md)

## The signed inequality

Use the conditioned singleton and the complete numerator $`N`$ of release 204990e. Its denominator is the mean assigned potential $`M>0`$ in ticks per core. The two extra payments, tabulated in normalized potential units per core, are as follows.

| New stationary floor | Lower witness | Upper witness |
|---|---:|---:|
| Old donor with degree above five | $`7.35584318670\,10^{-7}`$ | $`1.39160762527\,10^{-7}`$ |
| Three-core path with a high-degree mother | $`5.88485196511\,10^{-7}`$ | $`6.56967517127\,10^{-7}`$ |

For the old-donor payment at 746c825, the host is an exact two-vertex-table state outside the three-core transfer, its degree exceeds 18, the paying neighbour's degree is at most five, and the donor's degree exceeds five at reward time. Each payment recovers one ordered Fenchel–Young slot in the original edge row. The previous small-neighbour and older-payer masks cannot contain that slot. A pure degree-18 edge cannot contain the paying edge. The larger lower-witness run replaces its smaller nested run. [Ordered donor proof and receipts](../../edges/old-donor/README.md) · [Independent old-donor check](../../../../../../../reviews/cold/astra/checks/quintet-old-donor-sol.md)

The H3 payment at 2705dac recovers part of $`(-P)_+`$ for the exact outward charged row $`P`$ of one unordered physical three-core path, allocated to its middle. The author and reviewer prove that $`\max(P,0)`$ satisfies both underlying higher-host caps, including at mixed tail labels. Thus replacing $`P`$ by $`P+(-P)_+`$ preserves the common census allowance $`h`$ used in $`N`$. The new edge payment occurs on a different term of the signed decomposition. The path baseline uses the same upstream reciprocal with siblings priced adversely. Where the sibling increment pays, the two contributions sum to at most that path's single negative part. The one-child direct floor is disjoint from the sibling increment. [Exact H3 components and scopes](../../higher/full-h3/RELEASE.md) · [Physical cap and sibling accounting](../../../checks/quintet-h3-cohort-sol.md)

Write $`S=2^{40}`$ for the original scale, $`C=S(R_{\rm donor}+R_{\rm H3})`$ for the sum in ticks per physical core, $`\pi`$ for the stationary core law and $`\eta=\pi/\mu`$. For either original trial $`\lambda`$ and inward sign $`s`$ ($`+1`$ lower, $`-1`$ upper), every inward displacement $`a\ge0`$ satisfies
```math
 [s(\mu-\lambda)-a]M/\mu
 \ \ge\ \eta[N+C-a u],
```
where $`u`$ is the original census upper potential. Each additional floor is a bound on an already integrated occupation, so $`C`$ enters once on every unit of core mass. The newborn rate still multiplies only the original maternal lifetime value. Auxiliary edge and flux variables receive zero cost. The production drift derives the inequality before substituting a previous trial bound.

Both exact repaired dual margins on 613 equations and 5,554 columns are positive. They exceed $`0.6771303167`$ and $`0.5386302076`$ ticks in the measure $`\eta`$. A separate check discards the positive slack of the old dual strictly above its certified trial and divides $`C`$ by the original stationary upper bound on $`M`$. Using the upper values from joint24.log gives the outward check $`[0.565704609691780,0.565745628223974]`$. With only the old-donor payment the same check gives $`[0.565704134004711,0.565746158933064]`$. [Exact weights and earlier stationary upper bounds](../../stationary/README.md)

Run the following commands from the Astra root with NumPy, SciPy, and the original private gold1 copies named in the parent proof. The first command uses the two saved duals and needs no optimization. The optional flag --solve regenerates them.

~~~sh
w=round2/.work/replay-inputs
p=round2/occupation-quintet/conditioned-joint/composition
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -IBu "$p/check.py" "$w"
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -IBu round2/checks/quintet-composition-sol.py --side lower
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -IBu round2/checks/quintet-composition-sol.py --side upper
~~~
