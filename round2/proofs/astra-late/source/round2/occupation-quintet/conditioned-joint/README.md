> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# A stationary exponent bound using the mother's degree at promotion

The original Neon-2 five-vertex witnesses give the enclosure
```math
 \boxed{\frac{2828517697069}{5000000000000}\le\mu\le
             \frac{2828731356747}{5000000000000}}
 \quad=\quad[0.565703539413800,\ 0.565746271349400].
```
The lower endpoint rises by $`4.78758604\cdot10^{-7}`$ and the upper endpoint falls by $`6.2280579\cdot10^{-8}`$ relative to the outward cut-24 release at `d08d9f8`. The two displayed fractions are exact certified trials. [Original five-vertex review](../../../../../../reviews/cold/astra/attacks/neon2-quintet/README.md) · [Previous stationary bounds](../stationary/README.md) · [Independent check of both trials](../../../../../../reviews/cold/astra/checks/quintet-palm-joint-sol.md)

## Counting newborns at their actual mothers

Let $`\pi`$ be the stationary law viewed from a uniformly sampled core, $`\eta=\pi/\mu`$, and $`M>0`$ the mean assigned potential in the original integer ticks. Use $`\sigma=1`$ for the lower witness and $`-1`$ for the upper witness at its nominal trial $`\lambda`$. The potential scale is $`S=2^{40}`$. A living nonfounder of total degree $`D`$ and $`k`$ core children promotes a new core at rate $`b=(D-k-1)/D`$. Its own degree is the same immediately before and after that promotion. The newborn sees this vertex as its mother.

Let $`V_2,V_3`$ be the certified discounted lifetime values of the original singleton slack, conditional respectively on maternal birth degree 2 and at least 3. The lower values in both witnesses come from the frozen degree-360 conditioned Bellman computation at `fd63517`. That computation tracks the actual mother and latest child while integrating the focal receive. Its continuation covers all degrees. Stable age accounting gives $`\pi F\ge\pi[bV(D)]`$ for the singleton slack $`F`$, where $`V(D)=V_2`$ for $`D=2`$ and $`V_3`$ otherwise. Indeed the maternal-degree birth measure in $`\eta`$ units has mass $`\eta[b\,\mathbf1_{\{D\in A\}}]`$, and $`\eta b=1`$. Each newborn contributes its conditional lifetime value once. A degree-overflow column still pays $`bV_3`$. The existing virgin age balance implies the bounds on the degree-two maternal birth share used in the standalone calculation. [Conditioned lifetime proof](../lifetime/conditioned/README.md) · [Cold check of the four bounds](../../../../../../reviews/cold/astra/checks/quintet-conditioned-sol.md) · [Exact entrance values](roots.log)

Let $`E=C_f+C_g+H`$ be the original per-core outward drift cap. Its values are 19,186 and 19,783 ticks, respectively, with the signed upper singleton cap $`C_f=-1053`$. The host bounds $`h(D,k+1)`$ on the higher-family allowance and $`u(D,k+1)`$ on the upper assigned potential are the original stationary release's physical bounds. Their means satisfy $`\pi u\ge M`$. Write $`B`$ for the sum of the following stationary edge floors in ticks per core.

- Both witnesses use the original edge floor when the two endpoint degrees are at most 18, plus the frozen Fenchel–Young floor on larger paying edges with paying neighbour and donor degrees at most five. [Edge occurrence](../edges/README.md) · [Original Fenchel–Young floor](../edges/fy/README.md)
- The lower witness also uses the reviewed `3a06c62` floor with paying-neighbour degree above five, donor degree at most five, and at least one paying-edge endpoint above 18. All three tests use degrees at the instant of payment. The source receipt is in normalized potential units and is multiplied once by $`2^{40}`$ to obtain ticks. [Directed edge cut](../edges/fy/extended/README.md) · [Cold audit](../../../../../../reviews/cold/astra/checks/quintet-fy-extended-sol.md)

For a proposed inward displacement $`a\ge0`$, define
```math
 N=B-E+H-h+bV(D).
```
The signed drift and comparison bounds imply
```math
 [\sigma(\mu-\lambda)-a]\,M/\mu\ \ge\ \eta[N-a u].
```
The constant $`B-E`$ occurs once per unit of core mass. Only the newborn lifetime term is multiplied by the promotion rate. Division by $`\mu`$ puts all terms on the same occupation measure. This calculation spends the conditioned singleton source once, with zero further singleton residual. The higher-family allowance and the potential comparison are evaluated on each original census class. The complete five-vertex drift certificate and its proof at every degree are inherited from the independent original audit.

## Signed trials and reproduction

The inputs are the two original binary copies and `gold.json` described in the earlier stationary release. The checker authenticates these files and the original birth law, host code, entrance values, and receipts. The 613 equality rows and 5,554 columns are unchanged. A fractional linear program suggests a trial by scaling the same vector $`\eta`$ by the positive integral of its upper assigned potential. The numerical optimum is used only to choose a trial. At each fixed trial, the stored rational multipliers are priced against **every** original column of $`N-au`$, including zero-cost auxiliary columns. For each column group the largest positive defect is charged against its known upper mass $`1/L`$, where $`L`$ is the original exact lower exponent bound. A signed lower bound survives this repair in its original units. For the printed trials the exact positive margins in $`\eta`$-ticks exceed $`0.750302733676`$ and $`0.559140454944`$, respectively.

From the Astra repository root use a Python environment with NumPy and SciPy. The default check reads the stored dual certificates and performs no optimization.

```sh
w=round2/.work/replay-inputs
p=round2/occupation-quintet/conditioned-joint
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -IBu "$p/check.py" "$w"
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -IBu round2/checks/quintet-palm-joint-sol.py --side lower --stored
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -IBu round2/checks/quintet-palm-joint-sol.py --side upper --stored
```

To regenerate the two entrance rows, run `python3 -IBu round2/occupation-quintet/lifetime/conditioned/replay.py --bound lower` and compare stdout with `roots.log`. Running the main checker with `--solve` writes fresh rational dual certificates. Exact verification applies to every resulting certificate regardless of which numerical basis HiGHS selects. The independent checker rebuilds the physical birth kernel and census bounds and reads both saved multipliers. [Independent verification](../../../../../../reviews/cold/astra/checks/quintet-palm-joint-sol.md)
