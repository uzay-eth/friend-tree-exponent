> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Stationary comparison and joint pricing of the Neon-2 quintet

The joint degree-24 calculation gives **`[0.565703060655196, 0.565746333629979]`** for the core growth exponent. Relative to the exact endpoints of Neon-2 Claim 6 it moves the lower end inward by more than `0.882741282e−6` and the upper end by more than `2.073807611e−6`. Astra's original independent replay of the underlying whole-tree certificate is at `2d76a12`. Its separate singleton lifetime calculation is at `32986a8`. The exact fractions for the joint result and the positive trial margins are in the [production receipt](joint24.log) and [independent audit](../../../../../../reviews/cold/astra/checks/quintet-joint-sol.md).

| Independent replay, outward endpoints | Lower | Upper |
|---|---:|---:|
| Original five-vertex certificate | 0.565702177913 | 0.565748407438 |
| Host comparison and higher-family cap | 0.565702368490575 | 0.565748281318382 |
| Joint calculation with the singleton reduced reward | 0.565702661777197 | 0.565747732733322 |
| Joint calculation with both edge payments | **0.565703060655196** | **0.565746333629979** |

## Physical accounting

Write all rewards and potential values in the original ticks, with scale $`S=2^{40}`$. Let $`\pi`$ denote the distribution seen from a uniformly sampled core, $`M=\pi\Phi`$ the mean potential assigned to it, and $`d=+1`$ on the lower witness and $`-1`$ on the upper. At its original nominal trial $`\lambda`$, the sign of interest is $`d(\mu-\lambda)`$. The original outward cap per core is $`E=C_f+C_g+H`$, equal to 19,186 and 19,783 ticks respectively. This includes the signed upper singleton cap $`C_f=-1053`$. An unborn edge enters the singleton once at promotion. Thereafter its contribution belongs to the one existing unoriented edge.

The host allocation puts the singleton at its vertex and an edge at its younger endpoint. H3, S3, S4, P5 and B5 go to their middle, hub or center. A P4 splits equally between its two internal vertices. Multiplicities count distinct neighbours even when the census labels coincide. For each host census $`(D,c)`$, `weights.py` bounds the allocated potential by a lower value $`l`$ after the original recharge and an upper value $`u`$ before recharge. Its higher-family cap $`h`$ is the smaller of $`H`$ and the sum of the positive family caps times the actual eligible host multiplicities. Endpoint interpolation, the degree-1,999/2,000 seam and all larger degrees enter these bounds. The identity $`\sum(c-2)=-2`$ on each finite core tree removes the per-core recharge in the stationary mean. Thus $`\pi l\le M\le\pi u`$. The finite reduction visits all 1,998,999 nonfounder censuses per sign and matches the original finite comparison extrema. The [host audit](../../../../../../reviews/cold/astra/checks/quintet-stationary-sol.md) gives the individual host inequalities and the founder convention.

Let $`F`$ be the original unused singleton drift, $`V`$ its nonnegative lifetime subsolution, $`J=(2,1)`$ the actual newborn state, and $`[L,U]`$ the original certified interval. With $`A`$ the physical generator on the living focal core, the full residual is $`R=F+(A-U)V\ge0`$. Stationary balance including the newborn source gives
```math
\pi F=\mu V(J)+(U-\mu)\pi V+\pi R\ \ge\ L V(J)+\pi R_{24}.
```
The masked residual uses its exact affine intercept and slope at each admitted $`(D,c)`$. Its argument is the sum of the physical receives from children and mother. The same two receives appear in the gain rate and its edge marginals. The LP degree-overflow label 25 represents every $`D\ge25`$ and receives zero masked reward. Unresolved donors still enter the receive constraints. The original subsolution covers the unlimited degree continuation beyond its degree-12,000 table. [Singleton export and identity](../lifetime/dual/README.md) · [Separate exact audit](../../../../../../reviews/cold/astra/checks/quintet-dual-sol.md)

Let $`B_g`$ be the certified unused charge on edges with both endpoint degrees at most 18. Let $`B_{FY}`$ be the unused Fenchel–Young term on an ordered pair of neighbours for which the paying edge has at least one endpoint of degree greater than 18. This second term compares the maximum over the original knots to the value at the donor's actual reciprocal. Both registers start at the physical promotion that completes the incidence. Only the edge terms own these two credits. The higher-family bound retains its original transfers, and the present calculation takes no separate reward for a negative higher-motif occurrence. [Existing-edge proof](../edges/README.md) · [Outside-edge proof and masks](../edges/fy/README.md) · [Independent audit](../../../../../../reviews/cold/astra/checks/quintet-fy-sol.md)

Putting $`B=L V(J)+B_g+B_{FY}`$, these identities yield the single stationary inequality
```math
 d(\mu-\lambda)M\ge B-E+\pi\{H-h+R_{24}\}.
```
The original positive global comparison bounds enclose $`M`$. The sided LP also prices $`l-M_{\min}`$ and $`M_{\max}-u`$, giving refined bounds $`M_-\le M\le M_+`$. Bounds in the LP are per unit of $`\pi/\mu`$. Each cost is nonnegative, so multiplying a nonnegative repaired bound by $`L`$ gives a per-core lower bound. The twelve additional rows use physical birth probabilities computed with an upper discount at least $`U`$. Both the birth cutoff and its degree/arity saturation are kept in the replay.

For a proposed inward displacement $`a\ge0`$, use one column objective
```math
 H-h+R_{24}+a(M_{\max}-u).
```
Adding its lower bound to $`B-E-aM_{\max}`$ gives a signed margin $`m\le[d(\mu-\lambda)-a]M`$. For negative $`a`$, the comparison uses $`(-a)(l-M_{\min})`$ and the constant $`-aM_{\min}`$. Once $`a`$ is certified, its improvement is at least $`\max(0,m)/M_+`$. In the full run both proposals were positive. The independent check also evaluated each printed outward endpoint directly and obtained margins exceeding 327 and 364 ticks respectively. Both denominators and all repaired bounds are rational. HiGHS supplies candidate multipliers. Repricing and checking every original column yields the certificate.

## Replay

Start at the Astra project root. Place the authenticated copies `quintet_lower.bin`, `quintet_upper.bin` and `gold.json` from the [original certificate replay](../../../../../../reviews/cold/astra/attacks/neon2-quintet/README.md) in `round2/.work/replay-inputs`. The two original binary digests and the golden receipt digest are fixed in `weights.py`. Use Python with NumPy and SciPy. Both raw payloads are opened read-only and authenticated. `joint.py` also authenticates the singleton exports, edge receipts and physical LP inputs. The enclosed cutoff uses all 5,554 columns and 613 equality rows. The host-only calculation uses 5,542 columns and 601 rows.

```sh
p=round2/occupation-quintet/stationary
w=round2/.work/replay-inputs
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 16 python3 -IBu "$p/joint.py" "$w" --limit 24
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 16 python3 -IBu "$p/joint.py" "$w" --limit 24 --without-edges
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 16 python3 -IBu "$p/check.py" "$w" --limit 24
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 16 python3 -IBu round2/checks/quintet-stationary-sol.py --log "$p/host24.log"
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 16 python3 -IBu round2/checks/quintet-joint-sol.py
```

The [host receipt](host24.log) records the zero-extra calculation and its 388 labelled host comparisons per sign. The [joint receipt](joint24.log) records both exact endpoint fractions and each intermediate margin. The joint cold review uses the two frozen edge receipts. The separate degree-18 and Fenchel–Young solvers and their clean replay are documented in the edge package.
