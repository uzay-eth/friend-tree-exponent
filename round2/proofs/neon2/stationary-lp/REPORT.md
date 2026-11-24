# Physical stationary occupation LP pilot

**Pinned claim commit:** `859148a2d84895eb8ec142a6d7d73865c4c1653b`.  
**Review and dependencies:** No separate rival verdict for this pilot. Uses the physical stationary balance; the exact small witnesses are rational LP controls, not bounds on the exponent.

**Author:** Neon-2 lane (Gemini 4 Neon and GPT-6.2 Sol Ultrafast).

**No improvement on the supplied prior** `[0.5650, 0.5665]`. Both endpoints are numerically feasible for Level 1 at cutoff 8 and Level 2 at `(N,D3)=(5,4)`. A feasible endpoint already prevents excluding any interval touching it: the same occupation vector remains feasible on that larger interval. No in-prior Farkas certificate was produced. Larger solves below are incomplete where marked. The two shipped certificates are out-of-prior positive controls, not new bounds within the prior.

## Model and scope

- `N` is the largest resolved total degree; the tail starts at `N+1`. Level 1 has singleton masses and parent-to-child edge counts. A low-degree endpoint has its exact reciprocal. Each tail endpoint has a Bernstein mixture over `0` and `U=1/(N+1)`; the product lift on edges preserves both reciprocal means and their product. This is a lift of physical observations to nonnegative weights, not an assertion that the process has degree infinity or integer-reciprocal extrema. Core degree is exact for `c<=N`; the overflow `c>=N+1` is a single bin.
- Aggregating incoming edges gives one per singleton. Outgoing gives `(c-1)` per singleton except in overflow. There, the mass, reciprocal and complementary-reciprocal inequalities follow by multiplying `k>=N` respectively by `1`, `x` and `U-x`. The sum of outgoing masses is one. Promotion on tail bins is represented exactly in aggregate as `p-E[cx]`, including in overflow: sum the reciprocal of the promoter over its `k` outgoing and one incoming edges. Promotion at an overflow endpoint does not change its bin; a new edge still enters from that promotion. There is no discarded high-degree/high-arity mass.
- Edge generators use gains and promotions at both old endpoints, and add the new edge with the **post-promotion** parent state. Internal donor reciprocals are exact expectations. For a finite-degree endpoint the omitted donor totals are variables bounded by the number of omitted neighbors divided by two, and their sums over edge types match the corresponding counted incidence transport. For fixed parent `c`, the total parent-side ancestor credit is `(c-1)` times its incoming donor credit; sibling credit is `(c-2)` times its outgoing donor credit; for a fixed child state the child-side external total is its full outgoing donor credit. A tail degree gain cannot change any bin.
- For singleton tail reciprocal `H=1[D>N]/D`, the gain decrement is `(x-cx²+ρx²)/(1+x)`. For every tail `D,c`, it lies between `(x-Ucx)/(1+U)` and `x`, using `0<=ρ<=c/2` and `cx<=1`. Add the exact degree-`N` entry at reciprocal `U`. This supplies valid two-sided generator rows without assigning an asymptotic tail constant. Singleton balances, `π(b)`, and all resolved crossings `k=2,...,N` also have their separate rows.
- Level 2 has chains and unordered sibling pairs (the two child incidences recover the ordered-wedge marginal). For each resolved edge, exact motif incidence plus an explicit omitted count equals the corresponding multiple of edge mass. In these same identities each resolved third vertex contributes exactly `1/D`; the omitted credit lies between zero and its count divided by `D3+1`. The motif generator uses all member changes and the new motif from promotion on the corresponding old edge. External credit at a motif member is a nonnegative moment bounded by half the number of available slots times motif mass. This is a genuine four-vertex relaxation; no four-vertex stationarity or closure is imposed. Nor does Level 2 impose independent second reciprocal moments.
- The directional identities require a parent at each core and a child count of `c-1`. Here that is the supplied directional premise. `operator.md` constructs an unoriented annealed eigenmeasure and does not by itself assert uniqueness or fixed-time convergence. The occupation derivation only needs a stationary unimodular lift with such an orientation, not uniqueness. For the finite-seed construction in that document, retaining the finite birth orientation in its local weak limit is a natural route: the root-uniform seed frequency vanishes, and by marked-ball tightness so does hitting that seed within any fixed radius. The LP is not an independent proof of convergence or of the one-ended claim.

For integrability of the stationarity rows, an edge bin with a high-degree endpoint uses a test constant in that endpoint's actual degree. Its unbounded gain rate does not appear in the drift, and promotions are bounded by one. Every gain that changes an occupation indicator acts at a bounded-degree member and has rate at most that degree; the counts rooted at a parent can alternatively be rooted at the child for edges and at the last vertex for chains by mass transport. A wedge rooted at its parent is bounded by the motif cutoff and gains/promotions produce only bounded crossing counts. The reciprocal tail test has a bounded decrement intensity, as shown above. These observations allow the occupation drifts to be obtained by compensating crossings, the same argument used for the degree crossings in `degree-law.md`.

## Empirical diagnostic

`sim_core.cpp` uses the physical process on the full tree: total event rate equals the number of vertices, each event selects a uniform ringing vertex and adds a leaf at a uniformly chosen neighbor. `sim_pi_check.py` weights trees by their core sizes, as in the annealed law, and reports tree-cluster errors.

2,000 trees at `t=12`, seed `38129`: **2,279,602 cores**. Counts of nonzero singleton/edge/chain/ordered-wedge keys: `36,763 / 525,253 / 1,464,207 / 3,193,474`. `π̂(b)=0.565703716` (SE `0.000175`), `π̂(c)=1.998245308`, `π̂(1/D)=0.295127823` (SE `0.0000935`), `π̂(1/D²)=0.121896593` (SE `0.0000567`). Nine crossing residuals (`k=2..10`) have maximum absolute value `0.001870049` (at `k=10`, SE `0.000685`). Wall `163.29s`, peak `794,652 KiB`.

`check_sample.py` projects the same trees onto **every** one of the `6,849` Level-1 rows at `N=8` and `39,074` Level-2 rows at `(8,6)` (wall `63.22s`, peak `234,048 KiB`); see `audit.log` or replay. Physical pointwise donor caps, all motif marginal/credit rows, both tail reciprocal inequalities and the resolved crossings behave as expected. Finite-seed deviations are visible: total edge mass is exactly `1-2000/2279602`; the seed also causes incoming/outgoing bin defects. One edge stationarity row and 11 motif upper stationarity rows exceed 5 batch SE, largely sparse cells for which a normal approximation is unsuitable. This diagnostic is not a claim that every finite-time stationarity residual vanishes or that fixed-time frequencies are proven to converge.

## Solver results

| Level | `(N,D3)` | Variables / inequalities | Endpoint result | Wall s / RSS KiB |
|---|---:|---:|---|---:|
| 1 | `(8,0)` | `5,540 / 6,849` | Both feasible; row violation at most `9.77e-9` | `18.76 / 128,220` |
| 1 | `(12,0)` | `25,398 / 31,139` | Timed out before first endpoint | `70.08 / 297,292` |
| 1 | `(16,0)` | `76,360 / 93,661` | Timed out before first endpoint | `70.11 / 488,772` |
| 2 | `(5,4)` | `2,603 / 3,708` | Both feasible; row violation at most `7.6e-15` | `14.52 / 101,764` |
| 2 | `(8,6)` | `28,990 / 39,074` | Extended IPM run and PDLP both timed out before first endpoint | `220.03 / 321,680`; PDLP `45.03 / 186,496` |
| 2 | `(12,8)` | `172,937 / 226,845` | Timed out before first endpoint | `70.13 / 448,920` |
| 2 | `(16,10)` | — | Timed out during model construction | `70.19 / 932,304` |

In the completed Level-2 comparison the two endpoint primal vectors differ by `0.0215603` in singleton L1, `0.0120511` in tail singleton L1, `0.584402` in Bernstein edge L1, `0.0754782` in external edge donor L1, and `0.0102383` in external motif donor L1 (sum of coordinatewise maxima `0.0173995`). The biggest singleton shift is `0.00624018` in tail `c=3`. These are arbitrary feasible vectors, so the differences are evidence that the unresolved tail and external donor allocations have freedom, not a ranking of which missing constraints would be binding. No in-prior infeasible dual was available to make such a ranking. A concrete next tightening is to match the external motif-donor sums to existing edge credit sums after removing the explicitly resolved motif credits; this would share budgets across motifs instead of letting each motif spend the full slot cap independently.

## Exact checker and replay

For `Ax<=r`, `0<=x<=U`, a nonnegative rational vector `y` certifies infeasibility if `yr < Σ min(0,(yA)_j) U_j`. `verify_farkas.py` checks this with `Fraction` and no numerical solver or simulator, rebuilding the exact matrix using the LP's model-definition code. It is an independent arithmetic implementation, not a second implementation of the physical row formulas. Bounds used to absorb rounding error are analytic: singleton and total edge masses are at most 1; a chain is a subcount of edges (one possible ancestor per edge); a resolved sibling motif is bounded by `(D3-2)` times the total edge mass; donor and omitted counts follow from their corresponding finite-slot inequalities. Thus all columns admit finite caps even with unbounded total degree or arity.

Both control files exclude `[1/10,1/5]` with 25 multipliers, `yr - lower_bound ≈ -0.0453535894351` (use `--exact` to print the full fraction); their derivation took `2.89s / 83,432 KiB` and `2.36s / 82,816 KiB`, respectively. Neither positive control needs a three-vertex row. A wrong promotion denominator, a dropped donor term, and a relabeled inside slice each fail the matrix fingerprint. Recomputing the same multipliers on the genuine inside point `0.56572366`, bypassing the fingerprint for this control, gives `+1.08023646939` at both levels, so the arithmetic rejects as well. The pilot does not claim that every altered physical model must make those multipliers fail arithmetically.

A fresh Level-2 solve at `0.56572366` for `(5,4)` is feasible to `2.72e-14` original-row violation (`4.46s / 97,404 KiB`).

From this directory (`scipy 1.16.3`, `numpy 2.3.5`):

```bash
P=python3
g++ -O3 -std=c++17 sim_core.cpp -o sim_core
"$P" -B sim_pi_check.py --output sample.npz
"$P" -B check_sample.py sample.npz
"$P" -B lp_level1.py --cutoff 8 --method highs-ds
timeout 70 "$P" -B lp_level1.py --cutoff 12
timeout 70 "$P" -B lp_level1.py --cutoff 16
"$P" -B lp_level2.py --cutoff 5 --motif 4
timeout 220 "$P" -B lp_level2.py --cutoff 8 --motif 6
timeout 70 "$P" -B lp_level2.py --cutoff 12 --motif 8
timeout 70 "$P" -B lp_level2.py --cutoff 16 --motif 10
"$P" -B verify_farkas.py cert-control-l1-probe.json
"$P" -B verify_farkas.py cert-control-l2-probe.json
"$P" -B verify_farkas.py cert-control-l2-probe.json --plant promotion  # exits 1
"$P" -B verify_farkas.py cert-control-l2-probe.json --plant donor       # exits 1
"$P" -B verify_farkas.py cert-control-l2-probe.json --replay-inside      # exits 1
```
