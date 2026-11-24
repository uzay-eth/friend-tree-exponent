> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# GOLD2 original H3: independent all-age search

The independent physical-chain solve exceeds the requested `5e−7` lower-witness search target: **`6.700007633081629e−7` per core** for the companion threshold support. It is an alternative floor on that same support. A separate search mask above it gives **`1.4978705475163613e−8` per core**. Neither result is inserted into the frozen threshold audit or any release coordinate.

## Ownership and rate model

`cohort.py` loads only the authenticated original GOLD2 payload and receipt via the pinned Gold2 host reader. The loaded lower seams are H3 middle `D<80,c≤10`, H3 endpoint knot `D=88` with stock clamp at `11`, P4 internal `D<24,c≤7` and degree sum `≤32`, and S3 hub `D<20,c≤7`. There are no GOLD1 input weights, caps, or prior endpoints.

For the one real nonfounder mother–middle–child path `h→m→y`, calculate the original outward signed row `P`: base H3 occupation/discount, actual gains, promotions, newborn next motif at post-promotion arity, the two original `g→H3` transfers, and centering only for covered `H3→P4/S3` sends. The complement is computed on the original physical affine row, not on a projection gap. The verified production cap applies after Gold2 interpolation: the two incoming envelope terms are dominated by the independent corner expectation by the triangle inequality; tail gain/promotion coefficients remain in the original intervals and each tail's reciprocal receive is independently bounded by the larger middle jump at the knot and INF. Thus `P+A≤max(0,P)≤max(0,K)` for `0≤A≤[-P]_+`, including the clipped row in a local-cap aggregate. The labelled control checks the compact row formula by actual gain/promotion and deletion on a branching tree (42 paths on the lower payload), including mixed-tail corners.

The physical floor retains each actual exterior rate in both the reward and the corresponding gain event. At each stopped state it verifies a downward-rounded Bellman inequality separately with reward zero and with reward `−P`; the maximum of the two candidate values is valid for `[-P]_+`. Promotions continue the older cohort and create a longer one at the actual post-promotion arity. In stationarity, summation of the nested generator inequalities cancels the registration terms and only the one birth source `μV(J)` survives. The check uses the original `HI` as discount and multiplies by the original `LO` once; dividing raw values by `2^40` reports per-core normalized units.

The separate sibling search instead starts from the adverse exterior box `B`. For each actual outside child and for the unique upstream mother of `h`, use `z_i=max(0,b/2)−b/D_i` for the mother slope `b`. Keep the rest at their adverse boxes. Across these distinct ports, `Σ_i ([z_i−B]_+−[-B]_+)≤[Σ_i z_i−B]_+−[-B]_+`, by convexity and the zero at the origin. Therefore the base box floor plus the increments is at most that path's physical recovery. The three orders `hmys`, `hmsy`, `hsmy` cover exactly one registration per distinguished sibling; `uhmy` covers the unique upstream. Overflow retains the original donor identity and its bounded reciprocal; no new registration at a seam. This decomposition overlaps the physical floor on the same mask. The registration fixture checks 10 chains, 26 distinguished sibling incidences and 2 upstream chains. No unowned raw gap is priced.

## Runs

Command prefix (from `.`):

```sh
taskset -c 11 nice -n 19 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./round2/.venv/bin/python -Bu round2/occupation-quintet-gold2/higher/cohort.py
```

| Lower physical search | Exact normalized floor | Terminal states / potentially paying rows |
|---|---|---:|
| `--kind physical --mother 79 --young 18 --mother-stock 8` | `5988193558027003984698925131/8937592143119350239749734400000000` | 236,250 / 180,803 |
| `--kind physical --mother 100 --min-mother 80 --young 5 --mother-stock 8 --outside-threshold` | `5820589581257439366906303/388590962744319575641292800000000` | 66,528 / 14,553 |

Both floor sweeps exited 0 on the final checker, each printing `literal labelled H3 rows 42`. The earlier small mother `19..37`, middle and young `D≤5` pilot produced `1.1363540411267029e−7` per core for box plus sibling. The larger sibling search with `--mother 100 --young 11 --mother-stock 10 --middle-stock 5 --sibling 3 --upstream 12 --exclude-strict` printed exact box and sibling floors (`2.4995208005564656e−8` and `3.007934913920325e−7`) before timing out in the upstream pass. None of these sibling-search figures is added to either physical result.

The companion frozen threshold floor is `1084423030952180046412413909/1704710415481443450880000000000000`; see [companion replay](../post-fy/README.md). Its independent audit uses that exact receipt with the specified two FY floors. The figures above remain search evidence outside that coordinate.
