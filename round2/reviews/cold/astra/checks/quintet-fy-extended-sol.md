# Cold audit: lower-witness paying-neighbour FY extension

**Pass.** On the lower witness, the actual-degree mask `175/36/5; 10/6` has the following exact stationary bounds, normalized per physical core. The upper bound is on the **same lower-witness mask**, not an upper-witness credit or a bound on all remaining FY slack.

| Quantity | Exact normalized value |
|---|---:|
| Additional FY floor | `2225011506750759438377/5497558138880000000000000000` |
| Matching FY ceiling | `41290159092445185199539/27487790694400000000000000000` |
| Frozen FY floor plus additional floor | `488690443656182809066571/549755813888000000000000000000` |

These are approximately `4.04727235354723e−7`, `1.50212723719761e−6`, and `8.88922738624720e−7`, respectively. The last row excludes the separately frozen D18 `g` credit. No endpoint or board is changed by this audit.

## Physical disjointness and reward

The frozen `d08d9f80188d6b07100863a1f47545cec9ec06a8` joint solver (Git blob `bcc3d580e8423b8fb7d2ef1afa5f10d36e6ffd63`) binds exactly the old D18 full `g` receipt and the FY `175/5/8` receipt. The former takes the full residual only with **both paying-edge endpoint degrees ≤18**. The latter takes individual directed FY atoms only with **both other vertices at a host of degrees ≤5**, and also excludes the D18 paying edge. The extension retains a particular directed atom `(host, payer, donor)` only with actual `D_payer>5`, `D_donor≤5`, and `max(D_host,D_payer)>18`. It therefore intersects neither frozen reward. Reverse directions are separate atoms, even if they occur on the same physical path, and the complementary older-donor mask has the opposite payer/donor cuts. The gates apply at payoff time; neither a newborn payer nor a bin straddling a boundary can mature early. The degree-24 LP cutoff is not an edge mask; its other costs are original `f`, higher-family slack, and comparison potential, none of which spends this directed `g` atom.

The extension reads the original lower witness and computes the signed knot maximum separately at every physical payer degree, then subtracts the same host/donor primal at the actual reciprocal. The unbinned reward is nonnegative and owned by its one donor slot in the paying `g` row. The min/max over eligible degrees bounds that reward from below/above in each bin; arity and donor degree are kept fixed. Transition hazards still use each exact degree. There is no Jensen difference, physical age reset, or replacement of the original `g` cap. Because these masks partition the spending of that cap, the new floor is an additive **rebate to the original edge charge**. Its ceiling cannot be treated as a rebate floor or used to bound unpriced atoms.

## Registration, finite measures, and the rebate

At the second birth, a completed directed path is registered against the unique existing edge that remains when the new leaf is deleted. A mother edge starts at the mother's actual promotion, at post-promotion mother arity with the host at `(2,1)`. A child edge starts at the host's actual promotion. The two directed child–child orders correspond to distinct old-receiver and old-donor stages; the eventual new payer begins at `(2,1)` and continues through all degrees ≤5 with zero instantaneous reward. No positive payment is generated at entry into a bin. All continuations are zero after the respective role exits a cap, and physical degree/arity monotonicity prevents re-entry. The initial founder/seed edges carry no stable-age birth flux; bounded degree cuts also bound the number of retained exceptional registrations per such landmark. Finite starting contributions do not supply an extra newborn source.

One can keep three separate finite occupation measures, for the singleton, oriented old edge, and completed path, using the same core normalization. At finite size each has mass at most a cut-dependent constant times the number of cores before normalization: there is one singleton per eligible mark, at most a bounded number of oriented old edges and ordered paths per retained host, and all visible rates and test values are bounded on the cut. Thus the stationary identities for these bounded local tests are legitimate; finite-size Dynkin has a terminal term and does not license dropping a positive term at an arbitrary finite time. If `V_i` denote the three test potentials and `A_i` the physical continuation generators, their stationary identities have `μ·π_i V_i` on the left and `π_i A_i V_i` plus the preceding stage's birth input on the right. The source at one stage is exactly the output attached to the preceding stage's real promotion. Summing the Bellman inequalities cancels those outputs, leaving the singleton newborn source `μ V_0(J)` and the discount slack. On the floor side using `HI` makes `(HI−μ)·Σπ_i V_i` nonnegative; it may be discarded before multiplying the newborn value by `LO`. On the ceiling side using `LO` makes `(LO−μ)·Σπ_i V_i` nonpositive; discarding it and using `HI` gives the upper bound. These are not independent probability laws normalized afresh at each stage. Converting to ticks for the frozen ledger multiplies the new normalized floor by `2^40` once; no comparison denominator, new `f` credit, or extra `g` reward is implicit.

## Evidence and pins

```sh
taskset -c 5 nice -n 19 timeout 245s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  python3 -Bu round2/checks/quintet-fy-extended-sol.py --wide
```

Output: `live mask controls / registered payments 24464`; `exact reward bands / rational rate boxes 456 54`; both wide totals reproduced exactly with component, state-count and band-count agreement; `PASS paying-neighbour FY lower witness`; `wall=140.20 rss=63448 exit=0`. The mask enumeration, selected primal interpolation/bin extrema and rate corners are independent calculations. The full wide reproduction runs the pinned owner solver serially, and the narrow receipts are authenticated and checked arithmetically. Partial logs `cut36.log` and `cut22-exact.log` are not evidence. The replay uses the pre-authenticated lower payload read-only, with file identity checked; there was no full verifier scan.

| Input | SHA-256 |
|---|---|
| `extended/README.md` | `557cc363ee69f15273a77d601825f378d0e62bfdff84e8f12dfb73831c10dd75` |
| `extended/replay.py` | `4000aa83a0395eeb2324ef9fb97d2a64e0795c9fa08ab58f4774a0cf9d2d3aee` |
| `cut36-floor.log` | `56b31c406cfe1d794a85f6dea700211eaedbf297e5e1cd00470896e0bd28a824` |
| `cut36-ceiling.log` | `3bbe03c7c94bd7baa09bc0d4de225262cb2dd3c0cd134320395f41745252da24` |
| `cut22.log` | `317155cb7c6b5c21ca1a6688f80761deddbd07abcdc8639fd6b54d6af3898e3a` |
| `cut22-singletons.log` | `7682cb4e5c1ac986eb465422d6c0bb62ada204d8b2f167843a7f8718a4c8a676` |
| Authenticated lower payload | `81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986` |

The script also checks the edge/census dependency and original gold receipt digests, the detached source HEAD `b04cc50c1e391e0420eb7b006ac781d484c68543`, and the exact two input digests named by the frozen joint solver.
