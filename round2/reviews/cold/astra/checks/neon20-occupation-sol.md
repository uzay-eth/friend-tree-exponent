# D20 stationary occupation price — independent cold audit

**HOLDS, conditional on the original D20 all-state signed-row certificate and the physical mean exponent.** Both directions reproduce the outward interval **[0.565681716887067, 0.565766354357632]**. This is an inward repricing of the same two D20 witnesses; it does not change their trial discounts or borrow any D18 numerical inputs. I did not repeat the exhaustive GMP row verifier (835c's scope).

## Byte and row provenance

I pinned `ft2-neon@51f4cfadac675fdf1cef84d44935eebb5456aabb`, including the proof, binary format, GMP source and both checked logs. The five Git blobs were checked in the executable audit. Both local originals are 11,190,881,056 bytes and independently hash to the digests in those logs:

| Payload | Mode / exact trial | SHA-256 |
|---|---|---|
| `b3f20_lo.bin` | 2 / 565680/1000000 | `eb93985b2c90a516b65e85934f4ac847779f42b3fdbfcd2af4dc07fd30c37469` |
| `b3f20_hi.bin` | 1 / 565768/1000000 | `55dd6015eca85ca2fba6267d9444ec8807554c456999803d8f9abcc695c090e2` |

The headers encode `D0=20, KT=19, LF=16000, KF=200`, `P4=(12,9), Yt=(11,8), Yb=(12,9), S3=(10,7)`, `S=10^11`, `Sr=10^9`; there are 1,398,838,433 potential entries and 1,398,860,120 total int64 entries. The D18 files have 3,495,091,968 bytes each and trials `.565673/.565775`; the old-flow D18 interval is only a comparison. The six local proof/export files are byte-pinned in the checker. The rows author (4e4f) confirmed those pins at the read-only freeze; I relayed the audit timing and pins to the parent. Replays check that the payload file identity, size and timestamps do not change during each run.

The executable uses the existing binary layout reader and higher-family affine-row evaluator. It independently rebuilds the singleton extraction from raw ticks, all-stock controls, grouped comparison maxima, physical transitions and Bellman recurrences. The shared higher-row evaluator is checked against an independently assembled whole-tree drift and against the pinned source's allocation, rather than presented as a second all-state enumerator.

## Family ownership, signs and units

The decomposition is of `GΦ−λΦ` on the lower witness and `λΦ−GΦ` on the upper. Each source term at an event uses the original degree of the sender; on promotion the mother keeps degree, loses a leaf, gains a core child, and the newborn is `(L,k,h)=(1,0,1)`. The newborn singleton and newborn-edge **value** belong to the mother singleton drift. The newborn higher occurrence belongs to the drift of its connected predecessor motif after deletion of the newborn tip. The `g` objective in the physical fringe solve is the original, separate, charged derivative row on an existing edge, not another payment of the singleton's newborn-edge term. Transfers pair exactly: an omitted one-neighbour receive puts its centered derivative in the receiving larger motif and removes it from the smaller. For tail `f/g` it is the specified transfer polynomial rather than the raw finite difference; degree- and arity-restricted higher transfers switch off precisely when the receiving row is inactive. Potential arithmetic and derivatives are in ticks divided by `S`; receive centers are divided by `Sr` first. Signed charges have units of potential per unit time, Bellman values are discounted potential per newborn, and the `LO` factor makes them stationary charge per core.

The all-family reward is a sum with single apex ownership, independent of which vertex originally carries a comparison-potential charge. `hC`, `Yt` and `S3` are unordered at their symmetric slots; `Yb` distinguishes the continuing child from its sibling. For a motif inside a complete planted fringe, omitted *internal* neighbours use their actual reciprocals; only the apex parent remains an exterior variable. The same exterior rate enters both charge and gain in the joint recurrence. The unit of normalization is one physical core, not one edge or a conditional selected motif. No separate diagnostic bins or separate higher optima are credited in the reproduced endpoints.

## Exits and physical recurrence

* **All stock:** Directly re-extracted every physical stock below `LF` for `k=0,…,4`, including bare legal states. At `L=LF` there is no gain difference in `f`, but promotion reads `LF−1` and its exact newborn-edge value. For `L≥LF+1` it reads the clamped successor `f`; the edge is enclosed by restricting its Bernstein polynomial to `[0,20/(LF+1+k+1)]`. Independent power-to-Bernstein conversion and signed receive/transfer corners match all ten exported seam and above-seam envelopes per side. The upper/lower signs are applied before minimization.
* **Degree exit:** Observe the exact `D=64→65` gain; at `D≥65`, `k≤4`, the continuation depends only on `k`. Degree gains then have zero effect, and promotion, which increases `k`, has rate in `[1−(k+1)/65,1]`. Set `k=5` and above to zero there. Already-exited states cannot return; the degree-65 floor does not duplicate finite reward. By contrast the deficit solve pays zero on **all** degree exits. These are lower bounds on positive rewards, so no omitted terminal debit is needed.
* **Fringe exits:** All ordinary internal arrows except sends to the inert parent slot contribute total rate `n−1/D_root`, including transitions that increase size or core count beyond the cap. The actual outside-to-root arrow adds `r` to the root-gain rate and `r` times the total affine reward slope; outgoing root-to-exterior does not change the fringe. Pay zero after size 16 or core count 4 is exceeded (the exit event still contributes its rate to the denominator). The shape is physical, with no artificial kill, copy, or newborn mark hazard. All moves add one ordinary vertex, so backward solution is acyclic. Cross-multiplied inequalities at both exterior corners and two interior rates were checked exactly; the two corners cover arbitrary adapted exterior inputs because that inequality is affine. The singleton solver uses the same comparison principle with its wider receive interval `[0,(k+1)/2]`.

The whole-tree test includes sender-by-sender ordinary-arrow enumeration on small fixtures, compressed physical clocks at large stock, all families, the `19/20/21` seam, separate Level-4 exits, `k=19/20`, `k=199/200/201`, and stock around 16,000. For each sign it equates signed whole-potential drift to the sum of separately charged occurrences. Replacing sender degrees by recipient degrees changes the resulting drift and is rejected.

## D20 denominator and physical transfer

The original pinned comparisons give `4Sc0=113472636536/115454594604` for lower/upper respectively, `SC0=184219674916/200854054323`, and the two positive `SΦ(T3)=146041684821/146044513514`. The tighter ceiling independently uses **all eight D20 potential arrays**. Charge `g`, `hV`, `P4` to their bottom; charge a `Yb` to its continuing tip (one per `k_apex−1` choices of the other branch). Charge `hC`/`S3` to their apex and `Yt` to its middle, with the corresponding binomial multiplicities. At each charged vertex hold its own exact/embedded state and maximize over the other slots; absent restricted branches must be covered by zero when appropriate. For the 404 stock groups, independently maximize the one-body stock block and the bounded tail control contributions. The upper overflow uses the `k=201` cap; because the recharge coefficient is nonpositive, it also bounds all actual `k>201`. This is a globally valid envelope. It is distinct from the signed generator-row charge decomposition above.

In a finite tree `Σ(k−h)=0` exactly, including the founder. Adding `t(k−h)` ticks with `t=-200000000` (lower) or `-1000000000` (upper) therefore costs no potential. The global cap is respectively `181661912793` or `183051705651` ticks. At a particular physical state use the sharper convex envelope before taking stock-group maxima. The nonnegative deficit between the global cap and this sharper recharged envelope is physically priced up to degree 64, with zero continuation outside. Independent Fraction replay of that deficit at HI, multiplied once by LO, gives these exact per-core ceilings:

| | Lower witness | Upper witness |
|---|---:|---:|
| Denominator | `16959459852924085304227707/12500000000000000000000000` | `17064269776474933055949939/12500000000000000000000000` |
| Rounded up (15 decimals) | `1.356756788233927` | `1.365141582117995` |

The transfer can also be justified entirely in finite-tree Laplace means, avoiding an implicit passage of an unbounded local drift through a spatial limit. Normalize by `M_s=∫e^(−st)E C_t dt` as `s↓μ`. The physical mean asymptotic gives `M_s→∞`; the discounted birth count is exactly `sM_s−C(T₀)`. For each newborn the robust Bellman bound is uniform over its actual complete birth environment. All rewards used in this step are bounded *masks* (including the fixed-arity stock floors), so a change of the discount from `μ` to `s↓μ` changes their uniform birth integral by `o(1)` by domination; this also covers `μ=HI`. Multiply the uniform birth values by the normalized birth count, whose limit is `μ≥LO`; the discounted integral at `μ` dominates the integral at HI for these nonnegative masks. The deficit gives an upper bound on normalized potential by subtracting its lower payment from the constant cap. The founder has no birth credit. Stopped integration by parts for the original potential, bounded by `C0 C_t`, gives normalized signed charge `sign·(s−λ) [Φ]_s` up to a vanishing initial term, with `sign=+1` for lower and `−1` for upper. This justifies division by the respective **upper** potential ceilings, whether or not one separately names a spatial potential expectation. The original positive `c0` protects the comparison. Add for lower, subtract for upper, and round outward only at the end.

## Exact replay

From `[local path]`, run `ionice -c3 nice -n 19 sha256sum b3f20_lo.bin b3f20_hi.bin`, then for each `side=lower,upper` run:

```sh
taskset -c 8 nice -n 17 round2/.venv/bin/python -Bu round2/checks/neon20-occupation-sol.py "$side"
```

Both payload hashes matched. Both completed exact replays printed 2,079 singleton rows, five all-stock floors, 404 cap groups, 25 whole-tree fixtures with 1,671 charged occurrences and wrong-recipient negative control, and 3,914 retained fringe states with 1,259 size exits and 9,060 core exits (distinct cached states, classified core-first). The exact reward ticks on the common `10^29` scale and final outward integer ticks on `10^15` were:

```text
lower 227220866275491127308960 5718952004892475098000 565681716887067
upper 222182631325898847818400 2470851317190905657040 565766354357632
```

An earlier combined run with per-side `--hash` hit the 180-second wrapper deadline during upper fringe replay; the standalone hashes and both final no-rehash replays completed. No partial output is used as endpoint evidence. No production files or board entries were changed.
