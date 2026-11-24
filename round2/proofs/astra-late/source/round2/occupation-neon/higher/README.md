> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Physical prices for the six higher-order Neon D18 charges

The size-12 physical replay, combined with the disjoint first-body numerator in [the renewal calculation](../renewal/README.md), gives the outward bracket **[0.565675597658, 0.565772677869]**. Relative to that calculation's `[0.565674722414, 0.565773666188]`, this moves the lower endpoint inward by `8.75244e-7` and the upper endpoint by `9.88319e-7`. Size means the total number of ordinary vertices in the complete descendant fringe (cores plus their leaf stock, excluding the parent). The separately assigned `size > 12, cores ≤ 4` project uses the same top-core owner and excludes `f,g`; its support is disjoint from this replay's `size ≤ 12` at every moment. No credit from it is included here.

| Certificate | Higher-order stationary reward lower bound (approximate display) | Improvement from this reward alone using universal `C0` (approximate display) |
|---|---:|---:|
| Lower | `1.2769028496806224e-6` | `6.934983888926748e-7` |
| Upper | `1.5010294586749667e-6` | `7.66532674488982e-7` |

## Row ownership and GMP correspondence

Use sign `+1` for the lower blob and `−1` for the upper blob, applied to `ℒΦ−λΦ`. Price only `hV,hC,P4,Yt,Yb,S3`, once per embedding at its top core. In particular a pair of equal-state siblings still represents two distinct tip deletions when there are two corresponding transferred incidences. The replay assigns zero reward to the founder, to fringes outside the size cutoff, and (conservatively) to any fringe with a vertex of degree `D` such that `D+1≥18`. The original nonnegative GMP decomposition permits all these omissions. With size 12 all the visited degrees and their gain successors are strictly interior to `D0=18`, although the narrower motif cutoff tests still apply. Membership in those cutoff tables is checked separately at the old, gain and post-promotion states. An unsupported successor has weight zero.

The exact row consists of `−λq`, leaf-stock and within-motif send gains, the promotion at each vertex, center terms for routed outgoing sends, actual-rate terms for unrouted sends, and incoming centered transfers from connected tip deletions. The newborn part is owned by the old embedding left after removing the newborn. The incoming part is owned by the enlarged embedding once per deleted tip and equals `(1/D_tip−r̄) Δq_old` in its proper receive direction. These assignments follow §3 of `round2/proof-b3f-motif4.md` at the Neon commit and blob pinned below. In the tracked copy of the pinned GMP code (`round2/attacks/neon1/round2/cert/cert_b3f_verify_gmp.c`, historical reference), the interior loops use these connected deletions and newborn assignments:

| Replayed row | Incoming after deleting a tip | Newborn extension owned by this old row | GMP source |
|---|---|---|---|
| `hV` | Two `g` rows (upper/lower endpoint deleted) | `Yb` at top, `Yt` at middle, `P4` at bottom | 1375–1435 |
| `hC` | Two `g` rows (either child deleted) | `S3` at apex, `Yb` at either child | 1528–1584 |
| `P4` | Two `hV` rows (either endpoint deleted) | None | 1642–1684 |
| `Yt` | `hC` (top deleted), two `hV` (either fork tip deleted) | None | 1729–1774 |
| `Yb` | `hV` (short branch deleted), `hC` (deep tip deleted) | None | 1826–1868 |
| `S3` | Three `hC` rows (one per leaf deleted) | None | 1912–1956 |

Every unrouted directional coefficient `d` in a signed row has the GMP floor `min(0,d/2)` per donor. The extractor first checks the sum of this floor and the other row terms is nonnegative for each replayed embedding. It replaces only these unrouted terms by their true rates for in-fringe donors. All remaining exterior occurrences have the same physical variable `r=1/D_parent`; the top-owned sum is `a(z)+b(z)r`. A routed exterior send stays centered; its correction is in a larger embedding at the exterior ancestor, which this root does not own. Summing all eight reconstructed charges on complete trees, with no exterior edge, agrees exactly with direct `Cert.potential` drift on 114 controls per blob.

The first-body reward in the renewal calculation prices only `f`, including the newborn and centered terms belonging to `f`. A centered transfer into a higher row is paired with its subtraction in the smaller row; neither the smaller row floor nor its surplus is booked here. Thus the higher-order reward can be added to that renewal numerator. The older whole-drift estimates and the residual exploratory prices are not added.

**Suffix overlap:** `arms/` and `renewal/next/` both price the six top-owned higher families on `size > 12, cores ≤ 4`. The latter's `higher` run is the interior overlap, and its `both` run also contains that overlap alongside a `g` edge component. Do not add the two higher suffix prices, or add `arms/` to `both`. For the same certificate and drift numerator, take the stronger bound on the common higher objective, or first split into verified disjoint rewards. For example, with a separately certified `g`-only price, compare `arms + g-only` against the joint `both` bound and keep the larger; the neutral combination is at the reward level, before dividing by a common validated potential ceiling. The `size ≤ 12` price above and the `f`-only renewal numerator remain disjoint from either suffix. A split into three-core versus four-core charges is only additive if the `arms/` credit is split on the same partition before comparison.

## Uniform physical occupation

The projection at birth is always one core with one leaf and one parent slot. Birth labels can all be stored as `2`, since neither the score nor any intrinsic physical rate reads those labels. Keep exact physical gain and promotion moves on the entire descendant fringe. The only unknown projected jump is a gain at its root, of adapted rate `0≤r≤1/2`; changes elsewhere only change that rate. Every projected jump adds one ordinary vertex, so the backward solve at cutoff `n` is finite. Put `V=0` above cutoff, let `(q_i,z_i)` be intrinsic jumps, and let `z+` be a gain at the tracked root. At discount `HI`, compute and round down

```
V(z) = min_{r∈{0,1/2}} (a(z)+b(z)r + Σ q_i V(z_i) + r V(z+))
                       / (HI + Σ q_i + r).
```

The code checks both endpoint inequalities after rounding. The cross-multiplied inequality is affine in `r` and hence holds pointwise for all intermediate rates, irrespective of their history dependence. The reward itself is nonnegative there. Discounted Dynkin from the newborn gives `V(J)≤E∫ exp(−HI·t) (a+b r)_t dt` in every physical birth environment. Apply the stable-age identity with `μ∈[LO,HI]` to get `π(a+b r)≥LO·V(J)`, with no birth-class weights needed. The cutoff reward is a bounded local observable: a bounded-radius exploration either lists the whole fringe or sees more than the size cutoff, and reads the parent degree. The finite-tree-to-stationary transfer uses the same normalized mean drift argument as the renewal calculation. Add numerators at the original certificate `λ`; the combined code conservatively starts from the already outward-rounded first-body endpoints and applies the added reward divided by the same outward denominator ceilings. Final rounding is again outward.

`--ledger` gives a weaker, separately additive breakdown from the same physical process, discarding any path after its first exterior send and minimizing the reward over the shared exterior variable. The intrinsic rates are unchanged before that send and the extra killing hazard is at most `1/2`; this is a uniform lower occupation. The cap-9 upper breakdown is approximately `9.55624e-11` from row floors, `5.15988e-7` from actual in-fringe sends, and `1.99768e-8` from jointly minimizing the exterior send. The lower counterparts are `9.58220e-11`, `2.07288e-7`, and `1.75319e-8`. The shared component is the gap between a joint minimum and separate minima, not a new motif owner. `--dump PATH` records the rational mass and `(family, embedding IDs, row floor, intercept, exterior slope)` for each paid context of this suboccupation. It is not added to the stronger Bellman number.

## Portable replay

Run with Python 3.11+ with NumPy, assertions enabled, and the two immutable Neon D18 payloads. Set `ASTRA` to an Astra checkout root, `PYTHON` to the Python executable, and `NEON_UPPER` and `NEON_LOWER` to the payload paths. The entry point can run from any working directory; it uses no checkout-relative discovery in CLI mode. You can override the two helper paths if copying the inputs to a standalone directory.

```sh
AUDIT_SOURCE="$ASTRA/round2/attacks/neon1/audit.py"
PHYSICAL_SOURCE="$ASTRA/round2/finite/residual_entries.py"
REPLAY="$ASTRA/round2/occupation-neon/higher/extract.py"
"$PYTHON" -Bu "$REPLAY" "$NEON_UPPER" --audit-source "$AUDIT_SOURCE" --physical-source "$PHYSICAL_SOURCE" --cap 12
"$PYTHON" -Bu "$REPLAY" "$NEON_LOWER" --audit-source "$AUDIT_SOURCE" --physical-source "$PHYSICAL_SOURCE" --cap 12
```

Expected cap-12 output: 114 exact whole-tree identities and 20,298 Bellman contexts per sign; exact higher-order rewards `1876286823343708393724246511/1250000000000000000000000000000000` (upper) and `798064281050388986950514709/625000000000000000000000000000000` (lower). The CLI verifies payload size and full SHA-256 before mapping the tables. It rejects other payloads even when the header matches. The helper files used here have SHA-256 `d548a9aafa899ecab2b00a7f6ad0a4e6f828600088780f0a15f08927f17f0bfe` (`audit.py`) and `59f4ef6f915d6d7d932076dfe9170d7e6e5fa5d819cf3d0c11baf7932baa1277` (`residual_entries.py`).

Both CLI commands were replayed from outside the checkout with explicit inputs: upper emitted `0.565772677869` (7.50 s, exit 0) and lower emitted `0.565675597658` (7.54 s, exit 0). The size, optimized-Python, and cap-13 negative controls each exited 2 before reading the payload tables.

## Production obligations and provenance

1. **Pin the base proof and full payloads.** Neon source commit `cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`, tree path `round2/proof-b3f-motif4.md` (Git blob `3bbda54c3cbdbd235a6cf4af8ecf80ca75d16ae7`). Given `NEON_SOURCE`, the proof is available via `git -C "$NEON_SOURCE" show cec3855cc52ab4265ce3e36b8a7a38f1ae257f77:round2/proof-b3f-motif4.md`. The public GitHub remote did not resolve this commit when checked; use a checkout containing that pin, not a link to an untracked scratch clone. Astra tracks byte-identical copies of the GMP verifier (`round2/attacks/neon1/round2/cert/cert_b3f_verify_gmp.c`, historical reference) (SHA-256 `f2a1b91575a82e0039b40ac1eb247e86291c803219eb6176306393c70475796b`) and format header (`round2/attacks/neon1/round2/cert/cert_b3f_format.h`, historical reference) (SHA-256 `05055d3cb789c3ae83ca20321b7e7def5f76b30bf2960f7f60f03c3947cd3729`). Both payloads are 3,495,091,968 bytes. Lower SHA-256 `9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec`; upper `1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41`.
2. **Retain both full GMP passes.** The tracked production rehearsal archive (`round2/attacks/neon1/outputs/neon1-production-rehearsal.tar.gz`, historical reference) and [rehearsal instructions](../../../../../../reviews/cold/astra/attacks/neon1/REHEARSAL.md) preserve the pinned verifier passes: eight row families with zero failures per sign, 34 exact concrete/boundary trees per sign with zero failures, positive `4S c0` of `113465364604` lower and `115462837236` upper, and `S C0` of `184124847315` lower and `195820675182` upper. The full proof is needed to discard nonpriced rows, including cutoff and tail cases outside this replay. The sparse replay does not re-run those large checks. For a fresh full rehearsal, supply `NEON_SOURCE`, `AUDIT_PYTHON`, `REPLAY_OUT` (fresh path), and `CPU` explicitly to the production script (`round2/attacks/neon1/replay-production.sh`, historical reference), passing `"$NEON_LOWER" "$NEON_UPPER"` in that order.
3. **Retain the first-body dependency for the combined endpoints.** The output labeled `additional lambda improvement lower` uses only the pinned Neon `C0`. The output labeled `outward endpoint with first-body renewal` additionally depends on the f-only endpoint and denominator ceiling from [renewal](../renewal/README.md). Package `renewal/check.py`, `renewal/denominator.py`, and `rows/charged-f.json` with the replay. Run both renewal commands with the same explicit Python and blobs; verify the seven ceilings on each side against `YOUNG_CEILINGS`, and verify the outward f-only endpoints and denominator ceilings against `F_EDGE` and `F_DEN`. The hashes of the reviewed `check.py`, `denominator.py`, and `charged-f.json`, respectively, are `423cfdddf34559c705b55b7f50f8c4a57f655915226d23cf7170028a0ff81388`, `eb9662f8d91d104fa55c00bf10375463a014cd9f3b84ef150c6cdeb8e5510143`, and `415a712878e03b71c431a36bf4cf8eb11a047c2a47cafb5d9298fe23bc1561da`. These renewal inputs must travel alongside this directory for a cold replay of the combined number; they are not part of the base Neon GMP archive.

   ```sh
   RENEWAL="$ASTRA/round2/occupation-neon/renewal"
   "$PYTHON" -Bu "$RENEWAL/denominator.py" "$NEON_LOWER" "$NEON_UPPER"
   "$PYTHON" -Bu "$RENEWAL/check.py"
   ```

   The existing renewal scripts resolve their data and helpers relative to their source files, so retain the indicated directory layout. Both were replayed; all 14 reported integer ceilings matched, and the f-only bracket and its two denominator ceilings reproduced.
4. **Retain the prior and stationary transfer.** `HI=.565775` is the pinned Neon upper endpoint. `LO=.565673570322` comes from the [Astra joint physical enclosure](../../../../../astra/source/round2/joint-enclosure/README.md); it is used only as the lower multiplier for a nonnegative occupation, never in place of the Neon lower trial `λ=.565673`. A cold check of the numerical composition should distinguish reproducing that previously supplied prior from re-verifying its full certificate. The [spatial limit](../../../../../astra/source/round2/fringe/spatial-limit.md) and [birth identity](../../../../../astra/source/round2/degree-law.md) supply the analytic bridge from finite-tree drift to the stationary reward. Independently check the row attribution and this bridge; the 114 whole-tree checks are a telescoping cross-check, not a substitute for either proof.
