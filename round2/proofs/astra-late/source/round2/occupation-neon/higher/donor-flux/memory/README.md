> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Physical limits on four exact high-slope donor apertures

**The strongest sampled `Yt` aperture cannot deliver material stationary slack on that exact support.** Its apex is `(D,k)=(4820,8)`, and even allowing the one donor to send at full `1/2` or lose the full `1/2` gives total reciprocal or loss mass at most **`5.157e−7` per core**. That is about 700 times below the mass needed for a `10⁻⁶` upper-endpoint improvement with the existing denominator ceiling. Even using the *lower* GMP potential comparison (which gives a true ceiling on the normalized contribution), the single `Yt` component contributes less than `6.80e−9`. The large sampled lower `P4` aperture at apex `(1270,9)` also cannot reach `10⁻⁶` even at full send/loss. These are exact-descriptor results, not a ceiling for degree bands or the aggregated old-motif program.

The more accessible `P4` apex `(21,2)` admits an exact **positive** physical occupation on the same support: lower-witness send excess ≥`3.191241684301e−18` per core and upper-witness loss excess ≥`9.43974064652e−19` per core. These are intentionally small lower bounds; the exact conditional isolated-donor calculation is far sharper than the chronological entry bound. No endpoint increment is booked.

## Donor memory at a fixed old receiver

Let `B = n/2 − Σ_y 1/D_y` on the receiver's outside child set (`n` excludes children already in the old motif). Follow the same actual children. Each new promoted child adds one to `n` and contributes `1/2`, so it adds **zero** to `B`. Promotions at existing children do not change degree and neither reset nor remove their prior loss. Only their degree gains contribute, and the generator contribution is exactly `Σ_y g_y/[D_y(D_y+1)]`, with `g_y` the real degree-gain intensity, not a promotion clock. For an existing donor at degree 2, 3, or 4 the loss contributions are `0,1/6,1/4`; gain increments are `1/6,1/12,1/20`. These are pointwise statements before any aperture exit.

For the fourth descriptor in the table below, construct one actual six-core crown. Birth order is the apex, the two middle vertices, the terminal receiver, an extra child of the apex, and the receiver's sole outside donor. In that prefix, divide each wanted jump's intrinsic rate by the sum of all intrinsic jump rates, `HI`, and the exterior parent-port bound `1/2`. This is a lower bound on its discounted first-jump factor. All gains on the host and the fixed donor then remain live simultaneously. The donor starts at degree two with one leaf; retain its gains through degree four, with no restart. Host degree can rise to 21, and all other crown changes, including donor promotions, kill. At `(host degree, donor degree)=(21,2)` the conditional discounted episode has these brackets, allowing any adapted exterior rate `0…1/2` and any discount `μ∈[LO,HI]`:

| Episode integral | Lower | Upper |
|---|---:|---:|
| `∫e^(−μt) /D_donor(t) dt` | `.01842408830` | `.01876980866` |
| `∫e^(−μt) (1/2−1/D_donor(t)) dt` | `.0003051927659` | `.0003167478052` |

The reward stays zero until the apex reaches degree 21; after that, the next apex gain exits. The donor gain rate while resident is `D_donor−1+1/3`, since its receiver has degree three and it has not promoted; the reciprocal donor-to-receiver send also enters the receiver's exact kill rate. Solve backward in both degrees, minimizing/maximizing the fractional-linear first-step expression at the parent-rate endpoints. This honors the actual duration of a degree-two state and the monotone increases to degrees three and four. Multiply the lower entry/continuation value by `LO` using the physical stable-age law. The resulting weighted occupation lower bounds are `1.541077599076e−15` for send and `5.147117156370e−16` for loss. One selected start at each host is enough for this lower bound. The old `P4` occurs with exactly the original tuple of slot states throughout payment; the extra apex child and donor give at least five cores, and no old motif slot is virgin. This books only `c/D` for the lower witness or `|c|B` for the upper, after the inherited directional floor, never the same signed affine term twice.

## Full send/loss ceilings at identical support

For any exact apex `(D,k)`, the same physical stable-age identity gives a **uniform upper** mass: follow a newborn apex with root stock `l=d−j−1`, promotion rate `b=l/d`, and gain in `[l,l+(j+1)/2]`. Above either coordinate cut set continuation zero because both coordinates are nondecreasing. At the target pay one, and elsewhere zero. Solve at discount `LO`, maximizing over the two gain endpoints and rounding upward on `10⁻²⁴`. The cross-multiplied inequality is linear in the gain and is checked at both endpoints, so this covers adaptive exterior input. `HI` times the newborn value is an upper stationary mass. Multiply by the exact tuple's allowed embeddings per apex, and by `1/2` for the single external child donor. For `P4` the embedding cap is the product of child counts along its first three slots; for the `Yt` tuple, each allowed middle has exactly two tips with distinct states and at most one matching fork. The size and core-count restrictions are preserved in the lower construction and only dropped in the upper envelope while bounding a subset on the identical target support.

| Exact aperture (old states `(L,k,h)`, receiver index from zero) | Max weighted mass/core, either direction | Upper-witness contribution divided by ceiling | Lower-witness contribution divided by ceiling |
|---|---:|---:|---:|
| `Yt`: `(4811,8,1),(0,2,1),(0,1,1),(1,1,1)`, slot 2 | `5.157e−7` | `1.435e−9` | `1.396e−9` |
| `P4`: `(78,2,1),(0,2,1),(0,1,1),(7,0,1)`, slot 1 | `2.841e−4` | `6.636e−7` | `4.946e−7` |
| `P4`: `(1260,9,1),(0,1,1),(0,2,1),(1,1,1)`, slot 3 | `4.085e−6` | `4.562e−9` | `4.261e−9` |
| `P4`: `(18,2,1),(0,1,1),(0,1,1),(1,1,1)`, slot 3 | `.001644` | `2.205e−6` | `2.509e−6` |

Use `c>0` with reciprocal send and `c<0` with loss to interpret each cell. The above ratios cap what this method could certify **at the supplied denominator ceiling**. To cap the actual normalized contribution, divide its reward ceiling by the pinned GMP lower comparison `c0` instead. For example that also refutes a `10⁻⁶` contribution from either exact high-degree tuple (first and third rows), but the second and fourth tuples do not have this stronger impossibility at `10⁻⁶`. The code prints both quantities.

These masks may overlap full-row prices, including aggregated `pressure-memory/old-motifs` and high-apex `Yt/P4`. Compare on shared support, or add only a proved disjoint excess split. The fourth tuple is specifically apex `D=21,k=2`, unlike the old-flow `k≥5` band, but do not infer disjointness from a folder name. Spectacular estimates for a single extremely high-degree configuration in four finite samples are not estimated stationary masses for that exact integer descriptor.

## Replay

For each original blob, using the existing Astra helper layout and Python with NumPy:

```sh
taskset -c "$CPU" nice -n 17 timeout 45s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  "$PYTHON" -Bu "$SCRIPT" "$BLOB"
```

`SCRIPT` is this directory's `replay.py`. Both signs were run on CPU 3 at nice 17 with the original `b3f18_{hi,lo}.bin` payloads: `wall=30.07 rss=45076 exit=0` (upper), `wall=22.76 rss=44968 exit=0` (lower). Both verified the pinned payload hashes and passed the exact add-zero/no-reset controls and the physical recurrences. The original GMP verification and potential comparisons remain the inputs documented in [higher](../../README.md).
