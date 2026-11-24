# Cold audit: wider native complementary P4 replacement

**Correct in scope.** For degree caps `(16,9,9,16)`, arity caps `(4,4,4,3)`, side memory three, and mother memory twelve **only in the mother component**, the exact physical-Palm lower floors are:

| Witness | Exact floor per core | Approximate |
|---|---:|---:|
| Lower | `23575773396079784114802186917/25500000000000000000000000000000000` | `9.24540133e−7` |
| Upper | `43746697941322995419/56100000000000000000000000` | `7.79798537e−7` |

Replace the medium complementary P4 credit inherited from `e911a1c` in `f7cc2e1c77192ecea8f231519f4afd2bf50aab9b`. The new outward interval is **`[0.565679672983674, 0.565768735586042]`**. The old root P4 band, S3, Yb and hV credits stay booked. The result concerns only this low-apex complementary mask; no high-apex calculation enters it.

## Native row and lifetime arithmetic

- The input path caps put old descriptors at `D≤16` and every old gain lookup at `D≤17`, all within the original inner tables. The two intermediate P4 membership checks precede the P4 access, and parent/child orientation and dtype are carried through the authenticated maps. The native row matches the original P4 discount, own gains and promotions, and the two centered incoming hV transfers from deleting a connected tip. No size-five weight belongs in the P4 promotion term. This use of hV transfer coefficients does not spend the differently allocated residual hV row.
- Rates and the discount are integer on `lcm(10^6,2,…,Dmax)` (wide, including mother12: `9009000000`); charge rates and centers use `153153000000000`. The 64-bit rate construction fits before runtime checks: documented validated CLI caps are at most 16 on tags and 23 in memory, and the global `Dmax≤24` unit is `66927861000000`. The checked arrow bound `32*unit`, at most 20 arrows, plus discount keep the summed denominator below `2^58`. A transition tick is narrowed after each widened accumulation; tick differences are less than `2^61`; all Bellman products and residual sums fit signed 128 bits, conservatively below `2^120`. Signed division truncates inward in different directions, so the remainder test subtracts one for negative lower and adds one for positive upper. Unsigned wrapping is confined to the checksum. Checks are active in the optimized C++ build.
- The charge's final multiply by `10^7` must also fit **before** it can be narrowed. For these pinned payloads the independent binary scan found max absolute P4/hV ticks `2377223860/10037313027` (lower) and `2019021826/10079031623` (upper). Counting four gain rates each at most 16, four promotions and two centered transfers each at most one, at most sixteen directional floor copies each at rate `1/2`, and `|Δ|≤2·max_tick`, gives a 111-bit absolute bound including the final stretch. This bound is on the authenticated release, not arbitrary substituted tables. No certificate is claimed for hand-built inputs to the raw C ABI.
- Reachable exits are tested before ranks or table access. The four coordinates never decrease. The short tip owes one spine child, and a not-yet-registered donor owes another at its slot; when the corresponding child is born `k` rises and that demand disappears together. Neither a killed degree nor a killed arity-demand registration can reenter. A live registration can be outside the **reward** mask and later grow into it. Mixed-radix order is increasing in finite degree/arity, donor continuation/overflow, and mother continuation/overflow; donor attachment and path extension enter later blocks. Distinct simultaneous continuation alternatives are summed at one promotion rate. The actual child's reciprocal multiplies both its slope and the corresponding receive update; only arrows with identical unpriced successor are merged. The direct adverse-corner drift check verifies all four rounded value claims at each live row.

## Shared mother birth law and exact LP

The conditional ticks are expectations of discounted lifetime reward, **before** multiplication by `μ`, from the actual birth, with mother degree measured before promotion (unchanged by it). Write `x=π/μ` and `β(f)=x(bf)`, `b=(D−k−1)/D`. This is one normalized physical birth measure because `πb=μ`. For either sign, integrate the sum of the six matching conditional tick bounds against the **same** `β`. Exact bins 2 and 3 keep their old entries for the five other components; all classes `4…12` and `13+` still use those components' cap-three overflow bound. The mother component alone uses its fresh 2–12/13+ partition. The referee checked all 48 fine mother ticks against their old-class enclosing intervals and checked the replacement identity separately in every joined bin/lane. In particular degree 12 stays exact; the degree-13-and-up objective is constant even across the LP's several coarse cells.

The sided LP is a relaxation of the physical law: finite degree vertices lift their two receive rates into the exact affine boxes, as well as the true `μ` into the discount endpoints; the census test uses the same actual rates. The same edge incidence measure supplies both opposite-end reciprocals. In the degree tail `p=1/D` is represented by endpoint mixtures; at overflow arity `≥7`, `0≤b≤1−8p`, the excess stub count and the full child receive use separate nonnegative columns. The physical identity `kp=1−b−p` still fixes the reciprocal-weighted parent count. The LP's birth objective uses `b` and a value constant within each overflow degree class, so no untracked product needs to be preserved there. The additional `β≥βA` constraint is a lifetime direct-birth bound from each start class, with nonnegative reward dropped beyond degree 34 and a monotone incoming reciprocal bound. It does not impose a later-age pressure conditional on birth degree.

The companion reconstructs all 144 Palm kernel entries, then checks explicit exact lifts of nine physical censuses and one mixed/overflow edge, including degree 12/13, degree 18/19, arity saturation, and a large-degree high-arity case. It uses the pinned topology builder but independently constructs the objective and repairs every rational constraint: HiGHS proposes the dual, multipliers are rounded to `10^−12`, and the largest positive defect of each class is subtracted with mass bound `1/LO`. The class bounds follow from `x(1)=x(k)=1/μ`, the same total edge mass, `xρ=(1−μ)/μ`, and total Palm slack at most one. All **3,820 columns**, including zero-cost overflow, edge and slack columns, were checked for each of the eight sign/lane/cap solves; there are 462 rows. Multiplication by `μ` is rounded out over the incoming bracket, and the result is intersected with the shared-class uniform bound. The independent rational repair reproduces the producer's exact refined intervals. The cap-three shared floors are about `8.84750046e−7` / `6.67298549e−7`; they are controls, not additional payments.

## Exact replacement and disjointness

The current-descriptor old-band exclusion contains the **entire** original root P4 band. The size/core lower-count gate excludes both frozen higher-body reserves, just as in the medium audit. The broader run replaces the whole previous complementary six-part price; it cannot be added to it. hV, S3 and Yb have different allocated rows, even if an arrow or geometric support intersects. The corrected potential stays fixed, so use the original full denominators `1357493815758092773/10^18` and `273436004214812983/(2·10^17)`.

Subtract the old complementary floors `359297001126132137/10^24` and `6963149521191559/25000000000000000000000` in signed numerator units, then add the new floors above. The increments are about `5.65243132e−7` and `5.01272556e−7` per core. Using the exact `f7cc2e1` endpoints, the final exact fractions before directed rounding are:

```text
6853566920645426224860150084809360393/12115632305640977999025000000000000000
10848445658578037268151132921697974717/19174699795563760432875000000000000000
```

## Evidence and pins

The companion's SHA manifest pins the native source, Python bridge, pool, input JSONs, archived medium Python receipt and transitive row/LP sources. In particular the reviewed `solver.cc` SHA is `af0df8e91f6ea7575ef77a050d8bf55be9745277105c3f75b19227d19ae8abfb`, and `charges.hh` is `f9edec256b6ef4cd07d4a222110265facfa7c8580b49b3b4922dbaa041b281e0`. The tick fixtures additionally SHA-authenticate both original 3.49 GB witnesses. No producer source was changed.

- Native medium was freshly compiled and run on CPU 1, nice 15. All 24 archived Python component endpoints and 18 component counts agree: 6,017,604 live states and 24,070,416 drift checks. The complete conditional JSON was byte-identical (`30c8466ce23478b0a135b4838ab098bc36dbb85b66a84ed382b64acd55a6f76a`). The UBSan build on CPU 13 with `--row-parity` independently compared **98,280** terminal requests (91,794 supported) to the rational extractor and returned all the same entries without a UBSan diagnostic. This uses the archived all-six Python solve; it does not rerun that slow solver.
- A separate UBSan bridge against the pinned C++ definitions passed **2,768 accepted and 116 rejected** wide seam and random row fixtures, comparing all 44 accepted row ticks against rational reconstruction. It also checked **23,072** degree/arity exits, **58** short-prefix demand probes and **72** signed rounding controls. The wide fixture is a sampled independent rational check, in addition to each native solver's full drift checks.
- The reviewer's three fresh native wide invocations on CPU 1, nice 15, independently reproduced **every conditional JSON tick**, byte for byte. Head, tail and mother12 ran respectively **76,893,297**, **64,145,862**, and **25,476,804** live states and four signed drift checks per state; allocated slots equaled reachable states throughout. The last replay command below verified these files and counts together. The producer pool was also run separately on CPU 2, nice 15, and its two exact refined intervals matched the independent rational calculation.

| Fresh receipt under `OUT` | SHA-256 |
|---|---|
| `medium-ubsan.log` | `6cb04c650082ceafd714e0bb89beb15f0030c37b16abbf780fb9f1439e73ebe0` |
| `wide-head.log` | `a7be74f203e2a8a67f51392ec9d6ffcc8c4adc414d5a5265c32d749ca6c540fb` |
| `wide-tail.log` | `7ccd86037a602ef87e81a914f5079c00a2cec44b6b50cf3a39043f17ea74cb3e` |
| `wide-mother12.log` | `5663a40b5d7e854b0c605aeea408a72b22b129dc3d9941f7ad3d868d6fae7cb6` |

Replay local accounting from `[local path]` with the reviewer's fresh native receipts:

```sh
P=round2/.venv/bin/python
OUT=[local path]
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15 taskset -c 6 timeout 110s \
  "$P" -Bu round2/checks/path-fast-sol.py --wire b3f18_lo.bin b3f18_hi.bin \
  --scratch "$OUT" --fresh "$OUT"
```
