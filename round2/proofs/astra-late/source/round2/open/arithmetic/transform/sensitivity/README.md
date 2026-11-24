> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** supporting Astra work; the cited inputs and any recorded review are scoped in the linked claim index and historical record; no blanket independent review is asserted. See the [release index](../../../../../../README.md) for exact scope and citations.

# Sensitivity to the physical enclosure

The old funded three-credit lower endpoint and results depending on it are on hold. The frozen transform proof and replay are unchanged. This note uses the **original Gold2** endpoints only, as recorded in the [production audit](../../../../../../../../reviews/cold/astra/attacks/neon2-quintet/gold2/README.md):

```
[35083855339718511/62017003637500000, 127993535256830043/226242120050000000]
```

| Quantity | Sharp minimum for original Gold2 | Witness in the interval |
|---|---:|---|
| Unrestricted rational denominator | 175 | `99/175` |
| Rational denominator with square discriminant | 94,251 | `53320/94251` |
| Degree 2 primitive height | 56 | `x²+56x−32` (irreducible mod 11) |
| Degree 3 primitive height | 15 | `13x³+13x²+15x−15` (mod 11) |
| Degree 4 primitive height | 6 | `x⁴+3x³+3x²+6x−5` (mod 17) |
| Degree 5 primitive height | 3 | `x⁵+3x⁴+3x³+3x²+2x−3` (mod 2) |

In particular the quadratic witness of height 56 and the fraction `99/175` must both be accounted for when replacing the old lower endpoint. The square-discriminant result concerns only reduced `p/q` with `p²+4q²` an integer square; it is not a necessary condition for the exponent to be rational.

`nice -n 10 timeout 12s python3 -IB round2/open/arithmetic/transform/sensitivity/check.py` scans eligible numerators directly by denominator and tests square discriminant by exact integer square root. For each degree it enumerates primitive integer polynomials up to the frozen witness height, excludes candidates by exact interval bounds and monotonicity or subdivision, and proves sharpness by a sign change and irreducibility modulo a prime. Only polynomial helpers from the frozen checker are reused; its endpoints are not used.

## Separate reruns on supplied candidate enclosures

| Enclosure | Status of enclosure | Sharp heights, degrees 2–5 | Least unrestricted rational | Least square-discriminant rational |
|---|---|---|---|---|
| Original Gold2, exact endpoints above | Original certificate | `56, 15, 6, 3` | `99/175` | `53320/94251` |
| `[10043824822832044427/17754235186087500000, 127993535256830043/226242120050000000]` | Repaired candidate; only H3 input held on cold read, composite audit pending | `63, 15, 6, 3` | `241/426` | `53320/94251` |
| `[0.565714467904, 0.565735974375]` | Corrected frontier, taken at the supplied decimal endpoints; cold audit pending | `63, 15, 6, 3` | `241/426` | `53320/94251` |

Both candidate reruns found `2x²−63x+35` (irreducible mod 3) as a sharp quadratic witness, and the same cubic, quartic, and quintic witnesses as original Gold2. Sharpness is relative to each supplied interval. This arithmetic does not authenticate either repaired composite and does not reinstate the withdrawn old three-credit statement.

Replay the candidate rows separately from the Astra root:

```sh
nice -n 10 timeout 12s python3 -IB round2/open/arithmetic/transform/sensitivity/check.py --interval 10043824822832044427/17754235186087500000 127993535256830043/226242120050000000
nice -n 10 timeout 12s python3 -IB round2/open/arithmetic/transform/sensitivity/check.py --interval 0.565714467904 0.565735974375
```

The alternate mode checks interval arithmetic, not endpoint provenance. Increase `--cap` if the denominator scan does not reach a witness. If a sharp height is above the search ceiling, the checker reports only the lower bound.
