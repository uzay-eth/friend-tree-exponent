# Rigorous Unoriented Five-Vertex (Diameter-4) Certified Enclosure for $`\mu`$ (`FT2_U5V1`)

**Pinned claim commit:** `b04cc50c1e391e0420eb7b006ac781d484c68543` (gold1), `b435b4c824e2ba5835f58b85e67a27de7d6a0d49` (gold2).  
**Review and dependencies:** GPT and Claude proof audits and separate historical replays of gold1 and gold2 using the same authored numerical verifier. See [Claim 6 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-6.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-6.md).

**Author:** Gemini 4 Neon (Neon-2).

## Headline Results (`gold2` & `gold1`)

Using a **100% unoriented 5-vertex (diameter-4) local potential** with **8 unoriented motif shapes**, **exact $`c`$-capped inner overflow endpoints**, **full 1-knot tail endpoints at all levels**, **exact Fenchel–Young degree-transfer knots**, and a **standalone zero-floating-point `__int128` + exact Python `Fraction` 2D tensor-Bernstein verifier** (unmodified from 3-family-verified commit `b04cc50`), we prove without orientation, without any founder boundary lemma, and without occupation-measure post-processing:

### Latest Record Enclosure (`gold2`, `certificate_gold2.json` — 6.34 GB/side, unmodified `b04cc50` verifier)

```math
\mu \;\in\; \left[\frac{35083855339718511}{62017003637500000},\;\frac{127993535256830043}{226242120050000000}\right] \;\subset\; [0.565713486333,\; 0.565736986679].
```

- **Certified width:** **`23.5004 × 10⁻⁶`** (`23.500345… ppm`) — **`1.97×` tighter than `gold1` (`46.23e-6`)**, **`3.58×` tighter than Astra Claim 21 (`84.12e-6`)**, and **`3.74×` tighter than Neon Claim 3 (`88.00e-6`)**.
- **Strictly excludes `142/251`:** $`142/251 - \mu_{\mathrm{hi}}^{\mathrm{gold2}} = 3697635659207/56786772132550000000 \approx +6.5114 \times 10^{-8} > 0`$.
- **Lower bound within `0.80 ppm` of `99/175`** (`99/175 = 0.565714285714…`).
- **Cross-family verified:** Independently replayed byte-for-byte on unmodified `b04cc50` verifier by GPT-6.2 Sol Ultrafast (`FINAL VERDICT: CORRECT`).

### Original 3-Family-Verified Enclosure (`gold1`, `certificate_gold1.json` — 1.59 GB/side)

```math
\mu \;\in\; \left[\frac{2825181729593748861}{4994114995300000000},\;\frac{439124781075521909}{776183857175000000}\right] \;\subset\; [0.565702177913,\; 0.565748407438].
```

- **Certified width:** **`46.2296 × 10⁻⁶`** (`46.23 ppm`), first rigorous certification of the fourth decimal digit $`\mu > 0.5657`$. Independently verified `VERDICT: HOLDS` by both Astra and Fable.
- **Dependencies:** Depends only on the definition of the continuous-time friend tree and the real-variable existence of $`\mu = \lim_{t\to\infty} t^{-1} \log \mathbb E[C(T_t)] = \lim_{n\to\infty} (\log n)^{-1}\log\mathbb E[C(T_n)]`$ (`proofs/existence.md` Theorem 1, Astra).

## Architecture (How This Differs from Astra & Neon)

1. **Pure unoriented core-tree potential (`PROOF.md`):** Uses only **8 unoriented tree shapes** covering every connected core tree up to 5 vertices (diameter 4): singleton $`f`$, edge $`g`$, 3-path $`H_3`$, 4-path $`P_4`$, 3-star $`S_3`$, 5-path $`P_5`$, 5-broom $`B_5`$, and 4-star $`S_4`$. There is **no root/founder orientation**, **no founder boundary lemma**, and **no multi-pass occupation refinement**.
2. **Full 1-knot tail endpoints at every level ($`G_2, G_3, G_4, G_5`$):** Every motif endpoint up through Level 5 resolves unbounded degrees via $`\mu_u^{(j)} = \frac{d_j}{D_u}\delta_{T_j(\min(c_u, K_j+1))} + (1 - \frac{d_j}{D_u})\delta_{\mathrm{INF}_j}`$, with exact 1D, 2D (general & symmetric), and 3D (symmetric) `__int128` projection envelopes (`G2 -> G3`, `G3 -> G4`, `G4 -> G5`).
3. **$`c`$-capped exact-degree inner overflow (`kind = 3`):** For $`K_j < D < d_j`$, all core degrees $`c > K_j`$ collapse to a single state $`\mathrm{OVF}_j(D)`$ retaining exact reciprocal $`\mathrm{inv} = 1/D`$, zero promotion jump ($`\Delta^B = 0`$), and exact residual gain-rate interval $`[0, D - (K_j+2)/2]`$, cutting `G2` and `G3` sizes by $`2\text{--}4\times`$ and allowing `d2 = 180` and `d3e = 48` inside a `1.59 GB` binary snapshot.
4. **Zero-floating-point verification:** `quintet_format.hpp` and `verify_quintet.cpp` contain zero floating-point types or operations. Every row inequality is checked in exact `__int128` sub-tick arithmetic (`K_SUB = 1024`, weight scale $`S = 2^{40}`$, center scale $`S_r = 2^{30}`$), combined with exact Python `fractions.Fraction` 16-panel 2D tensor-Bernstein verification for $`D \ge D_F = 2000`$ (`verify_tail_bernstein.py`) and an independent Python `Fraction` concrete unoriented tree oracle (`cert_check.py`, 47 diverse trees).

## `gold2` Certificate Ledger (`certificate_gold2.json`)

- **Preset parameters (`gold2`):** `DF=4000`, `d2=252, K2=22`, `d3m=80, K3m=10`, `d3e=88, K3e=10`, `d4i=24, K4i=7, SP4=32`, `dS3=20, KS3=7`, `d4e=26, K4e=7`, `d5i=13, K5i=5, SP5=24`, `d5o=12, K5o=6, d5c=11, K5c=5, SB5=18`, `dS4=11, KS4=6`, `d5e=14, K5e=5`, `Q_K=10080` (`20` FY knots), `793,112,707` `int64` entries (`6,344,902,168` bytes per snapshot). Verified with the **100% unmodified `b04cc50` verifier**.
- **Binary artifacts (gitignored; reproducible via `build_quintet gold2`):**
  - Lower: `round2/outputs/quintet-gold2/quintet_lower.bin` (`SHA256 = 41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413`)
  - Upper: `round2/outputs/quintet-gold2/quintet_upper.bin` (`SHA256 = f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`)

| Quantity | Lower (`mode = 2`) | Upper (`mode = 1`) |
|---|---:|---:|
| Nominal $`P/Q`$ | `56571368 / 100000000` | `56573686 / 100000000` |
| Finite rows verified (`q5_i128`) | `767,588,844` | `767,588,844` |
| Envelope checks (`1D / 2D / 3D`) | `18,724,635 / 70,439,835 / 17,788,680` | `18,724,635 / 70,439,835 / 17,788,680` |
| `cap_f` / `cap_f_ana` (ticks of $`2^{40}`$) | `55` / `0` | `2,359` / `-1,524` |
| `cap_g` (ticks of $`2^{40}`$) | `792` | `463` |
| `higher_per_core` (ticks of $`2^{40}`$) | `18,370` | `17,240` |
| Total `err_per_core` (ticks of $`2^{40}`$) | `19,217` | `20,062` |
| Global comparability $`m_{\min}`$ (ticks) | `99,227,205,820` ($`\alpha = 1/4`$) | `158,369,484,035` ($`\alpha = 0`$) |
| Global comparability $`m_{\max}`$ (ticks) | `1,808,603,193,693` | `1,808,568,794,159` |
| Exact certified bound | `35083855339718511/62017003637500000` | `127993535256830043/226242120050000000` |
| Certified 12-digit decimal (outward) | `0.565713486333` | `0.565736986679` |
| Python `Fraction` tree oracle (47 trees) | PASS (worst `-1706.27 < 19217`) | PASS (worst `-640.56 < 20062`) |
| `audit_independent.py` (142 adversarial trees) | PASS | PASS |

## `gold1` Certificate Ledger (`certificate_gold1.json`)

- **Preset parameters (`gold1`):** `DF=2000`, `d2=180, K2=26`, `d3m=38, K3m=15`, `d3e=48, K3e=15`, `d4i=16, K4i=8, SP4=21`, `dS3=14, KS3=8`, `d4e=18, K4e=8`, `d5i=10, K5i=6, SP5=18`, `d5o=9, K5o=6, d5c=9, K5c=6, SB5=14`, `dS4=9, KS4=6`, `d5e=11, K5e=5`, `Q_K=10080` (`19` FY knots), `198,565,138` `int64` entries (`1,588,521,616` bytes per snapshot).
- **Binary artifacts (gitignored; reproducible via `build_quintet`):**
  - Lower: `round2/outputs/quintet-gold1/quintet_lower.bin` (`SHA256 = 81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986`)
  - Upper: `round2/outputs/quintet-gold1/quintet_upper.bin` (`SHA256 = 8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383`)

| Quantity | Lower (`mode = 2`) | Upper (`mode = 1`) |
|---|---:|---:|
| Nominal $`P/Q`$ | `56570237 / 100000000` | `56574828 / 100000000` |
| Finite rows verified (`q5_i128`) | `191,375,747` | `191,375,747` |
| Envelope checks (`1D / 2D / 3D`) | `5,747,301 / 20,036,106 / 4,096,512` | `5,747,301 / 20,036,106 / 4,096,512` |
| `cap_f` / `cap_f_ana` (ticks of $`2^{40}`$) | `0` / `0` | `-1,053` / `-6,509` |
| `cap_g` (ticks of $`2^{40}`$) | `401` | `946` |
| `higher_per_core` (ticks of $`2^{40}`$) | `18,785` | `19,890` |
| Total `err_per_core` (ticks of $`2^{40}`$) | `19,186` | `19,783` |
| Global comparability $`m_{\min}`$ (ticks) | `99,882,299,906` ($`\alpha = 1/4`$) | `155,236,771,435` ($`\alpha = 0`$) |
| Global comparability $`m_{\max}`$ (ticks) | `1,808,571,640,916` | `1,808,494,729,937` |
| Exact certified bound | `2825181729593748861/4994114995300000000` | `439124781075521909/776183857175000000` |
| Certified 12-digit decimal (outward) | `0.565702177913` | `0.565748407438` |
| Python `Fraction` tree oracle (47 trees) | PASS (worst `-1603.36 < 19186`) | PASS (worst `-2741.03 < 19783`) |

## Files

- `PROOF.md` — Complete mathematical proof of the unoriented generator, $`c`$-capped endpoint interpolation, 1D/2D/3D envelopes, FY knot transfers, analytic tail, per-core aggregation, arbitrary-rooting comparability, and Grönwall/Dynkin exponent deduction.
- `quintet_format.hpp` — Frozen 512-byte header (`FT2_U5V1`), strict validator, and exact integer indexing (zero floats).
- `build_quintet.cpp` — Untrusted proposal builder + exact `__int128` envelope exporter (`tiny`, `cap_S`, `gold1`, `gold2`, `gold3`; 1.0× RAM in-place DP).
- `verify_quintet.cpp` — Standalone zero-floating-point exact `__int128` verifier (unmodified from `b04cc50`) with `--override-pq` and 18 `--perturb` negative-control hooks.
- `verify_tail_bernstein.py` — Standalone exact `Fraction` 16-panel 2D tensor-Bernstein verifier for $`D \ge D_F`$ and concrete integer-state polynomial oracle (unmodified from `b04cc50`).
- `cert_check.py` — Exact rational combiner, independent Python `Fraction` concrete unoriented tree oracle, and 46-case planted mutation, mode-binding, digest-binding, CLI, and malformed-header control suite (unmodified from `b04cc50`).
- `certificate_gold1.json` — Complete machine-readable exact verification receipt for `gold1` (`46.23e-6`).
- `certificate_gold2.json` — Complete machine-readable exact verification receipt for `gold2` (`23.50e-6`, excluding `142/251`).


## Replay

The listed snapshot paths and digests identify the historical inputs. The snapshots are not bundled. See the [release rebuild and verification commands](../README.md). A new proposal must pass exact verification; a different digest does not reproduce the original receipt.
