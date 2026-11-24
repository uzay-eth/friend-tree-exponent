# Certified Stationary Marginals of the Friend Tree (`FT2_M5V1`)

**Pinned claim commit:** `6c8afec5cd1b7f98bbcd7c6e7848c0d631a50f7d`.  
**Review and dependencies:** GPT-family full audits and historical replay; Claude supplied independent simulation consistency checks and checked the degree-two sum. See [Claim 10 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-10.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-10.md).

**Author:** Neon-2 lead (Gemini 4 Neon).  
**Scope & dependencies:** Certifies exact rational 4- to 5-digit enclosures for **eight stationary core observables** of the continuous-time friend tree — $`p_2, p_3, p_4, p_5, p_6, \pi(J), \pi(2, 2), \pi(L_o = 0)`$ — plus exact rational corollaries (strict monotonicity $`p_2 > p_3 > p_4 > p_5 > p_6`$, tail mass $`\pi(D_o \ge 7)`$, bare degree-$`\ge 3`$ mass $`\pi(L_o = 0, D_o \ge 3)`$, and conditional neighbor-reciprocal expectations $`\mathbb E_\pi[1/D_{\mathrm{parent}} \mid J]`$ and $`\mathbb E_\pi[1/D_{v_1} + 1/D_{v_2} \mid (2, 2)]`$).
- **Dependencies:**
  1. Certified rational bracket $`\mu \in [\frac{565713486333}{10^{12}}, \frac{565736986679}{10^{12}}]`$ from **Neon-2 Claim 6 (`gold2`, commit `b435b4c`)**.
  2. Explicit contractive stationary fringe measure $`\mathcal M_\mu`$ and exact free-tree coboundary $`R_{\psi, \mu}`$ from **Neon-2 Claim 3 (`claim14-proofs/README.md`, commit `ea6a093`, contour-free)** — which Theorem 1(2) identifies with the downward-fringe marginal $`\pi_{\mathrm{fringe}}`$ of Astra Claim 4 / Claim 10 (`operator.md`, `spatial-limit.md`). (Even without $`R_{\psi, \mu}`$, `PROOF.md` Theorem 2(1) yields bounds on the discounted Cesàro occupation limits from $`T_3`$ using only positive mean growth $`e^{-\mu t}\mathbb E_{T_3}[C(T_t)] \to c_* \in (0, \infty)`$.)

---

## Certified Enclosure Table (`high` preset, 16 × 413.1 MB `FT2_M5V1` snapshots)

Each of the 16 snapshots (`obs{1..8}_{lower,upper}.bin`, stored gitignored in `round2/outputs/marginals/`) is verified by the standalone zero-floating-point 128-bit verifier `verify_marginal.cpp` (`49,204,104` finite motif rows checked at **both** endpoints $`\mu_{\mathrm{lo}}`$ and $`\mu_{\mathrm{hi}}`$, plus `8,035,467` $`\lambda`$-independent cross-level envelope checks applying to both endpoints, with in-memory SHA-256 binding), the exact rational 16-panel 2D tensor-Bernstein tail verifier `verify_marginal_tail.py` ($`D \ge D_F = 1600`$), and 35 exact `Fraction` concrete-tree oracle checks per side.

| ID | Observable | Certified Outward Decimal Enclosure | Width | Exact Rational Enclosure (`certificates_all.json`) |
|---:|---|---:|---:|---|
| 1 | $`p_2 = \pi(D_o = 2)`$ | **`[0.38337348, 0.38351583]`** | `1.4234e-4` | $`[\frac{164657657862788533}{429496729600000000},\; \frac{164718791390725953}{429496729600000000}]`$ |
| 2 | $`p_3 = \pi(D_o = 3)`$ | **`[0.14517013, 0.14528215]`** | `1.1202e-4` | $`[\frac{62350096588093383}{429496729600000000},\; \frac{6239820708693547}{42949672960000000}]`$ |
| 3 | $`p_4 = \pi(D_o = 4)`$ | **`[0.08033089, 0.08041288]`** | `8.198e-5` | $`[\frac{34501858047442641}{429496729600000000},\; \frac{34537068106288267}{429496729600000000}]`$ |
| 4 | $`p_5 = \pi(D_o = 5)`$ | **`[0.05222046, 0.05228311]`** | `6.265e-5` | $`[\frac{7008911917707}{134217728000000},\; \frac{11227711645728793}{214748364800000000}]`$ |
| 5 | $`p_6 = \pi(D_o = 6)`$ | **`[0.03715899, 0.03720440]`** | `4.540e-5` | $`[\frac{15959666843328983}{429496729600000000},\; \frac{15979165909025803}{429496729600000000}]`$ |
| 6 | $`\pi(J) = \pi(D_o = 2, c_o = 1)`$ | **`[0.25227667, 0.25230053]`** | **`2.386e-5`** | $`[\frac{108352006706978037}{429496729600000000},\; \frac{5418112564593513}{21474836480000000}]`$ |
| 7 | $`\pi(2, 2) = \pi(D_o = 2, c_o = 2)`$ | **`[0.13108737, 0.13122678]`** | `1.3941e-4` | $`[\frac{28150798688597501}{214748364800000000},\; \frac{56361470507031517}{429496729600000000}]`$ |
| 8 | $`\pi(L_o = 0) = \pi(D_o = c_o)`$ (all bare cores) | **`[0.14520643, 0.14540145]`** | `1.9502e-4` | $`[\frac{7795710920099093}{53687091200000000},\; \frac{62449445791953999}{429496729600000000} + \frac{40}{\mu_{\mathrm{lo}}}2^{-1600}]`$ |

*(Note: for ID 8, the certified finite-range upper bound on $`\pi(L_o = 0, D_o < 1600)`$ is $`\frac{62449445791953999}{429496729600000000}`$, and Theorem 3 of `PROOF.md` adds the exact analytic tail $`\frac{40}{\mu_{\mathrm{lo}}} 2^{-1600} < 10^{-479}`$, recorded in full exact rational form in `certificates_all.json`.)*

---

## Certified Corollaries (`PROOF.md` §§4–5 & `certificates_all.json`)

1\. **Strict monotonicity of the stationary core degree law ($`2 \le k \le 6`$):**

```math
p_2 > p_3 > p_4 > p_5 > p_6,
```

with certified positive rational lower gaps $`p_2 - p_3 \ge 0.23809133`$, $`p_3 - p_4 \ge 0.06475725`$, $`p_4 - p_5 \ge 0.02804779`$, and $`p_5 - p_6 \ge 0.01501606`$.
2. **Independent 3-way degree-2 consistency:**
   $`\pi(J) + \pi(2, 2) \in [0.38336404, 0.38352731]`$ contains $`p_2 \in [0.38337348, 0.38351583]`$ (consistent with $`\pi(2, 0) = 0`$ from `PROOF.md` Theorem 4(1) across six independent `413 MB` certificates).

3\. **High-degree stationary core tail mass:**

```math
\pi(D_o \ge 7) \;=\; 1 - \sum_{k=2}^6 p_k \;\in\; [0.30130165,\; 0.30174603].
```

4\. **Dormant / bare cores of degree $`\ge 3`$:**

```math
\pi(L_o = 0,\, D_o \ge 3) \;=\; \pi(L_o = 0) - \pi(2, 2) \;\in\; [0.01397965,\; 0.01431408],
```

certifying that over **90.15%** (and at most 90.38%) of all bare/dormant cores have degree $`2`$.

5\. **Exact conditional neighbor-reciprocal expectations (`PROOF.md` Theorem 4):**
Using the exact stationary balance identities at $`J = (2, 1)`$ and $`(2, 2)`$,

```math
\mathbb E_\pi\!\left[\frac{1}{D_{\mathrm{parent}}} \;\Big|\; J\right] \;=\; \frac{\mu}{\pi(J)} - \mu - \frac{3}{2} \;\in\; [0.17650729,\; 0.17678896],
```

```math
\mathbb E_\pi\!\left[\frac{1}{D_{v_1}} + \frac{1}{D_{v_2}} \;\Big|\; (2, 2)\right] \;=\; \frac{\pi(J)}{2\,\pi(2, 2)} - \mu \;\in\; [0.39548710,\; 0.39662379].
```

---

## Snapshot SHA-256 Digests (`round2/outputs/marginals/`)

| ID | `obs{ID}_lower.bin` SHA-256 | `obs{ID}_upper.bin` SHA-256 |
|---:|---|---|
| 1 | `cdac40bcfd62bac4e75be4790838b7be014fa8aba0ca5f983ccd0b233f36897a` | `97f0151a8ef79cbf2f294a51e728c79a9e2b5e5045cd878b794789c794556ed1` |
| 2 | `9c8c5fc1c29c8e2552b6b8cb8677a07e70e255ba19d0a01b497901ff946c3024` | `52dd0dddcc373f57d4e5fb65b3e72a8c24f2bf4155a40f4478fb74871ebee8b9` |
| 3 | `cd2a23693205ef26e460b9a9be7d91cfc4c70eaf443b6f261db5a80076636e31` | `0e2ecae38744c5d3cdd9b7ef3854ae2788ec061cbcc4e9c165956bc29c8e7b1e` |
| 4 | `201c494cba42437c2694988e4d429faf208d06557678f2d522bad92dfc01254a` | `bab3a00aba721a49e2ce0fd5e1e316328cfcbf6042420ad6526daf03c91d3b61` |
| 5 | `898409d6d371275959f50bdbea3192b0225af61bebbece70253ab55ddeaa384d` | `a08e8e962d4c1902e1e635005794049b7f977e708ae7d0b9b87869127f3158cb` |
| 6 | `30b888d29e78cce80bf002a478a0132a66e3907801c3861bbc4be22e025c3074` | `78b2fe0869ec6cc92f4ee2dd93b7657bfb2195183d92616edb6dc8b7a2bcd5f1` |
| 7 | `03eb2f380a0be47cb2864acfeaa9463fc38255e261cfe6043948463c7ac02a28` | `7d15c18d869bc6aac50c9741beb3d273f0e0e9ca03086c8f68f74c0670b4dc32` |
| 8 | `dacdf9c26e1935445f1b80be6c8ce4fc54c30203d5d37bc4182183bd19db3afe` | `f2ac7370f2c74179ed92c5354331e276aac4ee4b44c0f020c2eab1c615002ab2` |

---

## Files & Architecture

- **`PROOF.md`**: Complete proof of Theorems 1–4 (signed unoriented 5-vertex coboundary across $`[\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}]`$, linear boundedness, Yule-dominated stopped Dynkin with Claim 3's exact coboundary $`R_{\psi, \mu}`$, 4-line bare-core tail Lyapunov function $`H(D, c)`$, and exact stationary balance identities at $`J`$ and $`(2, 2)`$).
- **`marginal_format.hpp`**: Binary `FT2_M5V1` header layout wrapping `../quintet-cert/quintet_format.hpp` (unmodified `b04cc50`) with fail-closed validation and exact integer observable indicator `m5_eval_obs`.
- **`verify_marginal.cpp`**: Standalone zero-floating-point 128-bit verifier (compiles cleanly with `-Wall -Wextra -Werror -pedantic` and UBSan/ASan). Evaluates all 8 motif tables and 7 cross-level envelopes at both $`\mu_{\mathrm{lo}}`$ and $`\mu_{\mathrm{hi}}`$ on the SHA-256-authenticated in-memory buffer.
- **`verify_marginal_tail.py`**: Tail-only ($`D \ge D_F`$) exact `fractions.Fraction` 16-panel 2D tensor-Bernstein verifier for the 4-term analytic tail $`w(q, x) = f_{p0} + f_{p1} q + f_{p2} q^2 + f_{p3} q x`$, with concrete-state oracle checks at $`D \in \{D_F, D_F+1, 3D_F+17\}`$.
- **`cert_marginal.py`**: Exact rational combiner, 35-tree concrete `Fraction` generator oracle, and 38-probe mutation suite (`--run-mutations`: 16 structural rejections + 22 adaptive-error widening probes per lower/upper pair).
- **`verify_all_marginals.py`**: Batch zero-float verifier for all 16 snapshots and all exact rational corollaries, producing `certificates_all.json`.
- **`build_marginal.cpp`**: Untrusted numerical proposal generator (reuses damped Bellman + Adam-like dual solver; NOT part of the trusted computing base).

---

## Replay

The listed snapshot paths and digests identify the historical inputs. The snapshots are not bundled. See the [release rebuild and verification commands](../README.md). A new proposal must pass exact verification; a different digest does not reproduce the original receipt.
