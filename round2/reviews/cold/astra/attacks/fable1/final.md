# Fable Claim 1 — independent official-final audit

**VERDICT: HOLDS** for the stated finite-window estimates and the approximate confidence bounds **conditional** on the stated assumptions. No unconditional limiting-exponent digits follow. Use `8.83e−11` for the outward-rounded (C+) bound in the final archive; the board shorthand `8.8e−11` is rounded down and should not be read as a strict outward limit.

## Freeze and sample separation

The official board labels its final posting 15:51 UTC on 24 November 2026 (board file mtime 15:51:41 UTC). Its primary number and 19,000-run count match the genuine final data at `ft2-fable@470ee9f45cbc679bd257cbb1599eee890be3166e`; the board claim still links the older `9828f1c` form pin. `final.json` (SHA-256 `8f752cbe407758fc5648648d9e760b15eb03620aa9d12b6be4e2ace3e18ad48a`) pins the final code and artifact hashes, including `moments.npz` at `94f550688381015973d10e083cedd3546476d21d824d4d9e9e23841a839e3109`, final output at `1815236800f217cd44defb80e781393a9c38dd4e76dae4a6fa246608d6a709f5`, the frozen executable and all five weights. The checker compares the live committed files byte-for-byte to the Git objects at this pin.

The bounded raw replay pins **567 files, 752,840 distinct stream/run keys, 321,245,866 bytes**. It compares first and second raw moments of each retained file with the corresponding archive entries. There is no complete eligible file missing at inspection. The sole incomplete late file had 18 of 60 rows; its six completed trajectories are not in the 2,840-run **complete-file** cohort. Phase 1 contributes 3,000 main, 75,000 old rung and 7,500 middle rung runs, all excluded from the nominated estimates. The three phase-2 families are chunks ≥500 on streams 9100/9200/9300; late and block replication use streams 9500/9210. The A1 subset (5,500 runs, main chunks ≥550) is disjoint from the exploratory main portion, but is part of the primary sample and must not be treated as independent of it. The posted launch chronology puts the independent A2 late cohort after the analysis rule and weights were frozen. The checker finds no duplicate derived RNG pre-seeds among the inspected runs; the file/key/weight checks do not themselves attest to historical invocations or stochastic independence. The training streams 811–814 and the earlier window choice are separated by design. The archived early/IR-0 controls were recomputed from moments only; no raw replay of their millions of trajectories or new simulation was needed.

## Raw reconstruction

Use the ratio of sums of direct physical endpoint differences. For per-run differences `(u_i,d_i)`, the centered influence is `(u_i−R̂ d_i)/Σd`; the run SE uses `sqrt[n/(n−1) Σ(ψ_i−mean ψ)²]`. Independently delete each **whole file** from the numerator and denominator and use the exact leave-one-block ratios for the cluster SE. Subtract paired influences/leave-outs before taking the error; use quadrature across independent samples.

| Quantity | Runs / completed files | Estimate | Run SE | Leave-file SE |
|---|---:|---:|---:|---:|
| Primary full 16→24 | 19,000 / 190 | 0.5657236605826 | 6.15109e−9 | 5.94571e−9 |
| Full 12→16 | 260,000 / 52 | 0.5657237030344 | 1.83637e−8 | 1.26751e−8 |
| Full 16→20 | 25,500 / 51 | 0.5657236423605 | 1.88141e−8 | 1.86106e−8 |
| Late full 24→28 | 2,840 / 142 | 0.5657236688045 | 5.98171e−9 | 6.26079e−9 |
| Weight `w` (main) | 19,000 / 190 | 0.0942442731119 | 1.56447e−6 | 1.65765e−6 |
| Fresh block full 12→16 | 360,000 / 72 | 0.5657237008025 | 1.56372e−8 | 1.63525e−8 |

Late minus primary is `+8.22191e−9 ± 8.58002e−9`; paired late 24→28 minus 20→24 is `−3.43559e−9 ± 2.08720e−8`. The inverse-variance summary is `0.5657236648083 ± 4.28834e−9`: it combines two different finite-window targets and its SE is not an error bar for `μ`. Registered full−ANPQ controls reproduce `−4.63340e−9 ± 2.90284e−8` (A1, −0.16 SE) and `+4.90767e−8 ± 4.17387e−8` (A2, +1.18 SE); their leave-file SEs are `2.66976e−8` and `4.70223e−8`. Independently replaying the registered residual block permutations on the fixed rows gave original `T=0.431880, p=.0010`, fresh `T=1.016630, p=.6040`. This reproduces non-replication on the fresh stream; it does not prove a cause for the original anomaly or certify independence. The run-level SE is used as specified and is not always larger than the leave-file SE (notably late).

## Bias algebra and limits of the verdict

Dynkin's formula gives `F′=E num` for `F=E den` by finite-time Yule moments and localization of the fixed feature map. It does not give a rate of cancellation for `F′−μF`. For positive **population** increments write `B=R(16,24)−μ=w b16+(1−w)b20`. If the two proposed contractions hold with factor `r`, then `|B|≤q|b12|` for `q=r[w+(1−w)r]`. Three contractions for the late bias give `q=r³`. From `|x|≤q|y|`, `q<1`, the triangle inequality gives `|x|≤q/(1−q)|y−x|`. The coefficients increase in `r`, and the primary also increases in `w`; use `μ≥.565` and the one-sided sample upper limit for `w` (`.094247339472`, below the claimed `.094249` cap). The old-minus-primary and old-minus-late contrasts are `+4.24517746e−8 ± 1.93665473e−8` and `+3.42298604e−8 ± 1.93134190e−8`; they share the old estimate (correlation `.901595`). The two-sided contrast multiplier `2.2414` and one-sided `w` multiplier `1.96` give the primary approximate coverage by Bonferroni for each line; late needs no `w` limit. These are separate, not simultaneous, bounds.

| Window | Premises | Coefficients (outward) | Conditional bounds (outward) |
|---|---|---|---|
| Primary | (C), (C′) | `.020094`, `.164649` | `1.73e−9`, `1.42e−8` |
| Late | (C+), (C′+) | `.001138`, `.043391` | `8.83e−11`, `3.37e−9` |

Positive sample increments do not prove positive population increments. Approximate normality, correct code and independent runs remain statistical assumptions. The early ladders, known-truth controls, kernel scan, and final-window agreement do not rule out a slow correction or cancellation in this observable; in particular an asymptotic spectral envelope does not imply the required finite-window contractions. The replay therefore confirms S1–S4 and the numerical reading of K1–K4, not the truth of (C), (C′), (C+), (C′+), nor an unconditional interval for `μ`.

Replay from Astra root (no Fable reducer or production executable):

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 round2/.venv/bin/python -IB round2/attacks/fable1/final_audit.py
```
