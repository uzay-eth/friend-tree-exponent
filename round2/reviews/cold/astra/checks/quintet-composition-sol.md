# GOLD1 Palm-joint composition: independent exact audit

**Pass on both original GOLD1 signs.** The saved multipliers certify **`[0.565704609764900, 0.565745628170000]`** on the original degree-24 Palm matrix. Separately, the old trials and audited stationary denominator ceilings give the slightly weaker outward interval **`[0.565704609691781, 0.565745628223974]`**.

## Ledger and normalization

Keep the exact `204990e` numerator `N₀ = bV(D) + G − cap_f − max(0,cap_g) − H`, including the signed upper `cap_f`. Keep its exact occurrence vector `U=W_hi` for the same original positive stationary denominator, and the same `η=π/μ` lift with birth right-hand side 1. Here `G` retains D18 `g`, stock8 FY and lower-only older-paying-neighbour FY. `H` retains the original host saving. Neither the physical matrix nor the denominator is reoptimized or replaced as an input.

| Witness | `746c825` old-donor FY, normalized per core | Audited H3, normalized per core |
|---|---|---|
| Lower | `80878351158765429162857/109951162777600000000000000000` (wide floor replaces the matched lower floor) | `677678669854000165705787542011/1151564514912178502778811842560000000` |
| Upper | `382522191322228678519/2748779069440000000000000000` (matched floor) | `236418900054226499063534772363/359863910910055782118378700800000000` |

The [old-donor audit](quintet-old-donor-sol.md) places its rebate on original `g` FY slots outside the D18 whole-edge mask and with donor degree `>5`, whereas the other paid FY slots require donor degree `≤5`. Ownership is by physical paying edge and donor slot. The [H3 audit](quintet-h3-cohort-sol.md), now frozen at `2705dac`, prices clipping on the same interpolated higher row while preserving its nonnegative cap, so the H3 recovery sits underneath the host cap already used in `204990e`. It contains the sibling-box baseline, gated direct component, and sibling increment; no unrestricted direct H3 floor is also paid. These two credits are on separate original family accounts even when the same physical vertices appear in both.

Let `C_d=2^40 F_d` and `C_h=2^40 F_h` in ticks per core. For positive trial `a` the exact cost is `N₀ + C_d + C_h − a U` on occurrence columns and zero on auxiliary columns. In the physical inequality the two new terms are per-core constants. In eta units they contribute `(C_d+C_h)η(1)=(C_d+C_h)/μ`. The stored floors already include their one birth conversion; neither receives another `LO`, `μ`, or local birth coordinate. The original singleton term still uses its actual mother's birth coordinate `b` once. This binds original GOLD1 only; the checker authenticates the existing GOLD1 receipt and payloads.

## Separate algebraic controls

Put `x=s(μ−λ)` with `s=+1` lower and `−1` upper. For the already certified trial `a₀`, the physical inequality gives `xM≥a₀ M+C`, where `C` selects the old-donor, H3, or combined credit, and `0<M≤M_hi` from the earlier [stationary audit](quintet-stationary-sol.md). Thus `x≥a₀+C/M_hi`. This discards the positive old eta margin and does not use the new saved vectors. Use the **exact** upper ticks from `host24.log`, not the printed scaled decimals. For each control the script reparses the old receipt, checks all its original columns, reproduces its exact signed margin, and checks `(a−a₀)M_hi=C` rationally.

| Credit added | Lower endpoint, rounded down | Upper endpoint, rounded up |
|---|---:|---:|
| Old-donor only | `0.565704134004711` | `0.565746158933064` |
| H3 only | `0.565704015100869` | `0.565745740640311` |
| Combined | `0.565704609691781` | `0.565745628223974` |

As a parity control, using the slightly larger exact ceilings in the separately audited `joint24.log` gives the producer's combined outward interval `[0.565704609691780, 0.565745628223974]`. The saved trials improve strictly on both controls, in exact arithmetic.

## Exact checker

The checker (`quintet-composition-sol.py`) rebuilds the existing independent Palm inputs, host envelopes and 144-entry kernel. It decodes each saved integer multiplier as `integer * 10^6 / 10^12` in ticks and computes `max(0, max(Aᵀy−cost))` separately for each of the five column kinds. The uniform repair subtracts the sum of those defects divided by original `LO` from `y_birth`; the Palm-specific repair instead bounds Palm slack mass by one. A strictly positive exact eta margin certifies the trial using the common positive denominator. It also checks the composed costs by the pinned Palm audit's separate cost builder and compares both residual calculations. No floating objective, fresh multiplier proposal, or refit is used.

Executed from `[local path]` on one niced core (exit 0):

```sh
taskset -c 11 nice -n 19 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  round2/.venv/bin/python -Bu round2/checks/quintet-composition-sol.py
```

For each side the replay checked 613 saved multipliers against all 5,554 original columns. Positive defects were confined to occurrence and edge columns; the two repairs coincided. Exact results:

| Witness | Inward trial | Exact uniform-repaired eta margin (ticks) |
|---|---|---|
| Lower | `22397649/10000000000000` | `380122186866156533791375684409394182641978321/561372275736681942195006189484800000000000000` (≈ `0.677130316718`) |
| Upper | `265183/100000000000` | `623642384946756036867499734589754495934821221/1157830318706906505777200265812400000000000000` (≈ `0.538630207614`) |

The saved `composition/lower.json` and `composition/upper.json` are pinned by their SHA-256s in the audit checker. The producer's `composition/check.py` was inspected for its input convention but not imported. The original `204990e` inputs and the relevant `746c825` receipts matched their committed files.
