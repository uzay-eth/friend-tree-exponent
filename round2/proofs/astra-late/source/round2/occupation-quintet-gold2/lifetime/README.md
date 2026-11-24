> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# Gold2 original singleton repricing

With the original Gold2 witness alone, the exact conditioned singleton lower prices give the outward interval **`[0.565713784375142, 0.565736707276694]`**, conditional on the accepted finite certificate and physical stable-age law. The lower end remains below `99/175`. No Gold1 binary, receipt, payment, denominator or exponent prior is used.

| Stationary original singleton reward, ticks/core | Lower witness | Upper witness |
|---|---:|---:|
| Unconditioned floor | 91786.442259914 | 187969.086680319 |
| Conditioned floor | **207990.215716345** | **296272.606830282** |
| Conditioned all-stock ceiling | **1797874.252038426** | **1904610.996806601** |
| Unconditioned all-stock ceiling | 2073286.323318540 | 2249760.082417987 |

The exact signed original singleton score is `cap_f - sigma*(drift_f-lambda*f)` in ticks, `sigma=-1` lower and `+1` upper, with Gold2 caps `55/2359`. Centered rows use the actual stored center in the **charge**, while focal gains always use the physical `L+rho`. A promotion has rate `L/D`, uses the post-promotion singleton and clamped newborn edge, and owns the newborn `f+g` once. Sparse row extraction uses the independently read Gold2 header, including `DF=4000`, `d2=252`, `K2=22`, and 20 FY knots. The signed affine reward and gain continuation share one receive variable. Cross-multiplied integer Bellman inequalities are checked at both receive endpoints for every row; direct rational quotients additionally test the key layout/analytic seams. Lower pricing uses the original exact Gold2 upper discount, times its original exact lower birth normalizer. Ceiling pricing reverses the discount and multiplier. Both use a true nonfounder newborn `(2,1)`.

The full lower solve uses `D<=7999,c<=30` and a lower analytic `h[c]/D` continuation; the upper solve uses `D<=3999,c<=70` and an upper `a[c]+h[c]/D` continuation. The exact cleared residual polynomials at fixed arity are bounded by rational Bernstein coefficients on the whole reciprocal tail. `P(0)>=0` permits the lower bound on `P/x` obtained by replacing `P(0)/x` by its cut value. Gain decreases `h[c]/D` at relative rate at most one, and at least `(Dcut-c)/(Dcut+1)` on the analytic suffix. The real promotion rate lies in `[1-c/Dcut,1]`. For the ceiling the skipped `c>70` part costs at most `130000*S*(1+LO)^(-70)`: 70 promotions are dominated by a rate-one Poisson process at exponential age; `130000*S` globally bounds the pointwise score since before the analytic seam gain is at most `3999` with weight difference at most `32*S`, all other weight terms have bounded magnitude, and the analytic gain increment has an extra reciprocal-degree factor. The ceiling therefore covers the whole physical score, not just the retained arity mask.

The conditioned comparison uses the same physical rate construction as the already audited live-neighbor algorithm, with Gold2 inputs only. For `D<=430,c<=18` it tracks the live mother at degree two until she grows (thereafter receive at most `1/3`) and the latest child at exact degree up to 3, including its core arity. The latest child is replaced only at an actual promotion, keeping the old child as a physical untracked donor. Both directional sends along a tracked edge use the actual live degrees. At the conditional exit, use the unconditioned value at the true census. The physical fraction of births from degree-two mothers is enclosed by `[(1/2)/(HI+2),(1/2)/(LO+3/2)]`, using the actual competing rate from the mother's virgin residence, and no founder-mother entrance credit. Multiply the unconditional mixed lifetime by the birth normalizer once. This refines the same singleton account rather than adding another account. Each conditioned solve used 86,742 states and 750 direct multi-port corner checks.

For endpoint conversion subtract original per-core errors `19217/20062` once from the certified stationary reward. The net is positive on both sides, so divide by Gold2's original `m_max=1808603193693/1808568794159`. The inward gains from original exact endpoints are at least `0.000000298041786/0.000000279401751`. Under this same original-denominator conversion, the conditioned ceiling bounds possible inward gains by `0.000001177109148/0.000001168689677`; the lower ceiling does not rule out `99/175`, but the computed floor does not certify it.

Replay from `.`, after both lines of the private Gold2 `hashes.log` say `OK`:

```sh
taskset -c 8 nice -n 17 timeout 100s python3 -Bu round2/occupation-quintet-gold2/lifetime/check.py
taskset -c 8 nice -n 17 timeout 85s python3 -Bu round2/occupation-quintet-gold2/lifetime/check.py --degree 3999 --arity 70 --upper-control
taskset -c 8 nice -n 17 timeout 120s python3 -Bu round2/occupation-quintet-gold2/lifetime/check.py --conditioned 430
taskset -c 8 nice -n 17 timeout 120s python3 -Bu round2/occupation-quintet-gold2/lifetime/check.py --degree 3999 --arity 70 --upper-control --conditioned 430
```

All four commands completed both witnesses. Broad solve counts are 239,534 and 277,514 per witness for the lower and upper price respectively. Exact virgin and mixed lifetime fractions are printed; the table uses directed decimals. No board entries or prior files were edited.
