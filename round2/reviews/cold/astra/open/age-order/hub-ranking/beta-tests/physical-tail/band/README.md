# Cold audit: continuous tail band

**HOLDS.** For every old tree and every older label `b≥20`, the two upcoming births have strict conditional tail order for all `x∈[1/3,1)`. Averaging and splicing the paired Fable receipts for `2≤b≤19` gives strict annealed order on that interval for every adjacent pair starting at `2→3`. The founders tie. At `x=1` the strict tails vanish. No defect found.

The audited writer blobs, under `third-physical/band/`, are `README.md` = `cb7af0fb9f25c802d7c240a8a9abd8dc6ae853e6` and `check.py` = `e2aa75a98bf62981a91e33c3b3664e2f55d248c6`. The replay (`check.py`) verifies these before loading code, as well as the frozen endpoint ledger, both first-family splice checkers and the prior Fable audit manifest. The [endpoint audit](../README.md) is a separate freeze.

## Every threshold at label 20

Common completion is size 21, from an old tree of size 19. In each closed box `[l,h]`, write `q=1−x`, `q_hi=1−l`, `q_lo=1−h`. The pending normalized beta cap increases with `q`, using the positive power `q^(1−r)` and negative integrand exponent `r−1`, for `r=1/k`. For promotion at the actual next size the discounted integrand is proportional to `q^(1−r)(1−qv)^r`, with derivative of sign `1−r−qv`. If `q_hi≤1−r`, use the high endpoint for both factors. Otherwise the power at `q_hi` times the integral at `q_lo` is an upper bound. This second branch is needed for `k=2` on the band's boxes; simply taking both factors at `q_hi` would not be justified. Pending truncation adds a geometric remainder, while all omitted nonconstant promoted terms are negative. Rounding is upward, including the power and the Bellman ceiling. The discounted promotion term enters once at rate `1/k`; the nonpromotion terms use `q_hi` once. Taking both gain-rate corners and the minimum with the pending cap preserves the uniform upper bound by backward induction.

On the serial side the cutoff atom in the scaled charged-family size is at `bθ/(1−θ)>t`, where `t=b(θ−x)/(1−θ)`. The uncapped law can therefore be used for this payoff. The conditional shifted survival of `t` and the relevant family law are independent of the tested threshold after the factor `q^b` is removed. The exact loss charge telescopes: `bΣ_(j>f)1/((b+j)(j+1)) = b/(b−1) Σ_(k=f+2)^(f+b)1/k`. Taking the right uniform-bin endpoint decreases `f`, increasing both the charge and the stochastically charged share. Left `t` endpoints, right family-size endpoints with a downward CDF and residual mass at infinity, and upward payoff charges give a genuine lower sum; the replay proves it exceeds `.2124` at `b=20`. At `b=21..39` the previously replayed uniform coefficient exceeds `.211`; this coefficient is also threshold-independent after normalization.

The signed ledger therefore compares the true tails for each threshold in each box. At `b=20` the exact replay gives chain reserve on `[1/3,41/120]`, `[41/120,7/20]`, `[7/20,3/8]`, and `[3/8,2/5]`; it gives sitewise reserve on `[2/5,1/2]`. The first four use the frozen retained-core checks, strict star checks and positive skeleton reserve. The last uses the positive sitewise surplus. Shared endpoints are included; nothing is interpolated from values of the actual gap at a finite set of thresholds.

## Splice and quantifiers

- **`b≥40`:** The old `first-sender/joint/` statement itself says `x≥1/3`, `x<1`. Its quantile envelope separately handles `b=40` and all `b≥41`, and the shifted endpoint survival is bounded below by the same exponential for all of them. Thus the lower coefficient `261/1250` is uniform in `b` and `x`. After dividing by the positive serial weight and `q^b`, the paired debit is at most `sqrt((1−x)/(π(b+1/2)x)) · (4/3+(1/2)((1−x)/(bx))^(1/6))`. Both factors decrease in either `x` or `b`; the strict rational corner test at `(40,1/3)` proves the band and beyond. This is monotonicity of the comparison bounds, not of the physical gap. The independently rerun coefficient at mesh `(256,24)` is `2138962243701651396480805099/10240000000000000000000000000`.
- **`20≤b≤39`:** The band replay covers `[1/3,1/2]`. The frozen `first-sender/refined/` result and its rerun `b=20, x=1/2` corner cover `[1/2,1)` for all `b≥20`, overlapping at `1/2`.
- **`2≤b≤19`:** The paired archived Fable receipts at `88e13f36b39162f6ee60561a2b465a982b34af02` give exact saturation and positive universal core-to-leaf mass on all 18 rows, ending at `19→20`. The prior audited equal-size physical ordered coupling gives the weak order; the finite residual-forcing cylinder on that positive mass gives strictness at each `0<x<1`. The receipts are annealed and share rooted-law propagation lineage. Audit files and their receipt manifest were authenticated again; terminal transport runs were not repeated. No use of the single `20→21` receipt.

## Replay evidence

Run from the Astra repo root:

```sh
nice -n 19 timeout 65s python3 -IB round2/open/age-order/hub-ranking/beta-tests/physical-tail/band/check.py
nice -n 19 timeout 55s python3 -IB round2/open/age-order/global-transport/annealed/extension/first-sender/refined/check.py
sha256sum --strict -c round2/attacks/fable5/release.sha256
```

The first command passed all five `b=20` boxes, every label `21..39`, and the `b≥40` corner. The second printed `PASS exact red-family conditional weights and cap atom`, `PASS b>=20, x>=1/2: refined coefficient > 41/200`, and `PASS b>=21, x>=1/2: refined coefficient > 507/2500`. All seven entries in the last command returned `OK`. The writer's own checker also passed in two invocations (`--start 20 --stop 21` and `--start 21 --stop 40`).
