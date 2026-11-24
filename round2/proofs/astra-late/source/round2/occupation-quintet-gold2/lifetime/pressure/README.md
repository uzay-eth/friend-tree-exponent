> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# Gold2 singleton pressure at `99/175`

The original **lower witness** singleton reward has certified stationary bounds **`[331898.049402897, 1518437.847419405]` ticks/core** using the entire Gold2 exponent prior. The joint marked birth and host trial leaves **at most `493630.005011907` ticks/core to certify from disjoint edges**, about 124,245 less than the ~617,875 starting gap (20%). The lower price replaces the old D430 full-singleton floor; the new joint gap is before adding any edge credit. At `mu=99/175` itself the corresponding conditional bounds are `[331951.248889791,1518368.820765867]`, with gap `493548.877732559`. Use the uniform-prior result for an unconditional Gold2 ledger.

These ceilings provide **no upper barrier to closure**: even the original all-stock threshold `1114713.791608332` is below the tighter ceiling. The D24 potential-only sufficient threshold is `842625.934035200`; falling short of either sufficient threshold is not an impossibility result. Truly ruling out the singleton alone would require an upper bound for `pi F - (99/175-lambda) pi M`, with the actual `M`, below the per-core error (or an upper bound and a valid lower bound on `pi M` giving the same conclusion). None of those comparisons rules out the singleton alone here. A fortiori this does not settle whether edges plus singleton can close the gap.

## Physical states and certificate

- At birth the focal vertex is `(D,c)=(2,1)`. A donor's tracked `(d,h)` includes its actual degree and core arity, hence leaf stock `d-h`. The mother is marked by her **post-promotion** `(m,h)`; the birth flux is `(m-h_pre)/m` at the pre-promotion census and assigns mark `(m,h_pre+1)`. The original singleton reward owns its own prescribed promotion terms once; a mother mark is information about the recipient's future, not another reward.
- Retain mother's census through degree 5 and the two most recent real children through degree 3 in the `D<=700,c<=12` focal window. New promotions insert a real `(2,1)` child and release the oldest of the tracked pair when full. Released or out-of-range donors remain physical neighbors. Overflow mother and child sends are bounded by `1/6` and `1/4`, respectively. Live sends to the focal have rates `1/m` and `1/d_i`; the focal feeds each live donor at `1/D`. For a live donor, gain has known rate `d-h+1/D`, with up to `(h-1)/2` from unobserved neighbors, and its promotion rate is `(d-h)/d`. No transition on an unobserved neighbor changes a retained mark. The finite recursion is acyclic in the underlying monotone censuses; exits use the previous all-degree comparison.
- `F=a+b*rho` is the full original signed **singleton** score at Gold2's nominal `lambda`, with the original cap. On centered rows the stored center is used only in the charge; all physical gains use `L+rho`. The same unobserved receive appears in the reward and focal transition. Integer Bellman inequalities are checked against the worst signed affine residual on every state; 17,262 independent explicit corner checks per 2,800,083-state solve. Upper and lower comparisons use directed grid rounding and the respective ends of the prior. The all-state upper uses the earlier analytic reciprocal-degree continuation and separate high-arity overflow.
- `pi F = mu beta V`, where `beta` is the true birth law and `V` the discounted future singleton reward. One independent mixture bounds each maternal mark's mass by its actual discounted birth flux from a virgin, using the stationary identity `beta K_mu=beta` and normalization. No assumption that every row of `K_mu` has mass one. The sharper sided D24 calculation attaches the mark value once to each actual promotion, then uses normalized exit/birth flows and reciprocal edge constraints. The D24 host lift and the exact saved host input are reused by hash; the prior D24 `joint_rebate` is **not** reused because its objective includes a different singleton residual.

## Trial objective and accounting

For the lower witness use `delta=99/175-lambda=53/87500000`. With positive comparison denominator `M0=original m_max`, original allowance `H` for higher rows, and independently valid local bounds `M<=U` and `higher<=h`, price one joint objective:

`mu E_beta[V] + pi(H-h + delta*(M0-U)) - err_per_core - delta*M0`.

This objective uses `F` exactly once, owns no old-edge reward and retains the original higher-host and potential comparison bounds. Disjoint edge credits may be added once. D430, the broad reduced-f witness in `host/joint24.json`, the marked bound on its own, and the alternative finite/Dynkin prices are replacements or diagnostics; do not sum them. The exact finite residual uses `F + L_phys V - sV + b V(2,1)` for a degree-truncated test with the actual same `rho=u+v` in `F` and `L_phys`; at `s=99/175` its correction to the finite stopping functional integrates to zero on the physical law. It was weaker (the short-window replay gave 714,446.992238821 remaining ticks for the residual and 698,017.060625294 for the finite stopping objective). Neither pricing supplies all exact higher-row drifts; those enter via the existing disjoint host bounds.

## Old/high arity diagnostic

`c` counts the mother and children, not degree. At the trial the first physical promotion into `c=4` has standing flux `[0.055713371,0.066271400]` per core. Following the exact degree transitions to that entrance gives **at least `87200.472712784` ticks/core** in the `c>=4` band, and an upper envelope `1996550.331727579`, weaker than the full marked ceiling. The D24 overflow-pooled calculation is weaker still: it can place the overflow flux at an entrance value over 52 million lifetime ticks. This demonstrates that the high-stock path contributes but neither its flux nor the present relaxed continuation proves enough to close the remaining gap. In particular, the `c>=4` contribution is already part of the singleton and is not additional credit. A tighter test should constrain actual receive/degree-age occupation in the high-degree continuation rather than increasing an independent birth mark table alone.

## Replay

From `.`, with the authenticated Gold2 inputs in place (the SciPy solver proposes a dual; exact rational column checks with outward repair certify it):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 200s round2/.work/venv/bin/python -Bu round2/occupation-quintet-gold2/lifetime/pressure/replay.py --mother 5 --tip 3 --live-cut 700 --joint --stock-flux
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 200s round2/.work/venv/bin/python -Bu round2/occupation-quintet-gold2/lifetime/pressure/replay.py --mother 5 --tip 3 --live-cut 700 --joint --stock-flux --upper-control
```

Use `--trial-discount` for fixed trial occupation, and `--mother 4 --live-cut 25 --stock-flux --trial-discount` for a short high-arity replay (with `--joint` for the residual; with `--upper-control` for the ceiling).
