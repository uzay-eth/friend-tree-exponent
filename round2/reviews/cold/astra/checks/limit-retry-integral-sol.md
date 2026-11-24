# Q4 cold audit: integrated cut retry

**Holds.** Reviewed final `retry-integral/README.md` Git blob `59237b4725f4cdcda363727683589e1c286f40aa` and `check.py` blob `fa14ebf85c2ffcebc00373629915ac89940e286b`, with the previously audited planted angular translation estimate. The full-law set and `L¹` translation exponent is `σ=747/11494`; the full density is weak `L^(11494/10747)` and strong below that endpoint. The fractional Sobolev range is `1≤p<11494/10747`, `0<s<1/p−10747/11494`. No bounded-density result is claimed.

## Clock and age integration

The rate-two first race pays for either winner; an anchor win ends its attempt immediately *after* the physical ring. After a leaf win the anchor has at least two leaves and at most one core neighbour. Only the anchor can promote any of its leaves. Its rings choosing the core neighbour do not decrease its stock or change its core degree. Fresh uniform marks can put the leaf choices first, so subsequent anchor rings with mark at most `2/3` give a genuine conditional `Exp(2/3)` upper bound on the wait for preparation. This is an upper bound on that wait, not a statement that its duration is independent of the setup state. At preparation split fresh independent lower `1,1/2` stocks; conditioning there gives their defective failure law, without conditioning the physical future on success.

For one family `q(u)=P(u<T_ext<∞)=w/[2(2−w)]`, `w=exp(−u/2)`. The minimum of two has defective tail `(1/2+q)^2−1/4=q+q²`, not the full survival probability of that minimum. Expand the coefficients `(j+1)/2^(j+2)`, sum the first eight and place the remaining positive coefficient mass at the ninth decay rate for an upper bound. The exact `ρ` calculation includes the first-race factor `1/(2−r)` on **both** failure branches, then the additional `(2/3)/(2/3−r)` and defective extinction factor on the leaf-win branch only. At `r=747/10000` the result is `ρ<99991/100000`. The completed-preparation exponential budget is bounded by the same two setup factors without the extinction factor.

Let `R_j` be reached restart times and `U_j` completed preparation times (possibly infinite). They are stopping times; `S`, the preparation of the trial that never fails, is not. By conditioning at the real restart `E[e^(rR_j); R_j<∞]≤ρ^(j−1)` and by conditioning at each reached start `Σ_j E[e^(rU_j); U_j<∞]≤M/(1−ρ)`. Infinitely many failures have probability zero by the same restart bounds at zero. Thus `P(S>t)≤e^(−rt)Σ_j E[e^(rU_j); U_j<∞]`. This inequality uses no optional stopping at `S`.

At deterministic `t` the physical selector maximizes the current minimum incident stock of a core edge, independent of unobserved graphical futures. Conditional crossing cost is at most `min(1,4/(1+B_t))`, including noneligibility when `B_t=0`. On `S≤t` the two witness stocks are positive and their minimum is at most `B_t`; `1/(1+B_t)` is bounded by the sum of their two reciprocals. Extend lower families privately after abandonment for the finite-age read; those extensions do not enter restart decisions or assert any physical domination for a failed trial. Sum over **all** reached setups before applying their fresh conditional laws. For the geometric parameter `q` of the positive population, the sharper exact optional identity is

```
E[1/(1+Z_u) | Z_u>0] = q[-log q−(1−q)]/(1−q)^2,
```

with limit `1/2` at `q=1`. Its elementary upper by `E[1/Z_u | Z_u>0]` yields `C(1+u)e^(−u/2)`. Since `r<1/2`, bound that by `C_r e^(−ru)` at age `u=t−U_j`. The random age is known when conditioning at the setup; the subsequent family law on the graphical space is the fresh law. The sum costs `C_r e^(−rt)Σ E[e^(rU_j); U_j<∞]`. This preserves the full exponent `r`: there is no deterministic stock cutoff, age split, or product of future crossing probabilities. The factor `1+u` would matter at `r=1/2`; that endpoint is not used.

## Shape weight, exponents and fold test

On the selected deterministic-time no-cross sublaw, keep the genuine terminal shape weight `L` in the fixed-history kernel average. Selection uses ordinary physical history, so neither the setup witness nor an earlier failure mask conditions the Beta clock. For `N=|T_t|` and `A=a+b`, the angular kernel costs `C√N A^(−1/2)√δ`; transfer this **with** `L` to the retained physical law using `X≤A`. Cauchy–Schwarz costs at most `C e^(t/2)` using the unconditional Yule mean and `E X^(−1)<∞`. Removed mass and its translate cost at most twice the removed mass, not a bound on its density. Balancing gives `r/(1+2r)=747/11494`; the reported gain `14561/798833` over `13/278` is exact. The all-set estimate uses the same physical transfer; interval-only bounds would not suffice to test density superlevel sets. Interpolate the translate **norm** and then take its `p`th power in the Sobolev integral, giving the stated strict range. No endpoint strong `L^p` or Sobolev bound follows.

For two positive coefficients, the weighted fold-potential condition is necessary as stated. A positive pathwise fold coefficient and a measurable one-sided radius give a nonnegative lower bound on that conditional density. Tonelli provides the lower potential bound on the subdensity almost everywhere. Each potential integrand extended by zero outside its open one-sided interval is lower semicontinuous (including at the excluded boundaries), so Fatou and openness upgrade any essential upper bound on an open interval to a pointwise upper bound there. Choose disjoint half-open dyadic shells inside the strict radius; their weighted sum is bounded by the potential. The shifted square-root location density gives a concrete logarithmic divergence even though its location measure satisfies a half-power interval bound. This only demonstrates a limitation of the inference, not a failure of boundedness of the physical law.

## Replay

`nice -n 10 timeout 12s python3 -IB round2/open/limit-law/smoothing/retry-integral/check.py` reproduced the writer's exact `ρ` envelope and exponent, 16 geometric envelopes and 32 critical shells.

```text
nice -n 10 timeout 12s python3 -IB round2/checks/limit-retry-integral-sol.py
PASS exact defective-clock and full-law exponent arithmetic
PASS 24 geometric reciprocal partial sums and 15 integrated-age controls
PASS 40 critical fold shells and planted nested-mask error
```
