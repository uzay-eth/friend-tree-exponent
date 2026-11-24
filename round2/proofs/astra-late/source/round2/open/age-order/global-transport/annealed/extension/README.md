> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../README.md#claim-26) for exact scope and citations.

# A paired-edge threshold bridge and an annealed occupation criterion

Write `b=n+1` for the older birth label and `m=b+1` for the common post-placement size. The paired-edge comparison below proves strict adjacent order on **every** pre-placement tree for `bx≥141`, improving the inherited [232 bound](../../../../../../../../../reviews/cold/astra/checks/preplacement-256-sol.md). Retaining the factor lost in the birth-scale cutoff also gives strict order for **every `b≥20` and `10/11≤x<1`**. Averaging gives unconditional comparisons on the entire finite-label gap in this fixed upper-threshold interval; combined with Fable's finite certificates, this covers all labels there (labels 1 and 2 have equal laws). The proofs of both new bands are conditional before the two births. A separate criterion below is genuinely annealed and could lower the band using a degree/stock occupation estimate. It does not yet certify that occupation estimate on the unresolved labels. In particular there is no new result at hubhood.

## 1. Split the adverse parallel weight by its actual parent

Use `r_u=Σ_(v~u) 1/D_v`, `w_u=r_u/(D_u+1)`, `s_u=w_u/D_u`, `A=Σw_u`. In the exact signed comparison, only the endpoint on the larger-`s` side of an edge needs an upper bound; keep the coefficient `|s_u−s_v|`, not `s_u+s_v`. Write `X(T)` for the sum of these edge coefficients on which the larger-`s` endpoint is an **old leaf**. Total coefficients over all edges are at most `A`. The parallel tip that attaches at an old core has a post-placement parent of degree at least three; it is bounded by `Q_{m,1/3}(x)` in the locked urn. A parallel tip attaching at an old leaf has cap `Q_{m,1/2}(x)`. Every serial tip has the latter cap. This is about the degree of the recipient, not the degree of the initiator, and remains true for a bare old core.

There is a deterministic cap `X≤A/3`. At an old core `o` of degree `D` with `L` old leaves, write `A_l=L/(2D)` for the serial mass of its old leaves. Its core serial mass is at least `L/(D+1)=2D A_l/(D+1)`. The old leaves' joint contribution to `X` is

```math
 X_o=L(1/(2D)-s_o)_+
 \le A_l(1-2L/(D+1))_+
 \le A_l(D-1)/(D+1).
```

The last inequality applies when `L≥1`; at `L=0` the contribution is zero. Thus `3X_o≤A_l+w_o`, by direct comparison for `D≥2`, and sum over cores. At size two the equal endpoint weights give `X=0`. A refined occupation majorant is `X≤Σ_o [L/(2D)](1-2L/(D+1))_+`.

Let `c_b` be the explicit serial-middle lower coefficient from the protected physical initiator-colour prefix in the [232 proof](../../../../../../../../../reviews/cold/astra/checks/preplacement-256-sol.md#adverse-weights-and-numerical-bound). That proof gives `h_M≥c_b(1-x)^b` for `b≥20` too: the prefix and Markov argument work unchanged, with smaller coefficient at small `b`. Explicitly put `t=5/6`, `l=Σ_(j=1)^11 (−1)^(j+1)(5/11)^j/j`, `e=Σ_(j=0)^11 (−t)^j/j!` and use

```math
 c_b=\frac{5(b-1)}{16b+5}(1-l/t)e\Big(1-\frac{t^2}{2(b-t)}\Big).
```

All three varying factors are nonnegative and the coefficient increases for `b≥20`. The inherited high-beta upper is `Q_{m,1/2}(x)≤U=(1-x)^b sqrt(1-x)/sqrt(π_math b x)`. Direct density comparison on `z≥x` gives

```math
 Q_{m,1/3}(x)\le Q_{m,1/2}(x)\min(1,C_b (bx)^{-1/6}),
 \quad C_b=\frac{\sqrt{\pi_{\rm math}}}{\Gamma(1/3)}\frac{b+1}{b+1/2}.
```

Indeed the density ratio has the factor `((1-z)/z)^(1/6)` and `Γ(m−1/2)/Γ(m−1/3)≤(m−1/3)^(5/6)/(m−1/2)` by Jensen on a gamma variable of shape `m−1/3`. The other gamma ratio displayed above is absolute: integrating the odd Taylor lower polynomial for `exp(-t)` over `[0,1]` with weight `t^(-2/3)` gives `Γ(1/3)≥2.4223`; `sqrt(π_math)<9/5`. The replay gives `C_b≤3/4` for `b≥120` and `C_b≤4/5` for `b≥20`.

Consequently, with `y=bx`, `θ=min(1,C_b y^{-1/6})`, the signed numerator is at least

```math
 A c_b(1-x)^b-U\{(1+\theta)A+(1-\theta)X\}.\qquad\text{(1)}
```

For `y≥141`, use `X/A≤1/3`, drop `sqrt(1-x)`, and test at `b=y=141` using monotonicity (increase in `c_b sqrt(y)` and decrease in `θ`). For `x≥10/11` with `b≥20`, use `y≥200/11>(3/2)^6`, so `θ≤(4/5)(2/3)`; keep `sqrt(1-x)` and test at `b=20,x=10/11`. Both tests are strict by exact arithmetic. This bridge uses the same physical serial lower as the 232 theorem; its gain is the paired adverse debit, not a second likelihood improvement.

## 2. Genuine annealed reduction and what one old jump buys

At deterministic size `n`, set `ρ_n=E X(T_n)/E A(T_n)`. Taking expectations in (1) shows that `1/3` can be replaced by `ρ_n` in the sign test. This is the ratio of expectations with the actual serial weight, not the ordinary expectation of `X/A`. For example the replay gives the following sufficient tests for the actual unconditional chain:

| Desired band `bx≥Y` | Sufficient bound on `ρ_n` at that label |
|---:|---:|
| 140 | `19/60` |
| 135 | `1/4` |
| 120 | `1/8` |

An exact stock-only sufficient observable is the refined occupation majorant for `X` above. Its simplest coarse proxy, the old leaf-target flux `F=Σ_core L/D`, has the following exact one-jump identity. For a core recipient with core degree `c=D-L`, `ΔF=c/[D(D+1)]`; for a leaf recipient at a core `o`, `ΔF=1/2−1/D_o`. The latter includes both the loss of a leaf at `o` and the appearance of a new degree-two core with one leaf. Thus at `n≥3`,

```math
 n(PF-F)=\sum_{o\text{ core}} r_o\frac{c_o}{D_o(D_o+1)}
       +\sum_{o\text{ core}}\frac{L_o}{D_o}\left(\frac12-\frac1{D_o}\right).\qquad\text{(2)}
```

Here `P` is one **actual** old-tree jump. The drift is nonnegative but does not on its own bound `ρ_n`. The degree/stock means alone also do not close (2); `r_o=L_o+Σ_{v~o,core}1/D_v` involves the incidence of unequal neighbouring degrees, and the positive-part majorant cannot be replaced by the positive part of an averaged expression. For example a pure star has large leaf-target flux and *zero* adverse old-leaf weight; doubling each arm of a large star instead gives `ρ` approaching `1/7` (compute the centre, degree-two arm, and leaf weights as in the hostile note). A degree sequence by itself does not determine `ρ` either: the replay compares two reachable size-seven trees with degrees `3,2,2,2,1,1,1`, but different `X/A`.

A usable nonstationary occupation certificate would propagate a signed correction over old jumps. For example, for target bound `r`, it suffices to find nonnegative test functions `V_k` on physical trees with `P V_{k+1}≤V_k`, and `X-r A≤V_k-c_k` pointwise for a known `c_k≥0`. Then `E[X-r A]≤E V_{k_0}-c_k` at every later size, so any deterministic lower bound for `c_k` above the initial budget closes that range. More generally the finite-horizon backward test is exact: for a terminal local upper charge `J_k≥X-rA`, recursively set `J_{i}=P J_{i+1}`; `J_{k_0}` can be upper-bounded without enumerating all future trees if one can find local super-solutions. Equation (2) identifies the extra old-edge incidence needed even for the simplest stock correction. No certificate of the required sign in this propagation has been found here.

The hubhood retry floor from [hostile](../../hostile/README.md) is a different endpoint estimate; it does not replace `c_b(1-x)^b` on positive thresholds. No bound on the annealed adverse hubhood debit over the whole finite-label gap is proved here.

## Replay

```
python3 -IB round2/open/age-order/global-transport/annealed/extension/check.py
PASS true old-jump stock flux identity and paired edge caps on 4 physical fixtures
PASS degree-marginal obstruction: equal degree sequences, X/A = 2/59 versus 0
PASS paired cutoff 141; annealed cutoff tests (140,19/60), (135,1/4), (120,1/8)
PASS fixed threshold x>=10/11 for every n+1>=20
```

The replay checks rational margins and local jump identities. It does not grow large trees, estimate occupation weights from sampled histories, or test the unknown small-threshold sign.
