> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; Claims 6, 10, 11, and 18. See the [release index](../../../../../README.md#claim-22) for exact scope and citations.

# All positive moments at late birth labels

Let `Z_v=lim D_v(n)/n` be the limiting degree fraction of birth label `v`, `X=lim C_n/n^μ`, and let `A,B` be the single and signed adjacent scaled tail profiles of [scaled](../scaled/README.md). For **every fixed real `α>0`**, define

```text
b_α = α ∫₀∞ y^(α−1) A(y) dy = E[X] ∫ u^α ν(du) ∈ (0,∞).
E[Z_v^α] ~ b_α v^(μ−1−α),
E[Z_v^α] − E[Z_(v+1)^α] ~ (1+α−μ) b_α v^(μ−2−α).
```

Here `ν((y,∞))=A(y)/E[X]` is the **finite positive-mark intensity**, not a probability. In particular the conditional moment of `v Z_v` given hubhood tends to `b_α/A(0)`, not `b_α`. The adjacent comparison is eventual for each fixed power, not uniform over all powers or a claim about early labels. At `α=1` it reproduces the eventual coefficient of [mean](../mean/README.md); at `α=μ` the spectrum normalization gives `b_μ=μ K E[X]` and the adjacent coefficient equals `b_μ`. This does not re-derive early hub or mean comparisons or the pressure identity.

## 1. Actual placement and stopping-size domination

For a mark in a finite state at size `m` with degree `d`, use the same uniform for the target-gain Bernoulli and a degree-count urn gain of probability `D_urn/j` at size `j`. The physical gain probability is `Σ_(w~mark) 1/D_w(j) / j ≤ D_mark(j)/j`. The common uniforms preserve the order from equality, so the future limiting fraction is stochastically bounded by `Beta(d,m−d)`. This also holds conditionally on the whole history at a finite size stopping time by strong Markov. For a true leaf, start instead at its **promotion size** `M`, where `d=2`; on no promotion its mark is zero. In particular

```text
P(Z>t | current state) ≤ P(M<∞ | current state) Q_m(t),
Q_m(t) = (1−t)^(m−2) [1+(m−2)t],  0≤t<1, m≥3.
```

This uses `Q_M(t)≤Q_m(t)`: the difference `Q_m−Q_(m+1)` is `(m−1)t²(1−t)^(m−2)`. The same bound without the promotion factor applies to a mark already at degree two at size `m`, such as the middle vertex in the length-two placement. Integration (or the urn law) gives, for real `α>0`, the upper factor `Γ(α+2)Γ(m)/Γ(m+α)` times the indicated promotion probability when applicable. There is no claim that promotion probability and its later mass are independent. Using the degree-one urn at birth alone cannot give the additional rare-event factor `v^(μ−1)`.

On integrating over the actual placement of label `v`, the [late-birth finite-promotion estimate](../late/README.md#2-the-tail-estimate-the-tagged-leaf-must-stay-in-the-parent-process) gives `P(M_v<∞)≤7 E C_(v−1)/(v−1)`. Thus, for `n≥3` and all `y≥0`, the single-birth integrand has envelope

```text
n^(1−μ) P(Z_(n+1)>y/n) ≤ C (1+y)e^(−y/2).
```

Indeed, if `y<n`, use `m=n+1` and bound `1−y/n≤exp(−y/n)`; if `y≥n` the probability is zero. For adjacent differences one cannot subtract two absolute envelope estimates at this order. Instead use the signed original two-birth formula below.

## 2. The signed double-placement mass and an integrable adjacent envelope

With an old size-`n` state `T`, put `r_u=Σ_(w~u)1/D_w`, `δ_u=1/[D_u(D_u+1)]`. Use the exact two-birth identity with any bounded terminal test of the endpoint's ultimate mass. In the length-two placement at `u` let `I(u)` be its conditional expectation at middle minus tip, in the same final tree; in a double placement at `u,w` let `H(u,w)` be the difference at the two new tips. The orientations and rates are

```text
S(u)= r_u/(D_u+1) I(u) − r_u δ_u Σ_(w~u) H(u,w),
n(n+1) E[test(Z_(n+1))−test(Z_(n+2)) | T_n=T] = Σ_u S(u).
```

Indeed the two orders are identical on nonedges; on an edge the second target rate is reduced by `δ` at the first recipient, and the orientation sum leaves `(r_w δ_w−r_u δ_u)H(u,w)`. At a repeated recipient the two tips are exchangeable. The residual event is the second birth attaching to the first newborn, at rate `1/(D_u+1)`. Thus all differences use their original joint future, not two separately sampled trees. Using the harmonic test `E[Z^α|final tree]` in this identity is valid for any fixed power; no signed stochastic domination of two urns is needed.

For the tail test at `y/n`, all marks at the common final size `m=n+2` have degree one or two. Each endpoint tail lies in `[0,Q_m(y/n)]`, so **each absolute conditional difference is at most `Q_m(y/n)`**, not an estimate at the single-birth scale. For each core `o`, charge its own oriented `S` and those of its adjacent old leaves. The total absolute coefficient is at most `2 r_o/(D_o+1)+Σ_leaf 1/D_o≤3`: the two types of oriented charge each cost at most `r_u/(D_u+1)`, and `r_ell=1/D_o` for a leaf. Consequently

```text
n^(2−μ) |P(Z_(n+1)>y/n)−P(Z_(n+2)>y/n)| ≤ C (1+y)e^(−y).
```

The right side bounds the absolute **signed** integrand uniformly over the entire mass scale, even if the finite-label sign is unknown at large `y`. Apply layer cake, then dominated convergence using the compact profile limits. The single integral gives `b_α`; the signed integral gives `a_α=α∫ y^(α−1) B(y)dy`. From the scaled telescoping relation, `B=(1−μ)A−yA'` on `(0,∞)`. Integrate first on a compact interval; `y^α A(y)` tends to zero at both ends by boundedness at zero and the urn envelope at infinity. Hence `a_α=(1+α−μ)b_α`. Equivalently, summing the adjacent expansion and comparing with the single expansion forces this same coefficient. The proof uses original-size probabilities throughout.

## 3. Cold audit of the two input claims

- **Spectrum/mark law.** The promotion envelope must use `P(M_v<∞)`, not the probability of eventual positive mass as a substitute prior to applying the urn. The late bound for finite promotion supplies exactly this. For the full expected count, compact-scale Riemann summation uses `y=xv`; the small-label contribution is `O(ε^μ)` and the large-label contribution is controlled by the exponential envelope. The `L¹` upgrade truncates first in labels and then in finite promotion time: discarded late finite promotions are controlled in `L¹`, and the retained count is bounded by the core count at a deterministic size and thus is uniformly integrable after scaling. Infinite promotion sizes are excluded. With `K E[X]=∫ y^(μ−1)A(y)dy`, the annealed joint intensity really is `K⁻¹ 1_(0<y<u)y^(μ−1)dy ν(du)`. The Jacobian for `r=u/y` gives `u^μ r^(−1−μ)`; the `μ` in the normalized Pareto density is essential. The count-weighted law is not a quenched marked law. The unmarked threshold ratio alone does transfer quenched.
- **Mean gap.** At a finite restart, the conditional physical mass satisfies `W=G_m Z`, where `G_m` has Gamma size law independent of the embedded path, **not independent of `W`**. Thus in the mean exact identity each terminal mass is divided by the same size `n+2`; more generally the divisor for the physical `α` moment is `Γ(n+2+α)/Γ(n+2)`, not `(n+2)^α` exactly. The old-edge negative rate change, the unmatched length-two placement, and the grouping bounds justify the eventual mean argument. The core donor upper bound pairs each future promoted neighbour's discounted incoming budget (at most half its birth discount) against the lost leaf. Before a leaf's promotion, only its parent can promote it, and its own rings still increase that parent's degree. Those facts justify the two competing-clock leaf estimate in the mean note. The all-power proof above replaces the separate terminal mean bounds for the late-label result; it does not certify the numerical early-label intervals anew.

## 4. The shared-future obstacle and its resolution

For a rectangle in the label and mass coordinates, the existing proof identifies the expected count on the `x^(−μ)` scale and establishes uniform integrability by domination by the total count. It does not establish concentration conditional on the growing history. In particular, replacing the terminal indicators by their conditional expectations at their respective birth times does **not** yield martingale differences: the centered residual for an older birth is still a function of the same future at later births. In the exact variance expansion, the off-diagonal two-vertex terms do not vanish for that reason. One must also allow a random first-order factor such as `X`; the total count already has the nonconstant limit `KX`, so unnormalized deterministic self-averaging is impossible for that marginal. For a marked *proportion*, the question is cancellation of the common factor, which is not settled by the variance observation.

Even an a.s. unmarked count limit and convergence of all annealed marked rectangles do not logically settle this. For example, take `N_j=j` deterministically and, on one common fair coin, mark all `j` points by either zero or one. The annealed proportions are always `1/2`, yet the empirical marked fraction is a nondegenerate coin for all `j`. This is an obstruction to that inference, not a counterexample in the tree.

The separate [quenched spectrum proof](../spectrum/quenched/README.md#2-cohort-fluctuation-without-independent-terminal-marks) resolves the issue. It compares two backward martingales for a finite-horizon cohort count and its integrated birth forecast, using one-leaf perturbation bounds on their jumps instead of independent terminal marks. The subsequent simultaneous terminal-mass comparison and clock transfer give the almost-sure quenched joint limit and `L¹` convergence for each fixed rectangle. The unnormalized limit retains the common factor `X`, which cancels in the realized-tree proportions. See the [independent audit](../../../../../../../reviews/cold/astra/checks/marked-spectrum-sol.md).

## Exact controls

```text
$ taskset -c 8 nice -n 10 timeout 8s python3 -I -B round2/open/age-order/moments/check.py
PASS: 96 exact signed moment comparisons; 60 nonzero edge terms
PASS: 156 exact Beta moment/tail rows at actual promotion sizes
```

The signed fixture compares full rational two-step propagation against edge-paired propagation with a one-step continuation of the endpoint power. The proof for general real powers and the infinite horizon is analytic.
