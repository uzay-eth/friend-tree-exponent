# Independent audit of the common-time skeleton renewal

The checker reconstructs the skeleton of the **full nineteen-label residual
process**. It imports no project solver. All finite arithmetic and all sign
predicates use exact rationals. At `.56572365`, a 57-state smoke candidate
gives the full radius interval `[.340723077317, 5.568046915439]`, so the
near-critical sign is unresolved. At `s=4` the checker proves `<1` with
upper `0.338949925887`, and rejects an inward perturbation of that certificate.

## Skeleton first event

Let a fringe have `k` cores, root mark `a=1/i`, and immediate child marks
`a_j`. Retain all physical core degree slots. The total ring/copy rate before
killing is `C=k+Σ_j a_j`; the rate including death is `q=a+C`. There is a
unit clock for every core, a clock of rate `a_j` for every child (including
each duplicate), and a death clock at `a`. Between these clocks the leaf
counts are independent rate-one-per-leaf Yule processes **conditional on
the common waiting time**. Conditional on the grown leaf vector `l'`, the
operator applies the following intensities under the factor `exp(−(s+q)t)`:

- At every core `x`, a promotion at rate `l'_x/D_x` and one send to each
  physical core neighbour at rate `1/D_x`. At the root the external parent
  choice is a null event of rate `1/D_root`, and its successor is the
  **grown** state. No outside parent ever sends in.
- The toll is `b*g[min(D_root,20)]` at an intrinsic root promotion, in
  addition to continuing with its promoted state. Its expected
  discounted value uses the same first-event denominator. The toll's paid
  newborn is counted and set aside, separate from the physical promotion.
- The copy successor is the plus child at rate `p_root`, or the plain
  child at rate `a_j−p_root`. The same clock additionally continues the
  grown original at rate `a_j`. Death produces neither a continuation
  nor a toll.

The last point follows directly from the generator. Put `L_Y` for pure
Yule gains, `B_g=b*g`, and `E` for the non-killing skeleton event matrix,
with the core null in its diagonal and copy links all positive. Then
`R=L_Y−qI+E`; writing `H=(s+q−L_Y)^(-1)` gives `U=H B_g+H E U`. The rate
`q` stays fixed along each pure Yule run. Including the copy rate in `q`
requires its original-parent successor in `E`. One can algebraically cancel
those two terms and instead use the first **intrinsic-or-death** event plus
the expected copies **until** that event. The checker uses explicit copy
clocks throughout and tests missing parent retention as a planted error.

## Beta integrals and infinite tail

For `u=e^(−t)`, an initially positive leaf count `l` gains `r` leaves with
probability `binom(l+r−1,r) u^l (1−u)^r`. An initially zero count stays
zero. Let `m=Σ_i l_i`, `r=Σ_i r_i` and `λ=s+q`. The occupation weight
before the common clock for a specified vector of increments is

```
    [∏_i binom(l_i+r_i−1,r_i)] B(λ+m,r+1),
    B(v,1)=1/v;   B(v,r+1)=r B(v,r)/(v+r).
```

The convention for the coefficient at `l_i=0` is `1_{r_i=0}`. The total
increment has coefficient `binom(m+r−1,r)` by Vandermonde's identity.
Consequently truncating at total increment `T` gives exact nonnegative
remaining occupation mass `1/λ − Σ_{r≤T} binom(m+r−1,r) B(λ+m,r+1)`.
Equivalently it is `λ^(-1) ∏_{j=0}^T (m+j)/(λ+m+j)`, the product of
discounted probabilities of `T+1` successive leaf births before the clock.
The implementation uses this product and tests it against beta subtraction.
The reference explicit-clock routine uses the competing total rate
`λ+Σ l'_i` at each individual leaf birth, then emits one terminal skeleton
event. This rebuilds the discounted rows via a different recursion.

For scores bounded by a multiple of the original capped weight `W`, the
unseen part is uniformly bounded by that mass times a shape-only toll
budget. The budget accounts for core null/other sends (`k W`), all copy
marks (`Σ a_j [W(parent)+W(child_j)]`), birth changes to `W` with the maximum reachable
path factor, and the maximum immediate toll. Leaf growth does not alter
`W`. The optional `tail_extra` refinement applies the *baseline* envelope
`αW` to additional increment totals without allocating a joint vector.
For a one-site term `F_i(l_i+j)` its exact occupation contribution at
fixed total `r` is

```
  B(λ+m,r+1) Σ_j binom(l_i+j−1,j) binom(m−l_i+r−j−1,r−j) F_i(l_i+j).
```

This reduces the refined baseline tail from all joint allocations to
`O(k r)` marginal terms per total. It is valid because under `W` the
copy payoffs are independent of the grown leaf vector and the core reward
and promotion-weight increments are sums of one-site terms. General
candidate values need not have this structure.

For a genuinely correlated example, with two initial counts one and
degrees two (no local count change except Yule growth before a clock),
the averaged product of the two reciprocals conditional on total increment
`r` has unnormalised allocation sum

```
Σ_{j=0}^r 1/[(j+2)(r−j+2)] = 2 (H_(r+2)−1)/(r+4).
```

This is a single rational Mellin-moment series after beta averaging, with
recursive harmonic and beta coefficients. At `λ=9/2`, totals through 20
plus the tail bound `1/[2(r+2)]` for each payoff when `r≥21` give
`[.047990694888,.047991444309]` for the unnormalised occupation integral.
Alternatively the conditional one-tip reciprocal is
`u(−log u−1+u)/(1−u)^2`; square it before the common Mellin integral.
This reduction preserves dependence, unlike integrating each tip
separately. In fact, normalized shared-time mixing gives `Cov(L_1,L_2)=36/245`
for `λ=9/2` and initial counts one. No finite B-spline representation of
arbitrary candidate continuation tails is used or established here.

## Certificate interface

`audit.py FILE` verifies JSON with:

- `schema: "skeleton-fringe-v1"`, `s` an exact rational `>0.51`, `g` 19
  strictly positive rationals in label order `2..20`, `band` a nonnegative
  total-increment cutoff, and optional `tail_extra` a nonnegative integer.
- Optional `lower` and `upper` lists of `[state, rational]` overrides.
  A state uses `[mark_label, leaves, sorted_children]` recursively. All
  copies and slots are reconstructed by the code. Rationals are strings or
  JSON integers, never floats. Unspecified values use `0` below and `αW`
  above, with `α=max(g)/(s−.51+.7)`.

Every provided override must satisfy its skeleton Bellman inequality, with
the truncated finite term used for the lower RHS and an outward
remainder added for the upper RHS. Bounds are checked on the entire
19-entry newborn action using Collatz inequalities. All lowered upper
overrides and raised lower overrides help the inequalities off the support
by monotonicity. Put `P=HE`. The verified weight drift implies
`W−PW = H (s−R)W ≥ (s−.51)W/λ`. A skeleton step increases core count by at
most one along any successor. After `j` steps, `λ ≤ s+1/2+(3/2)(k+j)` from
an initial `k` cores. Thus `P^n W/W` is bounded by the product of factors
`1−(s−.51)/(s+1/2+(3/2)(k+j))` for `0≤j<n`, which tends to zero by divergence
of the harmonic series. First-event decomposition identifies the positive
reward sum with the full residual process, which is nonexplosive by the
original weight bound. The finitely checked global sub/supersolutions bound
the intended kernel. This interface checks finite state overrides, not
arbitrary formula certificates or signed approximations to individual moments.

## Replay

```
taskset -c 8 timeout 13s python3 -uB round2/sol-connected/skeleton-audit/controls.py
taskset -c 8 timeout 13s bash -o pipefail -c 'python3 -uB round2/sol-connected/skeleton-audit/controls.py --emit smoke | python3 -uB round2/sol-connected/skeleton-audit/audit.py /dev/stdin'
```

The control suite compares 50 exact beta/explicit-clock rows and 20 exact
marginal-tail/full-allocation rows. Planted independent-time factorization,
deletion of the root null event, deletion of an identical child's copy,
deletion of the original successor at a copy clock, and an inward
certificate override are all rejected.
