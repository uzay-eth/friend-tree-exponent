# Independent capped fringe kernel checker

This directory reconstructs the 19 labels and generator from `../fringe/renewal.md`
and `../fringe/finite-kernel.md`. It uses Python 3's standard library only.
See the continuation backend (`CONTINUATION.md`) for the extended certificate
fields, independently checked unbounded proofs, and batched signed pilots.

## Certificate interchange (`sol-fringe-v1`)

* `schema`: `"sol-fringe-v1"`.
* `s`: exact rational greater than `0.51`. Exact rationals throughout are JSON
  integers or strings such as `"17679/31250"` or `"0.565728"`; never JSON floats.
* `cutoff`: integer `N >= 2`, cutting at total fringe size (core plus leaf
  vertices, excluding the external parent slot).
* `g`: 19 strictly positive rationals in increasing birth-label order `2..20`.
  Bounds apply to this particular vector, at this particular `s`.
* `scale`: optional positive integer, default `10^24`; score grid is `1/scale`.
* `lower`, `upper`: optional continuation objects. Each has `kind` (default
  `"default"`), `base` (default `"default"`), and `values` (default `[]`). The
  default lower is zero; the default upper is `max(g) W/(s−0.51)`. Optional
  upper `base: "drift-W"` uses `max(g) W/(s−0.51+0.7)`.
  A value entry is `[state, rational]`, where a state is
  `[birth_label, leaf_count, sorted_core_children]`, recursively. Child arrays
  are in lexicographic order; duplicate children retain their multiplicities.
  Example: `[2,0,[[3,1,[]],[3,1,[]]]]`.
  * `kind: "bellman"`: supply any finite list of nonnegative overrides. All
    overridden states must pass the full exact Bellman inequality. An upper
    override must not exceed the selected fallback. Missing states evaluate
    to the fallback, including when checking successors of overridden states.
    The verifier applies the continuation only at actual exits. Finite witness
    states beyond `N+1` are allowed.
  * `kind: "external"`: same values, with a nonempty `proof` string identifying
    a separately certified claim. This mode does not itself verify that proof
    and marks the corresponding result **conditional**. Do not use a
    conditional sign as a root certificate until the specific claim is checked
    against the same `s`, `g`, states and values. The killed sign always remains
    independent and unconditional.

`--emit-trial` prints a JSON certificate line with a positive trial `g`. A
colleague may use another basis or projection to derive the values, but the
exchange values at an exit must be bounds on the *full scalar continuation*
`U_s g(S)`, not corrections to zero or matrix entries. A tighter lower witness
can supply only a subset of actual exits; a tighter upper witness can likewise
fall back to a global analytic upper. No generator data are taken from a
certificate.

## Checks and justification

The generator enumerates directed arrows over vertex paths in each canonical
intrinsic fringe. It builds the copy rates separately, uses the uncapped
current degree for them, verifies nonnegativity, and checks the intrinsic
mass `n−1/D_root`. Initial-centre terms are checked using a different
representation, with physical parents and chronological vertex IDs.

For a continuation `v`, each equation has diagonal `s+a+n−1/D_root`; its right
side is `b*g[label]`, plus retained transitions applied to the unknown and
exiting transitions applied to `v`. There are cycles (growth upwards in size,
copies downwards); the checker does not assume a triangular system. A floating
iterate only proposes a candidate and a positive vector `g`. For the finite
system the verifier computes the exact drift of the integer column formed by
rounding `1024 W`, and checks that this drift is strictly positive. An exact
scalar multiple of that column then corrects the candidate in each direction.
A separate pass reconstructs every residual using state keys and the original
intrinsic/copy rate maps. The finite inverse is positive by the killed
semigroup argument in `finite-kernel.md`. Hence residual signs prove outward
scores; the floats cannot certify anything. The report gives maximum
per-newborn arithmetic enclosure width in ratio units, and the maximum
absolute finite Bellman residual in score-rate units. Printed decimals are
rounded outwards using integer division.

In the full process, the checked weight proof gives `W >= 1`, and the stronger
drift `(s-R)W >= (s−0.51)W + h(c_root)` implies the `drift-W` option because
`h >= 0.7` and `b <= 1`. For a finite Bellman override, rows off the override
set keep the inequality by monotonicity: lower replacements only increase the
zero fallback, upper ones only decrease the upper fallback. Every such score
is bounded by a finite multiple of `W`. Discounting its residual semigroup
contribution to infinity gives zero by `s > 0.51`, so the global lower and upper
inequalities bound `U_s g`. First exits transfer those bounds to the finite
problem. Finally `min_i lower(J_i)/g_i` is a certified lower bound on `r(K(s))`,
and `max_i upper(J_i)/g_i` is a certified upper bound; strict signs about 1
translate into strict bounds on `mu`. For the killed action, the reported
min/max ratios enclose the radius of **that killed matrix** itself; for the
upper action, the maximum is an upper bound on the full radius.

## Replay and pilot results

From the repo root (on this machine):

```sh
taskset -c 3 timeout 13s python3 -Bu round2/sol-kernel-verifier/physical.py
taskset -c 3 timeout 13s python3 -Bu round2/sol-kernel-verifier/verify.py --pilot 4 6 8
taskset -c 3 timeout 13s python3 -Bu round2/sol-kernel-verifier/verify.py --drift --pilot 4 6 8
taskset -c 3 timeout 13s python3 -Bu round2/sol-kernel-verifier/verify.py round2/sol-kernel-verifier/smoke.json
```

At `s=0.565728`:

| Cutoff | Retained states | Killed radius enclosure | Default full-radius upper | `drift-W` full-radius upper | `drift-W` max continuation gap in row-ratio units, including arithmetic |
|---:|---:|:---|---:|---:|---:|
| 4 | 152 | [0.329002481250, 0.329002481251] | 64.383162212546 | 5.052408539391 | ≤4.723406058140 |
| 6 | 1330 | [0.400769312094, 0.400769312095] | 57.831778384692 | 4.635778830861 | ≤4.235009518766 |
| 8 | 14687 | [0.449594247514, 0.449594247516] | 53.507529484738 | 4.362129729726 | ≤3.912535482211 |

At cutoff 8, the widest arithmetic interval of any same-test newborn ratio
was at most `2.028e-12` (default upper) or `1.468e-13` (`drift-W` upper).
The maximum absolute residuals there were `6.538e-14` (killed), `3.198e-11`
(default upper), and `2.235e-12` (`drift-W` upper). The approximation error
from truncation dominates. These pilots give no improvement to the existing
physical interval `[0.5656581227, 0.5657979941]`.

Physical enumeration tests the generator for every individual marked-fringe
indicator, including immigration, on 2,956 chronological marked trees of
sizes 3–8. Two added reachable fixtures exercise duplicate identical core
 children and actual inverse degree `1/22` at a capped mark `1/20`. Five actual
planted defects fail: activating the inert parent slot, deleting a copy for an
identical child, substituting `1/20` for `1/22` in the copy rate, perturbing a
newborn score inwards, and assigning an unjustified zero upper continuation.
The check also encloses the hand-derived killed `N=2` action on all 19 rows.
The JSON smoke witness replays a finite lower continuation and certifies
`mu < 2` with full-radius upper `0.702175857135`.
