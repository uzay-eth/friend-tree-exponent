# Explicit residual gap for the capped current-mark gauge

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 13: Claude for transfer bookkeeping; GPT for the uncapped bound and root transfer; inputs include physical jump geometry and the Round 1 exponent and prior. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


For the **cap-36** gauge `q(D)=1/min(D,36)` this gives a rational
pointwise inequality `R_c W ≤ (141/250) W` (`ρ=.564`) on every finite
active and dormant state, with a positive weight at most linear in
physical size. It gives explicit first-size-exit bounds for its 35-label
kernel. This is a larger gauge than the requested cap 20; it still does
not store permanent marks inside any tree. The same construction did
not certify cap 20 or a binary gauge. The executable rebuilds its
rational certificate and checks every inequality.

## 1. Weight and exact drift decomposition

All degrees in a planted tree include the root's inert parent slot. Let
`d=8/25`, `γ=ρ−d`, `z(D)=q(D)/(q(D)+d)`, and
`k(D)=(z(D)−z(D+1))/d`. Write `P(D)` for a nonincreasing **modified** edge
factor with `P≥z`, and `h(D,c)` for the local profile. Both tables are
constant in degree beyond the cap; the tail of `h` is also independent
of `c`. For a tree `S` set

```
V(S) = Σ_x w_x h(D_x,c_x),   w_x = ∏_{y strictly above x} P(D_y),
W_a(S) = V(S)/(a+d),
W_[l,u](S) = [1/(l+d)−1/(u+d)] V(S).
```

`h(2,1)=1`. Thus the increase of a `V` summand under a promotion has
local contribution `h(D,c+1)−h(D,c)+P(D)`. A promotion preserves total
degree, so it does not change other old vertices' path factors. A gain
can change the local `h` and all descendant path factors; every path
factor change is nonpositive. Let `u_a=a+d`. The outgoing child copies,
after multiplying active drift by `u_a`, split exactly into the positive
`u_a z(D) Σ_Y V(Y)` and

```
u_a (p/q) z(D) Σ_Y [V(Y⁺)−V(Y)].
```

The latter has two terms per child: a local `u_a (p/q) z(D) Δh_Y` and
a nonpositive term from its factor `P(D_Y+1)−P(D_Y)`. Drop the latter.
Charge the former to the intrinsic gain rate at the depth-one child as
an extra virtual rate `u_a (p/q) z(D)/P(D) ∈ [0,1/2+d]`.
The resulting local gain interval at **all** vertices is contained in
`[L, L+c/2+1/2+d]`; the full interval at depth one is essential if
`Δh_Y>0`. This includes every possible actual incoming clock.

The two dormant emission terms are `u_a r k(D) Σ_Y V(Y)` and
`u_a p Σ_Y k(D_Y) Σ_{Z child of Y} V(Z)`, where `r` is the true root gain
rate and `p=1/D`. The same `r` decreases root path factors, and the
specified true `p` decreases the first child's path factors. No bound
on total subtree mass or degree arity is used. Compare each coefficient
of a depth-one `h` summand with

```
(u_a+r) P(D) − r P(D+1) − u_a[z(D)+r k(D)] ≥ 0.            (1)
```

For every depth-two or deeper summand, the factors below its first child
cancel from the comparison. With `j=D_Y`, the required coefficient is

```
[u_a(P(D)−z(D)−r k(D)) + r(P(D)−P(D+1))] P(j)
        + p P(D) [P(j)−P(j+1)] − u_a p k(j) ≥ 0.          (2)
```

Here the gain at `Y` is charged only at its guaranteed input rate `p`;
the rest of its factor change and all deeper changes can be dropped.
In (1) and (2) the `u_a P(D)` payment is available because
`−a V + u_a(V−h_root) = d V−u_a h_root`. The `−u_a h_root` is not needed
by this proof. The remaining local drift is at most `γ V` by the checked
degree inequalities

```
r_loc [h(D+1,c)−h(D,c)]
 + b[h(D,c+1)−h(D,c)+P(D)] ≤ γ h(D,c).                    (3)
```

This proves the active bound. On a dormant row the consuming
activation and killing together act on the weight exactly as `d W`,
since `(u−l)/(l+d)−u(1/(l+d)−1/(u+d))=d(1/(l+d)−1/(u+d))`.
The intrinsic part alone is at most `γ W` by (3) and nonpositive path
changes, so the same bound holds on dormants. In particular, the
grandchild term in (2) must use the **sending root's** `p`.

## 2. All-degree rational verification

The only states needing (1)–(2) have an old child. For (1) the rate fits
`[0,D−3/2]`: the root has `c>=2`, input at most `(c−1)/2`, and `L=D−c`.
For (2), the specified child supplies `p_Y`; with the other `c−2`
children at most `1/2`, the root rate lies in `[p_Y,D−2+p_Y]`.
Both formulas are affine in `r` and `a`, so checking those endpoints
and `a∈{1/36,1/2}` covers all interior values. At `j>=36`, `P` is
constant and `k=0` and one interpolates in `p_Y∈[0,1/36]` in those
rate endpoints. At `D>=36` all terms involving `r` vanish and one
interpolates in `p∈[0,1/36]` in (2). The profile tail sets
`h=P(36)/γ`, so at `D>=36` (3) follows from `b≤1`. All other (3) rows
are affine in `r_loc`; `b=0` rows do not read an unavailable successor.

`cert.py` proposes the edge factors numerically, rounds them to rational
ticks and performs the actual exact endpoint proof. It constructs `h`
by backward rational degree/core-degree recursion, normalizes the
virgin to one, and rechecks every inequality after rounding. The
passing default run gives 5178 edge corners and 665 profile rows;
the smallest edge slack is strictly positive. It finds `h≤1`,
`h≥1578186331/5000000000`, and `0<P≤1`. These facts give a uniform
positive minimum for active weights and, for the finitely many
nontrivial mark drops, a positive minimum for dormant weights as well.
Both weights are bounded by explicit constants times vertex count.

## 3. Finite cutoffs and relation to the physical root

Only an intrinsic jump increases the physical size `n`, and every such
jump increases it by one. All emitted children and grandchildren and
activation successors have size no larger than the source. The intrinsic
event rate of a planted shape is `n−1/D≤n`. Both active and dormant
intrinsic drifts of `W` alone are at most `γW`. For `0<ε≤1`, apply concavity
`(n+1)^ε/n^ε−1≤ε/n` to `W_ε=n^εW`. The exact intrinsic split and the
nonincreasing size of every emission give

```
R_c W_ε ≤ [ρ+ε(1+γ/2)] W_ε,
```

since `n≥2`. Use any rational `ε>0` satisfying the right exponent `≤s`.
If killed on first intrinsic size `N+1` along each lineage, the total
discounted exit weight `W` is at most `W(S)(n(S)/(N+1))^ε` by putting
discounted `W_ε` of each stopped lineage in an absorbing account. The
drift proof holds for finite stopping and then by positivity. In
particular the full newborn action on `g≥0` minus the size-killed action
is at most

```
 ||g||∞ (1/2+d) W_a(J) / [h_min (s−ρ)] · (2/(N+1))^ε.
```

Indeed `B_g≤||g||∞(1/2+d)W/h_min` and the full positive continuation is
its Laplace drift bound. The killed responses are finite rational for
rational `s,g`; one can instead propagate the exact first-exit weights
and use this uniform formula only as a convergence guarantee. The
exponent is small at the target, so the size cutoff is not itself a
near-critical numerical sign certificate.

For completeness the centre immigration is bounded below the physical
critical exponent, as required by the [physical transfer](../physical/README.md).
Take the centre weight `V` with physical centre degree in its outgoing
factor and set its local score to one (no planted parent slot there).
On noncentre vertices (3) holds for all actual gains. At the centre
itself the factor decreases on gains, while promotions add at most one
because `P≤1` and `h(2,1)=1`. Thus `L V≤γ V+1`, giving
`E V_t=O(e^{γ t})` for the initial three-vertex seed. The immigration
row applied to `W` is bounded by `C V`: ordinary copies are positive
and bounded by a constant times the child `V` sum, including plus
copies since their factor decreases and `h` is positive and bounded;
dormant root emissions have `r` bounded by the cap whenever `k(D)>0`,
and grandchild emissions use rate `p≤1/2`. These sums at depth one and
two are bounded by fixed multiples of `V` since `inf P>0`. The direct
newborn row is bounded. The same centre argument works at physical
degree two, with its score one and its actual degree used. The positive
transfer and direct-generation argument therefore apply to the cap-36
gauge: its kernel has the physical critical root. This changes neither
the physical model nor the fixed-mark kernel's characterization, and
does **not** identify the two kernels entry by entry.

## Replay and scope

```
taskset -c 8 timeout 13s python3 -Bu round2/current-mark/weight/cert.py
taskset -c 8 timeout 13s python3 -Bu round2/current-mark/weight/physical_check.py
taskset -c 8 timeout 13s python3 -Bu round2/current-mark/weight/cert.py --cap 20 --rho .56499
```

The labelled-arrow replay checks both active mark endpoints and dormant
activation on 95 shapes, including broad roots and children at and beyond
the cap. It rejects two actual rate plants, using an incorrect send rate
and mixing cap-20 coefficients into the cap-36 proof. The conservative
cap-20 branch of this certificate fails normalization; this is a failed
candidate, not a no-go theorem. No binary bound or near-critical sign is
claimed.
