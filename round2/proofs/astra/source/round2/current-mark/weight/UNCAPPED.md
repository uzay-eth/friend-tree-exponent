# Uncapped current marks: explicit gap and compact input approximation

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 13: Claude for transfer bookkeeping; GPT for the uncapped bound and root transfer; inputs include physical jump geometry and the Round 1 exponent and prior. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


The true reciprocal gauge `q(D)=1/D` has `R_c W ≤ .53 W` for all finite
shapes and all active marks `0≤a≤1/2`. This uses neither a cap in the
generator nor historical marks inside a shape. An alternative **21-label**
coarsening gives `ρ=.56495`. This is the smallest successful grid in the
short searches reported here, not a global label-minimality theorem.

## Weight and degree-infinite proof

Use the weight, active edge inequalities (1)–(2), and local inequality (3)
in [the cap-36 proof](README.md), with `d=.365`, `γ=.53−d=.165`. For
uncapped marks set `z_i=1/(1+di)` and `k_i=z_i z_(i+1)`. There is no plain
copy; in the formula charging its local plus increment the multiplier
`p/q` is one. The same proof applies. In particular (2) retains the
true sender rate and the decrement at the first child. All deeper
factor decrements may be dropped, since neither dormant channel gives
a positive coefficient outside its corresponding arm. Endpoints now
include active mark **zero**.

Beyond the *proof* seam `M=400`, use

```
P_i=C/i, C=(d+1/2)/d²,       h(i,c)=P_M/γ .
```

Thus `P_i≥z_i`, `P_i−P_(i+1)≥(d+1/2)k_i`, and the root gain term in the
edge inequalities is nonnegative. For all child degrees `j`, put
`K=max(max_{j<M} k_j/P_j, M k_M/C)`. This bounds `k_j/P_j`: for `j≥M`,
the function `j/[(1+dj)(1+d(j+1))]` decreases because `dM≥1`.
For a source `i≥M`, `(P_i−z_i)≥(C−1/d)/i≥K/i`; this pays the entire
grandchild contribution even if the gain of the first child's factor is
discarded. For a source `i<M`, a tail child uses the interval
`p_j∈[0,1/M]`, and `k_j/P_j≤k_M/P_M`; dropping that child's nonpositive
change reduces the check to the rational tail corner in the script.
This last corner is implemented by dividing the row by `P_j`, then
multiplying by the positive surrogate `P_M`; it does not lower bound
the true tail factor by `P_M`. For the tail of the local profile, h is
constant and `bP_i≤P_M=γh`. The exact run checks 638392 edge corners and
159598 local rate endpoints. The required virgin normalization is
`999609761/10^9 < 1`; `h_min=86271979/10^9`. A smaller replay at seam 160
uses `d=.34` and proves `ρ=.553`.

Active `W` has positive uniform minimum `h_min/(d+1/2)`. Dormant
weights can tend to zero with the mark drop. The proof controls
discounted **active** response and scaled dormant response. More
precisely, dividing the dormant score by `u−l` gives a uniform positive
weight `V/[(l+d)(u+d)]`; in that equivalent mean equation the active
emission rate is multiplied by `u−l`. In addition to intrinsic evolution,
the scaled dormant mean row acts on a score `f` by
`f_active(l; S) − u f_scaled(l,u; S)`. This is a similarity of mean
generators, not a change of the raw particle law. A branching realization
of the scaled row emits an active at rate one **while retaining the
dormant**, with independent dormant killing at rate `u`. In the raw
realization, activation consumes the dormant at rate `u−l` and independent
killing occurs at rate `l`. It is not necessary to control
the *unscaled* number of dormant particles by a positive minimum of the
unscaled weight. For nonexplosion, use `n²` on actives and `(u−l)n²` on
raw dormants. Its drift is bounded above by a constant times itself:
intrinsic jumps have rate at most `n`; the sum of squared emitted subtree
sizes at either of the first two layers is at most `n²`; `r(u−l)≤1`;
plus-copies replace at most one vertex per copied branch, with sums of
squared plus sizes at most `(2n)²`. On a dormant row activation emits
an active with rate `u−l`, so this part is bounded by its own weight.
Finite stopping then gives a uniform bound on the expected active `n²`
account on bounded time intervals. The *raw* rate of dormant emissions
per active is at most `2n²`; hence the stopped count of raw emissions has
bounded mean on bounded time intervals. All other particle births are
bounded by the active size account. Removing the stops yields
nonexplosion, since finitely many emitted or activated particles and
finite intrinsic growth of each cannot produce infinitely many events
in bounded time. The scaled version avoids counting raw dormants in the
size-exit error.

The finite size-cutoff argument in [the cap-36 proof](README.md) works
unchanged in the scaled positive system, since every emission is a
subtree or plus subtree no larger than the source. Its rate is now
`(2/(N+1))^ε` whenever `.53+ε(1+.165/2)≤s`. At fixed `N` new marks
come from degrees ≤`N`; no infinite mark enumeration is needed to solve
its rational equations for a fixed input mark.

To bound centre immigration, use the same centre score as in the
cap-36 argument. The drift of `V` is ≤`γV+1`. The dormant founder
terms remain bounded by multiples of `V` despite `inf P=0`:
`r k_i/P_i≤1/d` and `p_i k_j/(P_iP_j)≤(1/(2d)+1)/d` for
`i,j≥2`, using `P≥z`, `r≤i`, `k_i/z_i=z_(i+1)` and `p_i=1/i`.
Copies and plus-copies are bounded using `p/q≤1`, `z/P≤1`, and
`V(Y⁺)≤V(Y)/h_min`. Thus the founder weight grows at rate at most
`γ`, below the physical exponent. This justifies the physical renewal
comparison for the uncapped gauge when the newborn kernel is viewed
on the compact mark space and its positive reproduction operator is
used. Each isolated finite label has positive direct production from
every row; strict positivity on nonzero continuous nonnegative scores
follows since zero is the limit of those isolated labels.

## Compact-mark rate for input interpolation

Let `G_s(a)=U_a(J)=∫ e^{−(s+a)t} E_J^0 F_g(S_t) dt` and
normalize `||g||∞≤1`. The common forcing `F_g` has no dependence on the
input mark. Write `A=1/2`, `h*=h_min`, and

```
L = (A+d)/(h*(s−ρ)),
M_F = 1/h* + L [1/h*+(2+A/d)/d].
```

Indeed `|U|≤LW` by the residual gap. The plus-copy forcing weight
is at most `(1/h*)V`, since `p/(p+d)=z≤P` and
`V(Y⁺)≤V(Y)/h*`. The root-drop forcing coefficient per first arm is
at most `r k_i/P_i≤1/d`. The grandchild coefficient per second arm
is at most `p_i k_j/(P_iP_j)≤(1+A/d)/d`. The direct toll is at most
`1≤V/h*`. Thus `|F_g|≤M_F V`. The exact intrinsic profile gives
`E_J^0 V(S_t)≤e^(γt)`, and for each nonnegative integer `j`,

```
|∂_a^j G_s(a)| ≤ M_F j!/(s+a−γ)^(j+1).
```

This proves a quadratic linear-interpolation bound on every interval
`[a,b]`: the error at `x` is `≤M_F (x−a)(b−x)/(s+a−γ)^3`.
Over `[0,1/M]`, the Taylor series at `1/M`, in the positive coordinate
`1/M−a`, has nonnegative coefficients for `g≥0`. Keeping `m` terms has
error at most

```
M_F/(s−γ) · [1/(M(s−γ))]^m .
```

This follows by integrating the exponential tail term by term and
summing a geometric series. For `s=.56572365`, `M=20`, `m=16`, the
uniform row error on tests of norm ≤1 is below `5.044e-11` **if the
coefficients are evaluated exactly**, or plus the propagated coefficient
uncertainty otherwise. This gives a finite-rank approximation in the
input row variable; it does not evaluate the coefficient operators on
the infinitely many output labels or control conditioning of a numerical
derivative fit. Compactness also follows from the corresponding uniform
first derivative bound and Arzelà–Ascoli. The independent interpolation
control checks the derivative/linear-remainder formula on rational
single-rate test transforms, not on the unknown physical forcing.

## Coarsening and layer limitations

The certified coarsened gauge is `q(D)=1/max{j≤D:j in I}` where

```
I = {2,3,...,15,20,25,30,35,45,60,80}.
```

Its 21 labels satisfy `q≥p` (plain copies are generally present). The
same drift proof works with the exact jumps in `k`, a constant tail now,
and `ρ=.56495`. The checks are 25272 edge corners and 6318 local endpoints;
virgin normalization is `999635653/10^9`. Short grid/dropout searches
produced 23 then 21 labels; no successful 20-label grid was found in the
tested three-deletion family of the 23-label candidate. These are search
results, not a lower bound on the number of labels any proof needs.

To test what dropping the plus-copy surcharge might buy in a layer
construction, fix the **seam-160 edge table** and reuse the
same shared local profile and constant h-tail. Even relaxing the
gain interval from `[L,L+c/2+A+d]` to `[L,L+c/2]` everywhere, the
backward recursion for `ρ=.53` needs a virgin `h(2,1)>1`. This is a
limitation of those fixed edge factors and independent local tests,
not a no-go theorem for distinct edge or h tables by layer. The
successful uncapped tail already lowers the certified exponent from
`.564` to `.53`; no separate layer-dependent certificate is claimed.

## Replay

```
taskset -c 8 timeout 13s python3 -Bu round2/current-mark/weight/uncapped.py
timeout -k 1s 14s taskset -c 8 python3 -Bu round2/checks/uncapped_gauge_sol.py
taskset -c 8 timeout 13s python3 -Bu round2/current-mark/weight/uncapped.py --coarse
taskset -c 8 timeout 13s python3 -Bu round2/current-mark/weight/controls_uncapped.py
taskset -c 8 timeout 13s python3 -Bu round2/current-mark/weight/interpolate.py
```

The second command is the independent exact check at `M=400`, `d=.365`, `rho=.53`; it uses integer cross-products for the edge corners. See the cold referee report. The separate two-level check tested 288 exact rows through and far
beyond the proof seam, as far as degree 10000. It rejected using a
wrong grandchild sender rate. No near-critical kernel sign is inferred.
