CORRECT

# Independent cold review: pointwise eigenfunction obstruction on vertices and oriented edges

Scope. Primary review: `edge-obstruction.md` (`closure/edge-obstruction.md`). I re-derived the mixed-branch
identity by expanding the physical friend-tree generator from the model (uniform vertex, then uniform neighbour,
then attach a leaf), checked it exactly on reachable rooted skeletons with arbitrary leaf counts and with core-degree
variation, tested the role-dependent first differences and their commutation, and reproduced the bare-path
contradiction. I tried to break the argument at `c=2` and with a distinct, arbitrary founder score. Secondary,
separate review (section 5): `single-vertex.md` (`closure/single-vertex.md`). No author file was changed. My
separate exact short checker is `edge_obstruction_hy4.py`.

Model used throughout (matching `bounds-hub.md` (`proofs/bounds-hub.md`):24): a core `x=(D,c)` has `L=D−c`
leaves; gain rate `r_x=L+Σ_{core w~x}1/D_w` (each leaf neighbour contributes 1, each core neighbour `1/D_w`), and
gain sends `(D,c)→(D+1,c)`; promotion rate `b_x=L/D=1−c/D` sends `(D,c)→(D,c+1)` and creates one core child in
state `(2,1)`. Score `Φ=f_o(x_o)+Σ_{v≠o}f(x_v)+Σ_{v≠o}g(x_par(v),x_v)`.

## 1. Mixed-branch identity (Section 1) — re-derived, confirmed

Let `x` be interior nonfounder with parent `y` (`c_y≥2`), child `z` (`c_z≥1`), other neighbours fixed. Expanding
the drift term by term: a promotion has rate depending only on `x`, and changes incident edge scores additively, so
it cannot couple `y` and `z`. A gain at a third site has the reciprocals `p(y),p(z)` in its rate but its increment
involves only its own edge. A gain at `y` or `z` alone does not touch the other branch. The only products carrying
both branches are the two cross terms of the gain at `x`:

```
p(y)·Δ₁g(x,z) + p(z)·Δ₂g(y,x),   B_x(z)=Δ₁g(x,z),  A_x(y)=Δ₂g(y,x).
```

This is a **structural** identity, true for arbitrary `f,g,f_o`. My checker expands the full generator on the
skeleton `o→(y,y')`, `y→(x,y'')`, `x→z` with a fixed stem above `y`, takes the mixed difference over the two
branches, and requires equality with `p(y)ΔB_x(z)+p(z)ΔA_x(y)`: 24 cases with independent leaf counts on `y,z`, and
4 cases where `y`'s core degree changes (an extra outward pendant core). All exact equalities held.

Because `sΦ` is a sum of single-site and single-edge terms, it has zero mixed branch difference; hence for an
eigenfunction the two cross terms vanish. Fixing two types with distinct reciprocals on one branch and varying the
other gives `B_x(z)=v_x−α_x p(z)` and `A_x(y)=u_x+α_x p(y)`. My checker recomputes the apparent `α` from each
branch; on arbitrary non-eigenfunction tables the two branches disagree (e.g. −18.7 vs −6.67 in the run), as they must.

## 2. Commutation, c=2, distinct founder

On an edge between two adjacent equal-type interior cores `(D,c)` with `c=2` (I used `D=2,c=2`, the bare-path
state), the two evaluations of the same four-corner difference give

```
Δ₁A_b(a)=α_b(p(a⁺)−p(a)),     Δ₂B_a(b)=−α_a(p(b⁺)−p(b)),
```

which are negatives of each other when the types agree (`α_a=α_b`, `p(a⁺)−p(a)=p(b⁺)−p(b)=1/(D+1)−1/D≠0`).
Since `g` is a function, the four-corner difference is path independent, so `2α(p⁺−p)=0`, i.e. `α=0`. This uses
only that the parent has `c≥2` and the child `c≥1`; with `c=2` the interior vertex has exactly one parent and one
child, so the equal-type adjacent pair exists. I verified both the path-independence identity and the coefficient
relation in exact arithmetic. A planted edge score with `α≠0`, namely
`g(y,x)=α[p(y)D_x−D_y p(x)]`, has a nonzero joint difference (`−11/400` in the run) and yields inconsistent
branch slopes (`−11` vs `−2`), so it cannot be an eigenfunction edge score. The founder score `f_o` was kept
arbitrary (random table) and never enters the local argument; a distinct founder is harmless.

Hence `Δ₂g(y,x)` is independent of `y` and `Δ₁g(x,z)` independent of `z`, which is exactly the integrability needed
for Section 2.

## 3. Moving differences off the edges and the bare path

Absorbing `h_out,h_in` into `f` leaves `Φ` unchanged and makes every edge increment zero at an interior endpoint.
For adjacent interior `u,v` with fixed `c` and independent `D`, the only joint degree term is
`p(v)Δf(u)+p(u)Δf(v)`; its mixed difference forces `Δf(D,c)` constant in `D`, and the upper comparison `Φ≤bC`
on a fixed skeleton with `D_u→∞` forces that constant to be zero. My checker constructs the normalized score and
confirms `Δf=0` at interior states, then builds bare paths of increasing half-length with the founder in the middle
and all interior cores `(2,2)` (no promotion, zero gain increment). The drift stays constant (`1` in the run, from
the two terminal promotions) while `Φ` grows (`4→12`). Since `sΦ≥saC→∞`, no `s>0` can satisfy the identity.

Planted growth: replacing `Δf=0` by `Δf=k≠0` makes `Φ/C` unbounded, so the comparability hypothesis fails and the
contradiction no longer follows automatically — this confirms boundedness is load-bearing. A planted omission of
the `p(z)A_x(y)` cross term is also rejected by the mixed-identity check.

## 4. Verdict on edge-obstruction.md

The three steps are locally valid on every reachable skeleton I tested, including `c=2` and a distinct founder.
I found no counterexample and no gap in the material under review. The conclusion is a **non-existence** statement
for a pointwise eigenfunction of this restricted form; it does not assert bounds on drift inequalities, and the
final paragraph's application to the conditional limiting-amplitude function depends on that function satisfying
the stated comparison `aC≤Φ≤bC` on all seeds, which is an external input I did not verify. Scope note: I did not
rerun any production sweep (none is required here) and did not attempt a general proof beyond the model expansion
above.

## 5. Separate audit: single-vertex.md

Separate scope: `single-vertex.md` (`closure/single-vertex.md`) states two obstructions — (i) an additive
single-vertex score `F=Σ_core f(D,c)` satisfying `LF=sF` is identically zero unless `s=0` (constant) or `s=1`
(multiple of the total vertex count); (ii) a degree-clamped score depending on the whole rooted core tree, if a
pointwise eigenfunction at `s≠0`, vanishes on every all-high seed with at least two cores.

I re-derived (i). Decomposing the gain rate into its local part `D_x−c_x` and the neighbour part
`Σ_{w~x}1/D_w`, the connection contribution is `Σ_{core {x,y}}(a_x p_y+a_y p_x)`. The four-vertex interchange
(`x` attached to `u`, `y` to `v`, then swap) gives `(p_u−p_v)(a_x−a_y)+(a_u−a_v)(p_x−p_y)=0`, which my checker
verified as an exact identity against a direct recomputation of the connection sum. Setting `x∼u,y∼v` forces
`a_u=a_v` for distinct reciprocals; a third degree handles equal reciprocals; fixing two interior states with
different reciprocals forces `a` constant. Then `f(D,c)=kD+g(c)`, `F=kn+Σ_x h(c_x)` with `h(c)=g(c)+k(c−1)`
and `n` the **total** vertex count (cores plus leaves), and `LF=kn+Σ_x(1−c_x/D_x)w(c_x)` with
`w(c)=h(c+1)−h(c)+h(1)`. I verified this decomposition exactly on a concrete seven-core tree. Growing leaves at
one vertex gives `k(1−s)=0`; varying a single `D_x` gives `w(c)=0` for all positive core degrees, hence
`h(c)=(2−c)h(1)` and `F=kn+2h(1)`, yielding the stated classification. Statement (ii) follows the same pattern:
on an all-high seed the clamped score is unchanged by gains, `j_x` is independent of old `D_y`, and varying one
old vertex's leaves forces `j_x=0` and then `F=0`. Verdict for this file: CORRECT within its stated scope. I did
not test the reachability-by-reversing-leaf-deletions remark beyond noting that every tree I used is so reachable.

## 6. Checker, commands and planted errors

From the repository root, one CPU, under one second:

```
python3 -B round2/checks/edge_obstruction_hy4.py
```

Observed, ending `INDEPENDENT REFEREE CHECKER: PASS`:

```
mixed-branch identity PASS leaf cases=24 skeleton cases=4 (c=2 edge, distinct founder f_o arbitrary)
role-dependent differences + commutation PASS alphaA=-18.7 alphaB=-6.67
bare-path contradiction PASS Phi=4..12 drift=1 constant
single-vertex interchange + classification PASS identity=0
planted alpha!=0 rejected: joint difference=-11/400, alpha from two branches -11 != -2
planted drop-branch / alpha!=0 / unbounded-growth errors rejected
```

Planted errors rejected: (a) dropping the `p(z)A_x(y)` cross term; (b) an `α≠0` edge score, rejected both by a
nonzero joint difference and by inconsistent commutation slopes; (c) allowing `Δf=k≠0`, which breaks comparability
`Φ≤bC` and so invalidates the boundedness step. Each shows the corresponding step is load-bearing.

OPEN: none.
