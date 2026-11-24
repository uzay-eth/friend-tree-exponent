CORRECT

# Independent audit: spatial eigenmeasure, one-ended local limit, eigenmeasure uniqueness

Scope. Primary target `spatial-limit.md` (`fringe/spatial-limit.md`). Its prerequisites
`renewal.md` (`fringe/renewal.md`), `almost-sure.md` (`fringe/almost-sure.md`) and
`operator.md` are taken as reviewed, as the brief directs; I did not re-prove their renewal
estimate, martingale positivity, or existence of the probability eigenmeasure. I audited only the new spatial
arguments: the old-cell boundary, the constants behind (1), uniformity in environments, empirical downward
embeddings, and ancestor unbiasing. My own exact fixtures are
`spatial_limit_hy4.py`. No author file was changed; no long computation was run.

## 1. Old-cell boundary with arbitrarily many old neighbours

The only possible leak is at the initial root `o`, whose core degree `c_o` counts old neighbours outside the cell sum.
I checked each step:

* The birth mark of a new child of `o` is `a_Y=1/m(D_o)` and the current reciprocal is `p_o=1/D_o`; with
  `m(D)=min(D,20)` we have `a_Y−p_o≥0` for every `D_o`, so the immigration terms of
  `spatial-limit.md:16-18` are nonnegative even when `D_o` is huge.
* Old neighbours have their own cells and their promotions create children of *them*, not of `o`. They therefore
  change neither `c_o` nor the child sum, and never enter `I` as new children. This is the point of
  `spatial-limit.md:25`; I found no route by which a high-degree boundary leaks into `Z`.
* `Z(0)=h_0(c_o)≤2` for every environment because `h_0(0)=2` is the global maximum of the renewal weight; this
  is independent of how many old neighbours `c_o` counts.
* The one-event inequality (renewal (4), adapted to the capped factor `1/(1+d·m(D))`, minimum `1/4` for `D≥20`)
  holds for every `0≤c≤D`: I verified all `2≤D≤100` exactly, and for `D≥20` the factor is `1/4` with left side at most
  `1/4<q·min h=(9/25)(7/10)=63/250`, so the check extends analytically to all degrees. Hence `E Z(t)≤2e^{qt}` and uniformly `I_t W≤1/d+2e^{qt}` as in (2).

This closes the boundary challenge: arbitrarily many old neighbours neither inflate `Z(0)` nor create spurious
new-child summands.

## 2. The `4^k` descendant bound and the cylinder estimate

With capped marks `a=1/m(D)≥1/20`, the edge factor `z(a)=1/(1+d m(D))` has minimum
`z(1/20)=1/(1+3)=1/4`. With `h_0≥7/10` and `a_S+d≤1/2+3/20=13/20`, factoring the path to every depth-`k`
descendant gives `W(S)≥(1/(a_S+d))(1/4)^k(7/10)Z_k`, so `Z_k≤(13/14)4^k W≤4^k W`. I verified this on 800 random
capped-mark planted fringes for `k=0,1,2,3` (worst ratio `Z_k/(4^k W)=0.591`). The constant `4` is therefore
correct and not an artifact of the uncapped scheme.

The remaining bookkeeping also checks out. For `j≤k`, `t^j/j!≤2^k t^j/(2^j j!)`, so
`Σ_{j≤k}t^j/j!≤2^k e^{t/2}` as claimed. A cell member at depth `>k` has a unique depth-`k` ancestor, so it is read
exactly once by `F_{k,f}` at that ancestor; members at depth `≤k` are absorbed in the `2^k e^{t/2}` remainder.
With `‖F_{k,f}‖_W≤4^k‖f‖_∞` and the uniform error constant of almost-sure.md (3), dividing by `e^{μt}` gives
exactly (1) with coefficient `ν(F_{k,f})`.

## 3. Uniformity in environments

`h(G,o)=d_2∫_0^∞ e^{-μu} I_u H du` is bounded above uniformly because `I_uH≤(1/d+2e^{qu})·sup H` and `q<σ<μ`.
It is strictly positive at every input by the finite scripts of `spatial-limit.md:35`: a leaf-bearing root promotes
in finite time; a bare root is first given a leaf by one of its finitely many (locally finite state space) old
neighbours ringing and selecting it before any change at that neighbour. Conditioning on finitely many clocks has
positive probability in every larger truncation, so the script is legitimate in infinite environments. The
convergence in (1) is therefore uniform in the starting graph with an absolute constant, regardless of `c_o`.
I note that pointwise positivity of `h` is all that is needed; no uniform lower bound on `h` is claimed or used.

## 4. Eigenmeasure identification, uniqueness, empirical embeddings

Applying the operator.md probability eigenmeasure to (1) gives `π(f)=π(h)ν(F_{k,f})`; from `f=1` and
`π Q_n=e^{μn}π` one gets `π(h)=1`, hence `ν(F_{k,f})=π(f)` for every radius, and the same argument for any other
probability eigenmeasure `π̃` forces `π̃(h)=1` and `π̃(f)=π(f)` on all cylinders — uniqueness follows. From `f=1`
and `h>0`, `ν Z_k=1` for all `k`.

I verified the combinatorial core of the empirical embedding step. For a planted type `S` and fringe `X`,
`E_S` (downward root-to-root embedding, cores to cores with matching marks, `S`-leaves to any child) is increasing,
as is `E_S∩{|X|≥|S|+1}`; an embedding forces `|X|≥|S|`, so `1_{E_S}−1_{E_S∩{|X|≥|S|+1}}=1_{|X|=|S|,E_S}`, and at
equal finite size the embedding is a rooted isomorphism, i.e. `X≅S`. I checked this identity exactly on 311
small rooted pairs with a proper rooted-isomorphism test. The `3(r+1)e^{2u}` one-leaf coupling bound follows the
almost-sure.md (6) infection calculation with each dirty core contributing at most `r+1` vertices on its path to
the founder; this is consistent with the coupling described there. The extension of the uniform finite-seed
expansion (4) is legitimate since both indicators are bounded by `1<W`.

Ancestor unbiasing: each noninitial core has exactly one depth-`k` ancestor, so the empirical mass of a
distinguished ancestor type is at most one, and the total limiting mass is `ν Z_k=1`; hence no mass is missing
except near the founder (`o(1/C_n)`). I verified the underlying counting identity on 80 random trees. One end
follows: the sampled vertex has an infinite ancestor chain, every finite ancestor has a finite fringe a.s. (the
residual process is subcritical, `ρ<μ`), and any ray leaving the chain lies in one of those finite fringes.

## 5. Verdict

The four challenged points all survive. I found no gap in the material under review. The statement that the
conditional law of the marked tree at a uniform core converges almost surely to `π` inherits the sample-path
machinery of almost-sure.md, which I take as reviewed.

## 6. Checker, commands, planted errors

From the repository root, one CPU:

```
python3 -B round2/checks/spatial_limit_hy4.py
```

Observed, ending `INDEPENDENT FIXTURE CHECKER: PASS`:

```
one-event capped-label inequality PASS pairs=10099
Z_k <= 4^k W PASS fringes=800 worst_ratio=0.5909
increasing-event embedding identity PASS cases=311
ancestor counting / unbiasing PASS trees=80
uniform Z(0)<=2 with arbitrarily many old neighbours PASS
planted uncapped-factor / 4^k-bypass / drop-size / old-as-new rejected (4/4)
```

Planted errors rejected: (a) using the uncapped minimum `z(1/2)=10/13` instead of `1/4`, which breaks the capped
one-event inequality and admits a `Z_1/W` ratio above `13/10`; (b) replacing `4^k` by the tighter `(13/10)^k`
bound, violated on a sampled fringe with mark `1/20`; (c) dropping the size condition in the embedding event,
which breaks the isomorphism identity; (d) counting old neighbours as new children, which destroys the uniform
`Z(0)≤2` bound for large `c_o`. Each confirms the corresponding step is load-bearing.

OPEN: none.
