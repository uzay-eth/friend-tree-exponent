CORRECT

# Independent audit: geometric core degree and ordinary fringe size

Scope. Target `other-marginals.md` (`fringe/other-marginals.md`). I take as reviewed its prerequisites
`degree-law.md` (stable-age formula `πf=μE_β∫e^{-μa}f da`, `β=πB/μ`),
`spatial-limit.md` (`fringe/spatial-limit.md`) (`π` is the unique eigenmeasure and is the almost-sure local
limit), and the operator/branch estimates. I audited only the new one-clock construction, the rejection
compensator, the fringe-size tail and the crossing law. My exact fixtures are
`fringe_marginals_hy4.py`. No author file changed; no long computation run.

## 1. One-clock change of measure and the geometric degree law

The proposal clock is rate 1 originally and rate `r=1+μ` under the changed law, with the same acceptance marks.
For a path with `P(t)` proposals the exact Radon–Nikodym derivative is `e^{(r−1)t}r^{-P(t)}`. Stopping at the
`j`-th acceptance `T_j` gives density `e^{μT_j}r^{-(j+R_j)}`; multiplying by the age discount `e^{−μT_j}` and by
`r^j` yields exactly `r^{-R_j}`, which is (3). I verified this stopped identity pathwise on 200 exact rational
paths. Dropping the rejection factor `r^{-R_j}` (a tempting error) breaks it.

Since `R_j` is nondecreasing on `{T_j<∞}` and `r^{−∞}=0`, the quantities in (3) decrease to
`B=E^{(r)}_β r^{-R_∞}`, giving `Pr_π(c≥j+1)~B r^{-j}`, i.e. `Pr(c≥k)~B(1+μ)^{-(k-1)}`, and the point mass
`B·(μ/(1+μ))(1+μ)^{-(k-1)}` by differencing. I checked the tail/point ratio `1−1/r=μ/(1+μ)` exactly.

**Compensated rejections and `B>0`.** I followed the argument in detail. With rate `r` proposals accepted with
probability `b=L/D`, the leaf drift is `L+ρ−rb`; `e^{-t}L(t)` decomposes into a martingale with finite expected
quadratic variation (jump rate at most `D+r`, and `D` is Yule-dominated) plus an absolutely integrable drift
(`|ρ−rb|≤c/2+r`, `E c(u)≤c(0)+ru`). Its limit `W_L` has mean
`m−∫e^{-u}E[rb−ρ]du≥m−r>0` because `rb−ρ≤r` pointwise — I verified this inequality on exact states. On
`W_L>0`, `L(t)≥(W_L/2)e^t` eventually while `c(t)=O(t)` a.s., so the rejection compensator
`∫ r(1−b)dt = ∫ r c(t)/D(t) dt` is finite; stopping at successive integer compensator levels gives finitely
many rejections a.s., hence `r^{-R_∞}>0` with positive probability and `B>0`. The identity
`1−b=c/D` underpinning the compensator was checked exact on a grid. This argument needs only `m>r` and local
finiteness, so it works from every newborn environment.

## 2. Ordinary fringe size

The exact arrow count for a branch of `N` ordinary vertices is `N+p_p−p_v`: every branch vertex sends to its
branch neighbours and its leaves at total rate `1` except the root, whose arrow to its external parent has rate
`p_v`; the only incoming external arrow is `p→v` at rate `p_p`. I verified this by explicit enumeration on two
small branches (including the `N=2` case `1+1/2+1/D_p`) and the sandwich
`n−1/2≤n+p_p−p_v≤n+1/2` over a grid, using `0<p_p,p_v≤1/2`.

Hence `e^{-t}(N+1/2)` is a nonnegative supermartingale and the process is sandwiched on whole paths by the
Yule chains of rates `n±1/2` started at 2. For rate `n+ε`, `e^{-t}(N_ε+ε)` is an `L²` martingale; the excess over 2
is `NB(2+ε,e^{-t})` with PGF `[e^{-t}/(1−(1−e^{-t})z)]^{2+ε}`. I verified the forward equation
`∂_tG=(z−1)(zG_z+εG)` symbolically for `ε∈{0,±1/2}` (exact rational evaluation of both closed forms),
so the limit is `Gamma(2+ε,1)`. With `W_N>0` a.s. from the lower chain, `ne^{-H_n}→W_N` and
`n^μ Pr_π(N≥n)=E[(ne^{-H_n})^μ]→E[W_N^μ]=A_N`; dominated convergence is supplied by the `L²` maximal bound.
The Gamma moment identity `E[Gamma(α)^μ]=Γ(α+μ)/Γ(α)` gives the bracket (4) with `α=3/2,5/2`; my fixture
evaluates it for the generic rational `μ=1131/2000`, giving
`A_N∈[1.1611,1.5981]`.

## 3. Crossing law

The stable-age crossing identity `μ Pr_π(N≥n+1)=E_π[(n+p_p−p_v)1_{N=n}]` follows because the compensator of the
single crossing counts gives `Pr_π(N≥n+1)=E_β∫e^{-μa}λ_n 1_{N=n}da` and `πf=μE_β∫e^{-μa}f da`. I verified
the equality exactly for all three toy chains `ε∈{0,±1/2}` and `3≤n≤11` by computing both products
`Π_{k=2}^n (k+ε)/(k+ε+μ)`. The factor lies in `[n−1/2,n+1/2]`, so the crossing law is consistent with
`Pr(N=n)~μA_N n^{-1−μ}` obtained by differencing (2). The same structure as degree-law.md:35 is intact; no
new closure assumption is introduced.

## 4. Verdict

Each challenged ingredient is correct: the RN factor and its stopped form, the rejection compensator `r c/D`,
the branch rate, the NB/Gamma sandwich and the crossing law. I found no gap in the material under review.
Scope limits: I did not recompute the reviewed renewal/age formula or the numerical value of `μ`, and I used a
generic rational `μ=1131/2000` for the moment fixture — the algebraic identities and asymptotic forms hold
for every `μ∈(0,1)`, so this is independent of any numerical enclosure for `μ`.

## 5. Checker, commands, planted errors

From the repository root, one CPU:

```
python3 -B round2/checks/fringe_marginals_hy4.py
```

Observed, ending `INDEPENDENT FIXTURE CHECKER: PASS`:

```
mu fixture = 1131/2000 (generic rational parameter in (0,1); identities hold for every mu)
one-clock RN stopped identity PASS paths=200
geometric tail / point-mass ratio PASS
rejection compensator r c/D and positivity bound PASS
downward branch rate N+p_p-p_v and sandwich PASS
Yule NB forward equation and gamma moments PASS A_N in [1.1616,1.5996]
crossing law identity PASS
planted change-of-measure / ratio / branch-sign / NB / gamma / crossing errors rejected (7/7)
```

Planted errors rejected: dropping the rejection count from the RN derivative, using the wrong geometric ratio
`μ` instead of `μ/(1+μ)`, flipping the boundary sign in the branch rate, reversing the sign in the NB forward
equation, using `Γ(α−μ)/Γ(α)` for the Gamma moment, and omitting the crossing intensity. Each confirms the
corresponding step is load-bearing.

OPEN: none.
