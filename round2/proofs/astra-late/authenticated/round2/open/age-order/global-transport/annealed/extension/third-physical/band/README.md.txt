# Continuous all-label tail order on [1/3,1)

For every old tree and every older label **`20≤b≤39`**, the two next births have strict conditional tail order throughout **`1/3≤x≤1/2`**. Together with the [half-threshold theorem for `b≥20`](../../first-sender/refined/README.md) and the [one-third bound for `b≥40`](../../first-sender/joint/README.md), this proves strict pre-placement order for all `b≥20` on `[1/3,1)`. Averaging and using only [Fable's paired receipts](../../../../../../../attacks/fable5/README.md) for `2≤b≤19` gives the all-label result. Founders are tied; the single `20→21` receipt is unnecessary. This proof is about threshold tails and does not use moment order or assume monotonicity of the gap in the threshold.

## 1. Enclose the physical continuation on whole threshold boxes

Keep the [frozen stock-sensitive Bellman comparison](../README.md), its 14 forward sizes, and its local signed ledger. For all thresholds `x` in a closed box `[l,h]`, put `q=1−x`, `q_+=1−l`, and `q_−=1−h`. All row values are normalized pointwise by `q^(n−1)` at size `n`. The pending locked cap has shape `a=r=1/k`, `B=n−a`, and equals

```
[Γ(n)/(Γ(a)Γ(B+1))] q^(1−r) B ∫₀¹ v^(B−1)(1−qv)^(r−1) dv .
```

This normalized bound is increasing in `q`: both the power and the integrand increase. Use `q_+`. For the promotion term from size `n` to `n+1` use `a=1+r`, `B=n+1−a`. The extra Bellman factor `q` combines with its normalized beta cap to give the power `q^(1−r)`, and now the integrand contains `(1−qv)^r`. If `q_+≤1−r`, their product increases for all `v∈[0,1]`, because its derivative has the sign of `1−r−qv`. In that case evaluate both at `q_+`; otherwise use the power at `q_+` and the integral at `q_−`. This is a box bound for the **discounted promotion term**, not a claim that its two factors increase separately. In the nonpromotion terms use `q_+` for the single Bellman factor. The two possible extremal physical parent-gain rates are unchanged. All probabilities in the recurrence are nonnegative, so taking the minimum with the box pending cap remains safe. Backward induction therefore gives one upper row for the *whole* box.

The integral series and tail estimates are exactly the ones in the frozen proof. The checker slightly tightens the reflection constant using the alternating degree-eight upper polynomial for `sin(z)/z` at `z=3.14159265r`. This polynomial decreases on `[0,π/2]`: its derivative divided by positive `z` is `−1/3+z²/30−z⁴/840+z⁶/45360<0` there (even dropping the negative fourth-power term and using `z≤2` suffices). Machin's formula verifies the rational lower bound on π. The factor `q_+^(1−r)` is rounded upward by integer power comparison, and beta and Bellman rows are rounded upward on the frozen grid.

## 2. A sharper fixed-label lower at b=20

The frozen one-family construction has no threshold-dependent constant once normalized by `q^b`. At each label `b`, conditional on the red endpoint exceeding `x`, let `t=b(θ−x)/(1−θ)`. It has survival `(1+t/b)^(-b)`, independent of `x`. The first marked blue sender index is `f=floor((b−1)(1−v)/v)` for independent uniform `v`. Conditional on `f,t`, the uncapped scaled charge `a` has survival `(1+a/b)^(-f−1)` on the range that affects the payoff; the physical atom at the cap is above `t`. These are the [frozen family laws](../../first-sender/refined/README.md). Use the exact ceiling on the subsequent loss coefficient rather than its logarithmic upper bound:

```
C_f = b Σ_(j>f) 1/[(b+j)(j+1)] = b/(b−1) Σ_(k=f+2)^(f+b) 1/k.
```

The normalized lower payoff remains `(1−C_f(1+a/b)/(t−a))_+` on `a<t`, and zero otherwise. At fixed `b=20` there is no need to bound `f` uniformly over larger labels. Take `f` at the right edge of each of 320 equal uniform bins: decreasing `f` increases both the loss distribution and its penalty, so this is conservative. Use exact rational bin weights for `t` on a `1/24` grid through 8, rounding down and using the left edge. For `a`, take the exact CDF at the `1/24` grid endpoints rounded down, place increments at the right edges, and leave the unretained mass at infinity. This stochastic discretization lowers the decreasing payoff. The calculation certifies **`M=.2124`** for label 20 throughout the band. For 21–39 use the already replayed **`.211`** from the frozen `third-physical` note. No credit for physical first misses is needed.

## 3. Continuous ledger and splice

At `b=20,21` the boxes have endpoints `1/3, 41/120, 7/20, 3/8, 2/5, 1/2`. For `22≤b≤39` one box `[1/3,1/2]` suffices per label. The frozen signed ledger is evaluated on the box upper rows, with the same chain constants and local degree profiles; it therefore pays every actual adverse debit throughout each box. All boxes are closed, with consecutive endpoints shared. In particular the endpoint `b=20,x=1/3` is proved on every old tree, rather than inferred from an annealed transport receipt. The splice at `1/2` is also inclusive. The conditional claim at `x=1` would be equality, since both strict tails vanish there.

```
nice -n 12 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/band/check.py
PASS b=20 serial coefficient > .2124
PASS b=20 over [1/3,1/2]
...
PASS b=39 over [1/3,1/2]
```

The replay authenticates the frozen inputs and uses only rational arithmetic, local degree profiles, and a scalar family quadrature. There is no old-tree enumeration.
