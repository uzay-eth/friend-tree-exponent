> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Deterministic validation of signed D18 forecasts

**Result.** A fixed-size forecast can be certified by transporting its *signed conditional defects*, with a separate continuation on every seam. An affine reciprocal-degree tail improves the upper forecast on a real donor-driven `P4` covering class. It resolves the sign even at some larger entry degrees, but does not close the all-motif ratio. A constant-error certificate based only on waiting for whole-crown size exit is too expensive at the existing occupation cuts: with a symmetric terminal half-width of `0.001`, the Yule lower obstruction alone requires a crown cap of at least **48,600** to attain `2.931916e−6` total width per core. These are deterministic statements; no empirical error appears in their hypotheses.

## 1. Root ratio: validate its numerator, not a pathwise quotient

At a fixed multiplier `θ`, pair the original charged residuals as `q=δ₋−θδ₊` before aggregating, and track the two physical potentials with the same ownership and rates. The positive production witnesses give

```
π(q) = (μ−a) πΦ₋ − θ (b−μ) πΦ₊ .
```

For the root occupation, the stable-age identity is `π(q)=E_β[μ∫ exp(−μt)q(X_t)dt]`, with the true physical birth law `β`. The corresponding lower and upper continuation bounds must therefore cover **every birth environment or a certified set of birth-bin masses**, not an empirical mixture. The argument applies to bounded local scores; for larger signed charges use the two integrable positive production sums before subtracting. At a trial `r`, positive rational bounds `P_lo≤πΦ₋≤P_hi` and `Q_lo≤πΦ₊≤Q_hi` give a separator. To exclude `μ≤r`, require a rigorous lower bound on `π(q)` greater than `(r−a)P_hi−θ(b−r)Q_lo`; to exclude `μ≥r`, require a rigorous upper bound below `(r−a)P_lo−θ(b−r)Q_hi`. Here `r∈[a,b]`, `θ>0`. A forecast for a selected support must be combined with rigorous disjoint complementary bounds before applying either test.

## 2. Chronological defect certificate and the covering requirement

Let `Y_j` be the projected state at the `j`th monitored jump, `t_j` its time, `R_j` the discounted reward until the next jump, and `H_j(Y_j)` any independently generated fixed-size forecast. Include killing as a jump with continuation zero. If `τ` is the chosen seam, the conditional signed forecast error is exactly the expected discounted sum of the predictable one-step defects, plus the terminal mismatch:

```
defect_j = E[R_j + exp(−s(t_{j+1}−t_j)) H_{j+1}(Y_{j+1}) | F_{t_j}] − H_j(Y_j)
V − H_0 = E[Σ_{t_j<τ} exp(−s t_j) defect_j + exp(−sτ)(V_exit−H_exit)].
```

Subtracting these conditional expectations produces a stopped martingale, not an estimate of the defects. Squared increments, small observed residuals, and independence of *computation* of forecasts do not bound its predictable part. If exterior rates can change between projected jumps, a frozen first-step row is not the true conditional row. The rigorous version uses the continuous-time residual `e=q+(A−s)H` and stopped Dynkin, or a Bellman inequality valid at **each instant** for all admissible predictably varying rates. For a covering class with affine generator and reward in shared rate variables, verify rational endpoint inequalities on the **same face**. Keep deaths and seams as transitions. For signed rewards check the discount endpoints in the cross-multiplied inequality, or include `s` directly in the reward and check any product dependence. Do not independently vary a donor's contribution in the reward and generator. An envelope for `V−H` follows by solving the same two Bellman inequalities with running score `e` and terminal mismatch bounds. Algebraically this is equivalent to validating `H+correction`, and explains why changing forecasts alone cannot improve the exact robust interval on a fixed covering class.

`check_cover.py` uses an independent fixed-size-8 forecast under a fixed **comparison** policy and transports its exact chronological defects forward to cap 12. It agrees exactly with the separate backwards solve: forecast `−.0011808979003`, defect `+.0000030587240`, seam `+.0000016890765`, result `−.0011761500998`. This identity checks the correction bookkeeping. The physical bound below never uses that fixed policy as a model of the real exterior.

## 3. Affine-tail pilot on a persistent donor-driven face

Use the top-owned `P4` occurrence of four consecutive bare degree-two nonfounders. The arm below the fourth vertex sends at `x=1/d`, the parent above the first sends at `p`; keep the same parent rate in the reward and in face exit. Use `θ=1.2415`. The exact normalized row is

```
q = A + Bp + Cx
A = 287589315471350942439/(5·10²²)
B = −.005553926199875,  C = −.02262119047951.
```

`check_cover.py` independently contracts both D18 payloads using `Local.charge`, checks the paired identity on four ports including degree 2,097,152, and touches only selected table pages. The analytic identity is affine in ports by the pinned four-body charge; the port checks are exact arithmetic controls. The face dies at the first gain at a protected vertex, at rate `3+p+x`. While it lives, track the lower arm root by `(d,k)`, whose gain rate is in `[l+1/2,l+(k+1)/2]`, with `l=d−k−1`, and whose promotion rate is `l/d`. This covers arbitrary descendant trees; their current reciprocal sums need not be frozen or independent. Discounts range over `[.565674722414,.565773666188]`, and `0≤p≤1/2` throughout.

At degree above `K`, write `U=u−v/d`. The promotion term cancels; the gain increment is `v/[d(d+1)]`. For the upper envelope take `s=LO,p=0` (the candidate is positive and `B<0`), and `u=A/(LO+3)`. The upper gain bound `g≤d` reduces the required surplus, divided by `x>0`, to

```
u−C−v(LO+3+x+1/(1+x)).
```

Choose `v` to make this zero at `x=1/(K+1)`; it is nonnegative at all smaller `x`, including the limiting test at zero. For the lower envelope the simplest better tail here is constant, the minimum of `(A+Bp+Cx)/(s+3+p+x)` on the tail. It occurs at `p=1/2,s=HI,x=1/(K+1)`; the exact numerator is positive on the tail. One can also use a valid affine lower envelope with `u=(A+B/2)/(HI+3.5)`, `v=(u−C)/(HI+3.5)`, and discard the nonnegative gain term. On this cut it is weaker than the constant lower envelope. Use the better certified side separately.

The finite triangular solve checks both directions at all rational corners, rounded outward. At cap 16 and arm entry `(d,k)=(16,3)`, the conditional contrast is in approximately **`[+0.00039352,+0.00127936]`**; the constant upper tail gave `+0.00152707` (the affine upper tail makes the interval about 22% narrower). At virgin entry `(2,0)`, the same cap gives **`[−0.00159080,−0.00105493]`**. Both signs are thus certified on different attainable entries of the same face. This is an all-time conditional **face episode**, not a bound for the entire crown after face death, nor a new stationary pricing credit for an already-paid support. The sign change illustrates why a pointwise ratio of even the paired rows is a poor continuation model for old donor contexts.

## 4. Lower cost bound for a Yule-only cubature error

For a full planted crown with `n` ordinary vertices, each vertex rings at rate one and only the focal-to-parent sends leave the crown. Exactly `n−1/D_focal+p` is its birth intensity, with `D_focal≥2` and `0≤p≤1/2`. Invisible exterior events may change `p`; the size holding time still has conditional survival squeezed between exponential clocks of rates `n−1/2` and `n+1/2`. Iteration at each size hit yields, from a virgin (`n=2`), for the time `τ_K` to first size `K+1`,

```
B₋(K)=∏_{n=2}^K (n−1/2)/(n−1/2+HI) ≤ E exp(−sτ_K)
    ≤ ∏_{n=2}^K (n+1/2)/(n+1/2+LO)=B₊(K).
```

At cap 20 the exit factor is in `[.211628,.283683]`; at cap 1,024 it is still in `[.0230098,.0316901]`. With no other information about a terminal error in `[-E,E]`, the two constant opposite terminal errors realize uncertainty of width at least `2E·LO·B₋(K)` per core, even with **perfect quadrature before exit**. An exact conservative lower bound at cap 10,000 uses `log(1+y)≤y`, the integral bound on the shifted harmonic sum, `HI<3/5` and `(10000−1/2)/(3/2)<6⁵`: `B₋>3/(5·6³)`. The sharper directed-integer product in the replay places the first crossing of `2E·LO·B₋(K)≤2.931916e−6` at **48,600** when `E=.001`. This is a necessary cost for the uniform-terminal strategy under that error envelope, not a lower bound on the actual D18 terminal error or on algorithms that exploit a better face cover. The product is of power order `K^(−s)`, not exponential in the size cap.

A second obstruction is that deterministic nodes without a verified interpolation/variation envelope give no uniform cubature gain: two bounded functions can agree on every queried cell and take opposite values everywhere else, leaving twice the bound times the unqueried mass as uncertainty. Here such a worst-case bound can be reduced using finite exact polynomial faces, signed defect bounds or certified reward-weighted occupation. Yule control on size alone gives none of these for the donor reciprocal sums. A plausible cheaper extension is to track the most important live donors jointly while contracting unused arms by their verified affine rate enclosures, with Bernstein guard cells at the D18 seams. The finite forecast pilot closes one such face, but no rigorous all-family complementary residual or improvement over the existing positive occupations was established.

## Replay

From this directory:

```sh
python3 -I -B forecast.py
./round2/.venv/bin/python -B check_cover.py
```

The first command uses only the standard library. The second uses the existing read-only D18 operands in `the scratch directory` and already-installed NumPy reader. Both ran on CPU3 at niceness 15. See [`forecast.log`](forecast.log) and [`cover.log`](cover.log).
