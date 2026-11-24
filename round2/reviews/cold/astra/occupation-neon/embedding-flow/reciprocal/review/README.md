# Independent review: reciprocal polynomial closure

**Verdict:** The proposed common-denominator representation closes for the full singleton and edge balance of degree-`m` polynomial tests, **per high coordinate**, under the conventions below. `N≥2m+2` is sufficient for these rows, including genuine reciprocal-flux projections tested at degree at most `m`. It is not a license to use every stored degree-`N` Bernstein coordinate as a test of every balance or projection. A one-sided overflow bound remains a relaxation, not a moment equality. This review does not run the LP or certify a D18 endpoint.

## 1. Measures and cancellation

On a census cell `C` let `I_C` be its exposed high slots; use `y_i=H/D_i`, `a_i=1+y_i/H`, and `W_C=∏_(i∈I_C) a_i^m`. Denote raw labelled-occurrence measure by `M_C`; a raw port measure `F_C^v` is the same occurrence measure weighted by the **sum of actual exterior inverse sender degrees into role `v`**. Store positive tensor Bernstein coordinates of `M_C/W_C` and `F_C^v/W_C`, using **the same** `W_C` for all ports on `C`. Then a raw polynomial integrand `G` is read by testing the stored measure with `W_C G`. No extra `W` appears in the physical rate, and storing a donor-count measure instead of a physical flux would require a separate reciprocal factor at that donor. The basis coordinate sum is the *weighted* mass, not raw mass. Read raw mass by `W_C`, whose per-axis degree is `m`.

At a gain of slot `i` within its high band, `y_i'=y_i/a_i`. For the univariate monomial `y^j`, `0≤j≤m`, the clearers are

```
U_j(y) = y^j(a^(m−j)−a^m),
V_j(y) = H U_j(y)/y,
a=1+y/H.
```

`U_0=V_0=0`. In fact `U_j` is divisible by `y^(j+1)` for `j≥1`; its degree is at most `2m`, and `V_j` has degree at most `2m−1`. For a tensor test, the untouched factors are `a_l^m y_l^(j_l)` and have degree at most `2m` in each other exposed coordinate. Decompose the receiver gain rate into `D_i−k_i−1`, its included-edge reciprocals and its exterior port. The first term uses `V_j−(k_i+1)U_j` on mass; the included reciprocals multiply `U_j` by `1/D_l`, i.e. `y_l/H` at a high included sender, or a real constant at a low one. The exterior term uses only `U_j` on the **stored port**, not `V_j` and not a product of means. Thus the maximum degree is at most `2m+1` per coordinate. With multiple high roles the denominator is a tensor product, not a total-degree budget. A promotion at fixed high degree multiplies a polynomial by `b=1−(k+1)y/H`, and the shifted arity does not change the weight factor in that source. This also costs at most `2m+1`. Discount costs at most `2m`. These degree counts use **differences** in the same band before integrating, never the separate formally nonintegrable terms `D·P(y')` and `D·P(y)`.

## 2. Full weak edge and singleton rows

For a test `f` on `(p,x)` with `p` the mother of `x`, read its value as zero off its specified cells. Let `u^+` increment degree and `u^b` increment `k` at fixed degree, and `J=(2,0)`. Write `g^0_p=L_p+1/D_x`, `g^0_x=L_x+1/D_p`, `b_i=L_i/D_i`. In the convention `η=π/s`, with one occurrence measure per edge, the row is

```
Σ_C ∫ [s f − g_p^0 Δ_p^+ f − g_x^0 Δ_x^+ f
             − b_p Δ_p^b f − b_x Δ_x^b f] dM_C
 − Σ_C ∫ Δ_p^+ f dF_C^p − Σ_C ∫ Δ_x^+ f dF_C^x
 = Σ_B ∫ b_p f(p^b,J) dM_B(root).
```

Each difference is successor minus source, so a promotion of the mother appears **both** as a transition of every existing edge and as exactly one newborn edge; the newborn edge sees the *postpromotion* census of the mother. The mother's exterior port includes its own mother and its other children; the child's exterior port consists only of its children. These can be split into distinct physical ports. No extra factor of `k` multiplies the newborn term. The singleton equation is the same without any included donor or second endpoint, with RHS `f(J)` since `η(b)=1` at the physical rate. For a band-specific test, incoming and outgoing terms are summed over **source** cells: multiply each by that source's weight. At `D=H−1` an incoming gain to the high band evaluates the successor test at exact `y'=1`; only the unchanged high roles contribute variable factors. The rate uses old stock and old donor reciprocals. With a low-only test, this same transition is a loss. There is no reverse predecessor from `H` to `H−1`; tests on both bands retain both evaluations if nonzero. If more internal degree bands are introduced, cross-band strips must be represented explicitly (an arbitrary indicator of `y'` landing in another band is not a polynomial); the `H−1→H` seam itself is a finite census. If `H=2`, the newborn high coordinate is also the evaluation `y=1`.

If the source `k≤K` promotion exits to `k>K`, a test extended by zero into overflow has an exact loss; no old overflow vertex can return to a retained arity. But an overflow **mother** can still produce a new retained singleton, or contribute an exterior reciprocal on some other retained pattern. Do not drop either contribution. Exact tests depending on overflow `k` or rates with a `k` factor require moments of it or separate bounds; `p≤1/H` alone does not price these events. A complete low sector includes all feasible `0≤k≤D−1` for `D<H`; if a low-to-high seam enters overflow directly, route it there rather than a nonexistent retained cell. For the stated balance on finitely many `k` cells the high gain integrand is bounded after cancellation, including as `y→0`; compact-degree cutoffs of a high polynomial that has a nonzero limit at `0` should not be used without estimating the boundary flux. One rigorous route is discounted Dynkin from a physical promotion along a tracked endpoint or edge, as in the birth-conditioned edge identity (`identification/time-cuts/README.md`), with stopping at growing sizes and bounded discounted tests. For a raw edge sum use one tracked child per edge by mass transport. Conditional first moments along finite histories and bounded rates for the *test increments* (after cancellation) justify the integral, and terminal discounted test vanishes. The unnormalized birth term above follows directly by averaging the promotion census under the physical birth-Palm law.

## 3. Projection and weight changes

For a base cell `B` and port with `c` eligible new roles (one for a mother, `k_v−c_v` for a side child), let `E` run over the resolved extension cells. Incidence multiplicity is labelled; no division by symmetric embeddings in these projections. For any **nonnegative** overlap polynomial `Q` of degree at most `m` per coordinate, define

```
M = ∫ Q dM_B,    F = ∫ Q dF_B,
S = Σ_E ∫ Q∘proj dM_E,    R = Σ_E ∫ (Q∘proj) p_new dM_E.
```

All integrals are raw; their LP coefficients are respectively `W_B Q` on base mass/flux, and `W_E (Q∘proj)` and `W_E (Q∘proj) p_new` on extension mass. If omitted donors have `0≤p_new≤u` (e.g. `u=1/H` on high-arity overflow), impose

```
0 ≤ cM−S,         R ≤ F ≤ R+u(cM−S).
```

If known `p_new≥l` on the omitted class, strengthen the lower to `R+l(cM−S)≤F`. For mother transport at `c=1` the row is apex extension. For child transport a new role fills one eligible slot. Separate classes with different bounds or exact cells; never use the bound for high overflow to also cover omitted low donors. If all classes are accounted for, the identities are exact. If `c` varies in the base cell, use `∫Q c dM_B` rather than a nominal constant. Likewise root and edge marginals multiply by their own weights; a weight contributed by an eliminated high role **does not disappear** just because the projected test is constant in that role. For example, projecting weighted base moment `∫Q d(M_B/W_B)` instead would read extension measure with `(W_E/W_B)(Q∘proj)`, polynomial only when the overlapping weight factors match or cancel in the required direction. Raw tests avoid this issue. If independently chosen denominator exponents differ across overlapping roles, check the per-source degree budget again: a shared raw overlap test of degree `d_test` requires at least `m_source+d_test` on that source axis before reciprocal or birth-rate factors. A transition **between** storage measures is not needed for the weak rows; each source is evaluated on its own measure.

## 4. Uniform rational repair of numerical duals

For `min c·x`, `Ax=b, Gx≤h, x≥0`, round arbitrary proposed duals to exact `y,z`, with `z≤0`. Set `r=c−Aᵀy−Gᵀz`. If each nonnegative column belongs to one group `g` with `Σ_(i∈g)x_i≤B_g`, then

```
min c·x ≥ b·y+h·z − Σ_g B_g max(0,max_(i∈g)(−r_i)).
```

Use separate `−c` for a maximum. Bernstein coordinates sum to one on the cube **for each tensor**, so summing coordinates costs weighted measure mass once, not a factor `(N+1)^r`. Weight `1/W≤1` gives a valid bound by raw mass. For a fixed topology grouped across all degree censuses, including overflow, explicit finite budgets follow from arity moments without a degree cutoff. In an embedding of `r` roles charge the occurrence to a role of maximal arity `k_*` (break ties by role). Reverse to that role via physical mass transport. Its ancestors for all up-steps are unique and the at-most `r−1` downward choices have at most `k_*` each. This gives mass ≤ `z r Eπ(1+k)^(r−1)`; a group summing all exterior ports per occurrence pays at most `r(1+k_*)/2` more, so a valid budget is `z r² Eπ(1+k)^r/2`. Multiple discount splits may share the budget only if nonnegative splits add to the original moment; duplicated measures need separate budgets or a multiplicity factor. By the pressure-memory comparison (`pressure-memory/README.md#an-effective-joint-degree-and-child-count-tail`), `Eπ(k)_j≤j!/s0^j`, so expand `(1+k)^r` with nonnegative coefficients in falling factorials and use `z≤1/s0`. These are uniform in the degree cutoff and basis size for fixed order. Auxiliary slacks/defects need their own proven budgets. A repaired objective bound is valid on the physical image even if a bare LP omits some budget constraints, but do not claim validity for arbitrary LP feasible points unless the constraints apply there. Optional extra controls/test rows must retain exact rational scaling in the dual repair.

## 5. Independent replay

From `reciprocal/review/`, run `round2/.venv/bin/python -Bu check.py`. SymPy checks symbolic clearing and tensor degrees for 90 two-slot monomials, then contracts exact Bernstein tensors against five physical edge fixtures (both directions of the seam and high/high states). Planted mutations change the seam value, use a wrong high-shift denominator, misweight the port, reuse the newborn target weight on the root source, omit the eliminated role's denominator factor in the projection, or equate all flux with just its resolved part. All are detected. Grouped dual repair is checked on a rational small instance. No project solver, optimization, or binary payload is used.
