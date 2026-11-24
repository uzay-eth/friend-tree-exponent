# Full-horizon stock lifetime: independent check

**Signed-leaf reversal audit:** The occupation lemma and all unit-polynomial benchmark numbers are unaffected; they use only one-tree clocks and no signed profile. Application to a concrete signed expansion is conditional on independently verified coefficient bounds. In particular, any `M=max|h−1|` used to instantiate `amax` or `C2` must be re-certified if it came from the withdrawn `zero/surrogate` enclosure. The helper does not certify such a coefficient or a resulting scalar interval.

The proposed stock-lifetime inequality is valid. The immortal skeleton is a Yule process at rate **`λ=1−1/d` per immortal lineage**, conditional on its initial count. No independence between that skeleton and the physical core-child count is needed. `lifetime.py` evaluates the proposed bound and two refinements using rational upper endpoints.

## Physical lower coupling and immortal rate

Start with degree `d≥2`, `l` leaves and `k=d−1−l` core children. Match one forest particle to each initial physical leaf. A matched leaf sends to the root at rate one; that arrow adds a physical leaf and a matched forest particle. Use a separate rate `1/d` death proposal clock for each matched particle. Kill the particle at every proposal, and promote its matched physical leaf with probability `d/D`; other physical leaves have promotion rate `1/D` each. This recovers the marginal physical selections (each leaf has promotion hazard `1/D`, and nonleaf root selections can be supplied separately at their true hazards). Extra incoming arrows, including any live parent input, only add ordinary leaves and can stay unmatched. A forest death not accepted as a physical promotion merely releases the match. New forest particles receive independent clocks, and the forest rates are constant regardless of acceptance decisions. Thus `D_t≥d` and `D_t≥L_t≥Z_t` under this coupling, where `Z` is a forest of `l` independent linear birth–death lineages with rates `1,1/d`.

Each initial lineage survives forever with probability `λ` by the branching extinction equation (the smaller root is `1/d`), so the count `R` of initially surviving lines is binomial. Keep only individuals whose descendant line survives forever. For one retained individual, a split within `dt` has **two** retained continuations at conditional probability `1·dt·λ²/λ + o(dt)=λ dt+o(dt)`. Splits with just one retained continuation do not change the count; deaths without a surviving continuation are excluded. Conditional on the retained branches, their future evolution factors by the branching property with survival conditioning on each branch. This identifies the retained count as a pure-birth Yule process given `R`, started from `R`, including the identically zero case. This conditioning is internal to the forest; it does not preserve the physical promotion law or imply any forest/physical independence.

Given `R=j>m`, the Yule distribution with `p=exp(−λt)` is `Pr(Y_t=n)=C(n−1,j−1)p^j(1−p)^(n−j)`. Dividing the binomial kernel by `(n−1)_fall(m)` and summing the remaining negative-binomial series gives **exactly**

```
E[1/(Y_t−1)_fall(m) | R=j] = exp(−m λt)/(j−1)_fall(m).
```

Since `Y_t^m≥(Y_t−1)_fall(m)`, this bounds its inverse power. On `R≤m` use `D_t≥d`. Physical root promotions have predictable total intensity `L_t/D_t≤1`; in the physical marginal they can be dominated by a unit Poisson clock via thinning (or added rejected proposals). This proves `E[G(k_t)^2]≤E[G(k+Pois(t))^2]` for nonnegative increasing `G`. It is a separate marginal bound; the proof does not put this Poisson count independent of `R` or `Y_t`.

For `m=2r`, split the expectation over `R>m` and its complement, and apply Cauchy–Schwarz at each time. This yields `sqrt(A)e^(−rλt) sqrt(E G²)` and `d^(−r)sqrt(Pr(R≤m)) sqrt(E G²)` respectively. Apply Cauchy again with the discounted time measure, using the corresponding rate `s+rλ` or `s`. This is the proposed inequality, with `P_z[G²]=∫₀∞ e^(−zt) E[G(k+Pois(t))²]dt` (unnormalized). The first marginal moment gives the independent crude fallback `d^(−r) P_s[G]`.

## Rational refinements and API

- `staged(s,d,l,r,g)` uses the full binomial distribution. For each `j` and time lane, optimize among the powers `h=0,…,min(2r,j−1)` (always including zero), since `D^(−2r)≤d^(−2r+h)Y^(−h)`. Take the minimum only after independently integrating each upper bound on that lane. Sum over `j`, then apply weighted Cauchy on the lane. Also take the lane's first-moment bound and the proposed full-horizon bound. For a constant `G`, use the corresponding direct inverse-`r` bound without Cauchy.
- `cutoff(s,d,l,r,g)` is the fast large-degree variant. For `R0=l−f+1>2r`, replace `A` by `1/(R0−1)_fall(2r)` and bound `Pr(R<R0)` by `C(l,f)/d^f` using the binomial failure factorial moment. Optimize up to 64 cutoffs and take the minimum with the crude first moment. No full binomial or large state table is allocated.
- `proposed` is the independently checked original expression. `g` is a nonempty tuple of nonnegative rational coefficients in ascending powers of `k`. Poisson moments use Stirling numbers; Laplace moments, exponential enclosures and square-root upper endpoints are rational. Output table decimals are diagnostics; the API returns exact rational upper bounds.

The coefficient checker's revised cubic splits as `amax*(1+3k)+C2*(k+2)(6k²+22k+12)`; use `(1,3)` and `(24,56,34,6)` respectively. Price the two terms separately or pass the combined coefficients; the smaller of the two certified answers is valid. A signed remainder bounded in magnitude by this stock envelope can be paid with the resulting radius, or separate one-sign coefficient envelopes can be passed through independently. These are full-horizon intrinsic occupations, not zero signed feedback bounds by themselves.

```sh
taskset -c 13 timeout 10s python3 -Bu round2/influence/zero/stock-lifetime/lifetime.py
```
