# Neon-2 Claim 6 — independent analytic tail audit

**Pass for the pinned tail and its finite/analytic join.** The signed singleton cap for every physical degree `D >= 2000` is `0` ticks on the lower witness and at most `-868735957253/133466700` ticks on the upper witness (ceiling `-6509`). Both caps agree exactly with the pinned receipt. The `1999 -> 2000` gain jump passes for every source arity on both sides. This audit is conditional on the finite motif row and finite comparison certificates; those belong to the separate production run.

## Inputs and reproduction

Read-only inputs: `round2/.work/neon2-claim6-7M6Ys3rJ/source`, detached at `b04cc50c1e391e0420eb7b006ac781d484c68543` with empty `git status --short`, and the two adjacent 1,588,521,616-byte copies. `sha256sum quintet_lower.bin quintet_upper.bin` in their parent gave:

```text
81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986  quintet_lower.bin
8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383  quintet_upper.bin
```

From `[local path]`, run `python3 -B round2/attacks/neon2-quintet/tail/check.py`. It uses read-only random access, Python `Fraction`, and its own compact layout reconstruction. It does not import or run the production verifier. It reads the pinned JSON receipt for result comparisons and recombination of the finite inputs. The separate SHA command above supplies binary authentication; the rational script does not rehash the payloads. It touches only the seam and relevant tail slices of the tables. Signed code inspected: `verify_quintet.cpp` lines 118–158, 491–555, 1160–1168, 1346–1355; `quintet_format.hpp` header constraints and endpoint layout; `verify_tail_bernstein.py` lines 152–217; `cert_check.py` lines 61–94; and `PROOF.md` §§1–4.

## Physical domain and assignment

Use ticks throughout, with `S=2^40`. Let `s=+1` for upper and `-1` for lower, so both tests are **upper bounds on** `s*(generator - lambda*weight)`. The signed trials are `56574828/100000000` upper and `56570237/100000000` lower. At a physical core let `D` be total degree, `c` its number of core neighbors, and `L=D-c`, including the single-core star `c=0` and the no-leaf hub `c=D`. Gain from each leaf is one, and from each core neighbor is its actual reciprocal degree. Promotion has rate `L/D`, leaving `D` fixed and increasing `c` by one. For `L=0` there is no bud evaluation at an invalid physical state. In the polynomial its formal extension at `q-x<0` is multiplied by zero when `q=0`.

The analytic singleton and seam are beyond `d2=180`, so there is no centered singleton gain and no outgoing `f -> g` correction. Every internal support cutoff for the higher shapes is at most `d2`. In a new connected motif of at least three vertices containing the newborn `J=(2,1)`, the promoting vertex must be internal, hence none is supported at the seam or in the tail. Existing higher motifs **can** contain this vertex as an endpoint; their increments still belong to their own rows. Likewise every existing edge increment stays in the edge row. Singleton promotion owns precisely the singleton change plus the newborn singleton and edge.

With `k=K2=26`, write `T_j=g(T2(j),J)` and `I=g(INF2,J)`. The actual newborn charge is

```
N(D,c) = f(J) + I + (d2/D) * (T_min(c+1,k+1) - I).
```

All 27 classes are possible at births already at `D=2000` (choose `c=j-1`, realized by attaching `c` disjoint degree-2 arms and `L` leaves). Class 1 is indispensable even though `c>=1` on every core-core edge **before** promotion on a nonsingleton core tree: the star can create its first edge. Class 27 is the saturated class for all `c>=26`, not an upper bound on physical `c`. Since `s*q*d2/D` selects a signed upper envelope, use the **unsigned minimum** `T_j` for lower and unsigned maximum for upper. `f(J)=1099511627776` both sides. Actual tick data:

| Witness | `I` | `min(T_j)` (class) | `max(T_j)` (class) | Chosen slope `d2*(extreme-I)` |
|---|---:|---:|---:|---:|
| Lower | −76398360618 | −75356776991 (16) | −75312767888 (1) | 187485052860 |
| Upper | −76358851515 | −75260533495 (1) | −73647237304 (27) | 488090557980 |

For existing motif tail endpoints the separate interpolation rule is also valid at unbounded degrees: at level `(d,K)` put `t=d/D`, `j=min(c,K+1)`, and interpolate `t*T(j)+(1-t)*INF`. Its mean reciprocal equals `1/D`. Degree and class projections compose by multiplication of `d/D` and nested `min`. At fixed `c` the gain difference is `t*(INF-T)/(D+1)`; promotion produces `t*(T(min(c+1,K+1))-T(j))` and is identically zero when already saturated. The gain slope multiplier `gamma=r/(D+1)` is at most 1. For an exact `j=c<=K`, with fixed adjacent inner donor rate `a=1/D_neighbor` (or `a=0` when no fixed donor is used), the lower slope `(d-c+a)/(d+1)` follows since `(D-c+a)/(D+1)` is increasing for `D>=d`. The birth slope lies in `[(d-c)/d,1]`. For overflow, only zero is a uniform lower bound: `c` can be `D`, and all donor degrees can grow without bound. At an exact-degree overflow below `d`, the external gain after one adjacent sponsor is removed is bounded by `L+(c-1)/2 <= D-(K+2)/2`. These are signed affine checks: choose the larger rate endpoint only for nonnegative signed increment. For an adjacent inner neighbor, the additional gain credit from this tail is `t/d*((1-t)*delta_inf+t*delta_T)`; bounding its signed bracket by the maximum of the two signed corner differences is valid for each coordinate. With two adjacent tails, both actual cross-donor reciprocals stay inside the respective `r` factors, and bilinear averaging gives the exact `t_u*t_v`, `t_u*(1-t_v)` etc. weights. The zero `g(INF,INF)` corner has no drift assigned to it. These statements transfer the slope domains; the audit does not recheck the full finite higher-motif row arrays.

## Rational enclosure and join

Write `w=A+Bq+Cq^2+Eqx` with `q=L/D`, `x=1/D`. The exact gain maps `(q,x)` to `((q+x)/(1+x),x/(1+x))`, and promotion maps it to `(q-x,x)`. The normalized gain rate `rho=r/D` ranges between `q` and `(1+q)/2`; the first is a limit for positive `c`, the second is realizable with all core neighbors of degree 2. One cannot select the same signed gain endpoint for both witnesses just because both unsigned tail weights increase in `q`.

With `G=(1+x)^2*(w_after_gain-w)/x` and `H=w_after_bud-w`, expanding the two differences gives

```
G = B(1-q)(1+x) + C(1-q)(2q+(1+q)x) + E((1-2q)x-q x^2)
H = -Bx + C(x^2-2qx) - Ex^2
P_rho = s*[rho*G + (q*(H+f(J)+I+slope*x) - lambda*w)*(1+x)^2].
```

For the actual arity class the numerator uses its own slope, so this is exact for the relaxed newborn envelope and an upper bound for the physical signed residual. The script obtains the polynomial **from direct increment evaluations**, not these expanded coefficients. The removable `x=0` gain quotient is `(B+2Cq)(1-q)`. After clearing the positive denominator, the bidegrees are at most `(2,3)` upper and `(3,4)` lower. On each of 16 closed equal `q` panels and the entire closed `x` interval `[0,1/2000]`, evaluate on a rational tensor grid, recover the Bernstein coefficients via exact inverse evaluation matrices, and check two off-grid identities per panel. The coefficient maximum bounds the polynomial by the nonnegative partition of unity. For maximum `M<0` divide by `(1+1/2000)^2`, since dividing a negative value by a larger positive denominator raises it. For `M>=0` retaining `M` is safe. Apply integer ceiling once at the end; zero is allowed. Affinity in `rho` then proves the bound for every intermediate rate. This covers the full continuous rectangle, not just integer test states, and in particular every unbounded real `D>=2000` under the rate, event and arity hypotheses. The `x=0` boundary is a compactification, not a physical vertex.

| Witness | Panels × gain endpoints | Exact worst panel bound | Tick ceiling | Worst seam residual and source arity | Worst rounded seam cap |
|---|---:|---:|---:|---|---:|
| Lower | 16 × 2 | `0` | `0` | `0`, `c=1999` | `0` |
| Upper | 16 × 2 | `-868735957253/133466700` | `-6509` | `-185516957025668413/49950012500000`, `c=1921` | `-3714` |

The lower tail zero is the relaxed no-leaf, low-receive corner; a finite no-leaf hub has strictly positive core receive, so this does not assert physical equality. At the join use the table for the pre-event value and for the bud value (degree does not change on bud); use the analytic value at `q=(2000-c)/2000`, `x=1/2000` only for the gain destination. No continuity or equality of the two definitions at the switch is needed. All 2,000 seam gain jumps per side are positive (smallest `111510597/2000` lower and `114022` upper), so the signed seam tests pick the **low** gain rate for lower and **high** for upper. The script computes exact terms, then reproduces upward rounding of each signed term at `1/1024` tick and upward rounding of their sum to ticks. The upper seam ceiling `-3714` fits the receipt's more permissive all-tabular cap `-3068`; the lower seam reaches its cap `0` in the relaxed rate box. At the cut, worst-case products under the validated parameter bounds are comfortably inside signed 128-bit (e.g. the gain numerator after sub-tick scaling is bounded by `1024*D_F*4*(16*S)*D_F^2 < 2^90`).

## Comparison bounds, negatives, and transfer

For both actual polynomials, `C,E<=0` and `B+2C+E/2000>0` verify monotonicity in `q` and nonincrease in `x`. The minimum for physical `c=0,1,2` is consequently at `D=2000`; the global rectangle extrema follow from its corners. For `D>=2000`, a vertex hosts no higher term under the chosen motif accounting, even when it appears as someone else's endpoint. The script checks `121743` tail/parent corner pairs per side to enclose the one allocated parent edge, interpolating between `t=0` and `t=180/2000` for each parent corner and including the root's zero edge. The edge bounds are `[−174995607038,0]` lower and `[−174987189676,0]` upper. The lower comparison must include the exact `alpha*S=S/4` arity shift: subtract it at `c=1`, add none at `c=2`, add at least one at `c>=3`. At `c=0` there is one core and no edge, so use the star value directly. The resulting tail minima are `99882299906` lower (indeed the receipt's **global** minimum, near `q=0`) and `980048237951` upper. Unshifted tail maxima are `1808571640916` and `1808494729937`.

Controls per side: 200 direct rational drift probes at `D=2000,2001,10007,10^30`, all domain edges and both sides of the class clamp, with relaxed rate extremes and realizable sponsor degrees `2,7,10^50`; 180 one-knot endpoint domain cases over all four levels. Planted errors rejected:

- **Use pre-bud `c` or the opposite signed class extreme.** Both disagree with the true newborn increments on actual classes; the opposite extrema underestimate the signed birth term at the physically realizable extremizing class.
- **Clamp the physical arity in gain, or give overflow the exact-class slope floors.** A near-leafless `D=2000,c=1999` hub detects the first; a `D=c=2000` hub disproves both false overflow floors, even with degree-2 neighbors.
- **Drop the lower receive endpoint.** For lower, all 16 high-endpoint bounds are strictly negative (worst ceiling `-3573`). A real no-leaf hub with total core receive `10^-20` (2,000 arms with neighbor degrees `2000*10^20`) exceeds that false cap.
- **Lose the `+1` in the seam leaf count.** The seam star's destination is analytically evaluated at `q=1`, not `1999/2000`.
- **Use the negative numerator unchanged, or floor the signed cap.** A synthetic constant cleared numerator equal to the actual worst upper numerator `-2606207871759/400000000` exceeds the false undivided bound after evaluation at `x=1/2000` and division; its residual is `-868735957253/133466700 > -6510`. This is a counterexample to the rounding/division rule, not a claimed violation of the platform's physical star. The actual upper star residual was `-5697190481931459/25000000`.
- **Drop the two tail cross-donor rates.** The pinned double-tail edge slices gave zero omitted contribution in all 729 class pairs per side in the probe. A synthetic symmetric edge with value one only at two finite tail corners (zero at infinity), on a realizable double-star of degrees `(2000,6001)`, yields a nonzero omitted term `-648324027/1441680680113340000`. Thus a control using only the pinned double-tail slices would not detect this error.

Recombination with the receipt-only finite caps reproduces `cap_f=0/-1053`, errors `19186/19783` ticks per core, and the receipt's two exact endpoints and outward decimals. The transfer is pointwise: sum signed row bounds, charge each core once, and use `C-1` edges with a nonnegative edge cap and the finite higher-host allocation. No degree truncation appears in the summed generator. On finite trees the shifted lower comparison sums using `sum(c-2)=-2`, so a positive lower potential comparison converts the per-core error to the signed exponent adjustment. Starting at `T3`, stop at size `M`; the Yule vertex count has a finite second moment on each finite time interval, and the absolute generator is bounded by a constant times `N(2N+1)` from linear comparison and unit population jumps. Dominated convergence removes the size stop; Grönwall then transfers the signed rates to the mean exponent. This argument uses receipt finite premises, not a new end-to-end verification.
