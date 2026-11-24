> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Gold quintet: signed physical raw drift

At the exact trial `s=141431/250000=.565724`, a radius-five marked-core transfer for the full raw score has a two-sided error **less than `3.569e−7` per core** at `K=106, J=1999`, provided the high-degree rows are retained as reciprocal moments. This is a bound on the integrated error of the **raw generator**: the proof does not set `s=μ` or use stationary cancellation to discard the error. The joint moments in the retained transfer have not been priced. There is no new exponent endpoint.

## Physical ledger and integration

Use `W=Φ_lower+Φ_upper`, with each on its own `2^40` scale. Everything below also bounds the error for `Φ_lower−Φ_upper` and for either endpoint separately, since the bounds use sums of the individual absolute norms. This convention is independent of the original certificate's endpoint signs. Write `w` for a normalized occurrence weight. For every distinct copy of each of `f,g,h3,p4,s3,p5,b5,s4`, the pre-event ledger is:

| Term | Rate times reward | Connected pre-event roles |
|---|---|---|
| Discount | `−s w` | occurrence |
| Leaf gain at `v` | `(D_v−c_v) Δ_D^v w` | occurrence |
| Core gain `u→v` | `(1/D_u) Δ_D^v w` | occurrence, plus the actual donor if external |
| Old promotion at `v` | `(1−c_v/D_v) Δ_c^v w` | occurrence |
| New motif at promotion of a leaf of `v` | `(1−c_v/D_v) w(post)` | actual connected predecessors, with promoter `v`; just `v` for `f` |

A gain increments the receiver's degree; a promotion increments the promoter's core degree at fixed total degree and creates a core of mark `(2,1)`. A new core can only be a tip (or a singleton). The post-promotion lookup includes the promoter's updated mark. Evaluate the true successor even when it crosses a cap. No old-old support can be created by either change: the topological minimum core arity is automatic, and interior cutoffs and degree-sum masks only become harder to meet under the old-role changes. Support can drop. A zero intrinsic table word on supported input does not suppress its successor. Equal-arm copies are distinct incidences.

Split each event summand evenly among its pre-event roles, including the external donor; call the resulting root reward `R_s`. It uses at most six connected cores, so the actual unoriented radius-five core ball with retained `(D,c)` boundary marks determines it. In particular a boundary mark is not the degree in the displayed ball. This split gives exactly `Σ R_s=(L−s)W` on an ordinary finite tree. In the uniform-core mean law `π`, the integrated signed transfer is

```
π R_s = μ E_β ∫₀^∞ exp(−μ a) R_s(X_a) da
      = π R̃_(s,K,J) + ε,                     |ε| < 3.569e−7 .
```

Here `β` is the actual physical newborn environment, `X_a` the actual aged neighborhood, and `R̃` the retained row transfer below. The first equality is the [physical stable-age identity](../../../../../astra/source/round2/degree-law.md#birth-times) applied to both signs. Reversal of split incidences is by mass transport on the [physical rooted law](../../../../../astra/source/round2/fringe/spatial-limit.md); it moves a weight along its actual edges, including either direction of a core gain. Absolute summability follows from the per-incidence constants and moments below. The separate compensated balance `π R_s=(μ−s)π W` follows by integrating the absolutely summable ledger over a finite time interval in the stable mean law (same argument as the [D18 integrated transfer](../../occupation-neon/embedding-flow/full-tail/README.md#4-stationary-integration-and-the-completed-hierarchy)). With positive `π W` it makes the sign of a certified raw average informative. It does not give the sign at this trial from the known bracket, and is not a premise of the error calculation.

## Retained high-degree rows

Discard a term if any of its pre-event roles has `c>K`. Keep the actual topology, low marks, and exact `c` on all roles. For `D>J=1999`, put `x=1/D` in the boundary slot. The exact old and promoted singleton is `a₀+a₁x+a₂x²`, obtained by expanding `q=1−cx` in the payload polynomial. Each high outer motif slot supplies the exact affine interpolation `w∞ + d(w_knot−w∞)x`, separately in each argument; admitted interior roles cannot be high. A high promoter still pays `1−cx`, and a high donor still pays `x`. Newborn and discount rows are exact polynomials in the high coordinates. So are gains at low receivers: evaluate both low old/new entries exactly and keep their full high-coordinate interpolation. This includes the low-to-analytic seam.

The only expansion is for a high receiver. Its precise changes in `x,x²` are

```
Δx  = −x²/(1+x),       Δx² = −x³(2+x)/(1+x)².
```

Use `1−x+x²` and `1−2x+3x²` respectively for the two inverse denominators. Their pointwise errors are at most `x³` and `4x³`, by the geometric remainder and Taylor's third derivative bound. Multiplication by the physical leaf gain `x⁻¹−c` leaves a polynomial because each derivative increment has a factor of `x²`. The same approximated increment is used for all core gain arrows. Thus for fixed `K,J` there is a finite list of role patterns, low marks, high arities and monomial powers. For each such entry the signed transfer is its row coefficient times the **joint** occurrence moment `π Σ ∏_(high i) D_i^(−m_i)` with the specified degree tests, topology, roles, and share. The moment of the constant monomial is required too. The ledger is a definition of all these coefficients without enumerating a giant table. The other high marks can be correlated with the receiver and donor, including in these moments. If only fixed finite census columns are available, each moment can first be approximated with an additional fixed degree cap `J'`; its missing count is at most `p K^(p−1) P_(J')` per labelled pattern before coefficients. Uniform intermediate slice bounds are needed if the support of a relaxation itself grows.

## Rational error ledger

All costs are degree/arity **pre-event** losses. For a labelled type with `p` pre-event roles, pay one of its maximal-`c` roles on arity overflow. Re-rooting at that vertex gives a count at most `c^(p−1)` for the type. Summing the possible choices gives

```
A_p = p E[(1+N)^(p−1); 1+N>K],
N = Pois(A), A ∼ Exp(s₀), s₀ = .565702 ≤ μ.
```

At the actual common observation age a new core starts from `(2,1)`. A rate-one promotion proposal clock dominates later core degree minus one. At the same age, the gain rate is bounded by `D−c/2≤D−1/2`, so the proposal chain from degree two with rate `d−1/2` bounds its degree. Open neither physical mark into independent neighbor marginals. Comparison uses the older exponential age at `s₀`; the arity budget and the degree-with-surviving-arity budget are disjoint. In `A_p`, use `P(N≥K)=(1+s₀)^−K`, the geometric memoryless property, and factorial moments `j!/s₀^j`. For degree use the rational product `P_J=∏_(d=2)^J (d−1/2)/(d−1/2+s₀)`, rounded up in the implementation. For `J>2000` the product beyond 2000 is at most `(2002/(J+2))^s₀`, with the exact dyadic relaxation in the replay.

For motif families let `M=32` be the sum of the two verified table word caps, `d` the family's endpoint knot, `e` its number of tips, and `r` its size. The convex interpolation bounds old/new weights by `M`, and a change by `2M`. Leaf gain has norm `G=2M max(38,d)`: an interior jet or crossing at an exact endpoint has `D≤max(38,d)`, while in the interpolated region use `|Δ_D w|≤2Md/[D(D+1)]`. Per family the arity budget is exactly the following outward bound:

```
(s M + r G + (r−1) 2M + r 2M) A_r
    + r M A_(r+1) + e M A_(r−1).
```

The last two summands are external core gain and newborn creation. Each oriented internal or external send is at most `1/2`; counting `2(r−1)` internal orientations yields `(r−1)2M`. The singleton uses the analogous one- and two-role counts, the header coefficient bound `F` for `|f|`, and `2 F·2000` for leaf gain. This last bound covers the analytic region too: the analytic jump times `D` is bounded by `|a₁|/(D+1)+3|a₂|/D²`, with `c≤D`, and the header-based definition of `F` in the checker bounds that expression by `2 F·2000`. Newborn singleton norm is evaluated directly.

At a high receiver the total true gain intensity is at most `D`. Thus a motif's integrated derivative remainder per such occurrence is at most `2Md/(J+1)^4` times its incidence mass. Re-root **at this receiver**, which must be one of the `e` tips. The first edge has at most `K` choices and the remaining `r−2` edges start at motif interiors with core degrees bounded by `m`. The physical count is at most `e K m^(r−2) P_J` (`m=1,15,8,8,6,6,6` in family order after `f`). For `f`, at surviving arity the summed bounds `b₁,b₂` on `|a₁|,|a₂|` come from the two headers. Its error is at most `P_J[b₁/(J+1)^4+12 b₂/(J+1)^5]`. This bounds the integral of the absolute errors, using the sum of the absolute errors of all terms if desired.

| At `K=106, J=1999` | Outward error per core |
|---|---:|
| Arity overflow | `1.633376e−7` |
| High receiver denominator remainder | `1.935214e−7` |
| Combined | `3.569e−7` |

For comparison, literal degree omission and replacing high marks by their limiting weight while dropping their reciprocal sends use much larger bounds in this calculation. Re-rooting the missing degree role in the arity norm replaces `A_p` by `p K^(p−1) P_J` for literal deletion. In the limit replacement the coefficient differences use tensor telescoping (`2Md/J` per high argument), promoter rate error at most `K/J`, and donor reciprocal error at most `1/J`; the checker prints that budget as `lump`. Those upper bounds are not impossibility results. There are genuine signed **local** deletion obstructions at the analytic seam: a homogeneous high-degree path, read at its interior, has omitted root share about `−1.120046269e−4` at `D=2000`, while the branched-tail fixture rooted at its high vertex has omitted root share `+1.808135050e−1`. Six extra path edges in the latter fixture keep the far boundary away from the root read. The marked balls are physically attainable by finite sequences of gains and promotions, with positive holding time; their scores persist for the displayed full degree/arity marks. These are not integrated lower bounds without quantitative cylinder masses, and they do not fix the sign of the total omitted expectation.

## Replay

```sh
taskset -c 5 nice -n 19 python3 -IBu round2/occupation-quintet/true-drift/check.py
```

The sparse reader uses the two original gold payloads from the audited snapshot and checks the detached source pin. No payload scan or rehash is required. The oracle sums direct directed ordinary-tree births with independent connected-subset potential enumeration; the ledger uses the route/center enumeration. The replay checks signed tail witnesses, arity/degree seams, repeated arms, joined hubs, and the marked-ball ownership identity. All calculations before displaying the results use exact rational arithmetic with outward rounding where specified.
