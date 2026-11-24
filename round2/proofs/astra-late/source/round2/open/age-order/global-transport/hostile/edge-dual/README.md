> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../README.md#claim-26) for exact scope and citations.

# Directed edge excess for the mixed hubhood interface

**At threshold zero, the gap for labels `n+1` and `n+2` is strictly positive on every finite old tree of size `n≥2`.** With at least two cores the conditional gap is bounded below by `(b+2a)/[n(n+1)]`, with `a=.00634` and `b=.017`. The seed and stars use the separate discharge bounds. Taking expectations proves strict adjacent order from label 3 onward; this placement comparison does not compare the seed labels with label 3.

## Edge ledger

Use the exact rounded comparisons of [discharge](../discharge/README.md) and the notation `C,s,U,L` there. Let `a=0.00634`, `b=0.017`. At an old core `u` of degree `d`, core degree `c`, and leaf stock `l=d−c`, bound its **directed** outgoing debit to a core neighbour `v` of total degree `e` by

```
J(u→v) = (s_u − 1/[d e(e+1)])_+ (U(d+1,l) − L(e+1))_+.
```

Indeed, the opposite `s_v` includes the actual old incoming reciprocal `1/d`, with nonnegative other contributions. This lower bound is valid at either endpoint even when its stock is zero. For the exact paired debit use

```
debit(u,v) ≤ max(J(u→v), J(v→u))
           ≤ b + (J(u→v)−b)_+ + (J(v→u)−b)_+.
```

Only one `b` is spent per edge. This decomposition is what lets a vertex budget only its own **excess**, even if its opposite neighbour has smaller degree and much higher outgoing exposure. The stock-independent ceiling `J(u→v)≤s_u U(d+1,l)` is sometimes enough.

Let `K` be the number of old cores, with `K≥2`. Define the exact residual of this ledger at `u` by

```
R(u) = C(u) + a(c−2) − b − Σ_(v core neighbour) (J(u→v)−b)_+ .
```

The old core tree has `K−1` edges and `Σ(c−2)=−2`; hence the whole pre-placement numerator satisfies

```
N(T) ≥ b + 2a + Σ_u R(u).
```

This holds for arbitrary finite old trees with at least two cores; stars and the edge seed already sign by discharge. Here the `−a` for a terminal in its `R` formula is paid using `C−max(b,sU)>a`, valid **regardless of the adjacent degree**. The previous terminal reasoning restricted to high neighbours is not needed: for `2≤d≤9` check receive `0` and `1/2`. The residual bound is concave in receive, since `C` is concave and the subtracted positive part is convex. At high terminal degree use the stocked bound in mixed-stock and `sU≤1/[2(d+1)]`.

## Certified profiles and exact obstruction to this ledger

For `d≥5` and `c≥2`, the uniform bound `C+a(c−2)−b−c(sU−b)_+` is positive. Through degree 112 it suffices to check total core receive `0,c/2`, with minimum `101/50000`. For `d≥113`, the analytic mixed-stock inequality gives `C+a(c−2)>1/22` without any restriction on the neighbouring degrees, and `sU≤1/[2(d+1)]<b`; thus the same minimum lower bound is valid. Stocked degree four (`c=2,3`) passes the two-endpoint test too with surplus greater than `101/50000`.

Below that, keep each low neighbour (`2≤e≤9`) exact in `J`. For a high neighbour use `sU`, since it upper bounds its directed debit independently of how large its degree is. For each local low multiset and number of high slots, the unknown contribution to the total receive ranges from zero to `#high/10`. Checking the two ends is valid by concavity, now of the full local residual in `s_u`. The finite local checks, not a sweep of trees, add these profiles:

| `d,c` | Certified provided | Minimum residual |
|---|---|---:|
| `3,2` | at least one low neighbour | `1173389118653/270000000000000` |
| `3,3` | all three neighbours low | `15612335509/4000000000000` |
| `4,4` | at least two neighbours low | `2503587899/750000000000` |

Thus every tree with no other internal profiles has conditional hub gap at least `(b+2a)/[n(n+1)]`. This includes trees with arbitrarily many genuine interfaces between small and large *internal* cores, e.g. one stocked degree-three junction with one low and one high core neighbour. The mixed subclass has no size or maximum-degree cutoff.

Before adding transfers, the local ledger can be negative at degree two with neighbours `(2,20)` and at both internal degree-three stocks facing only degree-111 neighbours. These are feasible local censuses, and the residual is a lower bound, so none gives a negative physical gap.

## High endpoint transfers and chain suppression

Each high endpoint can spend `.00356` on every **direct** degree-three or degree-four edge and `.00625` on every incident suppressed degree-two chain. These payments are simultaneous. For a high internal endpoint, let `k,j` be the number of direct paying ports and chain ports; they contribute at least `k/4+j/2` in reciprocal receive. The upper receive bound remains `c/2`. The uniform lower ledger is concave in receive and affine in the two payment counts after subtracting the transfers. Thus only the vertices of `k,j≥0`, `k+j≤c`, at both receive limits, need checking. For degrees 10–112 the exact replay verifies all five distinct corners. The same stock tail as in mixed-stock drops all outgoing excess at degree 113 and above. At a corner paying nothing or all direct ports, use its two credit bounds with `a` replaced by `a` or `a−.00356` in the `c` term only. At the all-chain corner the receive is `c/2`, and `w≥113/228`. If `s≤s_l`, then `s_l−s≤1/[2d(d+1)]`; if `s>s_l`, the leaf toll vanishes for `l≥4` and is bounded by `l(l−1)/[2d(d+1)]` otherwise. In both cases the leaf serial credit covers it, leaving `C≥g·113/228`. The terminal bound has enough surplus for either payment.

The direct payment is credited at the actual low recipient. The local multiset checks including these receipts give positive minima for **all** internal degree-three and degree-four profiles, including when the other low ports go to an internal degree-two core. The minimum for bare degree three is `1/50000`.

Suppress each admissible maximal chain to a single skeleton edge. Endpoints retain their **actual** old profiles, so a high endpoint sees reciprocal `1/2` on each chain port; the skeleton core degrees of retained endpoints are unchanged. On a bare degree-two chain, credit is `2sG(3)` per site, and consecutive internal chain sites of weights `s,t` have paired debit `(U(3,0)−L(3))_+|s−t|`. At a high endpoint of degree `e`, bound the chain-side debit by `s(U(3,0)−L(e+1))_+` and the other direction by `.012` using the actual-stock exposure check and the tail. Keep both `1/e` and `L(e+1)` exact for `10≤e≤30`; for `e≥31`, use `L=0` and the two reciprocal bounds `0,1/31`. Convexity of debit minus credit in the endpoint reciprocals reduces each case to corners. Lengths one, two and three suffice: further sites have identical `s=1/6`, contribute positive credit and introduce no extra pair debit. The worst net debit is `70627717243/2400000000000 < .0295 = b+2·.00625`, so the single skeleton baseline and its two payments cover the chain.

A chain joining a high endpoint to a low terminal is handled by checking the terminal and chain together rather than by charging the low endpoint. At the terminal the profile has its unique core neighbour of degree two, so the last pair can use the actual paired debit with that terminal profile. The joint credit exceeds `a+b`. To retain the residual formula, assign `a+b` of this credit to the terminal and zero debit to the skeleton edge: reserving the baseline on that edge is safe, and the terminal's full `R=C−a−b` is zero. No high payment is needed for this case. When both ends are low terminals, the whole core is a path; checking the two terminal pairs directly gives a joint credit greater than `b+2a`. There is no direct payment along any suppressed edge. Summing on the skeleton thus keeps the numerator floor `b+2a`. The edge seed and stars use the earlier checks.

## Closing the remaining low endpoints

The [low-end replay](../interfaces/low-end.py) supplies the remaining chain transfers and allows a single ledger for **all** endpoints, including terminals. A retained endpoint of degree 3–9 pays `.0037` per incident chain; a terminal of degree two pays `.021`; high endpoints still pay `.00625`. At a low endpoint each port of old degree two may be either direct or suppressed. Taking the maximum of its ordinary directed excess and the chain payment covers either choice. With the other ports unchanged and the `.00356` direct high receipt retained for degree-three and degree-four profiles that need it, the low residual is concave in the total receive. The verifier checks every local low multiset at both bounds on the high-slot receive. The lowest low-endpoint residual is positive even for degree-two terminals; the minimum for degree three is `1/50000`.

For a chain with a low endpoint, the actual endpoint stock determines the two paired slopes, and the known adjacent bare site contributes `1/2` to endpoint receive. Fix the endpoint degree and stock and range the other reciprocals between zero and `1/2`. The paired debit is convex in the endpoint coefficient; thus the two resulting endpoint weights suffice. For the high side keep the degree and receive ranges used above. For lengths one through three the worst low–high debit **after both payments** is about `.016904845`, and the worst low–low debit is about `.010693899`, both below `b`. Each later interior site adds positive credit and no new paired debit. The high–high calculation already supplies the same bound.

For a core tree of at least two vertices, suppression leaves a tree on the retained cores: the original core tree has at least two terminals, and neither is suppressed. This includes the two-vertex skeleton when the entire old core is a path. In the uniform ledger use the actual old endpoint profiles. Direct edges cost at most one baseline plus their two directed excesses. Suppressed edges cost at most one baseline plus the two chain payments, with the suppressed sites' own credits already netted out. The direct high-to-low transfers cancel in the sum. All retained residuals are nonnegative under these simultaneous charges and the two tree identities give `N≥b+2a`, with no omitted pendant, stock, or root terms.

## Earlier annealed sufficient inequality

Let `T_n` be the actual old tree at size `n`, and in the event `K≥2` let `M` count certified internal profiles above and `T` the number of terminal cores. Let `B2`, `B32`, `B33`, `B44` count just the four remaining census types **that fail the certified-profile gate**. The same local checks including every high-slot count bound their individual negative residuals by `.017`, `.00414`, `.01066`, `.00432`, respectively. Each certified internal profile contributes at least `.00202` and a terminal at least `16436047/3000000000000` (including unbounded terminal degree). The following is consequently sufficient for an unconditional strict order at the adjacent labels at that size, with star/seed contributions available separately:

```
E[1{K≥2}( .017 + 2*.00634 + .00202 M + (16436047/3000000000000) T
           − .017 B2 − .00414 B32 − .01066 B33 − .00432 B44 )] > 0.
```

These are means at the **same actual entrance size**. Exchanging the average of the signed ledger for a stationary limit does not settle finitely many labels. The bound reduces the annealed problem to four adverse censuses and two positive counts but does not certify their finite-size ratios. Counts of ordinary degree-two vertices without their core stock would not supply the missing `B2` estimate. The local charge identity itself needs neither independence nor a random-root stationarity assumption.

## Replay

From Astra root:

```
python3 -IB round2/open/age-order/global-transport/hostile/edge-dual/check.py
```

The run checks the inherited mixed-stock scalar bounds, 12,424 uniform degree/stock corners, the terminal vouchers, the local multiset inequalities and annealed adverse budgets, the high port transfers, and the chain and pendant bounds. It also runs the low-end verifier to cover every remaining endpoint type. The printed negative controls refer to the older ledger before transfers and suppression.
