> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited. See the [release index](../../../../../../README.md#claim-27) for exact scope and citations.

# Signed whole-tree first and second moments

**Both conditional gaps are strictly positive on every finite old tree**, given the uniform normalized serial-middle bounds `M_1 ≥ .56319` and `M_2 ≥ .50833` for an actual degree-two middle with a fresh leaf and arbitrary old parent. The analysis needs no host-specific increase of these middle bounds. The first moment suppresses only bare degree-two chains; the second has a simpler vertex-by-vertex payment.

Let `n≥2` be the old size, `m=n+2` the common completion size, and `p=1,2`. For a fixed completion, the physical rate-one tree has `W=S Z`, with the actual size amplitude `S~Gamma(m,1)` independent of the embedded tree path from that completion. Thus `E W^p/p! = (m)_p E Z^p/p!`. The [paired placement identity](../../adjacent/README.md) is linear in the continuation test, so the normalized moment numerator `N_p` satisfies

```
E[Z_(n+1)^p − Z_(n+2)^p | F_n] = p! N_p/[n(n+1)(m)_p].
```

This only uses size-clock independence from a fixed completion, not after promotion or after a clock race. The exact result below is `N_1≥1/40` and `N_2≥1/50` for all trees. More precisely `N_2≥Σw/50` for trees with at least two cores, and `N_2>.03` for the seed and stars. Indeed, `Σw≥Σr_u/n=1` since each `D_u+1≤n` and every initiator distributes total rate one among its neighbours.

## Endpoint moment inputs and old-leaf charge

Use the [exact shifted first-send upper](../tip-moments/README.md) `U(k)` at actual post-placement parent degree `k`, with cutoff 64 and the analytic tails `U(k)≤1/k` for `p=1`, `U(k)≤(k+1)/(2k²)` for `p=2`. Put `g=M_p−U(2)`: the serial adverse tip has parent degree two even when the older host degree is large.

For a favourable parallel leaf with parent degree `k`, use `L(32)=0`, `L(k)=0` for `k≥32` and the lower recurrence

```
L(k) = [M_p/k + k L(k+1)]/[p+k+1/k],   k<32.
```

Before promotion the marked leaf remains degree one, so the parent has total gain rate at most `k` and its send to the leaf is exactly `1/k`. The sequence is decreasing (checked), so the maximal competing parent gain rate is pessimistic in the discounted generator. At promotion restart with the same uniform middle bound, discounting by `exp(−pτ)`. At cutoff, discarding later promotions is safe. The proof uses the actual physical continuation and does not assert independence of the clock and its payoff. The replay uses the input `M_p` literally, with no fitted exterior data.

An old core has total degree `d`, `c` core neighbours, `l=d−c` old leaves, and core reciprocal receive `r`. Set `w=(l+r)/(d+1)`, `s=w/d` and `t=1/(2d)` (the coefficient of each of its old leaves). Pay its serial credit and all old-leaf edges in one profile:

```
C = g(w+l t) − l[(t−s)_+(U(2)−L(d+1))_+
                  +(s−t)_+(U(d+1)−L(2))_+].
```

Both directions use the same `s−t`, and each old leaf is paid only at its own parent. Core edges use the analogous same-difference maximum with their actual shifted parent degrees. In particular the outgoing debit of a core edge is at most `s U(d+1)` when that core is adverse.

## First moment: chain suppression and retained budget

Suppress every maximal path of internal cores with `(d,c)=(2,2)` to a skeleton edge. Retain actual old endpoint profiles, including reciprocal `1/2` per port into such a chain. A skeleton edge pays a baseline `b=.008`; a retained core of degree at least three pays `f=.004` per incident chain, and a terminal of degree two pays `f=.010`. There is no other transfer. With `a=.0085`, a vertex with `k` chain ports retains at least

```
C+a(c−2)−b−k f−c(s U(d+1)−b)_+ .
```

This counts its whole outgoing excess even if some ports are chains, which is deliberately conservative. For a direct edge the debit is at most `b` plus the two excesses: it is a maximum with one signed difference, so `b` is charged only once. With `k/2≤r≤c/2`, the displayed residual is concave in `(r,k)` on the triangle with corners `(0,0)`, `(c/2,0)`, `(c/2,c)` in `(r,k)` coordinates. Exact replay through degree 63 gives minimum `1/2000` even with joint payments. This uses no small-neighbour census or bound on maximum degree.

For `d≥64`, the direct excess vanishes (`sU≤1/(d+1)^2<b`), and the core-side old-leaf debit vanishes since `U(d+1)<L(2)`. Put `z=(64/65)g+(g−U(2))/2>0`. The remaining stock credit is at least `(l/d)z`. If `l/d≥1/2`, this alone pays `2a+b`, and otherwise `c> d/2` so `(a−f)c` alone pays `2a+b`. Thus the inequality holds through unbounded degrees, including terminal cores.

For a suppressed chain, length `k=1` has one coefficient `(r_L+r_R)/6`, with `r_L,r_R` the two *endpoint* reciprocals. For `k≥2`, the boundary coefficients are `(1/2+r_L)/6,(1/2+r_R)/6` and all interior coefficients are `1/6`. Bare chain vertices receive `2sg` each and an internal chain edge costs `(U(3)−L(3))_+|s−t|`. For an endpoint of degree `2..9`, keep the actual stock and evaluate the two corners of the remaining receive interval (adjacent chain port fixed). For degrees `10..30`, retain the exact reciprocal and `L(d+1)`, bounding the endpoint-side debit by `U(d+1)/(d+1)`; for degree `≥31`, use `L=0`, endpoint-side debit `≤1/32²`, and reciprocal endpoints `0,1/31`. Boundary debit minus credit is convex in the ranged coefficients. Testing lengths one, two and three gives strict net slack after both bills; each additional site adds `g/3>0` and no new internal debit. The smallest replayed chain slack is greater than `.000291477`.

The skeleton still has at least two cores (old terminals are never suppressed). Thus its edge count is `K−1` and its core-degree excess sums to `−2`. The residual bound and chain bills give `N_1≥b+2a=1/40`, including the two-terminal skeleton obtained from a whole path.

## Second moment: pure sitewise discharge

Reserve `1/50` of each serial weight, then charge every core its own `c s U(d+1)` for core edges, without testing which direction is adverse. For `d≤63` the residual after this payment is nonnegative at both receive endpoints `r=0,c/2`, with strictness unless the impossible actual profile `l=r=0` is used. Concavity gives the entire interval. For `d≥64`, write `h=g−1/50` and `u=66/(2·65²)`. All core-side debits (both old-leaf and core edges) sum to at most `d s u=w u`. The other old-leaf direction sums to at most `l t U(2)`. The retained residual after the reserve is therefore at least `w(h−u)+l t(h−U(2))≥0`, by `h>u`, `w≥(64/65)(l/d)` and the replayed strict stock bracket. No suppression or terminal identity is used for this moment.

For both moments, the seed has no edge debit and total serial weight one. Stars of centre degree `2..63` are checked directly with the same signed profile (`c=0`); in larger stars the centre-side debit is at most `w U(d+1)` and the centre serial credit alone suffices, even discarding the leaves' credit. This proves the stated seed and star floors.

## Replay

`python3 -IB round2/open/age-order/hub-ranking/mean-ledger/check.py` checks both local certificates in exact rational arithmetic, using 42b5's frozen shifted upper checker. No whole-tree enumeration or asymptotic law is needed.
