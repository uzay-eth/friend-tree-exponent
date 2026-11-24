> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../README.md#claim-26) for exact scope and citations.

# Pre-placement sign: a scalar bound and a grouped transfer criterion

**Result.** The expected two-birth hub gap is strictly positive on every old tree whose cores all have at least 1,024 leaves. There is also an explicit sufficient degree/leaf-profile test allowing some low-stock cores. The all-tree sign remains open. A same-parity leaf slide can increase the tip's hub probability by order one; the separate tip-monotonicity approach fails.

## 1. Coupling before placement and the leaf slide

In an old tree write `w_u=r_u/(D_u+1)` and `s_u=w_u/D_u`, where `r_u=Σ_(v~u)1/D_v`. After the first birth at `u`, a reference experiment chooses an old second target with its *old* distribution `r_v/n`. Relative to the unnormalized fresh target rates, every old adjacent `v` has lost `δ_u=1/[D_u(D_u+1)]`, and the new middle has gained `D_uδ_u`. The only other increment is the weight one at `u`, from the new leaf's own initiations. The repeated-target endpoint gap is zero. In the reference experiment all endpoint-gap terms cancel by exchanging the birth order. Thus the numerator of the actual conditional gap is

```
R(T) = Σ_u Σ_(v~u) s_u [ I(u) − H(u,v) ]
     = Σ_u w_u I(u) + Σ_{u,v} (s_v−s_u) H(u,v),
```

where the last sum is over unordered edges, with an arbitrary orientation used consistently in each factor. Here `I(u)=h_serial(middle)−h_serial(tip)`, and `H(u,v)=h_parallel(first)−h_parallel(second)`. The first expression is a literal transfer of second-target mass, from the old neighbour `v` to the new vertex `middle`; they have the same parity. The actual gap is `R(T)/[n(n+1)]`. A proof by improving both endpoints separately would require `h_serial(tip)≤h_parallel(second)`, which is false:

* Take `v` the centre of a star with `k` leaves and `u` one of its leaves. The parallel second tip's parent has at least `k−1` other leaves, so its promotion and hence hub probability is at most `2/k` by the [first-send barrier](../early/README.md).
* In the serial completion the tip's parent starts at degree two and has a core neighbour `u`. Its first selection of the tip has rate `1/2`, the tip's return rate is one, and outside input to the middle is at most `1/2`. Promotion before a degree gain of the middle has probability at least `1/4`; following promotion a tagged birth–death family survives with probability at least `1/2`. Thus `h_serial(tip)≥1/8`, disproving tip monotonicity for `k>16`. The same slide also defeats a parity-based rate order: the tip's initial incoming intensity increases from `1/(k+1)` to `1/2`. This says nothing against the comparison of the two *differences*.

## 2. Uniform serial gap, including arbitrary physical exterior input

Following serial placement, consider the tip's first promotion as a specified arrow from the middle. Before it occurs, at parent degree `d` let `j` be the number of other leaf neighbours. The specified send has rate `1/d`, other-leaf removals have rate `j/d`, additions of other leaves have rate at least `j+1` (including one independent unit of return immigration from the designated tip). The degree increases by one with each addition and does not decrease upon removal. This counts current stock, not the number of leaves ever acquired.

The middle retains its old core neighbour, including when the original target was a leaf (that target became a core on the first birth). Thus `d≥j+2` always. The same is true of the parent of the *parallel second tip* in a core–old-leaf edge. On this domain the first-send payoff `f(j)=1/(j+1)` is a supersolution, with send payoff one. For `j≥1` its drift is at most `−1/(j+2)+1/d≤0`, using the rates above; at `j=0` it is negative. This supplies an infinite-time upper certificate at any cutoff.

For cutoff `m=112`, define `P(d,j)` by terminal payoff `1/(j+1)` at `d=m` and the exact scalar recurrence

```
P(d,j) = [1 + d(j+1) P(d+1,j+1) + j P(d,j−1)] / [1+d(j+1)+j].
```

Terms with zero coefficient are omitted. Exact arithmetic verifies `P(d+1,j+1)≤P(d,j)` and `P(d,j)≤P(d,j−1)` where applicable on the whole finite domain `2≤d<m, 0≤j≤d−2`. In particular additional exterior input contributes nonpositive drift, so the scalar function bounds first send for **any** physically admissible exterior. A pending race exits almost surely (if it has not sent, the tip's unit clock alone forces degree to the cutoff). Optional stopping followed by the cutoff barrier gives `h(tip)≤P(2,0)<0.611`. The tracked-degree recurrence in [early](../early/README.md) with cutoff 16 gives `h(middle)≥Q(2,1)>0.613`. The exact difference is larger than `t=0.0029`, uniformly in the old location. Use `U=0.611` as the parallel-second upper when the opposite old target was a leaf. This proof ignores any beneficial exterior immigration into the serial middle and gives the tip full hub credit immediately on promotion, so the gap bound is conservative.

## 3. Grouping both endpoints of a core–leaf edge

For a core `o` of initial degree `D`, stock `L`, let `b_o=min(1,2/(L+1))`, a parallel leaf upper when the placement directly at `o` leaves `L` other old leaves available. If the opposite old target was an old leaf of `o`, it is promoted in the parallel completion and only `L−1` of the old leaves are available, so use `b'_o=min(1,2/L)` when `L≥1`. In a core–core edge the adverse edge term is at worst `−s_u b_u−s_v b_v`. On a core–leaf edge, write `H(o,ell)` with the core first. If `s_ell≥s_o`, the lower edge bound is `−(s_ell−s_o)h_parallel(second)`; the second tag has a degree-two parent with **no** other leaf. If `s_o>s_ell`, the lower edge bound is `−(s_o−s_ell)b'_o`. In particular the separate leaf charge should not be bounded as if its opposite core charge were absent.

Combining both endpoints of every old core–leaf edge and distributing the adverse parts of core–core edges gives the sufficient test

```
Σ_(o core) [ t (w_o+L/(2D)) − (D−L)s_o b_o
             − L (1/(2D)−s_o)_+ U − L (s_o−1/(2D))_+ b'_o ] ≥ 0.
```

Here `b'_o` is immaterial if `L=0`. Strict positivity of the real gap follows even at equality since each serial bound is strict. The test uses `w_ell=s_ell=1/(2D)` for an old leaf, and works also for mixed stocks. For example, if every core has at least `M=1024` leaves, combine the two losses involving `b_o,b'_o` and bound them by `w_o(2/M)`. Also `r_o≥L`; hence, putting `c_o=D−L`,

```
L (1/(2D)−s_o)_+ ≤ [L/(2D)] [(D+1−2L)/(D+1)]_+
                   ≤ c_o/[2(M+1)].
```

Using `Σ_o c_o=2(C−1)≤2 C` on the core tree and `Σ_o w_o≥Σ_o L/(L+c_o+1)≥C(1−3/(M+1))` gives the lower certificate

```
R(T)/C ≥ (t−2/M)(1−3/(M+1)) − U/(M+1) > 0.
```

Bare and short-stock cores remain the bottleneck: the parallel tip can have a degree-two parent even when the associated serial middle starts from the same degree and stock. The scalar bound signs the serial gap, but throws away the correlation with the parallel endpoints. A general proof still needs either a stronger signed two-environment estimate or cancellation exploiting the low-stock core geometry. The obstruction in §1 only rules out a one-sided comparison of tip marginals.

## Scalar certificate replay

From this directory, `taskset -c 13 nice -n 10 timeout 12s python3 -uB check.py`:

```
scalar drift signs checked at 6,105 (degree, stock) rows
middle hub lower: 0.613194634; serial-tip upper: 0.610259660
uniform serial gap lower: 0.002934975
1,024-leaf criterion: numerator per core > 0.000348006
```

The replay is exact rational arithmetic, not a tree enumeration. Printed decimals are rounded to nearest for display; the one-sided inequalities use the rational thresholds `t=29/10000`, `U=611/1000`, and the exact rational stock criterion.
