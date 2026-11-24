# Cold check: eventual adjacent hub-probability ordering

**Verdict: holds**, using the deterministic late-label asymptotic and the uniform-core local and mean limits. I found no sign, conditioning or continuity obstruction. The edge contribution need not have a fixed sign, even on exact finite-horizon tests. The sign of the limiting adjacent difference comes from telescoping the *total* score.

## Exact antisymmetry and bound

For an old tree of size `n`, let `r_a` be its unnormalized rate of attachment at `a`, and write `Δ(a,b)` for the success probability of label `n+1` minus that of label `n+2` in the tree after those two actual attachments. When `a,b` are distinct old targets, exchanging the attachments produces the same final graph with the two new tags exchanged. Thus `Δ(b,a)=−Δ(a,b)`, including when one or both old targets were leaves. Repeated targets have interchangeable tags and give zero. After growing at `a`, the second rate at a distinct old `b` is `r_b−1{a∼b}δ_a`. This is because only the sending degree of `a` changes among old vertices, and the newly born leaf sends only to `a`. Subtracting the reverse order gives `(r_b δ_b−r_a δ_a) Δ(a,b)` on each unordered old edge, and zero on nonedges. When the second target is `n+1`, its only sender is `a`, at rate `1/(D_a+1)`. This leaves the serial term `r_a Δ(a,n+1)/(D_a+1)` with a middle core and a tip.

Equivalently charge each old vertex `a` once with

```
R(a) = r_a [ Δ(a,n+1)/(D_a+1) − δ_a Σ_{b∼a} Δ(a,b) ].
```

Then `G(o)=R(o)+Σ_{a leaf of o} R(a)` and the exact conditional, scaled difference is `Σ_{o core} G(o)`. Since `|Δ|≤1`, `|R(a)|≤2r_a/(D_a+1)`. A core contributes less than two since `r_a≤D_a`; each old leaf attached to `o` contributes at most `1/D_o`, so their sum is at most `L_o/D_o≤1`. The uniform bound three is per **prebirth** core, even when either attachment promotes an old leaf. Grouping by the final core tree instead would need a different root law.

## The probability of future hubhood, including an already promoted middle

For each of the finitely many placements at a fixed rooted environment, `Δ` is a difference of conditional probabilities of *unbounded future degree*. In finite trees this is equivalent to positive limiting fraction by the existing tagged-family argument; in infinite local states no global fraction is needed. The finite changes at placement preserve local convergence as tagged marked-core states. Exchangeable choices of original leaves may be represented by any fixed leaf slots in the coupled neighbourhood; only finitely many slots occur at a given root. The old-cell cutoff gives continuity for bounded finite-time history probabilities after each such modification.

Use a lower test `(1−1/d) P(hit d by time h)` with **`d>2`**, for either final tag. Both tags start at degree one in the two-old-target sector; in the serial sector the middle starts at degree two. Taking `d≥3` forces a fresh upward hit for both, so a newly acquired leaf supports the usual birth–death survival probability, without treating the initial core as a future hit. If one wants `d=2` for the middle itself, its current tip is already the requisite single leaf and supplies the lower success probability `1/2`. Neither test treats an initially promoted core as already a hub.

For the upper test use the current-neighbour first-future-send bound after time `h`, capped at one. Leaf neighbours may be charged one; for a core neighbour use the stock bound excluding the tag. The next degree gain uses one of the existing incident edges, even for an unpromoted tag. This bounds the conditional hub probability. On a nonhub path the tag eventually has constant finite degree and each of its fixed neighbours has unbounded degree, since an edge with both degrees bounded would eventually transmit another accepted arrow. Their leaf stocks diverge by the first-hit tagged-family argument. The upper charge tends to zero on that event, so bounded convergence pays the **nonhub atom**. The lower and upper continuous history probabilities then squeeze the ultimate success probability at each input. No infinite-time local observation of the actual indicator, no uniform-in-environment history length, and no conditioning on eventual success are needed. Multiplying by locally fixed rates and summing the finitely many choices at the root preserves continuity; the independent bound three supplies global boundedness.

## Limit, uniform integrability and strict sign

Let `d_n=p_(n+1)−p_(n+2)`. The identity gives `n^(2−μ)d_n = (n/(n+1)) E[(C_n/n^μ) (Σ_o G/C_n)]`. The empirical average converges almost surely to `πG`, and the `L²` bound on `C_n/n^μ` supplies uniform integrability for the signed bounded average. Therefore the limit is a finite `A=E[X]πG`. Because `μ<1`, the differences are absolutely summable from this limit alone. Telescoping, using `p_v→0`, gives `p_(n+1)=Σ_(j≥n)d_j`. Inserting `d_j=(A+o(1))j^(μ−2)` is justified by bounding the weighted error by its tail supremum times the summable power tail, even if `A=0` or its sign is initially unknown. Hence `n^(1−μ)p_(n+1)→A/(1−μ)`. The independently established limit is `c>0`, forcing `A=(1−μ)c>0` and eventual strict adjacent ordering. This establishes neither a computable cutoff nor ordering at all smaller labels.

## Exact controls

`timeout 15s taskset -c 13 nice -n 10 python3 -Bu round2/checks/adjacent-age-check.py`

```
456 exact harmonic rows passed; nonzero edge signs=[False, True]
T3 zero-horizon threshold-2: total=serial=7/6, old-target sector=0
```

The independent checker enumerates actual labelled births, evaluates finite-horizon degree-hit probabilities by exact harmonic recursion, and compares the direct two-birth expectation to the unordered-edge and per-core forms. It includes truncations with zero future births and with up to two future births. The fixture rejects reversing the old-edge sign and dropping the serial term; the latter already affects the zero-horizon event in `T3`. The asymptotic continuity and sign argument is analytic.
