> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited. See the [release index](../../../../../../../README.md#claim-27) for exact scope and citations.

# The whole-tree signed moment ledger for every real exponent in [1,3]

**For every finite old tree of size `n≥2` and every real `1≤p≤3`, the normalized signed numerator is at least `1/40`.** If `m=n+2` and `F_n` is the common pre-birth history, this gives

```
E[Z_(n+1)^p − Z_(n+2)^p | F_n] ≥ Γ(p+1) / [40 n(n+1) (m)_p],
(m)_p = Γ(m+p)/Γ(m).
```

This uses the frozen [Mellin lower](../../mellin/README.md); it needs no joint H3 input. The checker covers the closed interval with 256 rational cells of width `1/128` and keeps the integer first-moment constant. The `p=3` endpoint also belongs to the separate `p≥3` certificate in `mellin/`, so the two bounds together sign the conditional gap for all real `p≥1`. These are comparisons of pairs born after the old tree, including `(3,4)` from the seed edge; no separate founder-pair moment claim is used here.

## Interval endpoint bounds

The endpoint quantity is `E W^p/Γ(p+1)`. The normalization and signed placement identity are as in the [integer ledger](../README.md): the true size amplitude `S` is independent of the embedded path from the fixed common completion and `W=SZ`. No size-amplitude independence is needed after a clock race.

Fix one cell `[q0,q1]`, and use the row's certified middle lower `M` in the same normalized units. For `0<r<1`, `A_p(r)=Γ(p+r)/(Γ(r)Γ(p+1))` is the `p`-th moment of `Beta(r,1−r)`. It is decreasing and convex in `p` by the pointwise properties on `[0,1]`. For `r=1` it is constant one. Evaluate the integer chord at `q0` to bound it above over the closed cell. For `A_p(1+r)`, multiply the bound for `A_p(r)` by `1+q1/r`, by the gamma recurrence. All integer moments are exact rationals.

In the [tip ring cap](../../tip-moments/README.md), the weight on `A_p(a)` increases with `p`, so use `q1` there; for the positive average of the other two cap terms use the weight at `q0`. Round the resulting cap upward. The backward first-send recursion uses `q0` in the discount and the finite terminal upper at parent degree 64. For all degrees there is also `U(k)≤1/k`: with `r=1/k`, `A_p(1+r)=(1+p/r)A_p(r)≤p+r` for `p≥1`, so the closed-form tip upper is at most `(p+r)/(1+pk)=1/k`. Taking the minimum with `1/k` at each finite row is legitimate. Applying the pointwise physical recurrence first, then taking upper enclosures of its positive terms, avoids requiring the interval enclosure itself to be a Markov test.

For the favourable tip use the actual return race with parent-degree increment rate at most `k`, send rate `1/k`, terminal cutoff 32 and the uniform postpromotion payoff `M`. Use `q1` in its discounted lower recurrence, rounding down. Monotonicity for the comparison is a property of the exact fixed-exponent recursion: its discounted value from `k+1` is at most `M/(1+p(k+1))`, which implies that the backward value at `k` is at least the value at `k+1`. Thus extra parent gains cannot invalidate this lower comparison. Both this lower and the adverse upper use actual shifted parent degrees in the signed old-leaf and chain coefficients.

## Two uniform discharges

Use the integer ledger's old-leaf-adjusted credit `C`, serial gap `g=M−U(2)`, and `w=d s` for each retained old core. Both routes certify the seed (total serial weight one and equal edge coefficients) and stars separately, testing finite centre degrees and using the upper tail for larger stars. For `d≥32` a star has only centre-side potential debit and its centre credit minus this debit already exceeds `1/40`.

**Sitewise.** Set aside `(w+l/(2d))/40` at each core and pay its own outgoing core debits by `c s U(d+1)`. Concavity in reciprocal receive verifies this at range endpoints for `d<32`. For `d≥32`, set `h=g−1/40`, `u=1/33`: all outgoing centre-side core and leaf debits together are at most `w u`, and the opposing old-leaf debits are at most `l U(2)/(2d)`. The checked inequalities `h>u` and `(32/33)(h−u)+(h−U(2))/2≥0` give the result. The reserved serial weight sums to at least `1/40` because `Σw≥Σr_u/n=1`.

**Suppression.** Suppress exactly the bare degree-two internal chains and keep their actual endpoint profiles. Keep `a=.0085`, `b=.008` and the ordinary chain bill `.004`; a terminal of degree two pays `.0131`. At an endpoint of degree `d` with `c` old core ports, `k` of which are chain ports, the retained lower ledger is

```
C + a(c−2) − b − k·fee − (c−k)·(s U(d+1)−b)_+.
```

The `k` chain ports have paid their complete chain-side bill, so consume no direct-edge excess. At fixed integer `k`, concavity in `r` reduces the range `k/2≤r≤c/2` to its endpoints. If `U(d+1)/(d+1)≤b`, the direct excess is identically zero: then the residual is jointly concave and affine in the bill count, and only `k=0,c` are needed. At `d≥32` no direct excess remains, and the core-side old-leaf debit vanishes because `L(2)>1/33`. Put `z=(32/33)g+(g−U(2))/2`. With stock fraction `l/d`, the credit plus available core allowances is at least `z(l/d)+32(a−.004)(1−l/d)`. The two endpoints of this affine bound pay `2a+b`.

Each chain pays one baseline after both end bills, by the paired boundary and interior debit calculation from the integer proof. Use the two degree/stock receive endpoints below degree 10, exact degree reciprocal and a safe outgoing envelope through degree 30, and reciprocal endpoints `0,1/31` beyond that. Compare chains of lengths one, two and three; each additional homogeneous interior site adds positive `g/3` without further debit. These are the same rational corner tests as in the integer proof, accelerated by clearing denominators for the Cartesian product. The skeleton remains a tree with at least two retained cores. Its identities give the floor `b+2a=1/40`.

## Frozen inputs and replay

Paths below are relative to `round2/open/age-order/hub-ranking/`. The lower inputs are `mellin/README.md`, `mellin/check.py`, and `mellin/lower-1-3.tsv`, Git blobs `9e75bd41e8a204b100a086d21e9fea76ad844fb0`, `9786bbdf7230ee86b9aa518d421d9fb471866130`, and `defa8793b129e4b3c3f5ef6f23fd4033fc4ff8fc`. The tip proof is `tip-moments/README.md` (blob `a2535cca95299aca1b8163197976eebbc0f03d97`); the integer proof remains at `mean-ledger/{README.md,check.py}`.

From repository root, run `python3 -IB round2/open/age-order/hub-ranking/mean-ledger/real/check.py`. Python 3.12, standard library only. The replay authenticates the TSV, verifies its exact mesh and endpoints, and discharges all cells (`87` sitewise, `169` suppressed).
