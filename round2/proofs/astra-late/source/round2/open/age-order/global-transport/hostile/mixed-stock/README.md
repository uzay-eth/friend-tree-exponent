> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../README.md#claim-26) for exact scope and citations.

# High-degree core skeletons with arbitrary live leaf stock

At threshold zero, the whole-tree pre-placement gap is strictly positive if **every core with at least two core neighbours has total degree at least ten**. Such an internal core can have any number of actual leaf neighbours, including none. Terminal cores can have any degree. This extends 42b5's bare-skeleton result to live-stock hubs. A small internal core adjacent to a high-degree core is not in this theorem.

## Endpoint charges

Use the [retry and paired-edge certificate](../README.md) and its [local discharge notation](../discharge/README.md). At a core of total degree `d`, core census `c`, and leaf stock `l=d-c`, write

```
w=(l+Σ_core-neighbours 1/D_v)/(d+1),  s=w/d,  s_l=1/(2d),
C=w*g_(d+1) + l*s_l*g_2
    − l*(s_l−s)_+*(U(2,0)−L(d+1))_+
    − l*(s−s_l)_+*(U(d+1,l−1)−L(2))_+.
```

Terms with `l=0` are omitted. Here `g_p` includes the old-neighbour first-gain credit, `L(p)` is the favourable leaf lower bound, and `U(p,j)` is the first-send upper with `j` other leaves. Thus `C` has already paid the serial-tip and paired old-leaf debits. For a core–core edge with actual endpoint parameters `(d,l,s)` and `(e,j,t)`, the debit goes in the direction of the **same signed difference** `s−t`:

```
max{0, (s−t)*(U(d+1,l)−L(e+1))_+,
       (t−s)*(U(e+1,j)−L(d+1))_+}.
```

In particular the first-send bound at each endpoint has its own actual stock; it does not use a bare continuation when the endpoint has leaves. We bound this maximum, not the sum of two independent edge costs. The elementary tail barrier `U(d+1,l)≤1/(l+1)` from the first-send argument is valid at every degree. It also bounds the scalar upper used in the replay: at the cutoff it is equality, and in the recursion it is a supersolution. Since `s≤(l+c/2)/[d(d+1)]`,

```
s U(d+1,l) ≤ (d+l)/[2d(d+1)(l+1)] ≤ 1/[2(d+1)].
```

The same calculation with the actual opposite endpoint bounds the other direction. Thus an internal–internal debit is at most `b=1/22`, regardless of either endpoint's stock.

For a terminal core next to an internal core, its credit after paying its one core edge exceeds `a=0.00634`. To see this, let its total degree be `e`. For `2≤e≤9`, the terminal coefficient depends affinely on the adjacent reciprocal in `[0,1/10]`. The credit `C` is concave in that reciprocal. Check `C−s U(e+1,e−1)>a` and `C−b>a` at both endpoints, then interpolate. These pay respectively the terminal-outgoing debit and the internal-outgoing debit. For `e≥10`, the stock is `e−1`, so `s>s_l` but the old-leaf edge charge vanishes: `U(e+1,e−2)≤1/(e−1)<L(2)`. Hence `C≥(9/11)·0.16694235+(9/20)g_2`, which pays even `b` with surplus `>a`. These are the two-endpoint terminal bounds from 42b5, retaining a slightly larger common surplus for the mixed-stock census.

## Internal credit and census

For every `d≥10`, `2≤c≤d`, prove `C+a(c−2)>b`. Through `d=112`, exact scalar verification uses the two endpoint profiles for `Σ_core 1/D_v∈[0,c/2]`. Concavity in `w` then pays every actual neighbour list between them. This is a census check, not enumeration of old trees or their futures.

For `d≥113` split at `l/d=1/2`. If `l≤d/2`, then `s≤(d+l)/[2d(d+1)]` is **not** sufficient on its own to imply `s≤s_l`; instead bound whichever leaf debit occurs. In the `s_l≥s` direction it is at most `(l/d)·0.611/2`. In the `s>s_l` direction it is at most `l s·U(d+1,l−1)≤s≤1/(d+1)`, but this direction is in fact already paid by the favourable bound whenever `l≥4`, because `U(d+1,l−1)≤1/l≤1/4<L(2)`. When `l≤3` and `d≥113`, `s>s_l` can occur, but its debit is still at most `(l/d)·0.611/2`: use `s−s_l≤(l−1)/[2d(d+1)]`, and the first-send upper is at most one. Consequently, putting `z=l/d≤1/2`,

```
C ≥ z*(0.16694235*113/114 + g_2/2 − .611/2).
```

The bracket is negative, and `c−2≥109/2`. Substituting `z=1/2` and `c−2=109/2` yields `C+a(c−2)>b`. If instead `l>d/2`, then `l≥57` and again the `s>s_l` direction has zero leaf debit. In the `s_l≥s` direction use `s≥l/[d(d+1)]`: the debit is at most `.611(l/d)(1/2−l/(d+1))≤(l/d)/(2(d+1))`. The old-leaf serial credit `(l/d)g_2/2` pays this. What remains is at least `0.16694235*113/(2·114)>b`.

Let `B≥1` be the number of internal cores. They induce a connected tree (deleting the terminal cores), with `B−1` internal edges. There are exactly `2+Σ_internal(c−2)` terminals. The combined numerator is therefore greater than

```
 B*b − a Σ_internal(c−2) − (B−1)*b
     + a*(2+Σ_internal(c−2)) = b+2*a.
```

For no internal cores the old core tree has at most two vertices. Stars and the edge seed have the previous positive certificates. For two cores, if both are at most nine use the bounded-degree discharge. Otherwise treat the smaller side as terminal next to a high-degree side and apply the terminal charge just proved (the proof used only the opposite endpoint's high total degree, not its being internal). The high-degree side itself has positive credit by the terminal computation above. This exhausts the class. Divide the numerator by `n(n+1)` to obtain the conditional gap.

## Replay and remaining interface

`python3 -IB round2/open/age-order/global-transport/hostile/mixed-stock/check.py`

The replay evaluates the finite scalar inequalities in rational arithmetic and computes the paired debits from both actual endpoint profiles on mixed physical trees. It includes a bare degree-ten core adjacent to live-stock internal cores, including one of degree 113, and a second fixture with several different internal and terminal degrees. The general mixed tree, with both low-degree and high-degree **internal** cores, still needs a matching boundary transfer. This does not give a negative physical example or a uniform window at positive thresholds.
