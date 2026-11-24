> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../README.md#claim-26) for exact scope and citations.

# A tree discharge for the pre-placement hub gap

At threshold zero the older newborn's conditional hub probability is strictly greater on **every finite old tree of maximum degree at most nine**, with no size cutoff. The comparison also holds on all stars, and on trees in which every core either has degree at most seven or has degree at least ten with exactly one core neighbour. The last class permits arbitrarily many high-degree vertices with arbitrarily large leaf stocks. No all-degree comparison or physical counterexample is certified.

## Local charge

Use the signed numerator and bounds from the [retry certificate](../README.md). Write `g_p` for the serial gap lower bound including the old-neighbour first-arrow credit, `L(p)` for the return-race lower bound on a parallel leaf with parent degree `p`, and `U(p,j)` for the first-send upper with `j` other leaves. All coefficients below are for the original two physical births; to get the probability gap divide by `n(n+1)`.

For an old core, its profile `a=(d; d_1,...,d_c)` records its total degree and the total degrees of its core neighbours. Put `l=d-c`, `w=(l+Σ_i 1/d_i)/(d+1)`, `s=w/d` and `s_l=1/(2d)`. Pay all old leaf edges within the profile. With `z_+=max(z,0)` the resulting vertex credit is

```
C(a) = w*g_(d+1) + l*g_2/(2d)
       - l*((s_l-s)_+*(U(2,0)-L(d+1))_+
            +(s-s_l)_+*(U(d+1,l-1)-L(2))_+).
```

Omit the last line when `l=0`. This uses both parallel tips: when the old leaf endpoint is adverse, it has a degree-two parent with no other leaves; its opposite has parent degree `d+1`. When the core endpoint is adverse, its parent retains only `l-1` other leaves, because the placement at the adjacent old leaf promotes that leaf. A core–core edge between matching profiles `a,b`, with degrees `d,e` and stocks `l,j`, has debit at most

```
D(a,b) = (s(a)-s(b))_+ * (U(d+1,l)-L(e+1))_+
       + (s(b)-s(a))_+ * (U(e+1,j)-L(d+1))_+.
```

Both edge signs depend on the **same** coefficient difference; independently charging two worst-case endpoints would discard the correlation. For a nonstar, each core has at least one core neighbour, so the formulas for `U` always retain its needed permanent nonleaf slot. A single-core tree is a star: at centre degree `d≥2`, its coefficient `(d+1)^(-1)` is greater than the old leaf's `(2d)^(-1)`. For `d≥4` the adverse centre-side leaf has send cap at most `1/d≤1/4<L(2)`, so the edge cannot be adverse. Degrees two and three use exact retry-replay rows. The size-two tree is the seed edge and its hub comparison is the [previous strict `p3>p4` bound](../../../early/README.md#signed-finite-clock-pullback).

## Directed subtree potential

For a parent profile `a` of degree `d`, and one of its core-neighbour degrees `e`, prescribe `M(a,e)` equal to the charge for the **terminal-core** child `b0=(e;d)`: one core slot, `e-1` old leaves. Specifically, `M(a,e)=C(b0)-D(a,b0)`, using the rounded lower and upper charges from the replay. This prediction is a lower bound for a whole oriented subtree if, for every compatible child profile `b=(e; d,e_2,...,e_k)`,

```
C(b) + Σ_(i=2..k) M(b,e_i) - D(a,b) >= M(a,e).
```

The parent slot `d` is removed **once**, even if several neighbours have that degree. Prove the bound by induction from the leaves of the old core tree. At the root all slots are child slots, and the lower bound is `C(a)+Σ_e M(a,e)`, with multiplicities. Individual messages may be negative; no positivity of post-placement differences or vertex rows is required. The certificate says the terminal-core child attains the minimum of the **rounded comparison recursion**, not that it minimizes the true physical hub probability.

With maximum degree nine there are 48,602 profiles and 194,472 messages. The minimum verified root charge is `27404273/10^8`, with positive sign at every root. The implementation reduces parent/child checks by enclosing the parent coefficient tick in its interval for the fixed `(d,e,l)` class. The debit is a maximum of two nonnegative linear terms and zero. For a single child coefficient tick `y`, on integer parent ticks it is affine to the left of `y−1` and right of `y+1`; check all three ticks in between. The difference of the raw debit for a test child and the raw debit for the terminal child therefore has its maximum at interval ends or those specified ticks for either child. If the rounded debits are ceilings of the raw debits divided by the common tick scale, their difference is at most the ceiling of the raw difference divided by that scale. This preserves the subtraction of the reference debit. `--pairs` verifies the same induction directly for every compatible profile pair instead.

## One unbounded-degree class

For `--stocked 7` the extra symbol represents degree `d≥10` with one core neighbour. Root at a small core whenever one exists, so large cores are terminal. A small core treats each such adjacent reciprocal as in `[0,1/10]`. The vertex credit is concave in `w` (a linear function minus positive parts), so its lower envelope is the minimum of the two endpoint values. The same endpoint profiles delimit `s` for debit rounding. For the large core, `w≥9/11`, `(d−1)/(2d)≥9/20`, and `s≤1/11`. Its paired leaf debit vanishes: `(d−1)/(d(d+1))>1/(2d)` and `U(d+1,d−2)≤1/(d−1)≤1/9<L(2)`. The first-send bound for its **core-edge** endpoint is `1/d≤1/10`, with `d−1` other leaves. The lower credit and interval debit tests verify all 6,426 profiles and 21,020 messages in this mode. If there is no small core, the only nonstar case is a two-core tree; each credit is bounded below by `(9/11)*0.16694235+(9/20)*g_2`, and the sole core-edge debit is at most `1/110`.

## Replay and unpaid cases

```
python3 -IB round2/open/age-order/global-transport/hostile/discharge/check.py 9
python3 -IB round2/open/age-order/global-transport/hostile/discharge/check.py --stocked 7
python3 -IB round2/open/age-order/global-transport/hostile/discharge/check.py --pairs 6
```

All three commands pass. Their minimum root ticks are respectively `27404273`, `22124577`, `27404273` over core trees with at least two vertices. The degree-ten all-profile sweep did not finish within 50 seconds; it gives neither sign. The method reduces the general case to an unbounded profile inequality. In particular, the extra degree symbol does not represent a crowded **bare** core or one with several core slots; those have different debit slopes and terminal-child reference credits. No conclusion for those profiles follows from the capped sweep. This proof is at hubhood; the earlier tree-dependent right-continuity consequence is not a uniform positive birth-scale result.
