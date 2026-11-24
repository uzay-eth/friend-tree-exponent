CORRECT

# Cold challenge: hub-spectrum.md

Scope. End-to-end challenge of `hub-spectrum.md` (`fringe/hub-spectrum.md`), conditional on the already reviewed
`almost-sure.md` (`fringe/almost-sure.md`), `degree-law.md` and `hub-mass.md`.
No source edits, no commits, no long enumeration. The model definition was cross-checked against
`existence.md` (`proofs/existence.md`) (uniform vertex, then uniform neighbour, attach a leaf).

## Mathematical gaps

**None found.** I initially misread U_k(t)=N_{≥k}(T_t) as nonincreasing in t; it is nondecreasing, and with the correct
direction the moving-window sandwich U_k(v)−U_{k+1}(w) ≤ N_{=k}(t) ≤ U_k(w)−U_{k+1}(v) for v≤t≤w is valid, as verified
on real growing trees. Every item in the brief is addressed below; no gap was identified in the audited argument.

## 1. Quantitative uniformity of (3) and the common conditioning time

The three bounds in (3) are exactly (4), (7) and the Yule size identity of `almost-sure.md`, applied to
f_k=1_{D≥k} with |f_k|≤1≤W; the constants there are uniform in |f|≤W, so they are independent of k. I verified the
exponent bookkeeping at the mesh u=(1−a)j covering [j,j+1]: with a=(2μ−1)/8, g=η(2μ−1)/100, b=2g, the Chebyshev sum
over ⌈e^{bj}⌉ mesh points and 2≤k≤e^{g(j+1)} with weight k^{2(1+μ)} has exponent
1+4a−2μ+b+(3+2μ)g = −(2μ−1)/2+(5+2μ)g, which is strictly negative since η<μ<1. The conditional-mean error needs
(1+μ)g<ηa, which reduces to 1+μ<12.5. Both hold with margin for every μ∈[0.565,0.5665] and every σ<μ, as my checker
confirms. Conditioning the whole block at the same u avoids differencing two martingale errors, and the extension to
intermediate times via U_k(v)≤U_k(t)≤U_k(w) costs k(w−v)Z*_{u}=o(1) uniformly because kΔ≤e^{g−g j}→0 (b=2g).

## 2. Gamma-minus-Exp coupling, untracked gains and the compensator

In this model a leaf has degree one, so its rate-one ring always attaches to its unique neighbour, the focal core:
tagged leaves at the core never lose membership by their own rings and form a birth process of rate equal to their
count; promotions move a tag to an auxiliary rate-one Yule family, so the total tagged population never loses members
and has the Gamma(L,1) limit G_L. I verified the leaf-degree claim on 1,900+ sampled leaves of small trees. The
dominating population for the missing mass has immigration rate one and unit birth rate per member; its size plus one
is a Yule started at one (generator n→n+1 at rate n+1 in both descriptions), so its normalized limit is E~Exp(1).
Incoming untagged leaves are omitted only in a lower bound, which is the favourable direction. For the birth-epoch
bound (6) I verified the antiderivative d/dh[−e^{−k e^h}/k]=e^{h−k e^h}, giving ∫_0^∞ e^h e^{−k e^h}dh=e^{−k}/k
exactly at k=1,3,12; the conditional birth intensity n(T_t)e^{r−t} is the Yule mean. No independence between vertices
is used.

## 3. Maximum c, size/time change and powers of Y

The maximum-core-degree argument is valid: a vertex gains core neighbours only on its own rate-one rings, so c_v(T)
is stochastically at most Poisson(T), E e^{c_v(T)}=e^{(e−1)T}, and summing over E n(T)=O(e^T) gives O(e^{eT}), hence
max c_v(T)=O(T) a.s. Since k=e^{T−t}=(t+2)^2 and T≈t, the condition c_v(T)≤εk/4 holds eventually. The change of
size/time is handled consistently: Z_v^{CT}=YZ_v, e^{μT}k^{−μ}=e^{μt}, and X=d_2 Z*_∞/Y^μ; substituting e^{−t}=Yx yields
x^μ#{Z_v≥x}→A(d_2 Z*_∞/Y^μ)=AX with the **same** A as the fixed-degree tail (degree-law.md via d_k/d_2=P_k~A k^{−μ}).
I verified the ratio algebra numerically. The randomness of Y is absorbed by intersecting all rational meshes and
using monotonicity, which bounds intervening times with factors e^{±μ/m}; this is the implicit uniformity needed for
deterministic x, and it is present.

## 4. Tail inversion and intermediate cumulative scales

Generalized inversion of the monotone count gives Z_(j)~(AX/j)^{1/μ} (checked: Z_(j)^μ j=AX). For s>μ,
∑Z_v^s=s∫x^{s−1}#{Z_v≥x}dx~sAX/(s−μ), finite exactly for s>μ with (s−μ)∑Z_v^s→μAX; for s≤μ the integral diverges,
so the threshold is correct. Section 4's two regimes cover every intermediate scale: Section 1 gives k≤n^γ (γ<g) and
Section 4 gives k≥2⌈log²n⌉, k/n→0, and these ranges overlap; the union bounds are summable because l_n~(log n)^2 gives
exp(−c log²n), and future births are covered by (6) with n e^{−l_n}/l_n. The event is common across meshes and n.

## 5. Controls

Commands, from the repository root on one core:

```
taskset -c 0 python3 -B round2/checks/hub_spectrum_hy4.py
```

Output: `HUB SPECTRUM COLD CHECK: PASS`. Planted errors rejected: wrong thinning (dropping a power of k in the weight
breaks summability), wrong normalization (X without Y^μ changes the (1) constant), reversed squeeze (empty interval),
and lost tags (a leaf with degree >1 would invalidate the pure-Yule tagged population). The moving-window inequalities
were checked on 60 growing friend trees.

## 6. Limitations

The inputs from almost-sure.md ((3)-(9), Z*_∞>0, d_k/d_2=P_k), degree-law.md (A, tail), and hub-mass.md (∑Z_v=1) were
taken as reviewed and were not reproved. I did not rerun any long enumeration. Within this scope the claims of
hub-spectrum.md — simultaneous degree counts, the a.s. hub limit with the common A, and all intermediate cumulative
scales — are supported.
