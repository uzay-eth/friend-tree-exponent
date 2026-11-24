# Cold audit: p-adic tests on successive physical generation polynomials

**Verdict.** The interval exclusions survive independent arithmetic checks, including on the narrower post-EDGE interval. Genuine new primes are visible in **finite, exact, fixed-discount generation increments**, and they strengthen the elementary small-value bounds for reduced finite generation differences. They do **not** rule out any additional candidate for the real critical value. Even p-adic convergence with a nonzero limiting residue does not on its own transfer a real critical equation to the p-adic limit.

## 1. Scope of the cold audit

- Using the exact [EDGE](../../../../checks/neon-edge-sol.md) and [post-EDGE](../../../../checks/edge-memory-sol.md) fractions as inputs, independent exhaustive integer tests give the same three minima in both intervals: `99/175`, height seven (`3s²−7s+3`), and `32361/57200` for the square-discriminant rationals. The Farey determinant is one and the mediant is inside both intervals. The square-case enumeration uses `u²−v²` over `uv`, already reduced when `gcd(u,v)=1`. These are tests of the arithmetic consequences of the enclosures, not re-verifications of the drift certificates.
- The scalar first-degree-prime countermodel in the prior note (`README.md`) is sound in its stated scope. The factorial shifts make the spikes summable on the positive axis, introduce no new denominator primes in their expansions, and leave the greedy correction integral coefficientwise. The first prime `p` contributes at the claimed order; the ordinary fixed-value holding factor has the opposite valuation from the degree reciprocal when `s` is `p`-integral. The flat exponential modification has the stated full-expansion agreement, and preserves positivity when added with positive coefficient. Neither construction preserves the joint physical generator.
- The resultant estimate and the raw-determinant divergence argument in the prior note have the stated direction. The size cap `N` excludes the gain path to degree `N+1`; at `0<μ<1` its gain factors dominate `(d−1)/(d+2)`. The next nontrigger birth and virgin terminal rate give the stated `1/[(N+1)^2(N+2)(N+3)]` lower bound. The inherited spectral abscissa bound gives `|D_N(μ)|≥(Hδ)^M`, whose `H^M` cannot be carried through to the reduced scalar numerator. The chronological-prefix proof is a different argument: its smallest holding pole survives reduction because maximal-multiplicity paths have strictly positive leading Laurent terms there.

## 2. Genuine successive-cap congruence

Start with the inert-slot virgin from the [physical generation identity](../../../../influence/generation/check/README.md). For `a>0`, let `P_N(a,z)` count all discounted core births through the first arrival at physical size `N`, weighted by `z^depth`, with the initial root counted once. This is a polynomial in `z` over `Q(a)` evaluated at a rational trial, or over the symbolic rational-function field in the discount. Increase the **whole-tree** size cap, not the independent sizes of marked fringe particles. The increment `Δ_p=P_(p+1)−P_p` at prime `p` counts just promotions on the size-`p` jump. All path factors are the actual rates divided by `a+n−1/D_root` for the visited size `n`. The checker enumerates and merges canonical rooted trees, constructs the coefficients as rational functions in the discount, and verifies all cap increments through size 7 at `a=99/175` against a separate, unmerged vertex-ID arrow walk. The `s·1` occupation convention gives the same discounted birth toll at a physical promotion; here the toll is booked directly and no residence after crossing the cap is charged.

Let `v_p(a)=−e<0`, with `e` a positive integer, and set `u≡p^e a (mod p)`. Through size `p−1`, all rates have denominators prime to `p`, and every inverse holding factor has valuation `e` and normalized residue `u⁻¹`. At size `p` all non-star paths still have that holding valuation and integral rates: only the root, with its silent slot, can have degree `p`. The unique all-gain path to that root star has product of gain rates `(p−2)!≡1 (mod p)` for odd primes. It emits only at depth one, at rate `(p−1)/p`. Therefore:

| Case | Exponent `B` used in `p^(−B) Δ_p(a,z)` | Reduction |
|---|---:|---|
| `e≥2` | `(p−1)e−1` | `−u^(-(p−1)) z = −z` |
| `e=1`, `u≠1` | `p−2` | `−u^(-(p−2))/(u−1) z = −u/(u−1) z` |
| `e=1`, `u=1`, `w=v_p(p(a+p)−1)<∞` | `p−2−w` | `−[p^(-w)(p(a+p)−1)]^(-1) z` |

Every other path is strictly higher valuation in these cases. The exceptional zero holding denominator is a pole of the star term and is not a permissible positive-real trial. For `a=99/175`, this gives `4z` at `p=5` and `5z` at `p=7`, with unscaled valuations 7 and 5 respectively. At integral discounts the cancellation of the star holding factor changes the comparison; the unique-minimum proof above does not apply. This is a fixed-value result for successive generation polynomials, rather than a calculation of one coefficient at large discount. It is an isolated first-prime increment for a given prime, not a uniform estimate for every future increment at that prime.

## 3. Reduced forms: a quantitative improvement short of exclusion

Set `G_N=[z²]P_N−[z]P_N`. Cancel the polynomial gcd over `Q[s]`, take the denominator monic, then extract the primitive integer numerator `g_N` with rational prefactor `c_N`: `G_N=c_N g_N/d_N`. For `a=99/175` the integer `I=175^degree(g_N) g_N(a)` has the following exact valuations. The degree and valuations are after cancellation.

| `N` | Degree | `v_5(I)` | `v_7(I)` | Factor gained in the integer lower bound |
|---:|---:|---:|---:|---:|
| 6 | 8 | 1 | 0 | 5 |
| 7 | 14 | 2 | 0 | 25 |
| 8 | 20 | 4 | 1 | 4375 |
| 9 | 27 | 6 | 2 | 765625 |

In particular `|G_N(a)|≥5^v5 7^v7 |c_N/d_N(a)| 175^(-degree)`. At `N=9` even this strengthened floor is more than `10^95` below the magnitude of the actual finite difference. The test also treats the height-seven quadratic directly: in `Q[a]/(3a²−7a+3)`, the exact nonzero integer resultants have 3-adic valuations `6,9,15,20` for the same four caps. The resultant floors thus improve by factors `3^6,3^9,3^15,3^20`. Simply detecting a common projective root at infinity modulo 3 only yields the first factor of 3. A finite difference of depths 2 and 1 is not forced to vanish at criticality: the physical depth result is that differences tend to zero at **increasing depth after the full size limit**. The tables improve arithmetic bounds on these finite forms, not rational or degree-height exclusions for `μ`.

## 4. Why a mod-p difference still does not force a non-arithmetic root

These two exact examples use a single scalar generation kernel, with rational holding rate 2 and positive rational path weights. Write `a=99/175`, `L(s)=(a+2)/(s+2)` and `K_N(s)=t_N L(s)`. Each increment `t_N−t_(N−1)>0` is another retained path; generation `m` has coefficient `K_N(s)^m`. All finite volume/depth truncations are rational, the real limit has `K(a)=1`, and the corresponding real crossing is simple (`L'(a)≠0`). They are models of the proposed inference, not models of the friend-tree law.

1. With `t_N=1−5^(−N)` the real limit is 1 while the finite successive-depth difference at `a` has valuation `−N(m+1)` and leading residue `(-1)^(m+1)` after scaling. There is no 5-adic size convergence. So even an indefinitely recurring unit witness after clearing denominators is compatible with the rational crossing.
2. Set `θ=2/3`, `t_0=0` and `t_N=θ+(1−θ)5^N/(5^N+1)` for `N≥1`. Now the same sequence converges to 1 over the reals **and** to `2/3` in `Q_5`. Every successive-depth difference at the 5-adic limit is a unit, with residue `3·4^m (mod 5)`. This survives even taking the size limit p-adically first. The real depth limit still has all differences zero exactly. Thus p-adic Cauchy bounds alone, even when accompanied by persistent mod-p information, do not identify the real critical equation at the putative algebraic value.

A sufficient arithmetic bridge would be a sequence of complete finite generation differences with a real error at a hypothetical critical candidate smaller than the strengthened integer/resultant budget. It must increase both size and depth, account for normalization and conjugates, and establish nonvanishing for those same forms. Alternatively a special-value identity would have to identify the two completions. The known real, full-size depth decay and the finite prime increments do not supply either bridge.

## Replay

From `[local path]`:

```sh
timeout 25s taskset -c 13 nice -n 10 /usr/bin/time -f 'elapsed=%es rss=%M' python3 -IBu round2/open/arithmetic/critical-value/review/check.py
timeout 8s taskset -c 8 nice -n 10 python3 -IBu round2/open/arithmetic/critical-value/check.py
timeout 8s taskset -c 8 nice -n 10 python3 -Bu round2/open/arithmetic/check.py
```

The independent replay exited 0 in `4.79s` (`rss=14464` KiB). Both earlier replays also exited 0. The earlier critical-value replay reports 18 first-prime coefficient congruences and the greedy prefix `100100000000000`; the parent replay confirms the exact 19-transition reset and the omitted physical paths. The independent replay uses only Python's standard library and no precomputed physical rows.
