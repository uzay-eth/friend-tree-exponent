# Independent cold review: 141 and the all-label upper interval

**Both claims hold.** Reviewed the paired-edge proof (`open/age-order/global-transport/annealed/extension/README.md`) at `a9995a53eee76426fd0a47b2888e1070e5316582`, using the [previously audited 232 argument](preplacement-256-sol.md) and the pinned Fable early-pair transport and strictness proofs. Put `b=n+1` (older label) and `m=b+1` (common final size). For each fixed old tree the conditional adjacent gap is strict for `0<x<1`, `bx≥141`; it is also strict for every `b≥20`, `10/11≤x<1`. Combined with the *unconditional* early certificates this proves the unconditional all-label order on `[10/11,1)`, with equal founder laws and otherwise strict pairwise inequalities for each such threshold.

## Signed edge split and the one-third cap

An old edge `uv` contributes `(s_v−s_u)(h_u−h_v)`. If `s_v>s_u`, dropping the favorable term leaves a debit on the newborn attached at **v**, the larger-`s` endpoint. Hence it is the recipient's old degree that determines whether the locked cap is `1/2` or `1/3`. Other siblings, including the case of zero other siblings, are irrelevant. The serial-tip cap is `1/2` (its parent, the older newborn, has degree two). All tests are at the same physical size and threshold. The sum of absolute edge coefficients is at most `Σ D_u s_u=A` and need not be equal to `A`.

For each core of degree `D≥2` with `L≥1` leaf neighbours, the old-leaf serial mass is `A_l=L/(2D)` and `s_o≥L/[D(D+1)]`. Therefore the expensive leaf-endpoint debit is at most `A_l(D−1)/(D+1)`, while `w_o≥2D A_l/(D+1)`. Indeed `3(D−1)≤3D+1`, so `3X_o≤A_l+w_o` and summing proves `X≤A/3`. Cores with zero leaves add no such debit; the size-two edge has equal endpoint values and no adverse weight. If `θ` is the relative small-beta cap, the edge debit is bounded by `U[X+θ(A−X)]`; adding serial tips yields `UA[(1+θ)+(1−θ)X/A]`. This is at most `UA(4/3+2θ/3)` since `θ≤1`.

## Beta constants and rigorous corners

The protected initiator-colour argument gives the whole lower coefficient `c_b` of the earlier audit for every `b≥20`; the former cutoff 232 was used only at its last sign test. Here the log-series upper, exponential-series lower and integer ceiling on the protected prefix are unchanged. In particular `c_b` is positive and increasing in `b` (both factors depending on `b` increase).

For the upper cap `U`, gamma Jensen gives `Γ(m)/Γ(m−1/2)≤√(m−1/2)`, with the beta integral dividing by `m−1/2`, so the proof may keep the full `sqrt(1−x)` factor. In the smaller-cap pointwise ratio the exponent is `((1−z)/z)^(1/6)`: with `y=bx`, gamma Jensen bounds its prefactor after factoring out `y^−1/6` by `sqrt(π)/Γ(1/3)·(b+1)/(b+1/2)`. Specifically `b^(1/6)(m−1/3)^(5/6)≤m`. The degree caps can also be ordered by their increasing beta likelihood ratio, giving the minimum with one. Integrating the degree-five Taylor lower polynomial on `[0,1]` against `t^(−2/3)` gives `141077/58240>2.4223`; with `sqrt(π)<9/5`, the prefactor is below `3/4` for `b≥120` and below `4/5` for `b≥20`.

For `bx≥141` drop `sqrt(1−x)` and use the algebraic corner `b=y=141` (no physical endpoint inference at `x=1`). Using `57/25<141^(1/6)` and `π>157/50`, the independent check gives strict `c_141² (157/50)·141 > (4/3+(1/2)/(57/25))²`. For `b≥20,x≥10/11`, `y≥200/11>(3/2)^6`, so `θ≤8/15`; use `c_20`, the lower `y=200/11`, and upper `1−x=1/11` simultaneously. Both are conservative corner comparisons, and the signs improve as `b` or `x` increases within their stated band.

## Reachability and the all-label splice

The degree-identical examples have chronological zero-based parent lists `(0,0,0,1,2,3)` and `(0,0,0,1,4,5)` (first edge `01`). Each listed recipient exists and has a neighbour at its placement, so each listed history has positive physical probability. Both degree multisets are `3,2,2,2,1,1,1`. Direct incidence gives `(A,X)=(59/24,1/12)` and `(167/72,0)`, respectively. Thus the degree sequence alone does not determine `X/A`, though these examples do not obstruct a sufficiently refined occupation estimate for the actual law. An annealed use of the signed lower bound involves `E X/E A`, not `E(X/A)`.

The [pinned Fable transport](../attacks/fable5/flow/README.md) at `88e13f36b39162f6ee60561a2b465a982b34af02` covers older labels `b=2,…,19`, and its [strictness argument](../attacks/fable5/endpoint/README.md) upgrades each to strict tail order. The last needed early pair is `(19,20)` at common size 20; the new bridge starts at `b=20`, which is `(20,21)` from old size 19. Its strict pointwise conditional inequalities remain strict after averaging. No single-program late Fable receipt is needed. Founder symmetry gives `Z_1` and `Z_2` **equal in distribution** (not pathwise equality); therefore for each `x` in this interval the sequence of survival probabilities is `p_1(x)=p_2(x)>p_3(x)>…`. Adjacent order implies the nonadjacent order by transitivity, and the all-pairs statement implies the adjacent one by restriction. This is an order of tails on the indicated interval, not stochastic domination on all thresholds, nor a claim of strictness at `x=1` where all strict survival probabilities vanish. No conditional early-pair splice is inferred from Fable's unconditional certificates.

## Bounded replay

```text
taskset -c 8 nice -n 10 timeout 8s python3 -IB round2/checks/preplacement-141-sol.py
PASS 873 reachable parent arrays; equal-degree fixtures have X/A=2/59 and 0
PASS gamma constants and strict rational corner margins (141, and b=20 at 10/11)
PASS boundary ledger: 1=2; 2→3 through 19→20; 20→21 onward
```

The producer checker also passed all four of its reported sections. The independent replay does not rerun the early Fable flow networks or infer any asymptotic from the finite fixtures.
