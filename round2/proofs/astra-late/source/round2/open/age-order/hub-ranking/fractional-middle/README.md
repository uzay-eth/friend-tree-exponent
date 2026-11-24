> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited. See the [release index](../../../../../../README.md#claim-27) for exact scope and citations.

# Fractional physical moments: exact stock boundary and signed interval

Uniformly on every physical serial-middle completion and every **real `0<p≤1`**,

```
E[W_middle^p] / Γ(1+p) > 0.55 .
```

The signed conditional moment gap is **strictly positive on every old tree for `473/512≤p≤1`**. On a tree with at least two old cores the normalized gap numerator is at least `1/40`; the seed and stars are strictly positive by separate rows. The serial middle alone beats the serial tip throughout `[1/8,1]`; that comparison does not pay the old-edge debits. This is an original-physical-clock comparison; use the [common completion factor](../mellin/README.md) `Γ(1+p)/[n(n+1)(m)_p]` to convert its normalized numerator to the conditional `Z^p` gap, with `m=n+2`. At random intermediate clock events use discounting, not independence of the size amplitude at those events.

## Exact stock boundary without a moving cutoff

For `p>0` set `P_p(j)=(j)_p/Γ(1+p)` at positive integers and define

```
φ_p(0)=0,
φ_p(j)=e^(-1) Σ_(k=0)^(j−1) P_p(j−k)/k!,   j≥1.
```

This solves the discounted harmonic equation **at every positive integer**, including the absorbing boundary at zero:

```
j φ_p(j+1)+φ_p(j−1) − (j+1+p) φ_p(j)=0.
```

Indeed, the generating series `H(z)=Σ_(j≥0)φ_p(j)z^j=e^(z−1) z (1−z)^(-p−1)` solves `z(1−z)H'+(z²−(1+p)z−1)H=0`. The coefficient formula is positive and increasing in `j` and in `p`, since `P_p(j)=∏_(i=1)^(j−1)(1+p/i)`. Also `φ_p(j)∼j^p/Γ(1+p)` by the Gamma ratio and dominated convergence in the factorial sum (each term divided by `P_p(j)` is bounded by `1/k!`). This is a positive boundary from stock one even as `p↓0`; it does not require discarding stock below a shift of order `1/p`.

In the real tree, track any subset of the marked vertex's leaves. Their own senders give stock birth rate at least `j`; promoting a tracked leaf occurs at total rate `j/d≤1` when the marked degree is `d`. Track all fresh arrivals at the mark. Increasingness and the harmonic equation give nonnegative discounted drift for `e^(-pt)φ_p(j(t))`; at stock zero the drift is also nonnegative. Localize with the total-size Yule clock. Its rescaled population has moments uniformly bounded through any fixed higher power, so the stopped comparison extends to fixed time and then to infinite time by uniform integrability. Since `j(t)≤D(t)` and `e^(-t)D(t)→W`, the boundary yields `E[W^p]/Γ(1+p)≥φ_p(j)` from any starting state with tracked stock `j`.

## Finite-degree amplification in the actual tree

Use `V_p(40,j)=φ_p(j)`. Below degree 40 use the [physical recurrences in the Mellin proof](../mellin/README.md#the-serial-middle), with genuine return attempted from the just-promoted degree-two neighbour when the stock hits zero. The return truncates at neighbour degree 96. Explicitly, with `B_d(96)=0`,

```
B_d(k) = [1/k+(k−1+1/d) B_d(k+1)]/[p+1/k+k−1+1/d],
V_p(d,0) = B_d(2) V_p(d+1,1),
V_p(d,j) = [j V_p(d+1,j+1)+(j/d) V_p(d,j−1)]/[p+j+j/d].
```

The return row is decreasing in neighbour degree by `B_d(k+1)≤1/(1+p(k+1))`. Extra physical gains at the marked vertex are helpful: `V_p(d,j)≤V_p(d+1,j+1)`, and `V_p(d,j)` increases in `j`. The proof is the backward induction in the Mellin note, which only requires `p>0` and an increasing terminal row. In particular at stock zero an earlier additional gain restarts with at least `V_p(d+1,1)`; at positive stock the additional gain contributes nonnegative discounted drift. Thus `V_p(2,1)` is a physical lower bound, not the moment of a process obtained by literally suppressing arrows on the real tree.

For an exponent box `[lo,hi]`, evaluate the boundary polynomials at `lo` and all positive discount denominators at `hi`. The last term of the alternating `e^-1` partial sum is the negative term of order 17, and each nonnegative recurrence result is rounded down. This is also valid on the first box `lo=0`: only the boundary polynomial is evaluated at zero, with no inference from a zero-threshold probability. The exact replay certifies `.55` uniformly over 256 closed boxes. For example the pointwise lower rows at `p=1/4,1/2,1` exceed `.6798,.6240,.5637` respectively.

## Fractional tip caps and whole-tree discharge

The [shifted two-race tip cap](../tip-moments/README.md) remains valid for `0<p≤1`: the comparison uses positive monotone physical moments, not convexity of the power function. Evaluate `A_p(a)=Γ(a+p)/(Γ(a)Γ(1+p))` using the Euler product

```
A_p(a) = Π_(k≥0) [1 + p(a−1)/((k+1)(k+a+p))].
```

For `a≤1` the factors decrease in `p`, are at most one, and the omitted tail may be bounded by one. For `a>1` use `hi` in the factors and bound the product from `k=48` onwards by `exp(hi(a−1)/48)≤1/(1−hi(a−1)/48)` (the sum of `1/(k+1)^2` is at most `1/48`). Use endpoint discount denominators in the positive mixture form of the shifted first-race cap. For larger parent degrees a convenient bound is `A_p(1+1/k)/(1+pk)≤1/(pk)`: the Euler product increases in `p` at this shape and `A_1(1+1/k)=1+1/k`, while `p≤1`. Favourable parallel tips use the [physical lower recurrence](../mean-ledger/README.md#endpoint-moment-inputs-and-old-leaf-charge) at the high end of the exponent box with terminal zero at parent degree 32.

The signed certificate uses the first-moment skeleton ledger with parameters `a=.0085, b=.008, fee=.004, terminal_fee=.0131`. At each retained vertex pay the exact leaf-direction profile using the favourable caps, then pay directed excess for non-chain ports and their fee for chain ports. At fixed port count, the residual is concave in receive, so check the two actual ends. When directed excess vanishes, concavity allows just the two extreme port counts. At old degree `d≥32`, use `u=1/(33 lo)` for adverse core-parent caps and `s U≤u/33`. The core-side leaf debit vanishes if `L(2)>u`, and stock credit is at least `(l/d) z`, where `z=(32/33)g+(g−U(2))/2` and `g=M−U(2)`. With `c=d−l` core ports the remaining credit before paying `2a+b` is at least `(l/d)z+(a−fee)c`, a convex combination in `l/d` of `z` and `(a−fee)d`. The verifier bounds both terms, even using `d/2` in place of `d` in the second. For the chain endpoints use the two exact receive corners below degree 10, the one-sided upper `U(d+1)/(d+1)` for degrees 10–30, and the analytic cap `1/(lo·32²)` beyond. Chain debit is convex in the endpoint reciprocals on their retained ranges. Lengths 1–3 suffice because each further bare interior site gives positive serial credit and an equal-coefficient edge. No old trees are enumerated.

These bounds leave the signed moderate range below `473/512` open here. In particular, the all-fractional serial lower does not by itself imply order: it must still pay the adverse edges. The independently handled tiny-power argument is not used in this certificate.

```
nice -n 10 python3 -IB round2/open/age-order/hub-ranking/fractional-middle/check.py
PASS 256 boxes: physical serial middle > .55 on 0<p<=1
PASS serial-only margin on 224 boxes covering [1/8,1]
PASS signed whole-tree margin on 39 boxes covering [473/512,1]
```
