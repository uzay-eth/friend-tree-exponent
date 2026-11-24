> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited. See the [release index](../../../../../../../README.md#claim-27) for exact scope and citations.

# Discounted physical continuation near zero

Two changes improve the endpoint bounds in the fractional Mellin ledger: an exact discounted stock hit for the serial middle and a joint stock-and-discount upper for the adverse tip. They use the original clock and the normalization `E W^p / Γ(1+p)`. At a random hit or promotion we apply the physical Markov property with discount `e^(−pT)`. The common gamma size factor enters only at the fixed two-birth completion when converting to share moments.

## Middle: use the true degree floor before cashing the stock moment

At a degree-`D` hit, with `j` tracked leaves, the degree never falls below `D`. The tracked stock's birth rate is at least `j` and death rate is at most `δ_j=min(j/D,1)`. Until the stock hits `K≥D`, use a birth–death comparison with these rates, paying zero at stock zero (even if there is a later physical return). At stock `K`, pay the [exact constant-death harmonic moment](../../fractional-middle/README.md#exact-stock-boundary-without-a-moving-cutoff) `φ_p(K)`. This terminal moment accounts for the entire remaining continuation.

On an exponent box `[lo,hi]`, let `F_j` have boundary values `F_0=0`, `F_K=φ_lo(K)` rounded down. Use `hi` as the discount on this box. Solve the finite equations without matrix inversion:

```
R_0 = 0
R_j = j / (j + δ_j + hi − δ_j R_(j−1)),   1 ≤ j < K
F_j = R_j F_(j+1).
```

Here `R_j` is the discounted probability of hitting `j+1` starting at `j` before zero in the comparison. It is in `(0,1]` for `j≥1`, so `F` increases with `j`. The finite equations give nonnegative discounted drift on the physical path: extra stock births and any reduction in stock deaths help. Stop at zero or `K`, first with a finite time horizon; the bounded time-cap remainder vanishes for `p>0`. This proves the lower without identifying the true stock with the comparison or assuming anything about the exterior at the hit. At `lo=hi=0` and `K=D` the ratios recover the gambler's-ruin probabilities with per-jump death/birth ratio `1/D`. The positive-power solution varies continuously up to this value.

Use this boundary in the existing physical return recursion below degree `D`, retaining the actual newly promoted neighbour on last tracked loss. For the 128 boxes `[j/512,(j+1)/512]` covering `(0,1/4]`, take `K=D`, with `D=8,12,16,20` on consecutive groups `j=0..23,24..53,54..91,92..127`. The new serial floors improve the old ones by more than `1/4000` on every box. Examples:

| Exponent box | Old floor | New floor |
|---|---:|---:|
| `[0,1/512]` | .773906869 | .776331282 |
| `[31/512,1/16]` | .742995165 | .744150331 |
| `[63/512,1/8]` | .716414109 | .717080577 |
| `[127/512,1/4]` | .674842491 | .675115827 |

The old degree-40 cutoff is not optimal for this boundary: at that cutoff, the chance of reaching stock 40 before zero can cost more than the improvement from the reduced death rate. The listed cutoffs are fixed by exponent box and the replay checks each box against the original floor.

## Tip: combine sibling stock and Laplace discount

Before promotion, let the tip have parent degree `k` and `j` tracked other leaf siblings; at least one permanent nonleaf neighbour of the parent ensures `k≥j+2`. Set

```
f(k,j) = 1/(j+1+pk).
```

This is an upper bound on the expected discounted indicator of eventual promotion, with value one at promotion. To verify it, write `z=j+pk` and let `p>0`. The tip promotion rate is `1/k`, the loss rate of tracked siblings is `j/k`, and their stock and parent degree increase together at rate at least `j+1`. Track any extra fresh leaves at the parent as well. These gains decrease `f`. The generator residual after discount, multiplied by `k`, is

```
1−f + j[f(k,j−1)−f] + k(j+1)[f(k+1,j+1)−f] − pk f
 = j/z − k(j+1)(1+p)/[(z+1)(z+2+p)].
```

The term involving `j−1` is absent at `j=0`. For `x=k−j−2≥0`, the numerator of the negative of the residual after multiplication by `z(z+1)(z+2+p)` is

```
j(j+1)x
+ p[(j+2)²x + (j+1)x² + 3j+4]
+ p² k(x+2).
```

Every coefficient is nonnegative. Thus optional stopping at promotion, with the running discount and a finite horizon first, gives the stated upper for any physical exterior. At `p=0` this reduces to the original stock first-send barrier. For `p>0` it combines the two separate bounds `1/(j+1)` and `1/(1+pk)` in a single denominator. The formula is also valid with `j=0`, where the stock-loss term is zero.

At the degree-`k` cutoff, multiply `f` by an upward bound on `A_p(1+1/k)` for the postpromotion normalized mass: parent degree at actual promotion is at least `k`. This uses a uniform continuation upper conditional on the promotion history, with no independence assumption. Below cutoff use the same shifted two-race promotion payoff and downward-in-degree supersolution as the [stock-sensitive ledger](../../mean-ledger/fractional/README.md), replacing its boundary and taking the smaller Euler-product cap when helpful. With cutoff 32 the seed tip `(2,0)` on `[0,1/512]` changes from `.611822835` to `.611327797`; on `[63/512,1/8]` it changes from `.552539093` to `.551267234`.

## Ledger audit and replay

The stock-sensitive ledger's original cutoff and large-degree debit bounds are valid too. If `l` is the old stock at degree `d`, the core-edge exposure using its stock-only cap is at most `(d+l)(1+1/(d+1))/(2d(d+1)(l+1))`, hence at most half the unstocked prefactor divided by `d+1`. This follows from `d+l≤d(l+1)`. In the opposite type of old-leaf edge the shift from `l` to `l−1` is needed before multiplying the adverse cap by `l`. The large-degree branch pays any outgoing transfer in addition to the chain charges, rather than assuming transfers vanish when local edge excess does. The exact cutoff-and-stock result above can therefore be used to improve endpoint inputs without requiring changes to those conservative tail payments.

Run `python3 -IB round2/open/age-order/hub-ranking/small-powers/continuation/check.py` from the repository root. It checks the terminal equations, the zero-power hit identity, 128 exact lower improvements and the finite tip comparisons on six representative boxes. These are independent endpoint improvements; the signed all-tree certificate is maintained in the stock-sensitive ledger.
