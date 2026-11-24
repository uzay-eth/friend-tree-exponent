> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Signed continuation on a realized bare P4

**Review candidate:** the broader entry calculation gives a strictly positive, disjoint upper-witness charge of at least `8.65873655175e−11` per core. With the same potential ceiling as the arms calculation, that would move its upper endpoint inward by at least `5.70115391669e−11`, to outward `0.565772451566988461`. This extends the earlier eleven-event calculation on the **same** objective; do not add the two. It remains outside the combined ledger pending independent review. This is a local paired sign, not a sign for the entire D18 difference.

## Physical face and Poisson correction

Fix the common multiplier `t=1.2415`. At the focal nonfounder, select the top-owned `P4` row with all four participating vertices bare with exactly one core child and degree two. There is exactly one such path there. The fifth vertex is the live arm rooted at the fourth vertex's core child. Let `p` be the actual reciprocal degree of the focal vertex's parent and `d,k` the degree and arity of the arm root. The exact contrast of the two signed charged rows is

```
q = delta−_P4 − t delta+_P4
  = A + B p + C/d
A = 287589315471350942439 / (5*10^22)
B = −.005553926199875
C = −.02262119047951.
```

The two exterior coefficients and constant are paired **before** optimizing. The script checks this identity on concrete parent degrees 2 and 7 and arm degrees 2, 3 and 8, and compares the complete physical drift on their corresponding whole trees. The arm and exterior need not be independent. The individual charged rows are nonnegative by the production D18 verification.

The face persists until a gain at any of the four protected vertices. This cannot recur for the same focal vertex once it exits because its degrees along the selected four-path are monotone and a formerly bare degree-two vertex never returns to degree two after a gain. The exact exit intensity is `3+1/d+p`: six directed internal path sends, plus the actual two outer sends. The arm is *not* killed on an arm event. For arm leaf stock `l=d-k-1`, it instead gains at `g∈[l+1/2,l+(k+1)/2]` and promotes at `b=l/d`. The sender from the lowest protected vertex has rate precisely `1/2` into the arm; sends out of the arm use its actual `1/d`, already included in the face exit rate. Omitted descendants only affect the adaptive gain rate of the tracked root.

For every `s∈[.565674722414,.565773666188]`, construct a bounded signed supersolution `U(d,k)` with

```
(s + 3 + 1/d + p + g + b) U(d,k)
    >= A+B p+C/d + g U(d+1,k) + b U(d,k+1).
```

Check the endpoints of `s,p,g`; the cross-multiplied inequality is affine. Above the finite degree table use the **positive nonzero** continuation `U_tail≥A/(LO+3)`. Promotions and gains then leave the constant unchanged, so this comparison covers arbitrary physical stock and arity. For degree ≤12 the table gives `U(2,0)<−.0010542375` at virgin arm and `U(2,1)<−.0011412348` at bare arm. All endpoint tests are exact; the tail deliberately pays the wrong-sign part. Along a physical discounted stopped trajectory the inequality gives the coboundary correction `q+(A_killed−s)U≤0`. Its integral uses only the initial `U`; the terminal at face death is zero and the bounded infinite-time term vanishes. A positive value of `−U/t` at entry is thus a valid lower price for **the upper witness's actual lifetime charge**: `−q/t=delta+_P4−delta−_P4/t ≤ delta+_P4`. This is one signed continuation, without two separately charged absolute histories or an assumed sign on each live row.

## Higher-weight entry and disjointness

From the focal newborn, four consecutive tip promotions create the four bare vertices and a live virgin arm. The discounted selected prefix has lower mass `Π_(n=2,…,5) [(.5)/(HI+n)]`; there is just one tip leaf after each successful promotion. Until the first gain at a protected vertex, admit **every** physical arm transition until the arm contains nine ordinary vertices. There are 286 geometric exits, of which 211 have positive `max(0,−U/t)`; admitted exits include descendant degree 8. The parent input to the protected crown root is allowed to vary arbitrarily up to `1/2`. The arm has exact adjacency transition intensities, including incoming `1/2` from the bare vertex above it. At arm size `n` they sum to `n−1/d+1/2`. Add the protected-vertex kill rate `3+1/d+p`; using `p≤1/2` gives the common discounted denominator `HI+n+4`. Its use here is valid because these pre-entry values and payoffs are nonnegative. The stopped Poisson value at the virgin is `5.69403029778807e−7`. Exact forward discounted flux and a backwards rounded solve agree within the explicit rounding allowance. Allowing the arm to reach size ten gives a different, overlapping and weaker bound; it is a comparison, not an additional credit.

The focal crown at payment has at least `4+9=13` ordinary vertices and at least five cores, throughout the remaining episode. It is therefore disjoint from the `≤12` higher ledger and the `>12, cores≤4` suffix of the **same owner**. Nor does it overlap first-body charges. The selected segment is dominated by the nonnegative stationary observable on this physical face and size/arity suffix, even though selecting an entry uses path history. Each focal newborn's charge is paid once. Bounded local truncations of the nonnegative charge justify the stationary transfer if needed; multiplication by `LO` comes from the true stable-age identity, with discount already covered by the robust comparison. Divide the new numerator by the same ceiling `1.518769126087491` and subtract from arms-only upper endpoint `.565772451624`. No additivity with an unreconciled competing combined ledger is asserted.

## Replay

From the Astra checkout:

```sh
taskset -c 13 nice -n 10 /usr/bin/time -f 'cpu=%U+%S wall=%e' timeout 12s \
  ./round2/.venv/bin/python -Bu \
  round2/occupation-neon/ratio-transfer/next/entry.py b3f18_lo.bin b3f18_hi.bin
```

Passed, `0.31+0.02` CPU seconds, `0.33` wall. It prints `pi delta+ lower 8.658736551752106e-11` and `outward endpoint ticks at 1e-18 565772451566988461` for arm threshold 9. `spine.py` contains the earlier fixed-path control, not needed for this broader entry bound. Both use the pinned certificate payloads without rehashing them. The forward/backward check is an independent arithmetic path within this replay; independent reviewer sign-off has not been received.
