# Cold audit: all-size hubhood discharge

**Holds in scope:** all old trees of maximum degree nine, all stars (including the seed edge), and the stated class with cores of degree at most seven or terminal cores of degree at least ten. This is conditional hubhood, not an all-threshold statement. Independently, the same scalar inputs certify the bare-skeleton class below.

## Pins and the physical charge

All pins are Git blob IDs in `round2/open/age-order/`:

| Input | README | Checker |
|---|---|---|
| `global-transport/hostile/discharge/` | `e76a611435ec1244e1601ba401589a7753cd7906` | `feccb888c61fbc58cb8fd93a65c5d7ed422e9575` |
| `global-transport/hostile/` | `4602cfd0688772210e90a562f01e05e7726f49c3` | `b5b63d953c0ae51a6d145fa3c00452c4493c1b95` |
| `preplacement-transfer/` | `2d6eae9e8f0fcb157ee270d94a3a8023cb3493bd` | `fa54080065cc4157f7043ab78f5a9448baccb202` |
| `early/` | `22e27fc0c7e4a4754217b66be5a184ee95282c93` | `5ca7998a11fc514efcd1e90e859a42d8b1bf2b8c` |
| `adjacent/` | `b25578d819568143f158c9c3411e0ef846c94d4a` | — |

The tracked stock rescue calls the return race only after the last tracked leaf's **real promotion**, when that neighbour is degree two. This is not a hub lower bound for an arbitrary bare core. In a serial placement the middle starts with one ordinary leaf and an old core parent. The fresh first-arrow credit is a race with the real incoming old-neighbour send and uses the old parent's **postbirth** degree. The baseline gain and loss payoffs balance via the rescue recursion; an old-parent degree gain has baseline value, and the pessimistic total competing intensity in the input safely bounds this race. The favourable parallel-tip lower bound instead uses a still-degree-one tag and its return immigration into its parent. The frozen send upper `1/(j+1)` requires a permanent nonleaf slot beside the tip in its parent; for an all-leaf star parent use the separate `B=1` barrier (the bound is still `1/(j+1)` when retaining the tip's rate-one return). The leaf degree, other-leaf count, and sender slot are all read in the **common final tree**.

If an old target has weight `w=r/(d+1)` and tick coefficient before rounding `s=w/d`, serial credit is `w*g_(d+1)`. On a parallel edge with `s_u>s_v`, only the new tip at `u` is adverse. For a parallel old-core/old-leaf pair, the leaf has itself been promoted by the other placement: the new tip rooted there has parent degree two with no other leaves, and the new tip at the old core has `l−1` spare leaves. On old-core/old-core edges neither old parent spends a leaf on the other placement: counts are `l,j`. The adverse coefficient depends on one signed difference; one must not charge both directions separately. The degree cap is on the **input** tree only, so e.g. a bare nine-core's attached tip uses `U(10,0)`. The checker explicitly creates final parallel trees for its leaf charges and evaluates rounded messages from alternate physical roots.

## Rounded recursion and bases

After all pendant leaves have been paid, orient the core tree at any core in cap mode. `C(b)+Σ M(b,e_i)−D(a,b)` subtracts the entering edge once; its sum excludes just one parent-degree slot, not every occurrence of that degree. The induction is valid even for a negative message. At each fixed argument the terminal reference child has exactly the incoming slot and hence zero onward messages. No claim that it minimizes the physical gap is required. Root sums are computed with all slot multiplicities, including the `c=d` bare profiles.

The integer parent coefficient is floored, and either signed difference between two original coefficients is at most the corresponding tick difference plus one tick. Both slopes are rounded upwards and the raw maximum is rounded upwards again to pay the debit. For cap mode's accelerated test, if the child tick is `y`, the positive affine directions switch within `[y−1,y+1]`; all integer ticks there are tested. Outside these neighbourhoods the difference of the two raw maxima is affine. The maxima of the difference are therefore controlled by the interval ends and the three candidate ticks per child. Subtracting the terminal debit is legitimate since `ceil(A/T)−ceil(B/T)≤ceil((A−B)/T)`. Bounding the difference this way is essential: subtracting an upper bound on a reference debit as though it were a lower bound would not be safe.

For stocked mode use the two endpoint *sums* of all reciprocals from high neighbours. Individual unknown reciprocals need not agree: the profile sees only their sum. The credit is concave in the sum, while the debit is convex in the two coefficient arguments, so interval endpoints cover all allowed interiors. The high terminal symbol is valid as a child; root at a small core if one exists. For a large terminal the pendant term is nonadverse since `s>s_l` and `U≤1/(d−1)<L(2)`; its sole core-edge adverse tip has `d−1` other leaves, so `U≤1/d`. If no small core exists the core tree (other than a single-core star) has exactly two vertices, with debit at most `1/110` as only one side is adverse. The pin explicitly tests stars of old centre degrees two and three. For `d≥4`, centre-side parallel stock is `d−1` other leaves: one old leaf was promoted and the new tip itself is excluded. Thus the cap is `1/d` and the debit vanishes. The edge base case uses the independent finite-clock pullback, checked at horizon four with lower gap `0.022645238637`.

## Extension: a connected high-degree bare skeleton

Let a nonempty connected set of `B` old cores be **bare**, each with degree at least ten, and let every other old core be terminal (one core neighbour), with arbitrary degree `e≥2`. The terminal cores are attached directly to this skeleton. Each bare core has no pendant debit and has nonnegative serial credit, which we may drop. On a bare–bare edge `s≤1/[2(d+1)]≤1/22`, because every sender has degree at least two. Since `U≤1`, the debit is at most `1/22` per such edge.

On a terminal edge, if the bare side is adverse its debit is actually at most `1/22` by the preceding bound on its coefficient; for the finer bound needed at small `e`, check the sign: the terminal has `s=(e−1+1/d)/[e(e+1)]≥(e−1)/[e(e+1)]≥1/22` for every `2≤e≤9`. Hence **only the terminal can be adverse** in that range. Its debit is then at most `s U(e+1,e−1)`. The rational endpoint controls give `C−s U>3/500` both at `1/d=0` and `1/d=1/10`; concavity (subtracting this linear expression) covers the continuum, and the same endpoint controls also give `C−1/110>3/500` as surplus checks. For `e≥10`, the terminal lower credit is `(9/11)·0.16694235+(9/20)g_2` and its side can cost at most `1/110`. If instead the bare side is adverse, debit at most `1/22` suffices against this credit. The ratio and phase of high degrees on the bare skeleton never enter this bound.

A tree with `B` bare vertices has `B−1` bare–bare edges and at least `10B−2(B−1)=8B+2` terminals by the incidence count. Therefore the charged numerator exceeds `(8B+2)(3/500)−(B−1)/22>0`, since `8·3/500>1/22`. This covers arbitrarily many adjacent crowded bare cores of unbounded degree, beyond both audited profile classes; no mixed positive-stock, high-arity internal class is asserted.

## Replay

From `[local path]`, `nice -n 10 timeout 55s python3 -IB round2/open/age-order/global-transport/hostile/discharge/check.py 9` passed 48,602 profiles and 194,472 messages, with minimum root tick 27,404,273. The same script with `--stocked 7` under a 35-second timeout passed 6,426 profiles and 21,020 messages, minimum 22,124,577. The direct `--pairs 6` replay passed 912 profiles and 2,305 messages, minimum 27,404,273. No cap-ten general sweep is claimed.

```text
nice -n 10 timeout 25s python3 -IB round2/checks/preplacement-discharge-sol.py
PASS 37 cap-nine rerootings; 432 adverse debits and 59 bare child checks
PASS five stars, sixteen terminal endpoint bounds and 23 additional bare-skeleton rerootings
```

The physical fixtures test the rounded induction on concrete trees; they do not replace the exhaustive profile sweeps. The scalar kernels are deliberately pinned imports rather than an independent computation of the true hub probabilities.
