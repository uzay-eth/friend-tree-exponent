> Produced by GPT-6.2 Astra.

# Transferring an observed neighbour's leaf-gain term

This separate certificate uses the process, degree cutoffs, pooled marks and hub potentials in [the baseline derivation](bounds-hub.md). It adds one conditional neighbour descriptor and redistributes terms of the exact drift along observed core edges. The baseline source and reproduction are frozen separately.

For every $`n\ge3`$ the resulting bounds on the mean core size are

```
(15430450607/17611596558) product_{m=3}^{n-1} (1+0.563/m) <= E[C_n]
E[C_n] <= (15361771189/640204837) product_{m=3}^{n-1} (1+0.568/m).
```

Thus the liminf and limsup of the mean log-growth exponent lie in **[0.563, 0.568]**. The previous width of 0.043 is reduced to 0.005.

Run `./reproduce-transfer.sh` from this directory, or add `--screen-17` to rebuild the proposals at 17/30 as well. Each side verifies 9,639,742 finite states, 6,445 ordered endpoint pairs and every degree at least 208 by the analytic tail. Both certificates use neighbour cutoffs 128, 26, 10, 6, 5, 4, 4, 3, 3, 3, then 2 thereafter. The script uses a single low-priority CPU and the undefined-behaviour sanitizer. See [the archived audit’s replay results](../reviews/review-bounds-transfer-gemini.md#21-live-ubsan-reproduction-of-reproduce-transfersh---screen-17).

| Exact weight data | Lower | Upper |
|---|---:|---:|
| Root weight, in units of $`10^{-10}`$ | 15430450607 | 15361771189 |
| Finite minimum, same units | 0 | 640204837 |
| Finite maximum, same units | 17611596552 | 17605547966 |
| Tail weight, same units | $`17611596558q`$ | $`11228082783(1+0.568q)`$ |
| Global lower and upper bounds, same units | [0, 17611596558] | [640204837, 17605633804] |

On states with $`c\le2`$ and $`L\ge1`$, together with $`c=2,L=0`$ and a tracked $`(2,2)`$ neighbour, the exact **finite lower weight minimum** is $`3877736170/10^{10}`$. The cutoff $`h_2=26`$ ensures both degree-two neighbour types are tracked. In the infinite tail on this set, $`q\ge206/208`$, so the same positive minimum applies to the full table.

The upper table itself has an exact obstruction at 17/30. At $`T_3`$, let $`f_0`$ be the root weight, $`f_S`$ the weight at the centre of the four-vertex star, and $`f_P`$ the weight of either core vertex after the root converts a leaf. The exact scaled drift is $`2(f_S-f_0)+(2f_P-f_0)`$. There is only one core at $`T_3`$, so any redistribution across core edges leaves this value unchanged. Its ratio to $`f_0`$ is $`8723588643/15361771189`$. In particular $`30\mathcal A f-17 f_0=557549077/10^{10}>0`$. A fresh 17/30 rebuild on the same geometry gives a newborn ratio of 1.00248783668, above the upper screening threshold of one. The exact starting-tree obstruction rules out certifying 17/30 with these same weights, while the model's true exponent remains inside the certified interval if it exists.

## A witness beyond a single core neighbour

When a central vertex $`u`$ has just one core neighbour $`w`$, track one other core neighbour $`x`$ of $`w`$ whenever $`2\le k_w\le4`$ and $`d_w\le14`$ and the degree of $`w`$ is below the central observation cutoff. Record the exact pair $`(e,j)`$ of degree and core degree of $`x`$ for $`e<6`$, and a single pooled label for $`e\ge6`$.

* At the birth of $`u`$, the old core neighbours of the creator $`w`$ are eligible witnesses. Take the first in its pre-conversion multiset, sorted by degree and then label. Ties between indistinguishable neighbours can be resolved by vertex ID. The relevant pre-conversion core degrees of the creator are 1, 2 and 3. Their observation cutoffs are at least 6. For degrees below 6 the capped core degree is exact, including the case of degree and cap both 5.
* Whenever $`w`$ converts a leaf, use that new degree-two, core-degree-one vertex as the next witness if the new state calls for one. A gain of $`w`$ preserves the witness if the new state calls for one. If the degree or core-degree range is exceeded, or the central vertex acquires a second core neighbour, discard the witness.
* The witness itself gains a leaf at rate in $`[e-j+1/d_w,e-j+1/d_w+(j-1)/2]`$ and converts at exact rate $`(e-j)/e`$. These events update its label. A pooled witness keeps its label. While a witness is recorded, the exact gain rate of $`w`$ is
  `d_w - k_w + 1/D_u + 1/e + sum of the other k_w-2 reciprocal degrees`.
  The last sum is between 0 and $`(k_w-2)/2`$. The pooled witness supplies $`1/e\in[0,1/6]`$. The conversion rate of $`w`$ remains $`(d_w-k_w)/d_w`$.

Since degrees and core degrees never decrease, a sole neighbour enters the witness class only at initialization or a conversion. The latter always supplies a fresh witness. This explains why the descriptor can be updated using its recorded data even after multiple gains or conversions. The table weight may depend on history, as already happens with the greatest pooled mark.

## An exact identity on each ordered edge

Let $`r_w=L_w+\rho_w`$ and let $`b_{u,w}`$ be the change of $`u`$'s own weight when its core neighbour $`w`$ gains a leaf. For each directed edge $`w\to u`$ choose a signed rational $`g(w,u)`$. Subtract $`g(w,u)`$ from $`b_{u,w}`$ in row $`u`$. Add the same $`g(w,u)`$ to the coefficient of $`w`$'s own gain in row $`w`$. For every tree and every history the changes sum to
$`-r_w g(w,u)+r_w g(w,u)=0`$.
This cancellation uses the actual conditional rate before taking any extrema. The conditional expected change of the total weight is therefore the sum of the amended row drifts. The verifier selects rate endpoints using the signs of the amended coefficients. The central gain rate multiplies the entire amended central coefficient. Witness gains at distance two and conversions retain their original coefficients.

The transfer depends only on the true degree and core degree of the ordered endpoints. Write $`(D,c)`$ for the receiving row $`u`$ and $`(d,k)`$ for its neighbour $`w`$. Its possible nonzero support satisfies

```
1 <= c <= 4,  2 <= k <= 6,
d < h_c,  D < h_k,  k < K_c,  c < K_k.
```

Both rows can therefore identify the ordered pair exactly. Pooled neighbours and ambiguous capped core degrees cannot hide a nonzero transfer. Also any edge incident to a central-degree tail has transfer zero, so the previous analytic tail proof applies verbatim. The verifier checks the integer transfer in every in/out lookup slot against the same dictionary of ordered pairs and enforces absolute transfer weight at most 1.

The proposal for $`g`$ averages two differences for a neighbour gaining a leaf. In one representative state the remaining neighbours are newborn pendants. In the other they have pooled marks. Witnessed types use a pooled witness in these representatives. Backward sweeps damp the proposed transfers, quantize them to signed multiples of $`10^{-10}`$, and then freeze them while iterating the newborn values. The verified inequalities require only the final integer transfers, their shared physical interpretation and bounded size.

The baseline integer bounds also cover the amended checker. The sum added to the central coefficient is at most $`512\cdot10^{10}`$ in integer units. Multiplying by a central gain rate at most 512 adds less than $`2.7\cdot10^{15}`$ to the drift magnitude. The subtractions from all neighbour coefficients together add at most $`512^2\cdot10^{10}`$ more. These terms fit within the existing 64-bit residual and repair margin. Rate products, accumulated drift and tail comparisons still use signed 128-bit intermediates.
