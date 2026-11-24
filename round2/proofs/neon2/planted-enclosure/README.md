# Planted-Tree Bellman Quotient & 8-Digit Barrier Diagnostics (Track 1B)

**Pinned claim commit:** `288355468d6f21e8c39b60d34e706c9164f8f68f`.  
**Review and dependencies:** No separate rival verdict for these diagnostics. Uses Astra Claim 14’s zero-gauge equation; numerical comparisons use the quintet enclosure.

**Author:** Neon-2 lane (Gemini 4 Neon and GPT-6.2 Sol Ultrafast).

This directory documents our quantitative probe of whether **unmarked planted-tree Bellman enclosures** on Astra Claim 14's zero-gauge operator equation $`(s - L^0 - C\Delta) w = b`$ can beat free-tree local motif potentials or reach 8 certified digits for $`\mu`$.

## Summary of Measured Enclosures

All three scripts below are diagnostic floating-point probes (not certified bounds) measuring the exact Bellman-quotient gap when lumping unmarked planted trees outside a finite core into coarse root summaries $`(D, c)`$:

1. **Exact $`|S| \le N_{\mathrm{ex}}`$ core + Depth-0 coarse tail (`test_depth2.cpp`):**
   - Keeps every unmarked planted tree of size $`|S| \le N_{\mathrm{ex}}`$ exact (zero internal uncertainty; $`N_{\mathrm{ex}} = 12`$ gives $`K = 4{,}636`$ exact trees), collapsing exits $`|S| > N_{\mathrm{ex}}`$ to coarse $`(D, c)`$ states with worst-case external child-rate / $`\Delta w`$ bounds.
   - Measured enclosure at $`N_{\mathrm{ex}} = 12`$: $`\mu \in [0.5487, 0.5900]`$ (width $`\approx 4.13 \times 10^{-2}`$).
2. **Depth-1 root+children quotient (`bellman_proto.cpp`):**
   - Resolves root $`(D, c)`$ together with the multiset of coarse child states $`(D_j, c_j)`$ for $`D \le 8, c \le 4`$ ($`12{,}778`$ states), bounding grandchild effects via per-fiber worst-case child jumps.
   - Measured enclosure: $`\mu \in [0.5550, 0.5763]`$ (width $`\approx 2.13 \times 10^{-2}`$).
3. **Closed 1-body planted LP with 1-way child$`\to`$parent telescoping transfer (`test_transfer_clean.py`):**
   - Adds a sign-constrained telescoping transfer $`A(S) - \sum_Y A(Y)`$ ($`A \le 0`$ for upper $`\mu`$, $`A \ge 0`$ for lower $`\mu`$) with a leak-free analytic tail for $`D > D_{\max} = 16`$.
   - Measured enclosure: $`\mu \in [0.4573, 0.5983]`$ (width $`\approx 1.41 \times 10^{-1}`$).

## Structural Conclusion on the 8-Digit Barrier

- Pure planted-tree quotient closures are **$`100\times`$–$`1000\times`$ looser** than unoriented free-tree motif potentials (`round2/proofs/neon2/quintet-cert/`, width $`< 5 \times 10^{-5}`$): without bidirectional cross-edge cancellation (paired transfers across every edge of the core tree), every worst-case sibling/grandchild perturbation is double-counted independently at each edge, and raw size cutoffs decay only polynomially ($`N^{-\mu}`$) while planted tree counts grow exponentially ($`\approx 3.7^N`$).
- Conversely, finite-diameter local potentials (`quintet-cert/`) achieve geometric improvement per diameter step (diameter 1: $`\sim 1.6\times 10^{-2}`$; diameter 2: $`\sim 1.4\times 10^{-3}`$; diameter 3: $`\sim 1.0\times 10^{-4}`$; diameter 4 / 5-vertex: $`\sim 4.5\times 10^{-5}`$), so pushing to $`10^{-8}`$ (8 digits) via local potentials alone would require diameter $`\ge 7`$–$`8`$ (thousands of tree shapes and $`> 10^{11}`$ states) or a genuinely new global spectral/contraction representation of the infinite-dimensional Cox interface beyond finite-diameter local sums.
