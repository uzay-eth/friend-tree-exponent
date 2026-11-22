VERDICT: CORRECT

# Cold Audit of `attempts/strogatz-astra/bounds/TRANSFER.md` and `transfer.cpp`

Independent cold audit of [TRANSFER.md](../attempts/strogatz-astra/bounds/TRANSFER.md), [BOUNDS.md](../attempts/strogatz-astra/bounds/BOUNDS.md), [transfer.cpp](../attempts/strogatz-astra/bounds/transfer.cpp), and [reproduce-transfer.sh](../attempts/strogatz-astra/bounds/reproduce-transfer.sh).

**Summary:** The certified bounds for every $n \ge 3$,
$$\frac{15430450607}{17611596558}\prod_{m=3}^{n-1}\!\left(1+\frac{0.563}{m}\right) \le \mathbb{E}[C_n] \le \frac{15361771189}{640204837}\prod_{m=3}^{n-1}\!\left(1+\frac{0.568}{m}\right),$$
giving $\liminf, \limsup \frac{\log \mathbb{E}[C_n]}{\log n} \in [0.563, 0.568]$ (along with the $c \le 2$ subset lower minimum $3877736170/10^{10}$ and the exact $T_3$ obstruction at $17/30$ for the upper weights), are **verified and mathematically sound**.

---

## 1. Hand Audit of Definitions, State Projection, Witness Rules, Edge Cancellation, Rounding, and Analytic Tail

### 1.1 Multiset projection, saturation, and pooling (`BOUNDS.md` + `transfer.cpp` lines 48–109, 177–211)
- Cutoffs $h_c$ are strictly checked nonincreasing (`c > 1 && cap > layers[c-1].cap` throws), with $K_1 = 20, K_2 = 7, K_c = 5$ ($c \ge 3$), also nonincreasing.
- Non-pooled types ($d < h_c$) record $(d, \min(k, K_c))$ with exact credit $1/d$. Pooled types ($d \ge h_c$) record a lower-bound mark $b \in \{h_c, \text{powers of 2 in }(h_c, H), H\}$ with credit in $[0, 1/b]$.
- In `lookup()` (lines 103–108), sorting puts pooled entries last and resets all pooled entries except possibly the last to `pool` ($h_c$), matching `enumerate()` (line 83) and preserving valid degree lower bounds.
- On central conversion (`r.convert`, lines 183–191): since both $h_c$ and $K_c$ are nonincreasing in $c$, every known degree lower bound $v.d$ projects validly to `next.label(v.d, v.k)` (and saturated $k = K_c$ remains saturated at $K_{c+1}$).
- On neighbour gain/conversion (`r.moves`, lines 192–210):
  - Saturated $k = K_c$: conversion leaves $\min(k_{\mathrm{true}}+1, K_c) = K_c$ unchanged, so skipping `!grow && k == g.k_limit` is exact; gain rate lower bound drops to $1/D$ (`else if (k != g.k_limit)` skipped in line 346) while upper bound remains $\max_{j \ge K_c}(d - j + (j-1)/2) + 1/D = d - (K_c+1)/2 + 1/D$, both exact.
  - Pooled neighbours ($t \ge \text{pool}$): neither gains nor conversions alter their mark, and degrees only grow, so skipping `r.moves` for $t \ge \text{pool}$ is exact.

### 1.2 Witness at $c = 1$ (`TRANSFER.md` lines 28–38; `transfer.cpp` lines 27–30, 54–58, 97–102, 196–207, 223–230, 334–342)
- Active iff $c = 1$, $2 \le k_w \le 4$, $d_w \le 14$, and $d_w < h_1$. Since $k_w \ge 2$, $w$ always has at least one other core neighbour $x \ne u$.
- **Indexing check:** `bridge()` offset `(d-2)*(d+1)/2 + k - 1` for $2 \le d \le 5, 1 \le k \le d$ gives offsets $0..13$, and `witness_extra = 14` matches pooled `{d, k, 6, 0}` in lines 55–57 identically.
- **Initialization at birth (`newborn()`, lines 223–230):** Triggered iff $D \le 14$ and $1 \le c \le 3$ (so post-conversion $k_w = c+1 \in \{2,3,4\}$). Picks `r.bag[0]` (the minimum-degree pre-conversion core neighbour of the creator; note `layers[3].cap >= 6` is enforced at line 212). If `r.bag[0]` is pooled or has $e \ge 6$, `w.d >= 6` maps to pooled witness `{D, c+1, 6, 0}`. If $e < 6$, then $e \le 5 \le K_c$ for $c \in \{1,2,3\}$ (including $e = j = 5$ at $c = 3, K_3 = 5$), so `(w.d, w.k)` is exact. When the creator itself had $c = 1$ with a witnessed sole neighbour, `w.d, w.k` read the creator's distance-1 neighbour (ignoring distance-2 fields), which is exact.
- **Update on $w$ conversion (`grow == false`, line 200):** `g.label(d, k+1)` points to `ids[d][k+1]`, whose first entry (offset 0) is `{d, k+1, 2, 1}` whenever $2 \le k+1 \le 4$ and $d \le 14$, deterministically switching the tracked witness to the newly converted $(2,1)$ leaf of $w$ (and discarding when $k+1 = 5$).
- **Update on $w$ gain (`grow == true`, line 200):** `g.bridge(d+1, k, wd, wk)` preserves `(wd, wk)` (including pooled `wd=6, wk=0`) for $d+1 \le 14$ and discards at $d+1 = 15$.
- **Witness events (`m.witness`, lines 203–207, 334–338):** For exact witness (`wk > 0`), $x$ is at distance 2 from $u$ in the tree (hence $x \not\sim u$), has exact gain rate in $[e-j+1/d, \, e-j+1/d+(j-1)/2]$ matching `(twice*d + 2)/(2*d)` identically, and conversion rate $(e-j)/e$. Gain from $e = 5$ pools to `{d, k, 6, 0}`. Witnessed $w$ gain rate in lines 339–342 equals $(d-k+1/D) + \{1/e \text{ or } [0,1/6]\} + [0, (k-2)/2]$, exact.
- All 6 event sets ($u$, leaves of $u$, $w$, leaves of $w$, $x$, leaves of $x$) are pairwise disjoint in any tree.

### 1.3 Directed edge transfer cancellation (`TRANSFER.md` lines 40–57; `transfer.cpp` lines 119–165, 213–220, 321–333, 354–374)
- Support `eligible(D, c, d, k)`: $1 \le c \le 4$, $2 \le k \le 6$, $d < h_c$, $D < h_k$, $k < K_c$, $c < K_k$.
- **Mutual observability:** Because $d < h_c$ and $k < K_c$, receiver row $u$ sees $w$ unpooled and unsaturated with exact $(d, k)$ (and every witness variant `{d, k, wd, wk}` gets the same `slot.in`). Because $D < h_k$ and $c < K_k$, sender row $w$ sees $u$ unpooled and unsaturated with exact $(D, c)$. Because $\min_{c\ge 1} K_c = 5 > 4 = \text{receiver\_limit}$, a saturated neighbour $w$ can never falsely match an outgoing transfer either, and `certify()` lines 357–373 explicitly verify `slot.in` and `slot.out` against `transfer` across every `(c, D, t)` including all witness variants, asserting `0` whenever `t >= pool || w.k == k_limit`.
- **Multiplicity & exact cancellation:** At receiver $u$, `v -= edges[m.type].in` is multiplied by `m.count * r_w` (line 347); at sender $w$, `for (int t : r.bag) diff += edges[t].out` sums over `r.bag` with multiplicity before multiplying by $w$'s central gain rate $r_w$. Before taking rate extrema, every directed core edge $w \to u$ contributes $-r_w g(w,u) + r_w g(w,u) = 0$ identically; any edge incident to $D \ge H$ has transfer $0$ because $h_c \le h_1 < H$.

### 1.4 Directed integer rounding, overflow safety, and analytic tail (`transfer.cpp` lines 302–352, 376–457)
- `divide(__int128 a, Int b, bool ceiling)` correctly implements exact floor/ceiling for all signed 128-bit $a$ and $b > 0$ (checked all 4 sign/rounding cases).
- `coeff(x, a, b)` rounds UP for upper and DOWN for lower on every drift term; `compare` rounds DOWN for upper and UP for lower; boundary `up` at $D = H-1$ rounds UP for upper and DOWN for lower (sound because central gain rate $L + \rho \ge 0$).
- Magnitude check: $|g| \le 10^{10}$, `weight_limit = 16 * 10^10` (note reported finite max $\approx 1.76 \times 10^{10} \ll 1.6 \times 10^{11}$), `|diff| <= 5.5e12 << 2^63`, all intermediate products/accumulators use `__int128` with checked narrowing, and builds run with `-fsanitize=undefined -fno-sanitize-recover=all`.
- **Analytic tail ($D \ge H$, zero incident transfers):** With $q = L/D$ and $0 \le L+\rho \le D$, $-q/H \le -q/D \le Gq = (L+\rho)(1-q)/(D+1) - q/D \le 1-q$.
  - Lower $f = Bq$: $\mathcal{A}f \ge q(s - B/H) \ge \delta Bq \iff B(\delta H + 1) \le sH$, checked exactly in `__int128` against `actual.back()` (pooled-$H$ newborn weight) at lines 429–431.
  - Upper $f = A(1+\lambda q)$: $\mathcal{A}f \le A\lambda(1-q) + qs \le \lambda A(1+\lambda q) \iff s \le A\lambda(1+\lambda)$, checked exactly in `__int128` at lines 429–431.
  - Global bounds `norm = upper ? min(min_z, tail_coeff) : max(max_z, tail_coeff)` and $c \le 2$ subset tail minimum ($q \ge (H-2)/H = 206/208$) are exact.

---

## 2. Independent Reproduction and Experimental Tests

### 2.1 Live UBSan reproduction of `./reproduce-transfer.sh --screen-17`
Executed `./reproduce-transfer.sh --screen-17` (compiled with `g++ -std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all`, zero compiler warnings, zero UBSan errors, exit code `0`). Every single output line matched `transfer.log` and `TRANSFER.md` bit-for-bit:
- **Lower 563/1000:** `9,639,742` finite states, `6,445` ordered pairs, minimum newborn consistency `6,104,023`, minimum finite slack `0`, global weights `[0, 17611596558]`, $c \le 2$ subset minimum `3877736170/10^10`, `PASS lower 563/1000 product prefactor 15430450607/17611596558`.
- **Upper 568/1000:** `9,639,742` finite states, `6,445` ordered pairs, minimum newborn consistency `1,714,162`, `T3 exact local drift / weight 8723588643/15361771189; 30*drift - 17*weight = 557549077`, minimum finite slack `0`, global weights `[640204837, 17605633804]`, `PASS upper 568/1000 product prefactor 15361771189/640204837`.
- **17/30 screen:** `newborn ratios [1.00248783668, 1.00248783668]`, `screen failed` (exit code 1 as expected).

Independent Python arithmetic recheck confirmed:
- Exact eligible ordered pair count for `H=208, caps=[128,26,10,6,5,4,4,3,3,3]` is `6445`.
- $30 \cdot 8723588643 - 17 \cdot 15361771189 = 557549077 > 0$ (ratio $0.5678765 > 17/30 = 0.5666667$).
- Lower tail subset minimum $17611596558 \cdot 206 / (208 \cdot 10^{10}) = 1.74423 > 0.387773617$.

### 2.2 Independent exact-rational brute-force tree & history audit (Python `fractions.Fraction`)
Wrote a standalone from-spec Python implementation of the full history-dependent state machine (`pool_lb`, witness tracking, witness swap on $w$ conversion, witness pooling at $e=6$, creator witness init for $c \in \{1,2,3\}$, saturated $k=K_c$, and signed rational directed edge transfers $g(D,c,d,k)$) on `H=20, caps=[16,10,6,5,4,3]` (420 nonzero random rational edge transfers, random rational $f$ with affine tail):
- **Exhaustive small histories ($n = 3..7$, all 436 histories):** Exact rational equality between brute-force $n\mathbb{E}[\Phi(T_{n+1})-\Phi(T_n)\mid\mathcal{H}_n]$ and summed amended row drifts, exact zero net transfer cancellation, and exact rate-interval containment on **all 436/436 histories**.
- **Random & adversarial corner-case trees (48 trees, $n \le 90$):** Explicitly engineered trees exercising $c=3$ creator with exact $(e=5,j=5)$ witness, all-pooled creator converting a leaf, $d_w = 14 \to 15$ witness discard, $k_w = 4 \to 5$ witness discard, $e_x = 5 \to 6$ witness pooling, saturated $k=K_2=7$ gaining/converting, multiplicity $\ge 2$ identical neighbours with nonzero transfers, and boundary $D = H-1 \to H \to H+2$. **All 48/48 passed identically.**

### 2.3 Independent exact-Fraction audit of actual `transfer.cpp` dumps (small geometry `H=24, caps=[18,10,6,5,4,3]`)
Dumped complete integer tables, transfer dictionaries, and tail coefficients from `transfer.cpp` for `lower 545/1000` and `upper 590/1000` (`100,303` finite states and `456` ordered edge pairs each) and checked them with a standalone exact `fractions.Fraction` verifier (zero rounding):
- **Lower `545/1000`:** `100,303/100,303` finite states PASS (`min_exact_slack = 0.0000`), `456` transfers eligible and mirrored, analytic tail PASS, and actual brute-force tree drift $n\mathbb{E}[\Delta\Phi] \ge 0.545\Phi$ PASS on all `484` histories + adversarial trees.
- **Upper `590/1000`:** `100,303/100,303` finite states PASS (`min_exact_slack = 0.2000`), `456` transfers eligible and mirrored, analytic tail PASS, and actual brute-force tree drift $n\mathbb{E}[\Delta\Phi] \le 0.590\Phi$ PASS on all `484` histories + adversarial trees.

### 2.4 Negative controls (all caught)
1. **Wrong witness gain rate (dropping $1/d_w$):** CAUGHT by exact rate-interval check (`r_x out of bounds: 2 <= 29/14 <= 2`).
2. **Wrong saturated lower rate (using $d - K_c + 1/D$ instead of $1/D$):** CAUGHT on a tree with $d_w = 9, k_{\text{true}} = 9 > K_2 = 7$ (`r_w out of bounds: 7/3 <= 11/15 <= 16/3`).
3. **Asymmetric edge transfer mismatch (tree audit):** CAUGHT (`Edge transfer failed to cancel: net = 4`).
4. **Illegal transfer on saturated $k = K_c$:** CAUGHT (`Edge transfer failed to cancel: net = -11/45`).
5. **Omitting witness swap to $(2,1)$ on $w$ conversion:** CAUGHT (`good=('FIN', 2, 1, (('E', 4, 3, 2, 1),))` vs `bad=('FIN', 2, 1, (('E', 4, 3, 3, 1),))`).
6. **C++ `certify()` slot corruption (`+1` to one `slot[t].out` in a temp copy):** CAUGHT immediately with `std::runtime_error: edge lookup mismatch`.
7. **Stricter exponents & asymmetric transfer perturbations on dumped `H=24` tables:** Both lower (`560/1000`, slack `-2.02e8`; asymmetric transfer, slack `-2.00e6`) and upper (`575/1000`, slack `-7.77e7`; asymmetric transfer, slack `-2.00e6`) CAUGHT by the independent exact-Fraction verifier.
