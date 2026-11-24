# Occupation of physical donors in the four-vertex witness

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 9: Claude and HY4 Flash on specified components; GPT for full replay; inputs include the Round 1 physical exponent and prior, signed raw witnesses and disjoint occupation credits. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


This calculation proposes the interval
```math
 0.5656523639\ \le\ \mu\ \le\ 0.5658045676 .
```
The two endpoints rely on the full four-vertex row test and the occupation calculation below. Cross-family review of this numerical result is pending.

## 1. Inputs from the row test

Use the continuous-time physical tree started at its first core, the older-to-younger orientation of core edges, and $`l=113/200\le\mu\le u=1133/2000`$. Let $`C`$ be its core count and $`B=\mathcal LC`$ its total promotion intensity. The unit for every value in the tables is $`1/S`$, where $`S=2^{40}`$. In this proof all $`f`$, $`\Phi`$, $`T`$, $`Y`$, $`e`$, $`v`$ and $`y`$ are in **integer table units**, with the rate of time unchanged.

Write $`\sigma=1`$ for upper and $`-1`$ for lower. The signed row inequality from the [four-vertex witness](../motif-cert/proof.md) and its [higher rows](../motif-high/higher-rows.md) is
```math
 \sigma(\lambda\Phi-\mathcal L\Phi)\ \ge\ F_4+F_1-e C-E.             \qquad\text{(1)}
```
The error near the initial centre satisfies $`[|E|]_s\to0`$ in the normalization of the [bulk drift test](physical-test.md). One can discard all rows, and all unmatched transfers, that meet the founder or one of a fixed number of its generations. Their motifs have bounded size. Signed weights, newborn jumps and transfer centers are bounded. At a high-degree endpoint, leaf gains multiplied by their jump are bounded since the endpoint interpolation is affine in $`1/D`$. At the singleton, the analytic tail is a polynomial in $`L/D,1/D`$ with bounded first derivatives on the unit square, and each coordinate changes by at most $`1/D`$ on gaining a leaf. Thus the absolute remainders grow polynomially in the core counts within a fixed distance of the founder. Section 1 of the bulk drift test bounds their expected values by polynomials in time. The raw lower comparison permits the addition of a fixed multiple of $`\sum(k_x-1)`$ to its pointwise lower estimate; this sum over the oriented full core tree is $`-1`$. Its normalized contribution vanishes.

The quantities $`F_4,F_1`$ are sums of nonnegative real terms. For an occurrence of a four-core potential $`H`$, let $`x`$ be one of its vertices and $`z`$ any adjacent core outside the four-set. Its contribution to $`F_4`$ is
```math
 J(d,i)=\tfrac12\max(d,0)-i d,\qquad d=\sigma[H(x^+)-H(x)],\qquad i=1/D_z.   \qquad\text{(2)}
```
Here $`x^+`$ means one additional leaf at $`x`$. We only use entire four-sets with ordinary noninitial vertices, degrees at most six and an ordinary noninitial donor. For $`F_1`$ the identical formula uses the singleton potential $`f(x)`$ at $`72\le D_x<1399`$ and its ordinary neighbouring core $`z`$. Each pair is directed towards $`x`$. Different four-sets, donor slots, singleton incidences and directions correspond to different terms in (1). The signed residual cap $`e`$ was computed **before** retaining either source of (2). Credit (2) is the physical gap between maximizing a donor input on $`[0,1/2]`$ and using its actual value. It does not require a strict row inequality at a corner.

## 2. Following a connected set through time

Use all connected subsets of between one and five noninitial cores. Order their members by promotion time, oldest first. For a connected subset the oldest is the unique vertex nearest the founder. Every other member has its parent in the subset. Thus its shape is a list of parent indices $`p_j<j`$. The current state of an observed vertex is $`(D,c)`$, with $`c`$ its total number of core neighbours, including its exterior parent if that parent lies outside the set. Let $`V`$ be a nonnegative integer table on these ordered subsets, and let $`T`$ be the sum of their values on the physical tree.

A gain at member $`j`$ changes just $`D_j`$ to $`D_j+1`$. A promotion at $`j`$ changes $`c_j`$ to $`c_j+1`$ and gives *two* terms: the new value of the old set, and the value of the enlarged set obtained by adding the newborn in state $`(2,1)`$ to this post-promotion configuration. Drop the latter term when at the maximum set size. Every newly formed nonsingleton set arises once this way: remove its youngest member to recover a unique old connected set. Every promotion on the whole tree also creates one singleton of value $`v=V((2,1))`$. This includes a promotion by the founder.

For a positive observation of size $`r`$, the code uses the degree and core-degree ceilings in the two columns below. The same ceilings apply at each member of that layer. An out-of-range member is recorded by a single pooled state if it has at most one selected child. An out-of-range member with two or more selected children makes the score zero. Write $`d_0=1+\min(D_{\max},c_{\max})`$ for the pooled state's lower degree. Degrees and core degrees are nondecreasing, so the pooled state is permanent within the old subset. Dropping the possible new-set score at an unobserved promotion has the favourable sign. An enlarged set uses weakly smaller cutoffs at the existing positions.

| Set size | Degree ceiling | Core-degree ceiling |
|---:|---:|---:|
|1|192|30|
|2|48|14|
|3|16|9|
|4|9|6|
|5|6|5|

Only size five receives a source in the four-vertex calculation. For each leaf of the *undirected* five-set, delete it to obtain the prospective four-set and its unique contact. This includes deleting the oldest member if it is an undirected leaf. Read the remaining quartet potential in its physical orientation. Set its source to zero if any quartet member is pooled. For each exact donor use the floor of (2); for a pooled donor use the floor of the minimum on $`[0,1/d_0]`$. The size-five source is the sum over its deletions. The map from each deletion to a physical (four-set, donor, contact) incidence is injective, even if the five-set admits several deletions.

Two separate tables stop at size two. One observes the high-degree singleton at the older end of the pair, the other at the younger end. Their old-singleton ceilings are $`(1399,48)`$. At size two the target ceilings are $`(1024,24)`$ and the donor ceilings $`(16,8)`$. If the selected target is pooled its source is zero. At an ordinary selected target the source is obtained from (2) as above whenever its degree is in the specified singleton range. These two tables book distinct directed incidences of $`F_1`$.

For a nonpooled member let $`a_j`$ and $`b_j=(D_j-c_j)/D_j`$ be its true gain and promotion rates. Add the known reciprocal of each internal nonpooled neighbour to $`D_j-c_j`$. Each internal pooled neighbour contributes an unknown in $`[0,1/d_0]`$, each exterior core slot an unknown in $`[0,1/2]`$. This gives a rational interval $`[a_j^-,a_j^+]`$ containing the actual $`a_j`$, regardless of the exterior. Write $`V_j^+,V_j^b,V_{j,J}`$ for the score after gain, old score after promotion and additional enlarged-set score. The solver works backwards first in the number of members and then in each $`(D,c)`$ row. For every retained row it explicitly checks
```math
 uV\ \le\ f_{\rm src}
 +\sum_{j\ {\rm exact}}\big\{\min_{a_j^-\le t\le a_j^+}t(V_j^+-V)
                      +b_j(V_j^b+V_{j,J}-V)\big\}.                 \qquad\text{(3)}
```
The sign of $`V_j^+-V`$ determines the endpoint. With $`q=2\prod_j D_j`$ using $`d_0`$ at a pool, $`2000q`$ times (3) is an integer inequality. The floating trial merely proposes an integer $`0\le V\le S`$. The last assertion of each row substitutes that integer in the full signed-128-bit residual. The configured grids have $`q,2000q(5D_{\max}+10)+1133q<2^{62}`$, and scores and successors are at most $`S`$ and $`2S`$ respectively, so the products fit. Every zero-score omitted configuration satisfies its inequality because future values and sources are nonnegative.

Upon summation of (3), changes of every old set and births of new sets occur with their physical rates. Adding the three tests, their newborn values add, and
```math
 F_4+F_1\ \ge\ (u-\mathcal L)T+v B.                                \qquad\text{(4)}
```
A bounded score per set has $`\mathbb E T(t)=O(e^{\mu t})`$ by the global core-degree moment estimate in [the bulk drift test](physical-test.md). This argument allows two unbounded-degree intermediate vertices in the same set.

## 3. Mean size of the witness

The raw comparison assigns an outward ceiling $`m(D,c)`$ to the unsigned contribution of each actual noninitial core, including its assigned incident motifs. Then $`\Phi\le\sum m(D,c)+E_1`$ with a founder remainder of polynomial expected size. Let $`M`$ be the global pointwise upper comparison in that witness. The [one-vertex ceiling](lifetime.py) constructs a nonnegative function $`Y(D,c)`$ for $`D\le1399,c\le60`$ and uses the constant $`\lceil M/l\rceil`$ beyond. In reverse $`D,c`$ order use the smaller of that constant and the upward-rounded maximum of
```math
 \frac{m(D,c)+aY(D+1,c)+bY(D,c+1)}{l+a+b},\qquad
 b=1-c/D,\quad a\in\{D-c,\ D-c+c/2\}.                              \qquad\text{(5)}
```
Both entries in that minimum are valid upper row values: all successors are at most the constant and $`m\le M`$. The same proves the outer rows; no return from the outer region is possible. The quotient is linear fractional with positive denominator, so endpoints suffice. For $`Y_{\rm tot}`$ summed over old noninitial vertices and $`y=Y(2,1)`$ this gives
```math
 \Phi\le(l-\mathcal L)Y_{\rm tot}+y B+E_1.                           \qquad\text{(6)}
```
Equations (5)–(6) are used on each unsigned witness separately.

## 4. Combining the signs

Use $`[Z]_s=\int_0^\infty e^{-st}\mathbb E Z(t)\,dt/\int_0^\infty e^{-st}\mathbb E C(t)\,dt`$. The denominator diverges at $`\mu`$ and $`[B]_s\to\mu`$. For $`s>\mu`$ the integrated endpoints for $`\Phi,T,Y_{\rm tot}`$ vanish. Along any sequence $`s\downarrow\mu`$ choose a convergent subsequence $`[\Phi]_s\to m`$; the raw positive comparison makes $`m>0`$. Equations (1), (4), (6), $`u\ge\mu\ge l`$ and nonnegativity give
```math
 \sigma(\lambda-\mu)m\ge\mu v-e,\qquad 0<m\le\mu y.
```
For $`v>e/l`$ it follows that $`\sigma(\lambda-\mu)\ge(v-e/l)/y`$.

| Quantity in table units | Lower witness | Upper witness |
|---|---:|---:|
| Trial $`\lambda`$ | $`5656453/10^7`$ | $`56581165/10^8`$ |
| Total signed violation $`e`$ | 115075 | 155608 |
| Quartet donors, $`v_4`$ | 7402100 | 7709975 |
| Singleton donors, older receiving | 11253928 | 8148713 |
| Singleton donors, younger receiving | 301604 | 3226994 |
| Ceiling $`y`$ | 2654889165321 | 2655887184018 |

Substituting these integers gives the inward test before rounding,
```math
 \frac{1696971095548923704669}{3000024756812730000000}
 \ \le\ \mu\ \le\ \frac{1698065802472818376961}{3001152517940340000000}.
```
Taking the floor on the left and the ceiling on the right to ten decimal places gives the opening interval. Neither calculation uses Sabba's endpoint assertions as inputs.

## Reproduction

From the repository root, the following command reconstructs the two witness files, checks their frozen digests, scans all local rows and envelopes, tests the analytic high-degree part, recomputes the three sets of lifetime lower rows and the lifetime ceiling, and prints both final fractions. It runs on one CPU. The recorded environment used CPython 3.12.3, GCC and G++ 13.3.0, GMP 6.3.0 and its C++ bindings on little-endian 64-bit Linux. The two 611,658,376-byte snapshots are generated in private scratch space. The verifier uses integer or rational acceptance tests on the resulting bytes.

~~~sh
python3 -B round2/occupation/replay.py
~~~

The complete rebuild printed the following final line. [Full stable output](../outputs/occupation-verified.txt) · [Witness digests and row counts](../motif-cert/certificate.json) · [Small-tree control](check_five.py)

~~~text
Outward physical interval [0.5656523639, 0.5658045676].
~~~
