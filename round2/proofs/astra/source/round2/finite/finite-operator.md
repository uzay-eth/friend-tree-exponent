# A sequence of finite branching matrices for the friend-tree exponent

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 1: Claude and Gemini; inputs include finite seam and uniform one-port comparisons. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## 1. The branching law

Fix integers $`M\ge3`$ and $`B\ge1`$. The degree at which a vertex becomes a hub is $`M`$. The maximum number of active low-degree cores in one component is $`B`$. Put $`s=\lceil B/(M-1)\rceil`$. The initial type is the three-vertex tree with one core, two leaves and no missing neighbor.

There is one hub type $`H`$. Every other type is an unrooted tree of $`1,\ldots,B`$ active cores. At vertex $`u`$ store its full degree $`2\le D_u<M`$ and its number $`L_u\ge0`$ of leaves. With $`a_u`$ active neighbors, $`D_u-L_u-a_u`$ is the nonnegative number of inactive core-neighbor slots. Old connections to hubs and edges deliberately cut between components both contribute slots. Each edge of the active tree is ordinary. All inactive slots keep their place in the denominator and send no leaves into the vertex.

A hub stays alive and creates one singleton component $`(D,L)=(2,1)`$ at rate one. This newborn has one inactive slot connecting it to the hub. The hub uses an unlimited leaf supply in this auxiliary law, including when the corresponding physical vertex would have no leaves.

In an active component the transitions are as follows.

1. A gain at $`u`$ has rate $`L_u+\sum_{w\sim_{\rm active}u}1/D_w`$. Increase $`D_u,L_u`$ by one. If $`D_u`$ has reached $`M`$, replace $`u`$ by a hub and release the connected components obtained by removing $`u`$. Every former neighbor retains its old full degree and leaf count, so the lost active edge becomes an inactive slot. Inactive edges already incident to $`u`$ never reconnect.
2. A promotion at $`u`$ has rate $`L_u/D_u`$. Decrease $`L_u`$ by one, leave $`D_u`$ fixed and attach a newborn active core in state $`(2,1)`$. At size at most $`B`$ this completes the transition. At size $`B+1`$, choose uniformly among the edges whose two sides have at least $`s`$ cores and release the two components. Both endpoints retain $`D,L`$, making the severed edge a slot on each side.

The oversized tree has maximum active degree at most $`M-1`$. Take a centroid, so each branch has size at most $`(B+1)/2`$. Its largest branch has at least $`B/(M-1)`$ cores, and the other side is at least as large. Thus a legal cut exists. All components produced by one transition continue independently. The auxiliary total core count, denoted $`\widetilde C`$, increases by one at promotions or hub births and is preserved at every other transition. Its birth rate is at most $`\widetilde C`$.

## 2. The finite rational operator and its accessible classes

Index types reachable from the initial singleton. Write $`\nu_j`$ for the number of outcome particles of type $`j`$ when one particle of type $`i`$ experiences event $`e`$. Define the row generator acting on column weights by
```math
 G_{ij}=\sum_e r_{i,e}\big(\mathbb E_e\nu_j-\mathbf1_{\{i=j\}}\big).
```
For a cut, the expectation averages over all allowed edges, counting distinct physical edges before identifying isomorphic outcomes. The entries are rational and the off-diagonal entries are nonnegative. With $`c_i`$ the number of cores in type $`i`$ and $`c_H=1`$, the expected count from one particle of type $`i`$ at time $`t`$ is $`(e^{tG}c)_i`$. The finite-type backward equation follows by conditioning on the first event. The total particle birth rate and total event rate are bounded by constants depending on $`M,B`$ times total cores, ensuring nonexplosion and finite first moments.

It follows that the auxiliary mean exponent is
```math
 \lambda_{M,B}=\max_{S\ {\rm reachable}}\;\max\{\Re z:z\in{\rm spec}(G_{SS})\},        \qquad\text{(13)}
```
where $`S`$ ranges over strongly connected classes in the off-diagonal graph. For the upper bound put $`G`$ in block triangular form and use its matrix exponential. For the lower bound any path of positive rates to a maximizing class gives a positive expected number of its particles at time one. Within that class, the irreducible Metzler matrix has a strictly positive Perron vector after adding a positive diagonal shift, and the strictly positive reward $`c`$ dominates a multiple of that vector. Singleton classes have the same assertion directly. This proves a limit, including when a maximizing class can feed other classes. Each $`\lambda_{M,B}`$ is real algebraic. Dropping unreachable classes is essential in (13).

## 3. Counting all additional cuts

For the proof retain the full finite core genealogy, including inactive edges. This annotation is unnecessary for the matrix transitions. Designate as *seams* only the edges severed at size $`B+1`$. A release at a hub adds no seam. Existing seams stay in the annotated genealogy, even if an endpoint becomes a hub. Let $`J`$ be the total number of seams and let
```math
 P=\sum_{\text{active components }A}(|A|-s)_+ .
```
A new core increases $`P`$ by at most one before any required size cut. Replacing a low vertex by a hub cannot increase $`P`$, since $`\sum_i(n_i-s)_+\le(\sum_i n_i-s)_+`$ and the retained fragments have total size one less. A size cut with side sizes $`n_1,n_2\ge s`$ decreases $`P`$ by exactly $`s`$. At the founding singleton $`J=P=0`$. Consequently $`sJ+P\le\widetilde C-1`$ for every history. This also holds pathwise throughout any continuation from a reachable annotated history.

## 4. Uniform comparison on one time interval

At any stopping time with a finite reachable annotated auxiliary state, take a true lift of its full genealogy. Give each nonhub its stored leaves and each hub of core degree $`c_u`$ exactly $`\max(0,M-c_u)`$ leaves. Thus its true degree is at least $`M`$. Restore all core edges, including seams and edges to old hubs, as ordinary edges. The resulting object is a finite genuine tree with $`C_0=\widetilde C(0)`$ cores. The coupling below is conditional on the starting state and time is reset to zero.

Initially all cores are paired and clean. A clean low pair has matching full degrees, leaf counts and paired core neighbors, including the labels of inactive slots in the auxiliary genealogy. A clean hub pair has the same paired list of core neighbors and true degree at least $`M`$; no equality of leaf counts is required. Dirty vertices in either process remain dirty. Newly created dirty cores may be unpaired. Write $`I`$ for the number of dirty cores in the two processes combined.

* At a clean low pair, synchronize the leaf births and the rate-one core ring, choosing corresponding neighbors. A promotion creates a new clean pair. An active send to a clean low neighbor updates both. A low-to-hub send can add a leaf to a true high core and leaves that pair clean. If a low sender's chosen slot is a seam and its true send reaches a clean low recipient, mark both copies of that recipient dirty. A send to a dirty recipient needs no new marks. An ordinary synchronised gain reaching degree $`M`$ turns the pair into a clean hub pair.
* At a clean hub use its common rate-one ring. If the true choice is a leaf, both processes create a paired newborn. If the true choice is a core neighbor, the auxiliary law creates an unmatched child and the true law sends one leaf to that neighbor. Mark the source pair, its newborn and the recipient pair if clean. This adds at most five dirty cores. True leaf births at a clean hub require no auxiliary action. The sum of these exceptional ring rates is at most $`\sum_{\rm clean\ high} c_u/D_u\le2C_{\rm true}/M`$ by the connected core-tree handshake.
* Use the correct individual rates at dirty vertices. A core ring in either copy adds at most two dirty cores: a new child is dirty, or a hit on a clean recipient marks that pair. Leaf births at a dirty core affect only its own stock. No clean vertex acquires an unmatched neighbor: births adding an unmatched child have a dirty parent after the event. Size cuts change only the labels of inactive slots. Their random choices may be sampled independently from unused uniforms.

These rules give the correct marginal rates, including discarded low-to-slot rings as null auxiliary events. Clean low seam sends from both orientations have total rate at most $`J`$ because each orientation has rate at most $`1/2`$ and the true and auxiliary rings at their clean source are common. Other low-to-inactive sends go to a hub in the auxiliary genealogy; a clean high receiver absorbs them in the comparison and a dirty receiver already carries its marks.

Both expected core counts at time $`w`$ are at most $`C_0e^w`$. Localization at bounded counts, integration of the dirty-count drift and then monotone localization bounds give
```math
 \mathbb E I(t)\le
 \int_0^t\left(2\mathbb E I(w)+10\mathbb E C_{\rm true}(w)/M
                         +2\mathbb E\widetilde C(w)/s\right)\,dw .
```
The seam term uses Section 3 for the whole history, including seams present at the start. Since $`I(0)=0`$ and $`|C_{\rm true}(t)-\widetilde C(t)|\le I(t)`$, Grönwall gives the uniform bound
```math
 \frac{|\mathbb E\widetilde C(t)-\mathbb E C_{\rm true}(t)|}{C_0}
 \le d_{M,B}(t):=(10/M+2/s)(e^{2t}-e^t).                                            \qquad\text{(14)}
```
This compares to the true lift chosen at the start of each individual window. The lift does not have to reproduce the true history before that window.

## 5. Rates over arbitrarily many intervals

Let $`m(t),u(t)`$ be the infimum and supremum of $`\mathbb E_TC(t)/C(T)`$ over finite true trees with at least one core. From (14) and the Markov property, for any $`t`$ with $`d_{M,B}(t)<m(t)`$,
```math
 \frac{\log[m(t)-d_{M,B}(t)]}{t}\le\lambda_{M,B}
 \le\frac{\log[u(t)+d_{M,B}(t)]}{t}.                                                  \qquad\text{(15)}
```
In fact the two one-step bounds hold after every block conditional on its entire reachable augmented state, and iteration from the singleton gives (15). The same iteration in the true tree gives $`\log m(t)/t\le\mu\le\log u(t)/t`$. The means from $`T_3`$ and the original edge have the same exponent because the latter can reach $`T_3`$ at time one with positive probability and its only initial delay is the wait for its first ring.

Here are completely finite replacements. Use $`U_j(t)`$ from (10) of the [one-edge upper argument](one-port-upper.md), evaluated with its integer scale $`j`$. If $`r=t-h\ge2jk`$, with $`k,h,\eta_j,p_j`$ from Section 9 of that argument, put
```math
 L_j(t)=\max\left(1,\frac{p_j}{8k}e^{-\eta_jr}F(r)\right);
```
otherwise put $`L_j(t)=1`$. The planting argument in Proposition 6.2 gives $`m(t)\ge L_j(t)`$: its conditional bound is uniform in the starting tree. The upper argument gives $`u(t)\le U_j(t)`$. Thus when $`d=d_{M,B}(t)<L_j(t)`$ the explicit error bound
```math
 |\lambda_{M,B}-\mu|
 \le \frac1t\log\frac{U_j(t)+d}{L_j(t)-d}                                             \qquad\text{(16)}
```
follows since the interval in (15) with these replacements contains both numbers. Different integer scales may be used in $`L`$ and $`U`$. Only finite-time $`F`$ enters and Section 9 there provides its truncation error.

Finally take $`B=M^2`$, $`t=(\log M)/3`$. Then $`s\ge M`$, $`d\le12M^{-1/3}\le 1/2`$ for $`M\ge24^3`$. For $`r=t-h\ge2jk`$, the true growth inequality $`F(t)\le e^h(F(r)+1)\le2e^hF(r)`$ gives
```math
 \frac1t\log L_j(t)\ge \frac{\log F(t)}t-\eta_j
                                   -\frac{h+\log(16k/p_j)}t .
```
Combine this with (10)–(12), using $`L_j\ge1`$ and $`U_j\ge1`$, so the extra logarithmic loss from $`d\le1/2`$ is at most $`\log 3/t`$. Fixing $`j`$ and then increasing it proves $`\lambda_{M,M^2}\to\mu`$. Taking $`j=\lfloor\log t/(4\log\log t)\rfloor`$ for large $`t`$ gives $`|\lambda_{M,M^2}-\mu|=O(\log\log t/\log t)=O(\log\log\log M/\log\log M)`$. Equation (16) supplies the effective bound. Small matrices may be evaluated separately from the asymptotic choice of $`B`$.

## 6. Six exact matrices

The table reports rational enclosing intervals for the auxiliary growth rates. The largest class size shown refers to the maximizing class.

| $`M`$ | $`B`$ | Reachable types | Classes | Maximizing class size | Lower endpoint | Upper endpoint |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 4 | 3 | 2 | 0.61803398874989 | 0.61803398874990 |
| 3 | 3 | 9 | 5 | 5 | 0.67833737318518 | 0.67833737318519 |
| 3 | 9 | 27 | 11 | 17 | 0.67928074796982 | 0.67928074796983 |
| 4 | 1 | 8 | 5 | 4 | 0.55612325604897 | 0.55612325604898 |
| 4 | 2 | 23 | 8 | 16 | 0.61485063481498 | 0.61485063481499 |
| 4 | 3 | 68 | 14 | 55 | 0.62755802142835 | 0.62755802142836 |

For $`M=3,B=1`$, label the types as the founder, the hub, the singleton with one leaf and one inactive slot, and the singleton with zero leaves and two inactive slots. In that order the full row generator is
```math
 \begin{pmatrix}
 -3&2&2&0\\ 0&0&1&0\\ 0&1&-1&1/2\\ 0&0&0&0
 \end{pmatrix}.
```
The hub and the one-leaf singleton form the maximizing class, giving $`\lambda_{3,1}=(\sqrt5-1)/2`$. For $`(M,B)=(3,3)`$ the listed root satisfies $`x^5+17x^4+94x^3+162x^2-132x-416=0`$ with $`x=2\lambda_{3,3}`$. For $`(4,1)`$ it satisfies $`x^4+15x^3+62x^2+12x-270=0`$ with $`x=3\lambda_{4,1}`$. The $`B=1`$ values at $`M=3`$ and $`M=4`$ fall on opposite sides of the existing certified interval for the true exponent. Individual small auxiliary matrices therefore supply neither a general upper bound nor a general lower bound by themselves. [True-exponent bounds](../../proofs/bracket-motifs.md)

Every row, type, strongly connected class, positive integer witness and coefficient of the dominant characteristic polynomial is stored in the [matrix certificate](../outputs/finite-matrices.json). A type is encoded as a rooted nested triple of full degree, leaf count and sorted child encodings, or the string “hub”. The builder chooses the lexicographically least encoding across roots. Rows have triples of column index, numerator and denominator. Class indices refer to rows in that file. For a winning class $`S`$, the recorded polynomial is $`\det(xI-L G_{SS})`$, where $`L`$ is the recorded positive integer scale.

The independent checker reconstructs each row by enumerating individual leaf and core rings. It adds slots explicitly at both ends of each cut and identifies isomorphic outcomes by pruning to the tree center. It checks core-count drift and the seam-budget drift in every row, and recomputes reachability and classes. For each class $`S`$ with integer witness $`w>0`$, it evaluates all rational ratios $`(G_{SS}w)_i/w_i`$. The gap between the extreme ratios is less than $`10^{-16}`$ in every class. The dominant class has a lower bound strictly above all other class upper bounds. The checker verifies each polynomial using exact determinants at $`\lvert S\rvert+1`$ integer points. It certifies opposite endpoint signs and a strictly positive derivative throughout the reported bracket using rational interval arithmetic. [Generator](finite_matrix.py) · [Independent verifier](verify_matrix.py) · Run log
