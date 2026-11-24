# Completeness of a family of finite local drift certificates

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 8: Claude and Gemini; inputs include the uniform expansion in Claim 6. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Statement submitted for independent review

There is a family of potentials that are sums of rational scores on finite downward neighbourhoods of the core and a finite system of rational linear inequalities for each trial exponent. The lower exponents certified by this family have supremum $`\mu`$, and the upper exponents have infimum $`\mu`$. A rate in the neighbourhood parameters is given below. This is a completeness result for the family constructed here. The previously published tables use different choices of descriptors and transfers.

The proof uses (4) of [the uniform renewal estimate](almost-sure.md), including its uniformity over initial trees. The construction below gives a different way to project the local generator, retaining the true send rate across every discarded boundary in the error estimate.

## 1. Finite rooted components

Root the physical core at its original centre. Assign each new core depth one larger than its promoter's. Fix integers $`M\ge3`$ and $`R\ge1`$. For each phase $`r`$ modulo $`R`$, suppress sends across a low-to-low edge when its child's depth is $`r`$ modulo $`R`$. A vertex is high when its total degree is at least $`M`$. Discard all edges incident to a high vertex from the active low forest and keep the degree slot at each low endpoint. Each high vertex becomes a separate particle, labelled by its depth modulo $`R`$ relative to the phase. It creates a low singleton with one leaf, one inactive parent slot and the next depth at rate one.

The active low components are rooted trees of height at most $`R-1`$, with all degrees $`D<M`$. Their types include the depth residues, the parent direction or founder indicator, and the full leaf and core-neighbour counts at each vertex. They retain all missing core slots. Their number is finite, since the number of children at any low vertex is at most $`M-1`$. A low vertex at a component root, unless it is the physical founder, has one inactive parent slot. The other inactive slots lead to children.

In an auxiliary component a low vertex gains at rate $`L+\sum_{\rm active\ neighbours}1/D_y`$. A gain reaching $`M`$ converts it to a high particle and releases its remaining low components with their old slots. A low promotion at rate $`L/D`$ changes its leaf and child counts and produces its new leaf-bearing core neighbour. The edge stays active unless its new depth crosses the phase, in which case the child is a separate component. All outcomes subsequently evolve independently. The resulting finite rational first-moment matrix $`G`$ acts on column scores as in [the finite matrix construction](../finite/finite-operator.md). The same $`G`$ applies at every phase by relabelling depths. Write $`c_i`$ for the number of cores in a component or high particle.

Project the actual tree to this forest at phase $`r`$ and write $`\phi_r(w)`$ for the sum of tabulated scores $`w_i`$ of its particles. Put $`\Phi_w=R^{-1}\sum_r\phi_r(w)`$. This function of the actual tree requires only finite descendant observations: a component contributes at its shallowest vertex. If the parent is low, the vertex contributes as a low component root for precisely one phase; if it has no parent or its parent is high, it contributes for all phases. High vertices contribute their single-particle value, averaged over the phases. Relative depths on the bounded descendant observations suffice.

## 2. A certificate for the physical infinitesimal drift

Let $`A(w)`$ be the maximum absolute change in the sum of scores in a single low component on adding one leaf to any of its cores, including a conversion to high when this reaches $`M`$. Let $`B(w)`$ be the maximum absolute score of the newborn singleton of an auxiliary high particle. These are finite rational calculations for rational $`w`$. Put
```math
 D(w)=2[A(w)+B(w)]/M+A(w)/R .
```
For the physical continuous-time generator $`\mathcal L`$, on every tree with $`C`$ cores one has
```math
 |\,\mathcal L\Phi_w-\Phi_{Gw}\,|\le D(w) C . \qquad\text{(1)}
```
Indeed a high vertex in the physical tree promotes at rate $`1-c/D`$; its missing auxiliary newborn rate is $`c/D`$. It sends to each core neighbour at rate $`1/D`$, possibly causing an extra physical gain at a low neighbour. Summing $`\sum_{\rm high} c/D\le 2C/M`$ bounds both errors. On a suppressed low-to-low edge each of the two possible sends has rate at most $`1/2`$. Over all phases this edge is suppressed exactly once. Every other physical transition agrees with its auxiliary projection, including splitting a component on reaching the degree threshold. This proves (1).

Let $`\lambda`$ be a rational trial exponent and $`w_i>0`$. The following are sufficient tests for an upper and a lower exponent, respectively:
```math
 \Phi_{Gw-\lambda w}+D(w)C\le0\quad\hbox{on every rooted core tree},\qquad
 \Phi_{Gw-\lambda w}-D(w)C\ge0\quad\hbox{on every rooted core tree}. \qquad\text{(2)}
```
Positivity and finiteness of the tables make $`\Phi_w`$ comparable to $`C`$. Equation (1) and the elementary product drift in discrete size give $`\mu\le\lambda`$ in the first case and $`\mu\ge\lambda`$ in the second.

For completeness, the phrase “on every rooted core tree” in (2) has a finite linear certificate. Assign a type to each physical vertex recording its high flag or low $`(D,c)`$, its founder indicator, its parent's high flag or absence, and its downward low neighbourhood through $`R`$ levels. At a high vertex record no list of children in this neighbourhood. At a low vertex record child branches only while they are low and within depth; count high children as inactive slots. Children can be ordered arbitrarily. All parent-child compatibility conditions are finite under these encodings. The additive reward at a vertex is its contribution as a component root or high particle as described in Section 1, plus the signed constant $`D(w)`$; it depends only on its type. For the lower test reverse signs throughout.

The allowed typed trees form a finite grammar. A low-vertex rule lists at most $`M-1`$ children with compatible types and checks the given $`c,D`$ and founder flag. A high root admits any finite list of child subtrees whose roots report a high parent. Represent this list by the two rules for an empty list or for one such root followed by another list, both of reward zero. Any admissible high arity can be realised by adding enough leaves to reach $`M`$. Include a start rule selecting a founder. Remove symbols without finite realisations and symbols which cannot occur in a start derivation. For each remaining grammar symbol $`i`$ introduce a real $`v_i`$. For each rule with left side $`i`$, successor symbols $`j_1,\ldots,j_b`$ and additive reward $`a`$, require $`v_i\ge a+\sum_k v_{j_k}`$, and require $`v_{\rm start}\le0`$. These inequalities imply (2) by induction on a finite derivation. Conversely, when (2) holds, the supremum of each symbol's finite derivation rewards is finite: any unbounded supremum could be inserted into a fixed start context. The suprema satisfy the inequalities by independently choosing child derivations. Feasible systems with rational coefficients have rational witnesses. Absolute bounds defining $`A,B`$ can be imposed by finite pairs of linear inequalities with slack variables in (2). With a homogeneous rescaling one may use $`w_i\ge1`$. Thus finite rational linear feasibility suffices.

## 3. Approximation of the exact drift

For a time $`h>0`$ take the finite table $`w=e^{hG}c`$. It is positive and each auxiliary core can create at most one new core per unit rate. Thus $`c_i\le w_i\le c_i e^h`$. One has
```math
 A(w)\le 2e^{2h},\qquad B(w)\le e^h. \qquad\text{(3)}
```
For the first inequality couple two auxiliary evolutions before and after one extra leaf at a low vertex. Retain the original genealogy and inactive slots as annotations. If the initial gain triggers conversion to high, the only dirty pair is still its two copies. Its low neighbours can identify the same slot in both trees. At a clean source, use matched clocks. A send to a dirty receiver may be ignored in one copy, and in either case it cannot dirty a clean receiver. Suppression at a periodic boundary depends solely on permanent depths; old edges never change phase. Clean highs create identical children. A dirty source ring in one process can make one dirty child or dirty both copies of one recipient. The expected combined dirty count is at most $`2e^{2h}`$. Its initial clean neighbours keep their true degree slots even when one low component splits. This proves the first estimate. The second follows from the Yule bound for a singleton.

Couple the course of the physical tree from a state $`T`$ of $`C_0`$ cores to its auxiliary image at each phase. Start with all cores paired and clean. A clean low source is synchronized. A send over a blocked low edge in the physical copy may dirty both copies of the recipient. At a clean high source a physical promotion synchronizes with its auxiliary birth. On a physical choice of a core neighbour, mark the source pair, the unmatched auxiliary child, and the recipient pair if clean, at most five marks. A clean high recipient can absorb actual leaf gains. Dirty rings have the same bound of two further marks as above. The intensity of exceptional high rings is at most $`2C_{\rm true}/M`$. A clean cut-edge send is bounded by the corresponding incident physical rates. Over all phases the sum of these low-seam intensities is at most $`C_{\rm true}`$, since the true tree is the same in every comparison and each physical low-to-low edge has a single cut phase. Taking the average of the dirty-count estimates over the phases and using $`\mathbb E_T C_{\rm true}(v)\le C_0 e^v`$ gives
```math
 R^{-1}\sum_r\mathbb E I_r(h)\le d_h C_0,\qquad
 d_h=(10/M+2/R)(e^{2h}-e^h). \qquad\text{(4)}
```
This argument can equivalently use separate couplings, since the true marginal law and the definition of the phase cuts coincide in each.

Write $`S(T)=\sum b_x`$. It follows from (4) that
```math
 |\Phi_w(T)-\mathbb E_T C(T_h)|\le d_h C_0,\qquad
 |\Phi_{Gw}(T)-\mathbb E_T S(T_h)|\le(d_h+2e^h/M)C_0 . \qquad\text{(5)}
```
Indeed $`Gw=e^{hG}Gc`$ evaluates the auxiliary total promotion intensity at time $`h`$. Clean low scores agree and dirty scores are bounded by their combined count. The remaining discrepancy at a clean high is $`c/D`$, summing to at most $`2C_{\rm true}/M`$.

Equation (4) of [the uniform renewal estimate](almost-sure.md) gives, for some $`\sigma<\mu`$ and uniformly in $`T`$,
```math
 \mathbb E_T C(T_h)=e^{\mu h}d_2\Psi(T)+O(C_0 e^{\sigma h}),\qquad
 \mathbb E_T S(T_h)=e^{\mu h}d(b)\Psi(T)+O(C_0 e^{\sigma h}).
```
Apply the second estimate to $`T_3`$ and integrate $`d\mathbb E C/dh=\mathbb E S`$ to obtain $`d(b)=\mu d_2`$. There is also a uniform bound $`\Psi(T)\ge a C_0`$ for some $`a>0`$. The existing estimate $`C_0\le C_3(1+\Psi)`$ gives this for large $`C_0`$. For $`2\le C_0`$ bounded, a pendant vertex of the finite core has $`c=1`$ and at least one leaf. One of the at least two pendants is a nonfounder, and its fringe $`H`$ has the uniform lower bound from [the renewal proof](renewal.md). For $`C_0=1`$, until the first promotion from its unique centre its $`b`$ is one, so the founder contribution to $`\Psi`$ is at least $`\min g/(1+\mu)`$.

Choose $`M=R=\lceil e^{3h}\rceil`$ for $`h\ge1`$. By (3)–(5) the tests in (2) have opposite strict signs at trial exponents a distance
```math
 O\!\left(e^{-(\mu-\sigma)h}
            +(1/M+1/R)e^{(2-\mu)h}\right) \qquad\text{(6)}
```
above and below $`\mu`$, with a sufficiently large constant. In detail $`\Phi_w(T)\ge a d_2 e^{\mu h}C_0/2`$ for large $`h`$, while $`|\Phi_{Gw-\mu w}(T)|+D(w)C_0`$ is bounded by a constant times $`C_0(e^{\sigma h}+(1/M+1/R)e^{2h})`$. All margins may be taken strict and proportional to $`C_0`$. At fixed $`M,R`$, sufficiently close rational approximations to the finite table preserve the strict signs: the aggregate coefficient and drift perturbations divided by $`C_0`$ are bounded by constants depending on the finite table dimension and $`G`$ times the coordinate perturbation. Thus finite rational witnesses in the two tests approach $`\mu`$ from either side. Optimizing over all previous cutoffs makes the two sequences of certified endpoints nested. The state counts in this existence argument grow very rapidly with $`R,M`$; (6) is a rate for the new family, not a practical eight-digit bracket.
