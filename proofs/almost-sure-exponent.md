> Produced by GPT-6.2 Astra.
>
> Complete write-up: Claude Opus 5.5 (`claude-opus-5-5`), 2026-11-23. The argument, the constants and the conclusion are Astra's original ones. Every step is written out: the reserved randomness of each site, the construction at random tagged clocks, the exact product law of the preparation flags, the comparison inside one block, the block-to-block argument, and the return to discrete size. Section 7 lists what changed against the short version. No change alters a constant or the structure of the argument.

# The sample-path exponent equals the exponent of the mean

We use the notation and results of [the mean-exponent proof](existence.md), cited as "E-Lemma 1.1" and so on: the sources and admissible constructions of E-Section 1.1, the elementary dynamics (D1)–(D4), the constants $`p_1,p_2,p_3`$, $`h=k+1`$, $`\eta(M,\epsilon)`$, $`s_0`$, $`\theta_k(s)`$, the reference process and the triggering rings of E-Section 4, and the quantity $`Q_{T,v}(s)`$ of E-Corollary 5.3. Let $`\mu`$ be the exponent of E-Theorem 1.

**Theorem 2.** Almost surely,

```math
\lim_{n\to\infty}\frac{\log C_n}{\log n}=\mu .
```

The proof is in continuous time. The lower bound compares one block of time with independent copies of one finite experiment, one copy at each of at least $`C/6`$ disjoint sites. The upper bound is Markov's inequality.

## 1. Disjoint sites

Let $`T_0`$ be a finite tree with $`C\ge1`$ core vertices. Write $`\mathcal A`$ for the set of core vertices with $`c\le2`$ and $`L\ge1`$, and $`A=|\mathcal A|`$. As in E-Section 2, a core vertex is *of type B* if $`c=2`$ and $`L=0`$.

**Lemma 1.1 (sites).** There is a list of *sites*, determined by $`T_0`$ through a fixed rule, with the following properties.

- Each site is either an *A-site* $`(u;\ell)`$, where $`u\in\mathcal A`$ and $`\ell`$ is the lowest-labelled leaf of $`u`$, or a *triple* $`(u,w,z)`$ of type-B vertices forming a path $`u\sim w\sim z`$, in which case the neighbours of $`w`$ are exactly $`u`$ and $`z`$.
- The *reserved sets* $`\{u,\ell\}`$ and $`\{u,w,z\}`$ of different sites are disjoint.
- The number $`m`$ of sites satisfies $`m\ge C/6`$.

In both cases we call $`u`$ the *base* of the site.

*Proof.* If $`C=1`$, the core vertex is the centre of a star, with $`c=0`$ and $`L\ge2`$. It is the base of one A-site, and $`1\ge1/6`$.

Let $`C\ge2`$. Every core vertex has $`c\ge1`$. A core vertex with $`c=1`$ has $`L\ge1`$, so it lies in $`\mathcal A`$, and summing $`c-2`$ over the core tree gives, as in E-Lemma 2.1, that the number of core vertices with $`c\ge3`$ is at most $`A-2`$. The non-B core vertices are those of $`\mathcal A`$ and those with $`c\ge3`$, so there are at most $`2A-2`$ of them, and at least two, since the core tree has at least two vertices of core degree one.

Every type-B vertex has exactly two core neighbours and no leaves. The type-B vertices therefore induce a disjoint union of paths in the core tree. Let $`P`$ be one of these maximal paths. Since each vertex of $`P`$ has exactly two core neighbours, $`P`$ has exactly two edges to core vertices outside $`P`$ (for a single-vertex path, both edges leave its one vertex), and the two outside endpoints are not of type B. They are distinct, since the core tree has no cycle. Replacing each maximal path by a single edge between these two non-B neighbours turns the core tree into a tree on the non-B vertices, in which different paths give different edges. So the number of maximal paths is at most $`(2A-2)-1=2A-3`$.

Cut each maximal path, from one end, into consecutive disjoint triples $`(u,w,z)`$, leaving at most two of its vertices unused. Let $`J`$ be the total number of triples. In a triple, $`w`$ has type B and its two core neighbours are $`u`$ and $`z`$, so its neighbours are exactly $`u`$ and $`z`$. Counting core vertices,

```math
C\le A+(A-2)+3J+2(2A-3)=6A+3J-8 .
```

Take one A-site for each vertex of $`\mathcal A`$ and one triple site for each triple. Then $`m=A+J\ge A+J/2\ge C/6`$.

The vertices of $`\mathcal A`$ are not of type B, the triples are disjoint, and a leaf has exactly one neighbour, so the leaves $`\ell`$ of different A-sites are different and are not core. Hence the reserved sets are disjoint. $`\square`$

## 2. One block

Fix $`M\ge2`$, $`0<\epsilon<1`$, an integer $`k\ge\max(2,6M/\epsilon)`$, and $`s\ge s_0(k,M,\epsilon)`$. Let $`T`$ be a $`k`$-vertex tree maximising $`\mathbb E C_T(s)`$, and $`v`$ the root given by E-Corollary 5.3. Put $`\Lambda=h+s`$.

Fix the planting script of E-Lemma 3.1, stage 4, for $`(T,v)`$: a neighbour $`v_2`$ of $`v`$, and an order $`y_3,\dots,y_k`$ of the other vertices with parents earlier. Number the vertices of $`T`$ by $`1=v`$, $`2=v_2`$, and $`j=y_j`$ for $`3\le j\le k`$. For the $`j`$-th addition ($`1\le j\le k-2`$) let $`a(j)`$ be the parent of vertex $`j+2`$, let $`z(j)`$ be the lowest-numbered neighbour of $`a(j)`$ among vertices $`1,\dots,j+1`$, and let $`d_j`$ be the degree of $`z(j)`$ in the tree induced on vertices $`1,\dots,j+1`$, plus one if $`z(j)=1`$. These are deterministic, and $`d_j\le k-1`$. Put

```math
\hat p_4=\prod_{j=1}^{k-2}\frac{e^{-(j+3)}}{d_j},\qquad \hat p=p_1p_2p_3\hat p_4 .
```

Since $`j\le k-2`$ and $`d_j\le k-1<k`$, $`\hat p_4\ge (e^{-(k+1)}/k)^{k-2}`$, so $`\hat p\ge p(M,k)`$ of E-Section 3.

### 2.1 The reserved randomness

Start at time $`0`$ from a deterministic labelled tree $`T_0`$ with $`C\ge1`$ core vertices, in which each vertex is driven by its own source. List its sites by Lemma 1.1. For site $`i`$ with base $`u_i`$ the following randomness is *reserved*:

- the sources of the vertices of its reserved set, on the time interval $`(0,\Lambda]`$;
- a family of unused sources $`\sigma^{(i)}_0,\sigma^{(i)}_1,\dots,\sigma^{(i)}_k`$ (the *tagged clocks*) and a family $`(\zeta^{(i)}_j)_{j\ge1}`$ of unused sources for reference newborns;
- a family of independent uniform *thinning coins* $`\chi^{(i)}_1,\chi^{(i)}_3,\chi^{(i)}_{4,1},\dots,\chi^{(i)}_{4,k-2}`$. The coins never influence the tree.

Different sites have disjoint reserved randomness, and all of it is independent of the sources of the unreserved vertices.

### 2.2 The construction at random tagged clocks

The construction is the standard one of E-Section 1.1 except for the following rules, applied at every site $`i`$. Units are $`(0,1]`$, $`(1,2]`$, $`(2,3]`$, and $`(2+j,3+j]`$ for $`1\le j\le k-2`$; the comparison window is $`(h,h+s]`$.

- (R1) *Unit 1, triple site.* At rings of $`w_i`$ rank $`u_i`$ first. If the first ring of $`w_i`$ in $`(0,1]`$ chooses $`u_i`$, its newborn is driven by $`\sigma^{(i)}_0`$.
- (R2) *Unit 3.* At rings of $`u_i`$ in $`(2,3]`$ rank the leaves of $`u_i`$ first, by label. At the first ring of $`u_i`$ in $`(2,3]`$, if its coordinate $`U`$ is at most $`M/(M+2)`$ **and** the chosen vertex is a leaf of $`u_i`$, call that leaf $`v_i^*`$. Immediately after that ring has been processed, re-drive $`v_i^*`$ by $`\sigma^{(i)}_1`$, and let its newborn $`v_i'`$ be driven by $`\sigma^{(i)}_2`$. The promoting ring is a ring of $`u_i`$, decided by the source that drives $`u_i`$, so the later re-driving of $`v_i^*`$ does not affect it.
- (R3) *Stage 4.* If $`v_i^*`$ exists, the $`r`$-th vertex born in $`(3,h]`$ on the far side of $`u_iv_i^*`$ beyond $`v_i^*`$ is driven by $`\sigma^{(i)}_{r+2}`$ for $`r\le k-2`$, and by a new unreserved source for $`r>k-2`$. At a ring of a far-side vertex in the $`j`$-th unit, rank first the vertex that corresponds to $`a(j)`$ under the identification of the next paragraph, if there is one.
- (R4) *Window.* Call site $`i`$ *coupled* if at time $`h`$ the vertex $`v_i^*`$ exists, $`D_{u_i}\ge M+1`$, and the far side of $`u_iv_i^*`$ beyond $`v_i^*`$ consists of $`k`$ vertices driven by $`\sigma^{(i)}_1,\dots,\sigma^{(i)}_k`$ such that the map sending the vertex driven by $`\sigma^{(i)}_j`$ to vertex $`j`$ of $`T`$ is an isomorphism onto $`(T,v)`$. This is decided by the history at time $`h`$, which is a stage boundary. At a coupled site, run the window for $`u_i`$ and for the far side of $`u_iv_i^*`$ by the rules of E-Section 4 with $`u=u_i`$, $`v^*=v_i^*`$, the reference $`R_i`$ defined below, and the reserved newborn sources $`\zeta^{(i)}_j`$. At a site that is not coupled, the window uses the standard rules. Vertices outside all coupled far sides and bases follow the standard rules, or the rules of their own site.

Every newborn not named in (R1)–(R4) is driven by a new unreserved source.

In (R3) the identification is the one in which the vertex driven by $`\sigma^{(i)}_j`$ corresponds to vertex $`j`$ of $`T`$; ranking first a vertex that is not a neighbour has no effect.

**The reference at site $`i`$.** $`R_i`$ is the friend-tree process from $`T`$ in which vertex $`j`$ is driven by the restriction of $`\sigma^{(i)}_j`$ to $`(h,\infty)`$ (a point at absolute time $`r>h`$ is a ring of $`R_i`$ at time $`r-h`$), the $`r`$-th reference newborn is driven by $`\zeta^{(i)}_r`$, and choices use the coordinate $`U`$ with the lowest-labelled ranking. It is defined at **every** site, coupled or not, as a function of site-$`i`$ reserved randomness after time $`h`$ alone. It never reads the actual process.

**Lemma 2.1 (admissibility).** This is an admissible construction of the friend-tree process from $`T_0`$. Hence (E-Lemma 1.1) the annotated path has the friend-tree law from $`T_0`$.

*Proof.* Every assignment of a source is made immediately after a processed ring (a birth, or the re-driving in (R2)) or at the stage boundary $`h`$, and is a function of the history at that moment. The stage boundaries are the finitely many deterministic times $`1,2,\dots,h`$. No source drives two current vertices. Every leaf chosen under (R2), whether or not earlier flags of its site hold, is a leaf of $`u_i`$ that was not core at time $`0`$, so the far side of $`u_iv_i^*`$ contains no vertex that was core at time $`0`$; by E-Lemma 3.2 the far sides used in (R3) and (R4) at different sites are disjoint at all times, so the rules of different sites never act on the same vertex, and the construction is well defined for all sites simultaneously, before any restriction to successful sites. The tagged clocks of site $`i`$ are assigned only by (R1)–(R3) of site $`i`$, each to one vertex. A vertex re-driven in (R2) is a leaf of $`u_i`$, which is not in the reserved set of any other site $`j`$: the leaf $`\ell_j`$ of an A-site is a leaf of $`u_j\ne u_i`$, and triple vertices are core. The sources $`\zeta^{(i)}_r`$ are assigned to actual newborns only at a coupled site and only by rings strictly before its triggering ring, where actual far-side newborns are the reference newborns; such a newborn keeps its source afterwards, and after the trigger no further $`\zeta^{(i)}_r`$ is assigned. Each trigger affects only its own site. Each choice is uniform given the history strictly before the ring: by predictable ranking everywhere outside coupled windows, and inside a coupled window as checked in E-Section 4. $`\square`$

### 2.3 The preparation flags

For site $`i`$ define the following events. Each is an explicit event about the reserved randomness of site $`i`$ alone; none mentions the tree.

- $`F^{(i)}_1`$. *A-site:* the source of $`u_i`$ has no point in $`(0,1]`$, and $`\chi^{(i)}_1\le p_1e`$. *Triple:* the sources of $`u_i`$ and $`z_i`$ have no point in $`(0,1]`$, the source of $`w_i`$ has a point in $`(0,1]`$, and the first such point has $`U\le1/2`$.
- $`F^{(i)}_2`$. The source of $`u_i`$ has no point in $`(1,2]`$, and the *tagged leaf source* has at least $`M`$ points in $`(1,2]`$. The tagged leaf source is the source of $`\ell_i`$ at an A-site and $`\sigma^{(i)}_0`$ at a triple.
- $`F^{(i)}_3`$. The source of $`u_i`$ has exactly one point $`\tau_i`$ in $`(2,3]`$, that point has $`U\le M/(M+2)`$, the sources $`\sigma^{(i)}_1`$ and $`\sigma^{(i)}_2`$ have no point in $`(\tau_i,3]`$, and $`\chi^{(i)}_3\le e^{-2}/e^{-2(3-\tau_i)}`$.
- $`F^{(i)}_4`$. The source of $`u_i`$ has no point in $`(3,h]`$, and for every $`1\le j\le k-2`$, in the unit $`(2+j,3+j]`$: $`\sigma^{(i)}_{z(j)}`$ has exactly one point, at a time $`\beta_j`$, and it has $`U\le 1/d_j`$; the sources $`\sigma^{(i)}_r`$ with $`1\le r\le j+1`$, $`r\ne z(j)`$, have no point (the auxiliary clock $`\sigma^{(i)}_0`$ is not involved); $`\sigma^{(i)}_{j+2}`$ has no point in $`(\beta_j,3+j]`$; and $`\chi^{(i)}_{4,j}\le e^{-1}/e^{-(3+j-\beta_j)}`$.

Let $`\mathrm{Prep}_i=F^{(i)}_1\cap F^{(i)}_2\cap F^{(i)}_3\cap F^{(i)}_4`$.

**Lemma 2.2 (product law).** The events $`\mathrm{Prep}_1,\dots,\mathrm{Prep}_m`$ are independent, and each has probability exactly $`\hat p`$.

*Proof.* $`\mathrm{Prep}_i`$ is a function of the reserved randomness of site $`i`$ on $`(0,h]`$, and these families are independent for different $`i`$. Within one site, compute successively. The source of $`u_i`$ appears in all four flags and the tagged clocks recur in successive units, but every repeated use concerns a disjoint time interval, and a source's points in disjoint intervals are independent; within a unit, condition on the designated point time.

- $`F_1`$: at an A-site, $`e^{-1}\cdot p_1e=p_1`$; at a triple, $`e^{-1}\cdot e^{-1}\cdot(1-e^{-1})\cdot\frac12=p_1`$.
- $`F_2`$: $`e^{-1}\Pr(\mathrm{Pois}(1)\ge M)=p_2`$, independently of $`F_1`$, since only points in $`(1,2]`$ are involved. At a triple, $`\sigma^{(i)}_0`$ is involved only through $`(1,2]`$.
- $`F_3`$: the point of $`u_i`$ is at $`2+r`$ with density $`e^{-r}`$ and there is no further point with probability $`e^{-(1-r)}`$; the mark condition has probability $`M/(M+2)`$; the two tagged clocks are silent on an interval of length $`1-r`$ with probability $`e^{-2(1-r)}`$; the coin has probability $`e^{-2}/e^{-2(1-r)}\le1`$. Integrating over $`r\in(0,1)`$ gives exactly $`e^{-1}\cdot\frac{M}{M+2}\cdot e^{-2}=p_3`$.
- $`F_4`$: the silence of $`u_i`$ on $`(3,h]`$ has probability $`e^{-(k-2)}`$. In unit $`j`$, exactly one point with the mark condition has probability $`e^{-1}/d_j`$; the other $`j`$ tagged clocks are silent with probability $`e^{-j}`$; given $`\beta_j`$, the newborn clock and the coin together have probability $`e^{-(3+j-\beta_j)}\cdot e^{-1}/e^{-(3+j-\beta_j)}=e^{-1}`$. So unit $`j`$ contributes $`e^{-(j+2)}/d_j`$, and with the factor $`e^{-1}`$ per unit for $`u_i`$ the total is $`\hat p_4`$. The intervals of different units are disjoint.

The four events concern disjoint time intervals of the sources of $`u_i`$, disjoint intervals or distinct tagged clocks otherwise, and distinct coins, so the probability is $`p_1p_2p_3\hat p_4`$. $`\square`$

**Lemma 2.3 (what the flags guarantee).** On $`\mathrm{Prep}_i`$, site $`i`$ is coupled. That is, at time $`h`$, $`v_i^*`$ exists, $`D_{u_i}\ge M+1`$, and the far side of $`u_iv_i^*`$ beyond $`v_i^*`$ is a copy of $`(T,v)`$ whose vertex $`j`$ is driven by $`\sigma^{(i)}_j`$. Moreover that far side contains no vertex that was core at time $`0`$. This holds whatever happens at the other sites and outside all sites.

*Proof.* We follow E-Lemma 3.1 and check that every conclusion uses only the flags of site $`i`$ and (D1)–(D4).

*Unit 1.* At an A-site, $`u_i`$ does not ring, so by (D3) $`\ell_i`$ stays a leaf of $`u_i`$ and $`c_{u_i}\le2`$. At a triple, $`u_i`$ and $`z_i`$ do not ring, so by (D2) $`D_{w_i}=2`$ until the first ring of $`w_i`$, and that ring, having $`U\le1/2`$ with $`u_i`$ ranked first, chooses $`u_i`$. By (R1) its newborn is a leaf of $`u_i`$ driven by $`\sigma^{(i)}_0`$, and it stays a leaf because $`u_i`$ is silent. By (D3), $`c_{u_i}=2`$ at time $`1`$. In both cases $`u_i`$ has, at time $`1`$, a leaf driven by the tagged leaf source and $`c_{u_i}\le2`$.

*Unit 2.* $`u_i`$ is silent, so the tagged leaf stays a leaf, and each of its at least $`M`$ rings chooses its only neighbour $`u_i`$. Other vertices can only add leaves at $`u_i`$. At time $`2`$, $`L_{u_i}\ge M+1`$ and $`c_{u_i}\le2`$.

*Unit 3.* Up to $`\tau_i`$ the base is silent, so $`L_{u_i}\ge M+1`$ and $`c_{u_i}\le2`$ just before $`\tau_i`$, and $`L/D\ge(M+1)/(M+3)>M/(M+2)\ge U`$. With leaves ranked first, $`UD<L`$ gives $`\lceil UD\rceil\le L`$, so the ring chooses a leaf of $`u_i`$, and by (R2) this is $`v_i^*`$, re-driven by $`\sigma^{(i)}_1`$, with newborn $`v_i'`$ driven by $`\sigma^{(i)}_2`$. Which leaf it is may depend on arrivals from outside, but that it is a leaf does not. On $`(\tau_i,3]`$ the vertices $`u_i`$, $`v_i^*`$, $`v_i'`$ are silent, so by (D4) the far side at time $`3`$ is the edge $`v_i^*v_i'`$, with the two vertices driven by $`\sigma^{(i)}_1,\sigma^{(i)}_2`$. Now $`L_{u_i}\ge M`$, so $`D_{u_i}\ge M+1`$, and by (D1) this persists.

*Stage 4.* The base is silent, so by (D4) the far side changes only through rings of its own vertices, which are driven by tagged clocks. Induction on $`j`$: at the start of unit $`j`$ the far side is the copy of $`T`$ restricted to vertices $`1,\dots,j+1`$, with vertex $`r`$ driven by $`\sigma^{(i)}_r`$. In unit $`j`$ the only far-side ring is the point $`\beta_j`$ of $`\sigma^{(i)}_{z(j)}`$. Up to it $`D_{z(j)}=d_j`$, since all its neighbours are silent, and with $`a(j)`$ ranked first (R3) and $`U\le1/d_j`$ it chooses $`a(j)`$. By (R3) the newborn is driven by $`\sigma^{(i)}_{j+2}`$, and it is silent until the end of the unit. So the induction hypothesis holds at the start of unit $`j+1`$. After the last unit the far side is the copy of $`(T,v)`$ described.

Finally $`v_i^*`$ was a leaf just before $`\tau_i>0`$, so it was not core at time $`0`$ (D1), and every other far-side vertex was born at or after $`\tau_i`$ (the first one at the promoting ring itself), hence strictly after time $`0`$. $`\square`$

### 2.4 The block contributions

For site $`i`$ let $`E_i`$ be the event that the source of $`u_i`$ has no point in $`(h,h+s]`$ with $`U\le1/M`$. Let $`\mathrm{Acc}_i`$ be the event that every ring of vertex $`1`$ (the root $`v`$) of $`R_i`$ in $`[0,s]`$ has $`A\le D/(D+1)`$, with $`D`$ the degree of the root in $`R_i`$ just before the ring. Put

```math
Z_i=\mathbf 1_{E_i}\mathbf 1_{\mathrm{Acc}_i}\,C_{R_i}(s),\qquad Y_i=\mathbf 1_{\mathrm{Prep}_i}Z_i .
```

**Lemma 2.4 (independent copies).** The pairs $`(\mathrm{Prep}_i,Z_i)`$, $`1\le i\le m`$, are independent, $`\mathrm{Prep}_i`$ is independent of $`Z_i`$, and the $`Z_i`$ have a common law, which does not depend on $`T_0`$. Hence $`Y_1,\dots,Y_m`$ are independent and identically distributed, with

```math
\mathbb E Y_i=\hat p\,Q_{T,v}(s)\;\ge\;p(M,k)\,e^{-\eta s}\,\frac{F(s)-k-\theta_k(s)}{k}.
```

*Proof.* $`Z_i`$ is a function of the points of the source of $`u_i`$ in $`(h,h+s]`$ and of $`\sigma^{(i)}_1,\dots,\sigma^{(i)}_k`$ and $`(\zeta^{(i)}_r)`$ after time $`h`$, including their mark coordinates $`U`$ and $`A`$. $`\mathrm{Prep}_i`$ is a function of reserved randomness of site $`i`$ on $`(0,h]`$ and of coins. These families are independent, and they are disjoint for different sites. The law of $`Z_i`$: the source of $`u_i`$ gives $`\Pr(E_i)=e^{-s/M}`$, independently of the tagged clocks; $`R_i`$ is a friend-tree process from $`T`$ (E-Lemma 1.1); and, given $`R_i`$, the coordinates $`A`$ at the points of $`\sigma^{(i)}_1`$ are independent uniforms, so $`\Pr(\mathrm{Acc}_i\mid R_i)=\prod D/(D+1)`$. Hence $`\mathbb E Z_i=e^{-s/M}\mathbb E[C_T(s)\prod D_v/(D_v+1)]=Q_{T,v}(s)`$. The lower bound follows from $`\hat p\ge p(M,k)`$, E-Corollary 5.3 and $`\mathbb E C_T(s)\ge F(s)-k`$ (E-Lemma 6.1). $`\square`$

**Lemma 2.5 (the contributions are realised).** Pointwise, $`\sum_{i=1}^mY_i\le C(\Lambda)`$.

*Proof.* Let $`Y_i>0`$. By Lemma 2.3 site $`i`$ is coupled, and the window at site $`i`$ runs by the rules of E-Section 4 with reference $`R_i`$, whose vertex $`j`$ is the actual far-side vertex driven by $`\sigma^{(i)}_j`$. On $`E_i\cap\mathrm{Acc}_i`$, E-Lemma 4.1(b) applies: since $`D_{u_i}\ge M+1`$, no ring of $`u_i`$ in the window chooses $`v_i^*`$, no ring of $`v_i^*`$ chooses $`u_i`$, and so no triggering ring occurs in the window. Its proof uses only (D4) and the rules at site $`i`$. Rings at other sites or outside all sites can attach leaves to $`u_i`$, which only increases $`D_{u_i}`$, but by (D4) they cannot reach the far side of $`u_iv_i^*`$. So at time $`\Lambda`$ the far side of $`u_iv_i^*`$ equals the reference tree of $`R_i`$ at time $`s`$, every reference core vertex is an actual core vertex there, and the far side contains at least $`Y_i`$ core vertices.

The bases $`u_i`$ are distinct vertices that are core at time $`0`$, and by Lemma 2.3 each far side contains no vertex that was core at time $`0`$. By E-Lemma 3.2 these far sides are pairwise disjoint at time $`\Lambda`$. Sum over $`i`$. $`\square`$

A failing site may push leaves onto neighbouring old vertices during the block. All conditions in the flags are stated about reserved randomness, and Lemma 2.3 shows that such arrivals cannot spoil a site whose flags hold.

### 2.5 The block estimate

$`\mathbb E Y_i\le\mathbb E N_T(s)<\infty`$. Choose $`R_s<\infty`$ with $`\mathbb E\min(Y_i,R_s)\ge\mathbb E Y_i/2`$, which is possible by monotone convergence. Put $`\beta_s=\mathbb E Y_i/24`$.

**Proposition 2.6 (one block).** For every finite tree $`T_0`$ with $`C(T_0)\ge1`$ core vertices, the friend-tree process from $`T_0`$ satisfies

```math
\Pr_{T_0}\big(C(\Lambda)<\beta_s\,C(T_0)\big)\le\frac{96R_s^2}{C(T_0)\,(\mathbb E Y_1)^2}.
```

*Proof.* Use the construction of Section 2.2. By Lemma 2.1 it has the friend-tree law, and the event on the left is an event about the tree path, so it suffices to bound it in this construction. Write $`\bar Y_i=\min(Y_i,R_s)`$ and $`S=\sum_{i\le m}\bar Y_i`$. By Lemma 2.4 the $`\bar Y_i`$ are i.i.d. with mean at least $`\mathbb E Y_1/2`$ and variance at most $`R_s^2`$. By Chebyshev,

```math
\Pr\Big(S<\frac{m\,\mathbb E Y_1}{4}\Big)\le\Pr\Big(|S-\mathbb E S|>\frac{m\,\mathbb E Y_1}{4}\Big)\le\frac{mR_s^2}{(m\,\mathbb E Y_1/4)^2}=\frac{16R_s^2}{m(\mathbb E Y_1)^2}.
```

By Lemma 1.1, $`m\ge C(T_0)/6`$, so $`m\,\mathbb E Y_1/4\ge\beta_sC(T_0)`$ and the right side is at most $`96R_s^2/(C(T_0)(\mathbb E Y_1)^2)`$. By Lemma 2.5, $`C(\Lambda)\ge S`$. $`\square`$

## 3. The core count tends to infinity

**Lemma 3.1.** Almost surely $`C_n\to\infty`$ and $`C(t)\to\infty`$.

*Proof.* At size $`n\ge3`$ the discrete chain promotes a leaf, and so increases $`C`$ by one, with conditional probability $`\lambda(T_n)/n\ge1/n`$, where $`\lambda\ge1`$ is the promotion rate of E-Lemma 7.1. Since $`\sum_n1/n=\infty`$, Lévy's conditional Borel–Cantelli lemma gives infinitely many promotions almost surely, so $`C_n\to\infty`$. In continuous time $`C(t)=C_{N(t)}`$ (E-Lemma 1.3) and $`N(t)\to\infty`$, because every ring time $`\sigma_m`$ is finite almost surely. $`\square`$

## 4. The lower bound

**Proposition 4.1.** For fixed $`M,\epsilon,k,s`$ as in Section 2 with $`\beta_s>1`$, almost surely

```math
\liminf_{t\to\infty}\frac{\log C(t)}{t}\ge\frac{\log\beta_s}{\Lambda}.
```

*Proof.* Work with the process from $`T_2`$ in the standard construction, with filtration $`(\mathcal F_t)`$. Let $`t_j=j\Lambda`$ and $`I_j=\mathbf 1\{C(t_{j+1})<\beta_sC(t_j)\}`$. Put $`g(c)=\min\big(1,96R_s^2/(c(\mathbb E Y_1)^2)\big)`$ for $`c\ge1`$ and $`g(0)=1`$. By E-Lemma 1.2 and Proposition 2.6 applied with $`T_0=T(t_j)`$,

```math
\mathbb E[I_j\mid\mathcal F_{t_j}]\le g(C(t_j)).
```

By Lemma 3.1, $`g(C(t_j))\to0`$ almost surely, so $`J^{-1}\sum_{j<J}\mathbb E[I_j\mid\mathcal F_{t_j}]\to0`$ almost surely.

The process $`\sum_{j<J}(I_j-\mathbb E[I_j\mid\mathcal F_{t_j}])/(j+1)`$ is a martingale with respect to $`(\mathcal F_{t_J})`$ with increments bounded by $`1/(j+1)`$ in absolute value. It is bounded in $`L^2`$, since $`\sum_j(j+1)^{-2}<\infty`$, so it converges almost surely. By Kronecker's lemma, $`J^{-1}\sum_{j<J}(I_j-\mathbb E[I_j\mid\mathcal F_{t_j}])\to0`$. Hence the fraction of failed blocks $`J^{-1}\sum_{j<J}I_j\to0`$ almost surely.

The first ring happens at an almost surely finite time $`\tau_1`$, after which $`C\ge1`$. Let $`j_0`$ be the first grid index with $`t_{j_0}\ge\tau_1`$; it is finite almost surely, and $`C(t_{j_0})\ge1`$. By (D1) a failed block does not decrease $`C`$, and a successful block multiplies it by at least $`\beta_s>1`$. So for $`J>j_0`$

```math
\log C(t_J)\ge\Big(J-j_0-\sum_{j<J}I_j\Big)\log\beta_s ,
```

and $`\liminf_J \log C(t_J)/(J\Lambda)\ge\log\beta_s/\Lambda`$. For $`t\in[t_J,t_{J+1})`$, $`C(t)\ge C(t_J)`$ and $`t<(J+1)\Lambda`$, which gives the claim. $`\square`$

**Corollary 4.2.** Almost surely $`\liminf_{t\to\infty}\log C(t)/t\ge\mu`$.

*Proof.* If $`\mu=0`$ this holds because $`C(t)\ge1`$ after the first ring. Let $`\mu>0`$, and take an integer $`M\ge2`$ and a rational $`\epsilon\in(0,1)`$ with $`\eta=\eta(M,\epsilon)<\mu`$, and the corresponding $`k`$. By Lemma 2.4 and E-Theorem 1 in its continuous-time form ($`\log F(s)/s\to\mu`$),

```math
\frac{\log\beta_s}{\Lambda}\;\ge\;\frac{\log\frac{p(M,k)}{24k}-\eta s+\log\big(F(s)-k-\theta_k(s)\big)}{h+s}\;\longrightarrow\;\mu-\eta>0\qquad(s\to\infty).
```

So $`\beta_s>1`$ for all large integers $`s`$, and Proposition 4.1 along those $`s`$, a countable family of almost sure events, gives $`\liminf\log C(t)/t\ge\mu-\eta`$ almost surely. Take the intersection over the countably many pairs of an integer $`M\ge2`$ and a rational $`\epsilon\in(0,1)`$ with $`\eta(M,\epsilon)<\mu`$; since $`\eta(M,\epsilon)\to0`$ as $`M\to\infty`$ and $`\epsilon\downarrow0`$, $`\liminf\log C(t)/t\ge\mu`$ almost surely. $`\square`$

## 5. The upper bound

**Lemma 5.1.** Almost surely $`\limsup_{t\to\infty}\log C(t)/t\le\mu`$.

*Proof.* Fix a rational $`\delta>0`$. Since $`\log F(j)/j\to\mu`$, $`F(j)\le e^{(\mu+\delta/2)j}`$ for large integers $`j`$. By Markov's inequality, $`\Pr(C(j)>e^{(\mu+\delta)j})\le F(j)e^{-(\mu+\delta)j}\le e^{-\delta j/2}`$ for large $`j`$, which is summable. By Borel–Cantelli, almost surely $`C(j)\le e^{(\mu+\delta)j}`$ for all large $`j`$. For $`t\in[j,j+1)`$, monotonicity gives $`C(t)\le C(j+1)\le e^{(\mu+\delta)(j+1)}`$, so $`\limsup\log C(t)/t\le\mu+\delta`$. Intersect over rational $`\delta`$. $`\square`$

So $`\log C(t)/t\to\mu`$ almost surely.

## 6. Return to discrete size

Let $`\tau_n`$ be the time at which the continuous-time process from $`T_2`$ reaches size $`n`$. By E-Lemma 1.1, $`\tau_n=\sum_{j=2}^{n-1}e_j/j`$ with independent standard exponential variables $`e_j`$. The centred sum $`\sum_{j\ge2}(e_j-1)/j`$ has independent mean-zero terms with summable variances $`1/j^2`$, so it converges almost surely (Kolmogorov). Since $`\sum_{j=2}^{n-1}1/j=\log n+O(1)`$, we get $`\tau_n/\log n\to1`$ almost surely, and in particular $`\tau_n\to\infty`$.

By E-Lemma 1.3, $`C(\tau_n)`$ is the core count of the jump chain at size $`n`$, and the jump chain has the law of $`(T_n)_{n\ge2}`$. On the almost sure event where $`\log C(t)/t\to\mu`$ and $`\tau_n/\log n\to1`$,

```math
\frac{\log C(\tau_n)}{\log n}=\frac{\log C(\tau_n)}{\tau_n}\cdot\frac{\tau_n}{\log n}\longrightarrow\mu .
```

The statement "$`\log C_n/\log n\to\mu`$" is an event about the law of the sequence $`(C_n)_{n\ge2}`$, so it holds almost surely for the discrete chain. This proves Theorem 2. $`\blacksquare`$

## 7. Audit record

**Where independence comes from.** Every preparation flag and every block contribution is defined as an event or a variable about the *reserved randomness* of one site: the sources of its reserved vertices, its tagged clocks $`\sigma^{(i)}_j`$, its reference newborn sources $`\zeta^{(i)}_r`$ and its thinning coins. These families are disjoint for different sites, so the product law (Lemma 2.2) and the i.i.d. property (Lemma 2.4) are statements about independent random variables, proved without reference to the tree. The tree enters only in Lemmas 2.3 and 2.5, which show that the flags force the intended configuration whatever the rest of the tree does.

**Random tagged clocks.** The leaf $`v_i^*`$ promoted in unit 3 is a random vertex; which leaf it is may depend on arrivals at $`u_i`$ from other sites. The construction re-drives it, at the promoting ring, by the tagged clock $`\sigma^{(i)}_1`$, and gives its newborn $`\sigma^{(i)}_2`$ (rule (R2)). Re-driving immediately after a processed ring is admissible (E-Section 1.1), and it never decides whether the promoting point itself is a ring. The re-driving happens only when the chosen vertex is a leaf of $`u_i`$, so it never touches the reserved vertices of another site. The flag $`F^{(i)}_3`$ is stated in terms of the source of $`u_i`$ and the tagged clocks, so it does not depend on which leaf was chosen. The same device is used for every vertex born on the branch in stage 4 (rule (R3)).

**Reference and switch.** As in E-Section 4, the reference $`R_i`$ keeps its own drivers on the whole window and never reads the actual process; it is defined at every site, coupled or not. Only the actual process switches, after the triggering ring has made its prescribed choice. On $`E_i\cap\mathrm{Acc}_i`$ no triggering ring occurs.

**Filtrations.** Proposition 2.6 is a statement about the friend-tree law from an arbitrary deterministic tree. Proposition 4.1 uses it through the Markov property at the deterministic times $`t_j`$ (E-Lemma 1.2), in the filtration of the standard construction. No block construction has to be carried over from one block to the next.

**Differences from the short version.** No change alters a constant or the structure of the argument.

- (1) The number of maximal type-B paths is bounded by $`2A-3`$; the short version says $`2A-1`$. This is harmless.
- (2) All sites plant the same rooted tree $`(T,v)`$, so the stage-4 flags have the same probability at every site without the minimum over rooted trees mentioned in the short version. Thinning coins are used in stage 1 at A-sites, in stage 3 and for the newborn silence in stage 4, which makes every flag probability exact.
- (3) The block contribution is $`Y_i=\mathbf 1_{\mathrm{Prep}_i}\mathbf 1_{E_i}\mathbf 1_{\mathrm{Acc}_i}C_{R_i}(s)`$. The short version also requires the reference history to be good and assigned to $`v`$. Dropping that restriction does not decrease $`Y_i`$, and the lower bound on $`\mathbb E Y_i`$ is the same.
- (4) The short version says "The existence and positivity of $`\mu`$ give the lower bound $`\mu`$". Positivity is not needed and is not proved; Corollary 4.2 treats $`\mu=0`$ separately.
- (5) The tagged clocks, the rules (R1)–(R4) and the flags are written out explicitly, and Lemma 2.1 checks admissibility.

**Corrections made to the expanded draft before release.** An independent cold read of the first expanded draft of this file, which was never released, found two misstatements. The "before" wording below quotes that draft. Both were in the new write-up, not in the short version, and both are corrected here.

- *Stage-4 silence in $`F^{(i)}_4`$.* Before: the draft required silence of "the sources $`\sigma^{(i)}_r`$ with $`r\le j+1`$, $`r\ne z(j)`$". Since the tagged clocks are indexed from $`0`$, that wording also required silence of the auxiliary clock $`\sigma^{(i)}_0`$, which would give $`j`$ silent clocks plus one per unit and an extra factor $`e^{-(k-2)}`$ in the probability, contradicting the exact value $`\hat p`$ of Lemma 2.2. After: the condition reads $`1\le r\le j+1`$, $`r\ne z(j)`$, the indices of the branch vertices present in unit $`j`$, which is what the script needs and what Lemma 2.2 computes. The formula for $`\hat p_4`$ and the value of $`\hat p`$ are unchanged. Nothing is weakened: removing the extra silence condition enlarges the flag event, and Lemma 2.3 still holds on it. In unit $`j`$ its proof uses only that the designated clock $`\sigma^{(i)}_{z(j)}`$ rings exactly once with the prescribed mark, that the base and the other existing branch clocks $`\sigma^{(i)}_r`$, $`1\le r\le j+1`$, $`r\ne z(j)`$, stay silent, and that the newborn clock $`\sigma^{(i)}_{j+2}`$ stays silent from the birth to the end of the unit. None of this involves $`\sigma^{(i)}_0`$. At a triple, $`\sigma^{(i)}_0`$ may still drive its leaf and ring after unit 2, but those points are irrelevant to the flags and to what they guarantee. The corrected $`F^{(i)}_4`$ has probability exactly $`\hat p_4`$, and the full preparation event $`\mathrm{Prep}_i`$ has probability exactly $`\hat p`$.
- *Start of the block count in Proposition 4.1.* Before: the draft set $`j_0=1+\lceil1/\Lambda\rceil`$ and called "the first ring happens before $`t_{j_0}`$" an almost sure event. That event has probability $`1-e^{-2t_{j_0}}<1`$. After: $`j_0`$ is the first grid index after the first ring time $`\tau_1`$, which is finite almost surely, so $`j_0`$ is an almost surely finite random index with $`C(t_{j_0})\ge1`$. No deterministic time is claimed to follow the first birth almost surely. Nothing is weakened: the displayed lower bound on $`\log C(t_J)`$ holds for $`J>j_0`$ with a random but finite $`j_0`$, and dividing by $`J\Lambda`$ removes it in the limit.

Neither correction changes a hypothesis, a constant or the statement of any lemma, proposition or theorem.

**Real gaps found.** None.
