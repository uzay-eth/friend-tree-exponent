# Exact audit: local eigenfunctions and response transcendence

**Pinned claim commit:** `dd711c233c59f86d9c7249a84b47fa1526ed7205`.  
**Review and dependencies:** The lead reported a Gemini proof audit without a separately archived verdict. Claude independently checked the scalar response, and GPT-family audits cover the physical trigger and kernel entries. See [Claim 4 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-4.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-4.md).

**Author:** GPT-6.2 Sol Ultrafast; proof and code audited by Gemini 4 Neon.

**Verdict.** The nineteen-mark physical trigger $`T(s)=F(s)`$, the unmarked zero-gauge output $`U_s(J)`$, and each of the 361 permanent-mark kernel entries $`K_{ij}(s)`$ are transcendental over $`\mathbb C(s)`$. The proposed all-diameter eigenfunction obstruction is **open in this audit**: its first-round generator identities and its last step are valid, but the descent and the proposed degree integration have gaps. The examples below refute intermediate inferences, not the theorem itself. Function transcendence in $`s`$ does not decide the arithmetic of the crossing $`\mu`$.

## A. What the local calculation establishes

Write $`A_x=\Delta_x`$ for a leaf-gain difference, $`B_x`$ for a promotion difference (defined when $`L_x>0`$), $`p_x=1/D_x`$, and $`b_x=L_x/D_x=1-c_xp_x`$. At zero leaves the product $`b_xB_x`$ is zero, without defining a nonexistent promotion. The generator is

```math
\mathcal L\Phi=\sum_x (L_x+\sum_{w\sim x}p_w)A_x\Phi+\sum_x b_x B_x\Phi.
```

An induction that covers all distinct vertices, including endpoints, avoids circular use of the core-degree restrictions. Gain is legal at every vertex of a connected core with at least two vertices: a bare core has a sending core neighbour. All the decorated examples here are attainable from $`T_3`$ with the specified founder of core degree at least two, by creating leaves and promoting the appropriate children in order. An unoriented example can be realized by forgetting such a founder.

### Top round and the missing lower-round term

If $`\Phi`$ is a diameter-$`D`$ motif sum, terms responding to two edits at old sites more than $`D`$ apart cannot occur in one motif. For a promotion this includes changed motifs incident to the promoter or containing its new child; the child's distance to the other site is larger by one. Thus $`AA`$, $`AB`$, and $`BB`$ vanish at those separations when the edits are defined. This works with oriented edges and even with a founder flag.

For $`d(u,v)=k\ge1`$, take $`z\sim v`$ away from $`u`$. In $`A_zB_u\mathcal L\Phi`$, promotion preserves $`D_u`$. The gain rate at $`u`$ changes by $`-1+1/2=-1/2`$, and the corresponding rate-change term vanishes by $`A_z A_u\Phi`$ on the promoted tree. The new child's own gains and promotions, and any changes in its rate, are independent of the $`z`$ edit in the needed differences: the child is at distance **$`k+2`$** from $`z`$. This requires allowing $`c_{\rm new}=1`$ and also allowing $`z`$ to be a core endpoint. For other old gains the mixed increment commutes with $`A_zB_u`$; the only rate hit by $`A_z`$ that can see $`B_u\Phi`$ is the send $`z\to v`$. The other neighbours of $`z`$ are at distance $`k+2`$ from $`u`$. For a promotion at $`w\ne z`$, the rate is unaffected by the $`z`$ gain and the remaining mixed increment vanishes by old $`AB`$ on the appropriate trees (including $`w=u`$, with a zero term if no second leaf remains). At $`w=z`$ the change in **its own** promotion rate contributes

```math
(A_z b_z) B_z B_u\Phi=-c_z(A_zp_z)B_zB_u\Phi.
```

If $`L_z=0`$, instead use $`b_z(T^{+z}) B_z B_u\Phi(T^{+z})`$; only the gained state can promote. At the top round this whole term vanishes by diameter locality. In a later round the stated $`P(k)`$ gives no $`BB`$ cancellation. A degree-free diameter-2 score $`C+\sum_{i\sim x\sim j,\,c_x=2}\min(c_i,3)\min(c_j,3)`$ (each unoriented two-edge path counted once) satisfies both parts of $`P(1)`$ identically and $`C\le\Phi\le10C`$ on every tree, but has a nonzero term $`1/6`$ in the lower-round ledger. No eigenfunction claim is made. Merely taking a bare $`z`$ does not fix the inference, because $`A_z`$ creates an eligible leaf. If $`z`$ has a leaf, the surviving equation under $`P(k)`$ ties $`A_vB_u\Phi`$ to $`c_zB_zB_u\Phi`$ rather than setting it to zero. If $`v`$ has $`c_v=1`$, there is no off-path $`z`$ in the first place.

### Gain–gain shift and symmetry

At separation $`k+1`$ between $`y`$ and $`v`$, old $`AA`$ and $`AB`$ suffice for the gain–gain calculation. Intrinsic leaf-rate differences multiply vanished mixed increments. Promotion rates only see their own site, and the crossed increment is old $`AB`$. The only remaining gain terms are the two neighbour sends towards the middle of the path. Even when their receiver is the same ($`k=1`$), the two rate contributions simply add. Multiplying the resulting equation by the appropriate degree factors gives the stated sign-reversing shift with

```math
\mathcal I_T(x_0,\ldots,x_k)=\left(\prod_{i=0}^k D_{x_i}(D_{x_i}+1)\right)A_{x_0}A_{x_k}\Phi(T).
```

This is an identity **within a fixed tree**. It gives no extra neighbour where a path ends; adding one changes decorations and can change the value even if the original path decorations are held fixed. For example, a diameter-2 fork contribution at centre $`x`$ can assign $`D_iD_jc_w`$ to each choice of arms $`i,j`$ and third arm $`w`$. Its $`A_iA_j`$ changes from 1 to 2 as $`c_w`$ changes on an otherwise identically decorated $`i`$–$`x`$–$`j`$ path. Vanishing longer-distance gain differences does not remove such topology or core-degree dependence.

In a genuinely symmetric corridor, orientation-preserving reflection about an LCA exchanges the two equal branches and fixes the stem to the founder. Reflection shows **equality** of the reflected mixed increments; zero requires an **odd** number of sign-reversing valid shifts to a reflected window. An odd-length path cannot itself have an endpoint-swapping founder-preserving reflection (its middle edge would be reversed), but windows can be exchanged across an LCA on an even corridor. E.g., on a symmetric corridor of length $`2k+1`$ the two length-$`k`$ end windows are shifted $`k+1`$ times, which is odd for even $`k`$; a founder-preserving reflection of that odd corridor is impossible. Conversely on an even corridor of length $`2k`$ the end windows require $`k`$ shifts, so this configuration proves zero when $`k`$ is odd. These are statements about the special corridors; they do not prove that an arbitrary motif coefficient can be transported to such a corridor, or that the required symmetries survive passing through or placing a stem near the founder.

### Terminal step, and two suggested repairs

If $`P(0)`$ were established **unconditionally**, both mixed conditions for every pair of distinct sites imply, on every edge $`u\sim v`$,

```math
H_u+H_v=0,\qquad H_x=D_x(D_x+1)A_x^2\Phi.
```

The two nonzero terms before rescaling are the changes in the reciprocal send rates times the receiving site's second gain difference. Terms from the change in $`L_x`$ vanish by $`AA`$, those from $`b_x`$ by $`AB`$ (at zero leaves use the gained tree as above). Since $`A_v^2\Phi`$ is independent of $`D_u`$ by $`AA`$, $`H_u`$ is constant as $`D_u`$ varies on the fixed skeleton. Thus $`A_u^2\Phi=h/[D_u(D_u+1)]`$. Boundedness on a fixed skeleton forces the limiting first difference to be zero, hence $`A_u\Phi=-h/D_u`$. The harmonic sum and boundedness force $`h=0`$. After all ordinary degree dependence is removed, each promotion increment is degree-independent, wherever the promotion is legal. Varying $`D_u`$ with at least one leaf, holding the other degrees fixed, isolates $`b_u`$, which varies for $`c_u>0`$; this forces every promotion increment to vanish. Hence $`\mathcal L\Phi=0`$, contradicting $`s\Phi>0`$. This conditional obstruction holds for any bounded-on-fixed-skeleton positive score with the two $`P(0)`$ properties, without a motif representation or a path-limit argument. It includes diameter zero as a direct consequence.

- **Fixed distance, higher powers.** $`A_yA_z\mathcal L^m\Phi=0`$ at $`d(y,z)=D+1`$ is a valid necessary condition. However one must calculate $`\mathcal L^m`$ on gain and promoted configurations, not extrapolate its far-separation identity to closer pairs. Already at $`m=2`$ the outer expansion contains $`b_y A`$-differences in $`B_y\mathcal L\Phi`$, with the inner promotions evaluated on the new skeleton. For example, the product rule generates $`(A_yb_y) A_z B_y\mathcal L\Phi(T^{+y})`$ (and a remainder $`b_y A_yA_zB_y\mathcal L\Phi`$). At distance $`D+1`$, the inner top-round identity replaces the first factor involving $`\mathcal L`$ by $`(A_zp_z)A_{z'}B_y\Phi`$ at the gained tree. Thus the two-denominator term can contain a *promotion/gain* difference, not just a distance-$`(D-1)`$ gain/gain difference. Also $`A_y b_y=-c_y A_y p_y`$ has the same degree shape as a send difference; distinct products of reciprocals alone are not a coefficient-separation proof when the other factors still depend on the varying degrees. The checker includes a diameter-2 exact iteration at fixed distance 3; it gives a nontrivial second condition, not a sufficient criterion. No closed higher-power descent was established.
- **Diameter-2 integration.** Even removing *all* mixed gain differences is not enough to integrate away diameter-2 topology. The bounded degree-free example above has zero gain differences but a nonzero double promotion difference at distance two. Any diameter-1 motif sum has zero double promotion difference there. A diameter-2 reduction therefore needs additional control of promotion increments, not merely integration in $`D`$ as in the edge argument. Amplitude uniqueness, if invoked with its hypotheses, reduces the obstruction to nonlocality of the unique amplitude; it does not supply a nonzero distant mixed difference. No such amplitude estimate was established here.

## B. Transcendence by infinitely many denominator primes

The earlier isolated-transfer proof detected a terminal star through the reward $`b=L/D`$. One cannot reuse that *terminal* calculation for the trigger: its terminal reward is zero on a large star. One further action (the independent nontriggering direct newborn production) detects the required prime. There is no need to locate meromorphic poles; poles of finite truncations alone would not justify conclusions about the infinite response.

### Formal-to-analytic bridge

Let $`z=1/s`$. For the positive marked generator $`R`$ or $`\widehat R`$, the coefficient of $`z^{n+1}`$ is a generator iterate on its reward. For example one may truncate the positive chronological expansion at $`n`$ actions: from a virgin each row in this fixed prefix has bounded size, all holding parameters are positive, and the weighted sum of one-step numerators is finite. For any fixed length $`m`$, successive first $`m`$ transitions give $`O(s^{-m})`$ from the fixed starting state as $`s\to+\infty`$; the continuation is $`O(s^{-1})`$ on the finite set reachable at that length, by the weighted positive bounds in the kernel/trigger source. Each retained factor has the ordinary expansion $`(s+\alpha)^{-1}=z-\alpha z^2+\cdots`$, equivalently insertion of generator diagonals. This justifies a complete Poincaré expansion, without asserting that it converges. For the intrinsic generator, polynomial size bounds on fixed iterates and Yule moment domination justify repeated stopped Dynkin, or the same argument using a continuation bounded by $`n(S)/(s-1)`$ for a linearly bounded reward. The zero-gauge identity below then justifies its full expansion.

If a rational-coefficient formal series were algebraic over $`\mathbb C(z)`$, finite rational elimination on the coefficient equations of its finitely many monomials descends a relation to $`\mathbb Q(z)`$. Equivalently, take a minimal polynomial $`P\in\mathbb Z[z,X]`$ after descent, and write $`r=\mathop{\mathrm{ord}}\nolimits _z P_X(z,f)`$ with nonzero leading coefficient $`c`$. In characteristic zero $`r`$ is finite. For $`n>r`$, coefficient comparison at $`z^{n+r}`$ gives $`c f_n=-[z^{n+r}]P(z,f_{<n})`$. The denominator primes are therefore confined to those of $`c`$ and the finite initial block. This is the finite-prime consequence of Eisenstein needed here. A polynomial relation satisfied by the analytic function gives the same formal relation on its full asymptotic expansion, so a contradiction proves transcendence of the function.

### Zero gauge

Let $`N`$ be planted size (parent slot excluded) and $`G(s)=(V_s\sum_x b_x)(J)`$ be the total intrinsic promotion-rate resolvent. The established whole-tree transfer in `influence/zero/critical.md` gives $`sM_s(J)(1-U_s(J))=1`$ for large positive $`s`$, and integration by parts gives $`G=sM_s(J)-1`$. Hence $`U_s(J)=G/(1+G)`$, the same Möbius transform as the zero-mark newborn response in `isolated-transfer.md`. This directly provides a second implementation of the signed coefficients.

For a prime $`p\ge5`$, all states before action $`p-2`$ have $`N<p`$. All degrees, including the counted top slot, are bounded by $`N`$; a **nonroot** vertex has degree at most $`N-1`$. Reaching degree $`p`$ in $`p-2`$ actions forces the top star of size $`p`$, using $`p-2`$ consecutive intrinsic top gains. Only that terminal state can have a reward denominator divisible by $`p`$. Its rate product is $`(p-2)!`$ and its reward is $`(p-1)/p`$. At that coefficient the Möbius correction uses only smaller orders and is $`p`$-integral. Therefore

```math
u_{p-1}-(p-2)!\,(p-1)/p\in\mathbb Z_{(p)},\qquad v_p(u_{p-1})=-1.
```

There is also a direct $`C\Delta`$ count: $`U_s= (s-L^0-C\Delta)^{-1}b`$ pointwise/formally, with $`C\Delta f(S)=p_S\sum_Y[f(Y^+)-f(Y)]`$. A child has size at most $`N-1`$, so $`Y^+`$ has size at most $`N`$; a $`C\Delta`$ operation never increases size. Any term with such an action cannot reach size $`p`$ in $`p-2`$ actions. All its source and reward denominators are prime to $`p`$. This handles all orders of feedback, including the shifted child, not just $`V_s C\Delta V_s b`$.

### Physical degree-two trigger

$`T(s)=F(s)=(s-\widehat R)^{-1}h(J_2)`$. The marked generator has intrinsic changes, diagonal $`-(N-1/D+a)`$, shifted child copies at $`1/D`$, unshifted child copies at $`a_Y-1/D`$, and independent direct births to $`J_{\min(D,20)}`$ at $`b`$ only when $`D>2`$. Copies do not increase size; direct births reset size to two. At order $`z^p`$ there are $`p-1`$ actions. A denominator divisible by $`p`$ first becomes possible at the source of the final action, after the unique pure star growth to size $`p`$, with rate product $`(p-2)!`$. For $`p>20`$ all mark denominators are automatically prime to $`p`$. For $`5\le p\le 19`$ none of the smaller source sizes could have generated mark $`1/p`$; a mark created on the final action contributes no denominator because $`h`$ is mark-independent. At this final large star the diagonal leaves $`h=0`$, there are no existing child copies, and intrinsic gains and promotions both leave a nonvirgin top. Only its independent direct birth pays $`h=1/2`$. Thus for **all primes $`p\ge5`$**,

```math
a_p-(p-2)!\,(p-1)/(2p)\in\mathbb Z_{(p)},\qquad v_p(a_p)=-1.
```

For $`p\ge23`$ this is the direct newborn $`J_{20}`$; the smaller test primes land on $`J_p`$ instead. Factorials in these residues are nonzero modulo $`p`$ simply because each factor is less than $`p`$.

### All entries of the permanent nineteen-mark kernel

Here $`K_{ij}=(s-R)^{-1}(b\mathbf1_{\{\min(D,20)=j\}})(J_i)`$; $`R`$ has no direct births. Fix $`p\ge23`$ so all marks are prime to $`p`$. For $`j=20`$ the terminal large star proves valuation $`-1`$ at $`z^{p-1}`$ as above, with arbitrary input mark $`i`$. For $`j<20`$ all nonzero terminal rewards are $`p`$-integral, so we count the largest possible number of $`p`$-singular **action** factors. Before the earliest possible action source with degree $`p`$, at action count $`p-2`$, there are no such factors. Equality in the initial size bound forces the same pure top-star growth with product $`(p-2)!`$.

For $`3\le j<20`$ take $`p+j-3`$ actions, giving at most $`j-1`$ singular factors. To attain that many, every action after the star must be singular. While the top stays degree $`p`$, a singular intrinsic action can promote at the top or send a leaf from top to a child; singular diagonal actions can also retain the old shape. A top gain has integral rate unless a child has degree divisible by $`p`$, which is impossible here: all children start after the first $`p-2`$ actions and there are at most $`j-1\le18`$ actions left. When a copy leaves the top for a child, the new subtree and all its possible degrees in the remaining actions are too small to admit another singular action (a child starts at size 2, and the shift and subsequent increases each cost an action). Thus the copy must be the last action. A shifted copy needs exactly one top promotion and $`j-3`$ leaf gains to the selected child, all via the $`1/p`$ parent send; an unshifted copy requires an additional gain. There is no room for a diagonal, second promotion, or another branch. The unique maximal term gives $`p^{j-1}`$ times the coefficient congruent to $`(p-2)!(p-1)(j-1)/j\not\equiv0\pmod p`$. Its valuation is therefore $`-(j-1)`$ at $`z^{p+j-2}`$. Mixed-rate transitions in this count can be split into their individual arrow rates before counting valuations.

For $`j=2`$ take $`p`$ actions, with at most two singular factors after the star. The only two-factor term paying a degree-two reward is promotion at the star and then its **unshifted** child copy; the shifted copy has degree three. The singular part of the unshifted rate is $`-1/p`$, even though the entire rate $`1/20-1/p`$ is positive. The residue after multiplying by $`p^2`$ is $`-(p-2)!(p-1)/2\ne0\pmod p`$. Thus this entry has valuation $`-2`$ at $`z^{p+1}`$. The finite-prime lemma applies separately to every entry. No full numerical kernel coefficient census at $`p\ge23`$ was attempted.

## Exact replay

Run from this directory. The obstruction checker uses a formal independent symbol for each isomorphism class of decorated motifs, so its assertions compare coefficient vectors, not samples of numeric weights. The transcendence checker uses `Fraction` for each rate and coefficient and accumulates rows in scaled integers over a common denominator to avoid repeated gcds.

`taskset -c 5 nice -n 10 python3 -B verify_obstruction.py`

```text
D=1,2,3: 24 decorated cases, all first-round AA and AB coefficients match term by term
Founder-oriented V symmetry: 3 explicit reflections; reflection alone does not force zero
Controls: size at s=1, constant at s=0, fake constant promotion at s=1, dropped send detected
Lower-round gap witness: surviving (Δ_z b_z) B_z B_u = 1/6
Diameter-2 fork: same decorated path, mixed increments (1, 2)
D=2 fixed distance 3: sample (AA Φ, AA LΦ, AA L²Φ) = ('0', '-1/6', '-1')
```

`taskset -c 5 nice -n 10 python3 -Bu verify_transcendence.py`

```text
Signed generator and independent total-promotion transfer: 14 exact coefficients agree
u_1 = 1/2
u_2 = -1/12
u_3 = 7/36
u_4 = -71/480
u_5 = 18023/259200
u_6 = 17989/54432000
u_7 = -2222211137/45722880000
u_8 = 7876800793/240045120000
u_9 = 427459808447407/4032758016000000
u_10 = -12129870400995004123/37262684067840000000
u_11 = 44851936025213328611837/129115200295065600000000
u_12 = 3929440236628179152093142121/13293701022379954176000000000
u_13 = -13854293232845618199530331495834659/7451917045105307112898560000000000
u_14 = 55079938779596593618586075305207429415467/16112236958244890827224750489600000000000
a_1 = 1/2
a_2 = -1
a_3 = 7/3
a_4 = -767/144
a_5 = 52673/4320
a_6 = -1603691/57600
a_7 = 4608541543/72576000
a_8 = -1653965961821/11430720000
a_9 = 3167626116913423/9601804800000
a_10 = -12167337367986712799/16131032064000000
a_11 = 12931007087589151113887/7452536813568000000
a_12 = -333263360241520696139096803/82633728188841984000000
a_13 = 2362577699901342749909885919176777/248149085751092477952000000000
p=5: v_p(u_4)=v_p(a_5)=-1; both star-subtracted remainders p-integral
p=7: v_p(u_6)=v_p(a_7)=-1; both star-subtracted remainders p-integral
p=11: v_p(u_10)=v_p(a_11)=-1; both star-subtracted remainders p-integral
p=13: v_p(u_12)=v_p(a_13)=-1; both star-subtracted remainders p-integral
u_1..u_14 denominator prime supports: [(2,), (2, 3), (2, 3), (2, 3, 5), (2, 3, 5), (2, 3, 5, 7), (2, 3, 5, 7), (2, 3, 5, 7), (2, 3, 5, 7), (2, 3, 5, 7, 11), (2, 3, 5, 7, 11), (2, 3, 5, 7, 11, 13), (2, 3, 5, 7, 11, 13), (2, 3, 5, 7, 11, 13)]
a_1..a_13 denominator prime supports: [(2,), (), (3,), (2, 3), (2, 3, 5), (2, 3, 5), (2, 3, 5, 7), (2, 3, 5, 7), (2, 3, 5, 7), (2, 3, 5, 7), (2, 3, 5, 7, 11), (2, 3, 5, 7, 11), (2, 3, 5, 7, 11, 13)]
```
