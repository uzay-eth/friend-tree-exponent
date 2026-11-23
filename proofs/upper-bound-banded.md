> Produced by GPT-6.2 Astra.

# Degree cutoffs depending on both ends of an edge

The target upper exponent is $`\lambda=5657/9983`$. The exact difference $`17/30-\lambda=1/299490`$ makes this a test of the proposed fraction. From a clean checkout run

```sh
bash verify/run_all.sh
```

The banded suite compiles its C++ verifier with undefined-behaviour checks, decompresses the integer recipe and constructs the weights. The SHA-256 digest of `recipe.txt.gz` is `584bb93c2797138ae3de54cb85c7d8e487bb7d7cc5320c021b5a64ec0cf0f60a`.

This is an upper drift certificate for the original random friend tree. Write $`D=L+c`$ for the degree of the focal core vertex, where $`L`$ and $`c`$ count its leaf and core neighbours. Use the tracked neighbours, one distinguished distance-two witness, saturated core-degree labels, and directed transfers from [the preceding certificate](bounds-transfer.md). The following changes concern which information is retained, and how terms in the generator are summed.

## States at each degree

A row of the cutoff file specifies an exclusive upper endpoint of a band of central degrees, followed by cutoffs $`h_1,h_2,\ldots`$ for observing neighbour degrees. Missing entries are 2. The first band begins at $`D=2`$, and the last ends at $`H`$. The cutoffs are nonincreasing in $`c`$ and on passing to the next band. The supplied cutoffs end at central degrees 16, 32, 80, 200, 500, 1,000 and 2,400. Every band uses a cutoff of 620 for a single core neighbour. [Cutoff table](../verify/banded-upper/geometry.txt) The core-degree caps are 20 for $`c=1`$, 7 for $`c=2`$, and 5 otherwise.

For $`d<h_c`$, the neighbour label gives the true degree $`d`$ and its core degree capped as above. At $`c=1`$, when $`d\le14`$ and the core degree is 2, 3 or 4, it also describes a distinguished other core neighbour of that vertex: its exact degree $`e`$ and core degree $`j`$ for $`e<6`$, or $`e\ge6`$. A gain of the near vertex preserves the identity of the distinguished neighbour as long as the new state uses a witness. A conversion of the near vertex makes the fresh core vertex $`(2,1)`$ the new witness when one is required. On the focal vertex's own birth choose the minimum recorded label from the creator's previous core neighbours when a witness is required. In that case the creator had between one and three core neighbours, and all corresponding cutoffs are at least 6, so this identifies the needed exact type or certifies $`e\ge6`$. No tie between indistinguishable labels affects the recorded data. The first core in the three-vertex tree has no core neighbours.

All labels are assigned along the growth history. For ties choose the earliest-born neighbour. An isolated centre starts with an empty label list. The rules below therefore define a potential on each realised history, and the conditional drift is taken in the filtration of that history.

A pooled label gives a lower bound $`b`$ on the neighbour's degree. The marks in a band are $`h_c`$, powers of two greater than $`h_c`$ and less than $`H`$, starting at 4, and $`H`$. Only the largest pooled mark is kept. Every other pooled mark resets to $`h_c`$. Thus an observed neighbour of degree $`d`$ supplies exactly $`1/d`$ to the focal gain rate, and a pooled one supplies an unknown amount in $`[0,1/b]`$. When a central conversion raises $`c`$, or a central gain crosses a band boundary, project known degrees and core degrees down to the new cutoffs. A previous pooled mark projects to the greatest permitted mark at most its old lower bound. For a gain crossing bands at $`c=1`$, preserve the witness. A conversion also inserts the new $`(2,1)`$ neighbour. Marks on an already pooled neighbour need no update when its actual degree grows.

The centre gains a leaf at rate $`L+\sum_{w}1/d_w`$ and converts a leaf at rate $`L/D`$. A conversion creates one new core vertex. Given an exact neighbour type $`(d,k)`$, the neighbour's conversion rate is $`(d-k)/d`$, and its gain rate lies in $`[d-k+1/D,d-k+1/D+(k-1)/2]`$. At a saturated $`k`$, the latter lower bound is $`1/D`$, and conversions preserve the label. For a near neighbour with a witness, replace one of the $`k-1`$ unknown terms in its gain rate by the known $`1/e`$ or the interval $`[0,1/6]`$. An exact witness gains at rate in $`[e-j+1/d,e-j+1/d+(j-1)/2]`$ and converts at $`(e-j)/e`$. All rates here are $`n`$ times the one-step probability. For a fixed centre, a birth moves either its own $`(L,c)`$, the type of one near neighbour, or the type of the distinguished far neighbour. A leaf promotion changes core degree at its creator and creates a new core vertex. A gain at an already pooled neighbour changes none of its stored marks. The tree has no cycle, so the far witness is distinct from the centre and its near neighbour. These cases exhaust the local transitions. Summing one conditional expectation for the change of each old weight and adding each newborn weight once gives the drift of their sum, including steps which change several old weights.

## The generator on a labelled history

For each old core vertex $`x`$ let $`a_x=L_x+\sum_{w\sim_K x}1/D_w`$ and $`b_x=L_x/D_x`$. A gain at $`x`$ attaches a new leaf directly to $`x`$. A conversion at $`x`$ attaches a leaf to one of its leaf neighbours, making that neighbour a new core vertex. Number the vertices in birth order. Ties between leaves in a conversion may be resolved by choosing a particular leaf: all yield the same new label multiset and newborn weight, and the conditional choice in the actual process is uniform. At size three the centre has $`c=0`$. The same definitions then give $`a_x=2`$ and $`b_x=1`$. At every subsequent step exactly one of these events occurs. Indeed, its attachment target is either a core vertex or a leaf with one core neighbour. The sum of their rates is $`n`$, since each ringing vertex sends total rate one.

Write $`z_u`$ for the old weight. Let $`\Delta_x^+ z_u`$ and $`\Delta_x^b z_u`$ denote its changes at gain and conversion at $`x`$. Let $`z_{\rm new}(x)`$ be the new weight at conversion. Use $`S_u`$ for the set of individually observed core neighbours of $`u`$. Use $`J_u`$ for the singleton containing its witness when the witness has exact degree below six, and the empty set otherwise. For a finite focal vertex set $`I_u=\{u\}\cup S_u\cup J_u`$. Above the central cutoff set $`I_u=\{u\}`$. For each fixed $`u`$ this is a disjoint union. The witness is at distance two and the tree has a unique path from it to $`u`$. With $`\Phi=\sum_u z_u`$ the exact identity is

```math
n\,\mathbb E[\Delta\Phi\mid\mathcal F_n]-\lambda\Phi
=\sum_{u\in K}\left\{\sum_{x\in I_u}
  \big(a_x\Delta_x^+z_u+b_x\Delta_x^b z_u\big)
  + b_u z_{\rm new}(u)-\lambda z_u\right\}.
```

To verify the identity, expand each event's change of the sum by the old vertices whose weights it changes, then interchange the two finite sums. An event away from the centre can change a stored exact label only at the near neighbour itself or its distinguished witness. A pooled mark records a past lower bound and remains fixed on a gain or conversion there. A witness already of degree at least six keeps that label. Renewing a witness on a near conversion is counted in the near transition, since that is the vertex at which the conversion occurs. Projecting a whole bag under a central transition is likewise one change of the focal weight, even when many labels are dropped. The new weight has no old summand and is counted in exactly the row of its creator $`u`$. This accounts for every event-weight pair, including simultaneous changes in different rows.

## Two exact cancellations

Let $`g(D,c;d,k)`$ be a signed integer in weight units. It can be nonzero only if both ends know the true degree and core degree of the opposite endpoint under the cutoffs for their own central degrees. On a focal gain include $`\sum_w g(d_w,k_w;D,c)`$ in the change of its weight; on a gain of a tracked near neighbour $`w`$ subtract $`g(D,c;d_w,k_w)`$ from the corresponding change of the focal weight. Each added term cancels on summation against its subtraction, with exactly the same gain rate. Values on witness variants of a given pair coincide. All transfer pairs are zero when either endpoint has degree at least $`H`$.

There is also a nonnegative table $`a(D,c)`$ of finite support, and a threshold $`T`$. For each ordered core edge from $`w`$ to $`u`$ with $`d_w\ge T`$, subtract $`a(D_u,c_u)/d_w`$ from the row for $`u`$ and add the same quantity to the row for $`w`$. The rate $`1/d_w`$ is the rate at which $`w`$ chooses $`u`$. This cancellation holds on every realised tree, independently of the local states. Let $`U`$ be the change on central gain after including the signed transfers. For an exact neighbour the adjusted contribution is $`(U-a(D,c)\mathbf1_{d_w\ge T})/d_w`$. For a pooled neighbour with mark $`b`$ it is bounded above by $`\max(0,U-a(D,c)\mathbf1_{b\ge T})/b`$. When $`b<T`$ omitting a possible subtraction is safe because $`a\ge0`$. The central term from leaf neighbours is $`LU`$. At a sender of degree $`D\ge T`$, add $`(1/D)\sum_{w}a(d_w,k_w)`$, using the maximum over all true types consistent with its label for each term. Set $`a=0`$ outside its tabulated support.

Here is the local inequality after these cancellations. Write $`\rho_u=\sum_{w\sim_K u}1/D_w`$, $`g_{uw}=g(D_u,c_u;D_w,c_w)`$ and $`a(u)=a(D_u,c_u)`$ for the arrow table. Write $`U_u=\Delta_u^+ z_u+\sum_{w\sim_K u}g_{wu}`$. At a finite vertex the exact amended row before subtracting $`\lambda z_u`$ is

```math
\begin{aligned}
& L_u U_u+
  \sum_{w\sim_K u}\frac{U_u-a(u){\bf1}_{D_w\ge T}}{D_w}
  +b_u\big(\Delta_u^b z_u+z_{\rm new}(u)\big)\\
&+\sum_{w\in S_u}\left[a_w(\Delta_w^+ z_u-g_{uw})
                              +b_w\Delta_w^b z_u\right]
  +\sum_{x\in J_u}(a_x\Delta_x^+z_u+b_x\Delta_x^b z_u)
  +\frac{{\bf1}_{D_u\ge T}}{D_u}\sum_{w\sim_K u}a(w).
\end{aligned}
```

The coefficient of $`g_{wu}`$ added at $`u`$ is $`L_u+\rho_u=a_u`$, exactly the coefficient subtracted at $`w`$. A nonzero $`g_{wu}`$ requires that both endpoints be observed without saturation, so both lookups exist, even when their labels retain different witnesses. The two occurrences of $`a(u)/D_w`$ have the same threshold on the same directed edge. These identities also hold when a recipient or sender lies above the central cutoff: the signed transfer then vanishes, the recipient has zero arrow value if it lies in the tail, and a tail sender still pays.

For the programme's finite bound replace the newborn by its trial upper reward. For every remaining random rate in an interval $`[l,h]`$ multiplying a signed difference $`x`$, use $`l x`$ when $`x<0`$ and $`h x`$ otherwise. This bounds each term separately, even for dependent rates. Equal labels give the same transition target; multiplicity multiplies both endpoints. The intervals for near gains follow from summing the contributions of leaves (one each) and those of core neighbours (at most one half each). With exact core degree $`k`$ the lower endpoint is $`d-k+1/D`$. With a capped label $`k`$ and true core degree $`k'\ge k`$, the valid endpoints are $`1/D`$ and $`d-(k+1)/2+1/D`$, since the upper expression decreases with $`k'`$. At saturation a conversion preserves the recorded core degree. With a witness, $`k`$ is exact: the $`k-2`$ other contributions lie in $`[0,(k-2)/2]`$. Use the exact contribution $`1/e`$ for a known witness or $`[0,1/6]`$ for a pooled witness. For an exact far witness, its near neighbour contributes $`1/d`$ and each of its $`j-1`$ other core neighbours contributes at most one half. Conversion at an exact vertex changes its core degree but leaves its degree fixed. Its rate is its leaf fraction; if it has no leaves its rate is zero. These cases give all the terms in the displayed row.

## Boundary at large degree

Beyond $`D\ge H`$, take $`f=A+Bq`$, where $`q=L/D`$ and $`A,B`$ are nonnegative in the same weight units. In a gain $`q`$ increases by $`(1-q)/(D+1)`$. In a conversion it decreases by $`1/D`$. Consequently
```math
Gq=(L+\rho_u)(1-q)/(D+1)-q/D\le1-q,
```
because $`L+\rho_u\le L+c/2\le D`$. The new core vertex sees a single pooled creator marked $`H`$. Let $`s`$ bound the weight of that newborn, and let $`a_*`$ be the largest arrow table entry. Signed transfers vanish here and the arrow charge at the sender is at most $`a_* c/D=a_*(1-q)`$. With $`a(D,c)=0`$ in this region, the adjusted drift including the newborn is bounded by $`(B+a_*)(1-q)+sq`$. Checking $`B+a_*\le\lambda A`$ and $`s\le\lambda(A+B)`$ is sufficient for all these vertices. A finite row which grows into the tail uses the upper rounding of $`A+B(L+1)/H`$.

## Checking the finite certificate

The recipe consists of integer trial newborn rewards in order of the $`c=1`$ rows, the signed edge table, and the nonnegative arrow table. The cutoff geometry and weight scale are specified. All other weights can be constructed in decreasing order of $`D`$, then $`c`$, then sorted local labels in reverse lexicographic order. In this order every non-newborn successor of a finite row is already known. For a fixed row propose a positive integer weight $`z`$. Take an upper bound on every signed rate times its change using the appropriate endpoint, round each term up, and verify that their sum is at most $`\lfloor\lambda z\rfloor`$. The conditional drift is nonincreasing in the trial row weight $`z`$ before subtracting $`\lambda z`$. Raising a failed proposal and checking it again therefore preserves the direction of the certificate. The trial reward is used only as a positive-rate upper bound for the newborn. The checker must verify that the actual newborn weights in the completed table are at most the trial rewards. Any floating-point proposal for $`z`$ can affect only how quickly an integer meeting the checks is found. The integer inequalities and positivity suffice for the result.

More explicitly, a transition term has target $`v`$ and rational rate interval $`[r_-,r_+]`$. Its contribution at row weight $`z`$ is bounded by $`\lceil\max\{r_-(v-z),r_+(v-z)\}\rceil`$. For the central promotion the target is the post-conversion weight plus the trial reward of the child. For a central gain the target includes the incoming gain transfers, split across the leaf and individual core-neighbour rates. For a gain of an observed neighbour the target includes the opposite signed transfer. The redirection charge at the sender is a separate constant. Summing these integer upper bounds and checking the result against $`\lfloor\lambda z\rfloor`$ is sufficient. For a fixed central degree and core degree, each stored near or far transition strictly advances the changed label: gain increases exact degree or enters the pool, conversion increases an uncapped core degree, and a far transition advances the witness label within its fixed near type. Sorting the bag therefore advances its lexicographic index. If an exact neighbour enters the pool it begins at the first pool mark, so this transition never displaces a larger pooled mark. Reprojection, which can lower labels, occurs only on a central gain or conversion, whose successor already precedes this row in the calculation. Unchanged labels have zero increment. The checker also checks the degree tag whenever it reads a previously computed weight.

With $`m>0`$ a global lower bound on the weights, $`mC_n\le\Phi_n`$. The two cancellations and the checked inequalities give $`\mathbb E[\Phi_{n+1}\mid\mathcal F_n]\le(1+\lambda/n)\Phi_n`$. Consequently $`\mathbb E C_n\le(\Phi_3/m)\prod_{j=3}^{n-1}(1+\lambda/j)=O(n^\lambda)`$.

For the supplied input the root weight is $`1539738114735`$ and $`m=3908512529`$ in the same integer units. The relevant lines of the production output are reproduced below.

```text
verified rows 226545617, signed transfers 63309, minimum integer drift slack 3
scale 1000000000000, minimum newborn reserve 4340525, weight range [3908512529, 1764716280715], centre 1539738114735
PASS upper 5657/9983
```

If the exact check at $`\lambda`$ leaves a common positive reserve $`e`$ between the trial and actual newborn weights, it also gives an optional strict improvement. Let $`M`$ bound every weight from above and $`S_n=\sum L/D`$. With actual newborns the drift is at most $`(\lambda\Phi_n-eS_n)/n`$. Since $`\mathbb E[\Delta C_n\mid\mathcal F_n]=S_n/n`$, the shifted potential $`\Phi_n+eC_n`$ has drift at most $`\lambda\Phi_n/n`$. The inequality $`\Phi_n\le MC_n`$ then gives upper exponent $`\lambda M/(M+e)`$. This uses the reserve at the large-degree boundary as well as in the finite rows.

## Correspondence with the integer verifier

Line references below are to the [production C++ source](../verify/banded-upper/check.cpp). Its integer weight unit is $`10^{-12}`$. Dividing every weight and transfer by the same scale leaves all inequalities homogeneous.

| Mathematical operation | Function and lines | What is checked or constructed |
|---|---|---|
| Input and support | [`Recipe`, 35–66](../verify/banded-upper/check.cpp#L35), [`Checker`, 169–206](../verify/banded-upper/check.cpp#L169) | Positive rewards, bounds on signed and arrow entries, zero off the dictionaries, observable signed support on both endpoints, arrow recipient below the cutoff. |
| Exact or pooled label, unordered multiset | [`Layer`, 89–145](../verify/banded-upper/check.cpp#L89) | Enumerate all multisets in sorted order; at most one pooled label can exceed the first pool mark. Canonicalization keeps the largest. The $`14`$ exact witness labels correspond to the $`2+3+4+5`$ possibilities; the last is the pooled witness. |
| Bands and projecting gains and conversions | [`Band`, 151–158](../verify/banded-upper/check.cpp#L151), [`transitions`, 222–239](../verify/banded-upper/check.cpp#L222) | Check decreasing cutoffs within and across bands. Re-encode each remembered neighbour, using its lower bound if pooled. An increase of core degree adds the newborn's label; an increase in degree alone preserves the witness when $`c=1`$. |
| Old near and far targets | [`transitions`, 240–260](../verify/banded-upper/check.cpp#L240) | Replace one representative of each equal-label group in the bag; save its multiplicity. No event at a pooled near neighbour changes the bag. A near conversion resets the witness to the newborn when required. |
| Newborn witness and reward | [`newborn`, 262–270](../verify/banded-upper/check.cpp#L262), [`solve`, 298–301](../verify/banded-upper/check.cpp#L298) | Choose the first old label in the creator's sorted bag if a witness is needed, then index the reward in the initial $`c=1`$ alphabet; this alphabet is common across bands. |
| Central gain including signed and arrow adjustments | [`solve`, 276–297](../verify/banded-upper/check.cpp#L276) | Look up both orientations of the transfers, add the incoming total to the central successor, subtract the receiver's arrow if the observed degree or lower mark reaches the threshold. Use rate $`L`$ for leaves, exact reciprocal degree for tracked neighbours and zero to reciprocal mark for pooled neighbours. |
| Near and far rates | [`solve`, 302–323](../verify/banded-upper/check.cpp#L302) | Subtract the outgoing signed transfer from the near-gain target. Encode the preceding intervals with common integer denominators, including the exact witness or the $`1/6`$ endpoint. |
| Sender's charge | [`receiver_max`, 212–221](../verify/banded-upper/check.cpp#L212), [`transitions`, 240–244](../verify/banded-upper/check.cpp#L240), [`solve`, 323](../verify/banded-upper/check.cpp#L323) | Sum maximum recipient charges consistent with the marks, grouped by multiplicity, divide by sender degree and round up. |
| Tail and boundary | [`Checker`, 199–206](../verify/banded-upper/check.cpp#L199), [`solve`, 285–289](../verify/banded-upper/check.cpp#L285) | Check the two affine endpoint inequalities by integer cross-multiplication and initialize the global extrema with $`A,A+B`$. On a boundary gain use $`A+\lceil B(L+1)/H\rceil`$. |
| Acceptance, minimum and initial value | [`up`, 26–28](../verify/banded-upper/check.cpp#L26), [`solve`, 326–358](../verify/banded-upper/check.cpp#L326), [`run`, 360–377](../verify/banded-upper/check.cpp#L360) | Each signed term is rounded up, even when negative. Require $`z>0`$ and nonnegative integer slack. Check every completed degree-two newborn against its trial reward, take the global minimum and read the degree-two root with its empty bag. |
