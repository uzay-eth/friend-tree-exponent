# Positive current-parent-mark renewal

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 13: Claude for transfer bookkeeping; GPT for the uncapped bound and root transfer; inputs include physical jump geometry and the Round 1 exponent and prior. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


The construction is valid for the 19-valued cap and for finite, nonincreasing binary gauges with `q(D)≥1/D`. For cap 20 the current and permanent kernels have the same physical critical root, but are not equal as functions of the discount. The proof uses the permanent exponential renewal remainder to obtain a strict weighted growth gap for the current residual.

## 1. Active and dormant rows

Each particle has an unmarked planted shape `S`. Its root has an inert parent slot. An active label `a` is the current parent gauge when a physical summand is read; within an individual branching particle it is held fixed. Intrinsic evolution `L⁰` on the shape stores no other marks. Write `D,L,p=1/D,b=L/D,d=q(D),d'=q(D+1)`, and `r=L+Σ_Y 1/D_Y` for the root's own total leaf-gain rate **with the slot silent**. `Y⁺` means a new leaf at `Y`'s root.

The active particle evolves intrinsically, dies at rate `a`, and emits the following independent particles, retaining itself:

| Particle | Rate |
|---|---:|
| active `Y` labelled `d`, for every old child | `d−p` |
| active `Y⁺` labelled `d`, for every old child | `p` |
| dormant `(d',d;Y)`, for every old child, if `d'<d` | `r` |
| dormant `(q(D_Y+1),q(D_Y);Z)`, for every old grandchild `Z` of `Y`, if the labels differ | `p` |
| direct active `J` labelled `d` | `b` |

The residual generator `R_c` removes only the last row. Equal children and grandchildren are counted with multiplicity. No old-child or grandchild term reads the active particle's own label `a`. A dormant `(l,u;S)`, `l<u`, evolves by `L⁰`, dies at rate `l`, and **switches** to active `(l;S)` at rate `u−l`; no dormant is retained on switching. Its generator is `L⁰−u` on its own component plus `(u−l)` times the active `l` component.

For a mark-independent active observable (with zero dormant reward), or a discounted running reward independent of the root label, the dormant response is exactly `f_l(S)−f_u(S)`. The active response has original-residence representation `f_a=V_{s+a}G` with `V_λ=(λ−L⁰)⁻¹` and common forcing `G` consisting of direct reward and residual emissions. The dormant response is `(u−l)V_{s+u} f_l`, equal to the difference by the intrinsic resolvent identity. The finite-time version requires no convergence at infinity: subtract the active first-residence formulas, using

```
e^(-lt)−e^(-ut) = ∫₀ᵗ (u−l)e^(-uv)e^(-l(t−v)) dv
```

and the Markov property of the shape at time `v`. This holds for the full and the residual process, and for the vector rate observable `B₀g=b g(q(D))`, dormant value initially zero. For finite-time mark-dependent terminal data it generally fails: already at `t=0` the dormant value is zero while the active difference need not be zero.

## 2. Physical generator, including the founder

For each ordinary nonfounder `v`, physical evaluation is `f_{q(D_parent)}(S_v)`. Relative to its own intrinsic generator, `v` receives a leaf from its parent at rate `p_parent` with its own mark unchanged. Charge the difference `p[f_d(Y⁺)−f_d(Y)]` to the parent's row. On any intrinsic root gain, its old children's shapes are unchanged while their labels change from `d` to `d'`. Charge `r[f_{d'}(Y)−f_d(Y)]` to that root. The missing gain at `Y` from its external parent occurs at rate `p` and changes all old grandchildren's labels; charge their differences to the sending root. There is no label change on a promotion, which preserves the promoter's degree; it creates the new evaluation `f_d(J)`. Every shape or label change is counted once, since a gain affects labels on just its outgoing edges to old children. In particular the arrow's source rate in the grandchild term is the **parent's** reciprocal, not `1/D_Y`.

The decomposition `(d−p)f_d(Y)+p f_d(Y⁺)=d f_d(Y)+pΔ f_d(Y)` accounts for the artificial killing of `Y` after summation. Thus, for active functions and their dormant **differences**, the sum of positive rows gives precisely the physical drift. The centre uses its ordinary degree **without** a parent slot, contributes no active evaluation and no artificial kill, and supplies the immigration row `I` consisting of all five emission lines in §1, with that degree. In particular it supplies both classes of dormant immigrants, even when just one old child receives the centre's send. There is no incoming forcing from an exterior parent to the centre. A planted isolated physical seed uses the same formula at its top with its inert slot retained; if the top is included as an active evaluation with arbitrary fixed `a`, there is instead the single uncancelled artificial kill there.

For a terminal observable `f` independent of the active root mark and zero while dormant, use `P^c_{t-v} f` as the score in this identity. The finite-time difference property makes the dormant score exactly the required difference. Starting at `T₃`, with no nonfounder cores, the result is

```
E Σ_{v≠o} f(S_v(t)) = ∫₀ᵗ E I(T_u) P^c_{t-u} f du.        (1)
```

A general physical seed adds its initial nonfounder evaluations. This proves equality between current and permanent schemes for such aggregate physical observables with **their own** immigration and seed terms; it is not an equality of their newborn matrices.

## 3. Nonexplosion and the precise scaling

Let `n` count real vertices in a shape, omitting the slot. Give an active particle size weight `n`, and a raw dormant `(l,u;S)` weight `(u−l)n`. For fixed finite labels this is bounded above and below by positive constants times `n`. The active copy payment is at most `d∑ n(Y)+p #Y≤n`, since child subtrees are disjoint and `#Y≤D`. The root-drop payment is at most `sup_D D[q(D)−q(D+1)]∑ n(Y)` since `r≤D`. The grandchild-drop payment is at most `p(1/2)∑ n(Z)≤n/4`; all grandchild subtrees are disjoint. Intrinsic growth of size costs at most `n`, and direct birth at most 2. A dormant's activation term in its drift is bounded by its own weight, and its intrinsic drift is also bounded by its weight. Thus total weight has drift `≤C` times itself uniformly. Stopping at total weight cutoffs gives the finite-time first-moment estimate and precludes explosion: a bounded total weight bounds the number and size of all particles, hence their total event rate. This holds for the raw construction even though its *number* of root-drop emissions before a specified row is bounded only by `O(D #Y)` without using the finite support of the drops. For either finite gauge, a nonzero drop also forces `D` into a fixed finite range.

Alternatively, multiply each emission of dormant `(l,u)` by `u−l` and use the **divided difference** score on that dormant. Its row becomes `L⁰−u` with unit active production at label `l`. For a genuine branching interpretation, kill the dormant at rate `u` and independently emit active `l` at rate one **retaining the dormant** until that death. This is not a probability-one switch competing with a death at rate `l`. It is a diagonal similarity of first-moment systems: active expectations and accumulated active rewards are identical, not sample counts or laws. In this version size weight `n` suffices for both kinds. In both versions the finite event estimates bound a backward derivative for a bounded terminal observable by a fixed polynomial in `n` on bounded time intervals. Physical size is Yule with finite moments of every order; these dominate the localized jump formula and allow passage to (1).

## 4. Strict residual gap on the kernel: transfer through the full time kernel

The [gap proof](gap.md) transfers the permanent expansion on individual planted shapes to the current full semigroup. Subtracting isolated child branches cancels the strict-descendant terms exactly. The resulting weighted remainder gives analytic continuation and, by positivity, genuine convergence of the current residual Laplace transform below `μ`. The [truncation bound](truncation.md) converts the weighted gap to a finite-size approximation with rational finite blocks and a qualitative error rate.

## Exact replay

```
taskset -c 8 timeout 10s python3 -Bu round2/current-mark/physical/check.py
```

145 exact planted/centre generator checks reject omission of the root-drop transfer (134 fixtures), omission of the grandchild transfer (93), and use of the receiver's send rate (66). Five derivative orders check raw and scaled dormant differences; using death rate `l` as the entire dormant diagonal fails three orders. An independent permanent-mark row builder gives the newborn unit-payoff derivative differences (orders 0–6)

```
0, 0, 0, 0, 1/96, −113/3456, 14227/124416.
```

The leading nonzero term implies `K_c(s)1(J_1/2)−K_perm(s)1(J_1/2) ∼ 1/(96 s⁵)` as `s→+∞`, by local Taylor expansion and the finite-time exponential size bounds. Thus identifying the matrices, even just their row sums for all discounts, is false. It does not decide entrywise ordering across discounts.
