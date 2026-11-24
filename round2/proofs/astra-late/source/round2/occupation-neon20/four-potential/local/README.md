> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Local count Poisson correction at the four-potential trial

At `s=.56572365` the degree-25, order-five local correction sourced from the 16 adverse controls does not give a tenfold improvement that survives new trees. The unrestricted cell fit essentially interpolates the controls using founder cells. The nonfounder fit with a certified absolute stationary potential budget of `1/50` changes the training half-range only from `.003438276734` to `.003430545134` and increases the held-out half-range from `.002048210817` to `.004849289294`. A separate exact dual certificate on untouched physical trees leaves a barrier of **`8.194807e-6` per unit of the original stationary denominator** for *every* signed four-way mixture and *every* added degree-25 local count correction, even after an arbitrary constant offset. The barrier applies to the specified indicator family, not to stationary cancellation on the true law.

## Actual cells and the physical operator

A cell is a connected oriented subtree of one to five core vertices, marked by the exact `(L,k,h)` at each retained role. Child permutations identify cells but never identify two physical branches. Keep the true degrees at all roles, with `D≤25` for the **weighted cell**, and let `N_j` count occurrences. The proposal corrector is `H=Σ z_j N_j`. Its columns are `(ℒ−s)N_j/C` from true leaf/core arrows, including entries whose **pre-event** cell carries weight zero. The fixed source is the rational rounded four-way mixture `(.00190648,.28049101,.03306576,.68453675)` from the earlier sum-one contrast fit; its coefficients do not change in the local fit.

The 16 trees are the union of the earlier fixed supports. Their old and one-jump cells give 1,023 active generator columns. The 135 held-out trees are asymmetric brooms, new chains and stars (including shifted stars), and overflow chains/stars. None appears even in the earlier 541-row discovery catalog. The optimizer proposes minimax and then min-`l1` coefficients. Every displayed fit uses rational-rounded coefficients and direct exact residual recomputation. The bounded nonfounder fit is checked by three routes on the 16 originals and seven held-out controls: compressed post-state counts, exact ordinary-arrow whole-tree subset enumeration, and the pre-event marked flow ledger in `marked_flow`.

| Allowed correction cells | Training half-range | Held-out half-range | Stationary `abs(πH)` ceiling |
|---|---:|---:|---:|
| None | `.003438276734` | `.002048210817` | `0` |
| Singleton, including founder | `1.015e-12` | `.005239381097` | `.004340373` |
| Up to three, including founder | `1.368e-12` | `.002048210817` | `.000386432` |
| Up to five, including founder | `9.921e-12` | `.002048210817` | `.000165277` |
| Up to five, nonfounder only, unrestricted | `.001849935915` | `6.107658466` | `80.709728841` |
| Up to five, nonfounder only, budget `1/50` | `.003430545134` | `.004849289294` | `.019999996` |

For the budget column, a cell with `k` physical children and `m_c` retained children of identical cell type `c` has at most `(k)_(Σm_c)/∏ m_c!` assignments, multiplied by its children's bounds. Sum `|z_j|` times the resulting product for nonfounder cells; this is a rigorous ceiling on `|πH|`. Founder-owned cells have bounded count and their stationary mass is zero. The fixed original mixture has `πΦ≥.20696242982264076` (rounded value of the exact lower endpoints), so the budgeted corrected denominator stays positive. The unrestricted live fit's enormous bound gives no such certification.

There is also a tight, exact obstruction to the **budgeted fit**, independent of LP tolerances. In training positions 12 and 13 (zero-based), the maximum absolute generator contrast of a nonfounder cell divided by its multiplicity bound is `.002028912653` rounded up. Thus with multiplicity-weighted absolute budget `1/50` the corrected half-range cannot fall below **`.003417987484`**, even with all 1,023 candidate columns. All unseen columns are zero on both rows. This is within 0.37% of the tested budgeted fit and rules out tenfold improvement already on the training set in that budget class.

## Exact adversarial certificate

On 111 held-out trees no active local correction column is visible except a possible newborn singleton. Entries from existing gains/promotions are evaluated, not assumed absent. A particularly transparent subfamily has all pre-existing degrees above 25; monotonicity keeps them out, and no single jump creates anything more than one newborn. Its only new count is the virgin singleton, with marked per-core birth rate `b(T)=Σ L_v/D_v / C(T)`. The dual includes this column and the constant column with target zero. Take the exact box-corner representations of `(d,0,0)` in the augmented rows `(r_T,b(T),1)` as in the earlier four-potential proof. Fourteen ordinary whole-tree controls independently validate the union of used supports. Thus for every `w` and every such `H`, the supremum of `|(ℒ−s)(Σw_iΦ_i+H)/C−a|` on these trees is at least `8.194807e-6 |d·w|`, with `d=π(Φ_i)_i`. In the sum-one gauge the corresponding lower bound is `3.1504411e-5`. Neither assertion replaces `d·w` by the corrected denominator when the latter is unrestricted. For the budgeted positive fixed mixture, denominator inflation is already bounded as above.

## Exact boundary ledger required for a global stationary sign

Use the same owner convention for every term. Each event that increases degree at a tagged core has gain increment `Δ_D z` of the same extant connected pattern. The neutral loss `−s z` belongs to the original pattern. A promotion contributes `(L/D)(z_post−z_pre)` to that pattern, and pays a new connected extension via `(L/D) z(post,virgin)` to its unique pre-existing pattern obtained by deleting the newborn tip. The newborn singleton is paid by the promoting singleton. For an exterior gain, append its **actual donor** on the correct labelled incident edge and pay `Δ_D z / D_donor` there; a five-core test thus requires a six-core incidence column, including an exterior mother where relevant. This rule is exactly `marked_flow`, independently checked against the complete physical jumps. In particular, all tests with an out-of-mask pre-state and an in-mask successor, as well as the `25→26` outgoing loss, have explicit coefficients. An exterior donor need not satisfy the degree-25 cap.

Let `E` be any certified finite retained list of those marked terms and the corresponding original raw-potential terms. Write `J_E` for their signed stationary sum, and `U_−,U_+` for exact outward bounds on the negative and positive omitted sums. Founder values have zero stationary count, but their outgoing promotions into nonfounder counts must be assigned before taking the limit; their direct intensity is negligible. For positive corrected denominator `d_*`, the exact test is

```
J_E - U_- <= (μ-s) d_* <= J_E + U_+.
```

For a localized count correction, this follows by the integrable marked balance `πℒH=μπH`, not by declaring `π(ℒ−s)H=0` at an unknown trial. The correction itself leaves the true stationary sign unchanged wherever the corrected denominator is positive. A useful finite separator must put the signed *retained plus boundary* expression on one side of zero. To obtain tenfold reduction in a *certified interval width*, the sum of retained occupation uncertainty and both omitted error budgets must shrink accordingly; no bound on pointwise training error substitutes for these budgets.

The exact missing observables are joint occupations of (i) marked pre-event neighbour cells, (ii) their exterior incoming-port reciprocal moments by receiver, including six-role incidences for five-core columns, and (iii) promotion/exit seams and a controlled overflow complement. For example the contribution of an omitted exterior port with coefficient `δ` is `δ π[Σ_occ Σ_donor 1/D_donor]` on its actual marked tuple. If an omitted port has `|δ|≤B` and its positive reciprocal incidence mass is at most `m`, the loss is at most `Bm`; the degree cap on the **receiver** gives no substitute for that mass. A bound on count-only occupancy or a product of independent donor and receiver marginals is insufficient for the signed expression. The finite fits give no numerical stationary enclosure on these omitted moments, hence no central sign here. The adversarial dual pinpoints the need for overflow observables or an occupation-weighted omission; merely adding more bounded-degree cells cannot remove that exact blind control.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time -p timeout 130s \
  ./round2/.venv/bin/python -Bu \
  round2/occupation-neon20/four-potential/local/replay_local.py
```

The replay exited zero in 84.69 s, including 16 rational denominator-box corners and the independent checks above. The payloads are the same four inherited pins as in the parent experiment. No full scan or new payload authentication is used. Outputs are finite-control diagnostics and rational obstructions; there is no fitted endpoint or Fable input.
