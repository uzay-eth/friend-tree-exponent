> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Old-flow: positive age-tilted adjacency moments on live old branches

The certified finite-band addition to the [high-parent + size-20 ledger](../../high-parent/README.md) is at least **`3.0951968e−9` lower / `4.0254818e−9` upper** per core, giving the cumulative outward bracket **`[0.565676636772, 0.565771871520]`** at the fixed denominator ceilings. This is a small recovery; dropping the old branches is unnecessary, but the remaining adversarial receive-rate corners are expensive. The separate unbounded trial beginning at degree 65 pays about `3.11e−13` per core per sign, and is **not** added to the cumulative bracket.

| Charge | Lower price (analytic replay) | Upper price (analytic replay) |
|---|---:|---:|
| `hV` | `3.06204334292295e−9` | `4.00077770896578e−9` |
| `Yt` | `1.49696294463021e−12` | `2.90279385539250e−12` |
| `P4` | `3.16565034087608e−11` | `2.18012978802676e−11` |

The table decimals are rounded down. See `controls.py` for the exact arithmetic on these declared floors. The endpoint increments before rounding are about `2.28119e−9` and `2.94439e−9`.

## Exact birth forcing and age-tilted stationary moments

Anchor each embedding at its top physical core and retain its directed ancestry. `hV` and `P4` are ordered spines of length two and three. `Yt` extends top→middle→older child by the birth of the younger child of that middle. A prefix ends in a list of actual vertices, with topology fixed and the per-slot census `(L,k,h=1)`. Edges not included in the tags stay in the tree; they enter the measured arity `k` and the receive interval at their slot. The process does not erase old branches, condition on no side births, or demand independence after witnessing a tagged birth. Leaf stock and arity of every tagged vertex update at **every** promotion of that vertex, including an off-spine promotion at a crowded apex.

For a prefix of length `j`, write `M_j` for the stationary intensity of its descriptors, counted per core with adjacency multiplicity. In finite time each prefix comes into being at the promotion creating its last slot. For `Yt`, each unordered pair has a unique younger sibling; the pilot conservatively retains only half of the resulting nonnegative price. There is no factor `k` in a birth transition: a fixed prefix produces one extended prefix per appropriate promotion. The exact physical forcing, for the admitted birth-ordered prefixes, is

```
F_1 = mu δ_(newborn)  (or mu β(d_birth) δ_(newborn) within a birth class),
F_(j+1)(φ) = M_j[ b_i φ(x_i^bud, newborn) ],
M_j[(mu − A_j) φ] = F_j(φ).
```

Here `i` is the prescribed parent slot for the new tag. `A_j` tracks already existing tag identities and does not remove them at birth. The aggregate transport proof is to compensate the physical last-tag births and count at their top, or transport each charge to the last tag and use that birth-Palm law. A stationary age tilt of a prefix uses the **joint post-promotion environment**, not a product of per-slot ages. These are identities before imposing calculational masks. For bounded finite descriptive tests they follow directly in finite-tree time means; the [oriented rank-one local limit](../../../../../../astra/source/round2/fringe/spatial-limit.md) passes to the limit after controlling multiplicities by the exponential core-child moments. For fixed lengths, products of consecutive arities can be controlled by Hölder and edge reversal, or one can first stop the finite test at the bounded tagged arities used in the pilot; there each participating apex has a bounded count of paid embeddings. No large-degree leaf moment is assumed. The founder and a bounded number of generations from it vanish in the mean normalization.

## Certified backward actions without killing old branches

Compute a positive value for each prefix of the stencil. A promotion at `i` contributes `b_i V_j(x^bud)` for the old prefix plus `b_i V_(j+1)(x^bud,newborn)` when `i` is the correct parent slot. A receive changes the relevant `(D,k)` to `(D+1,k)` and a promotion to `(D,k+1)`. Thus the finite recurrence is triangular in the number of tags, their degree and their arity. Known sends between tags have exact rates `1/D_sender`. For each untracked child group there are `k_i−(# tagged children)` slots, with aggregate receive in `[0,(k_i−# tagged children)/2]`; at the apex its own untracked parent has inverse in `[0,pmax]`. The exact charged reward for a fixed stencil is affine in these same rates. Both rates are shared between reward and generator; minimization is over endpoints of an affine inequality, checked afterwards in exact arithmetic. Per-state minimizing endpoints may differ and may adapt to the full history. Values include all admitted future branches; cuts kill only the particular tag state that exits its census bounds. In particular, there is no cap on physical fringe size. The finite calculation uses `HI≥mu`, with normalization `LO≤mu` at the top newborn, and cannot count an outgoing branch twice.

The physical Palm ladder uses the exact capped parent degree at birth, with classes `{2}`, `{3,4}`, `{5,…,9}`, `{10,…,19}`, `{20}`. Parent degree never decreases, so `pmax=1/min(class)` holds at every age. Sum the universal value plus the **excess only** with the directed floors from [the separate physical class mass solve](../next/README.md). The residual 19-mark Perron coefficients are a different law and are not used in this step. For lower `hV` the support `D=18..32`, arity ceiling 7, `--palm` raised the unstratified price `1.1016686896e−9` to `3.0620433429e−9` in matched support. Four other birth classes each improved the value. Upper `hV` uses a wider degree mask and the universal parent bound, because the full wide Palm pilot hit its time limit. These are admissible subsolutions, not an optimum over all adjacency-aware stationary measures.

## Exact limiting degree action

At `D=65`, retain the old physical stencil and its unknown receive controls. After this entry continue the apex gain forever; stop the face at any tagged non-apex gain or tracked promotion. **Any** off-tag changes are still allowed. At each corner of the uncaptured parent/child receive box, the exact charged reward multiplied by `(18+18/D)^3` is a polynomial of degree at most seven in `t=65/D`. Only one slot is a Bernstein tail; a gain at that slot replaces `18/D` by `18/(D+1)`. In the term with physical gain intensity `D−k−1+ρ`, the possible simple pole cancels against the gain difference, which vanishes at `t=0`. Center and routing flags stay fixed on the face, and promotion/reciprocal coefficients have degree at most one. Slot weights and tail transfers elsewhere are constant or cubic in this one varying coordinate. Exact interpolation plus two independent extra degrees verifies the evaluator; lower bounds are certified by nonnegative Bernstein coefficients after factoring the leading power. Taking the worst coefficient of every corner gives a uniform `c t^a` lower reward. Since apex gain is bounded by `D`, the relative loss of `D^(−a)` is at most `a`; all other tracked exits have rate at most `1+Σ_(other slots)(D_i+1)`. The scalar face continuation `c t^a/(HI+exit_bound+a)` is therefore a bounded exact subsolution for every adaptive outside process. It remains a valid gain-band boundary without a finite upper degree cap.

With `--tail-only` and the hV mask `5≤k≤8`, both witnesses give positive constant (`a=0`) continuation on all 48 tested faces. The full age-weighted prices are `3.11093710371861e−13` and `3.10492287419809e−13`, respectively. This is a lower bound from entering the tail with a small existing stencil; it does not price paths newly seeded after age enters the tail. Compare by **max**, not sum, with 2417's overlapping `D≥65` higher terms.

## Exact ownership masks and remaining term

All paid stencils have every non-apex tag nonvirgin and at degree `≤3`. Apex bounds are inclusive. Lower `hV`: `18≤D≤32, 5≤k≤7`. Upper `hV`: `18≤D≤64, 5≤k≤8`. Both `Yt/P4`: `18≤D≤36, 5≤k≤6`. `k≥5` implies at least six cores in the full apex fringe; `D≥18` implies physical fringe size greater than 12. Thus these cannot overlap the size-20 three/four-core higher rewards or the size-12 prefix. The high-parent ledger pays only mixed `g` and virgin `hC` here. 4e4f's additive older `hV` mask instead has apex degree 13–17 or arity 3,4 at degree 18–64. Other methods on any of the same embeddings should be combined by max unless charge/descriptor/occurrence ownership is explicitly partitioned. The rotated bare-P4 calculation with an apex of degree two is outside this mask and was not added here.

The unresolved term is the true **off-stencil reciprocal receive contribution** in `hV`, and then its shared `Yt/P4` continuation: the minimum still permits all surplus stubs to behave as degree-two donors concurrently, and independently at different receivers. Improving a tail interpolation alone leaves that gap. A scalable next constraint is to transport the adjacency-weighted receive excess back to the actual surplus edge at each receiver, jointly across the two ancestor generations: use the prefix transport forcing above with one additional edge census, and subtract the already booked corner floor rather than paying the receiver again. A certified shell that permits older inner slots of degree 4–8 and higher arity is another way to target this gap; tiny four-slot payments at degree ≤3 do not show that the entire four-body class is small.

## Replay

Set `PY=./round2/.venv/bin/python` and `DIR=round2/occupation-neon/renewal/old-flow` from repository root. Each command below was run pinned/niced using `taskset -c 8 nice -n 17 timeout 30s env OPENBLAS_NUM_THREADS=1` (shorter ceilings were used on the one-shot diagnostics):

```sh
"$PY" -Bu "$DIR/pilot.py" b3f18_lo.bin --degree 32 --arity 7 --tip 3 --palm
"$PY" -Bu "$DIR/pilot.py" b3f18_hi.bin --degree 64 --arity 8 --tip 3
"$PY" -Bu "$DIR/pilot.py" b3f18_lo.bin --degree 36 --arity 6 --tip 3 --family Yt --family P4
"$PY" -Bu "$DIR/pilot.py" b3f18_hi.bin --degree 36 --arity 6 --tip 3 --family Yt --family P4
"$PY" -Bu "$DIR/pilot.py" b3f18_lo.bin --degree 64 --arity 8 --tip 3 --tail-only
"$PY" -Bu "$DIR/pilot.py" b3f18_hi.bin --degree 64 --arity 8 --tip 3 --tail-only
"$PY" -Bu "$DIR/controls.py" b3f18_lo.bin b3f18_hi.bin
```

The last command passed 18 exact tail signed-corner fixtures plus three inner extractor identities per sign and printed `combined outward 0.565676636772 0.565771871520`. The raw pinned payloads and previously verified physical decomposition are inherited inputs. No payload hash rerun or comparison scan is claimed.
