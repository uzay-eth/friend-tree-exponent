> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Persistent side-donor flow for older motifs

**One persistent donor gives a strict improvement over independent receive boxes on matched stopped tagged flows.** In the two-extension `P4` test the lower signed price crosses from negative to positive at all three tested apex degrees. In the corresponding `Yt` test the upper signed price at apex degree 18 nearly halves. These are conditional entry prices, not new exponent or whole-complement bounds.

## Physical information added

Fix the tagged topology. At one receiver `x` distinguish one actual side child from the other `k_out−1`. The state stores no donor, exact donor degree `m=2,…,H`, or an irreversible tail bin `m≥H+1`. For the distinguished donor, `B_known=1/2−1/m`; this nondecreasing share of the total deficit never resets. The other donors remain a separate box `[0,(k_out−1)/2]`. A receiver with no side child starts without a donor; its first surplus-child promotion installs the newborn with `m=2`. A receiver with an existing side child chooses one by a fixed birth-order rule at entry and retains it. The selected slot is never an extension parent on these masks, so the stored donor is never reclassified as a spine tag. At apex entry all allowed initial donor bins are included in the extrema; no degree-two condition is assumed there.

The exact bin `m` contributes reciprocal `1/m` to the receiver's gain, in both its generator and signed reward. After the donor's own degree gain its deficit jumps by `1/[m(m+1)]`. The rate of that jump is in `[1/D_x,1/D_x+m−1]`: the tracked parent contributes exactly `1/D_x` and each of the `m−1` other neighbours contributes at most one. The donor's own promotions leave its degree and its parent's tracked census unchanged; their future effects are covered by the same gain-rate box. Once the donor crosses `H`, the tail contribution stays in `[0,1/(H+1)]` and never returns to a young bin. For the one known donor this uses the same first-moment drift bounds as the proposed `T`, but retains the finite upward transition and the birth source as well. No independence between donor and receiver histories is assumed. Simultaneous independent boxes are a relaxation for all unobserved rates, not a claim about realizability of every corner.

## Matched signed integration

For `c=c_lo−1.2415 c_hi`, use the original pinned charged row at every full tag inside the mask and zero running reward for proper prefixes. Every tag promotion updates the old tag, and a promotion at the required parent additionally emits one new tagged prefix. The `P4` and `Yt` two-tag starts traverse two possible continuation births (`2→3→4`), so their prices include how often the final tag can be formed before exit. The `Yt` extension selects only the younger of the two siblings at birth and counts each unordered pair once. The older side child used for the apex variant stays attached across both extension births. In the `through-terminal` variant the terminal begins at its actual birth state and its own later promotion first introduces the distinguished donor.

At exit from a tag's degree/arity cap that tag's continuation is zero; the incipient birth and the gain or promotion action at the last interior state are still evaluated. This is a stopped occupation observable: there is no claim about the sign of the discarded complement, no reuse of a gauge exit value and no re-entry into that cohort. The baseline and memory runs use the identical caps, rows, signed discount, prefix forcing, rule for the younger `Yt` sibling, and zero exit continuation. In the baseline only, the chosen receiver's total side pressure may independently take any value in `[0,k_out/2]` on every action. The exit and the source both appear in the exact quotient recursion; no floor at zero is applied to signed values. The quotient's extrema follow by sorting the extra rates by their action values and adding each endpoint while it improves the quotient. The script cross-checks that optimizer against all corners of a signed fixture.

Default discount range is the entire pinned `μ` bracket `[.565673570322,.565775]`. The optional extra kill rate in the quotient lets the discount vary across states in that range; it is a relaxation enclosing any fixed physical discount there and avoids assuming monotonicity of signed rewards.

## Exact replay

Shown are conditional intervals rounded for display, with `H=5`. The script uses rational arithmetic throughout. The apex starts at `D∈{18,24,64}`, `k=5`; `cases` in [replay.py](replay.py) specifies the exact starts, inclusive caps and physical entry variants. Each interval is in units of integrated signed charge **per entry**, rather than per core.

| Case | D | Independent boxes | Persistent donor |
|---|---:|---:|---:|
| `P4` through terminal | 18 | `[-9.77732e−9, 7.29118e−9]` | `[1.32343e−9, 7.16522e−9]` |
| `P4` through terminal | 24 | `[-2.90474e−9, 2.21452e−9]` | `[8.82816e−10, 2.18190e−9]` |
| `P4` through terminal | 64 | `[-5.30008e−11, 4.34447e−11]` | `[3.46385e−11, 4.31645e−11]` |
| `Yt` through terminal | 18 | `[-2.87473e−8, 2.69281e−8]` | `[-2.85266e−8, 1.34757e−8]` |
| `Yt` full-tag entry, terminal with no child yet | 18 | `[-3.48824e−5, 2.98170e−5]` | `[-3.48077e−5, 1.69736e−5]` |
| `hV` apex, unknown initial donor degree | 18 | `[-4.62972e−7, 3.82351e−7]` | `[-4.38988e−7, 3.33137e−7]` |

All 54 exact comparisons are weakly inward and 51 strictly inward; the three ties are the `Yt` upper bounds with apex donor only. The cheap two-bin ablation `H=2` still flips the three `P4` through-terminal lower prices and gives a `Yt` degree-18 through-terminal upper of `1.43520e−8`, so most of the first improvement comes from preserving the virgin bin and the irreversible `p≤1/3` seam. On `hV`, the calculation works directly with the charged row, including the routing responsible for nonzero curl; no scalar gauge is needed.

For scale, a unique first-entry cohort on a fixed two-tag start can have normalized entry intensity at most `μ`: each nonfounder birth creates one physical parent–child prefix, and age discount from that birth to its first matching entry is at most one. Even using that coarse maximal intensity, the degree-18 through-terminal changes in the reported `P4` lower and `Yt` upper certified accounts are below `6.3e−9` and `7.7e−9` per core, respectively. Those are improvements to the stated relaxation for that cohort, not upper bounds on what more physical information could recover. The full-tag conditional gain is larger, but no corresponding entry intensity was bounded here. No sum over overlapping apex/terminal variants was taken and no charge is booked toward either complementary target. Next use a disjoint physical entry mask with an intensity bound, or keep two known donors while retaining the residual bin; for a signed complement, carry the exit flow into the adjacent cohort.

```sh
taskset -c 8 nice -n 17 timeout 30s env OPENBLAS_NUM_THREADS=1 ./round2/.venv/bin/python -Bu round2/occupation-neon/pressure-memory/old-motifs/replay.py b3f18_lo.bin b3f18_hi.bin
```
