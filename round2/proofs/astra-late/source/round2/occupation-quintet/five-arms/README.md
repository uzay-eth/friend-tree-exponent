> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Gold1 P5/B5/S4 physical seams and strict birth-mask ceiling

**Result:** For the explicit mask below, the original negative charged rows contribute between **66.9235674 and 129.9087021 ticks/core** on the lower witness and between **49.8483147 and 104.2611800 ticks/core** on the upper witness (decimal displays; exact values below). The ceilings are below **1.182e−10 / 9.483e−11 potential units/core**. The `1e−6/core` target is `1,099,511.627776` ticks/core, more than 8,400 times either ceiling. This does not bound unrestricted five-core rows or other birth masks. No exponent endpoint is updated.

## Physical row and support

Input is the previously authenticated private Gold1 pair in `round2/.work/replay-inputs/`, appraised by the detached clean source `b04cc50c1e391e0420eb7b006ac781d484c68543`. Replay maps the pair read-only using `round2/attacks/neon2-quintet/transfer/check.py`, without rehashing or scanning the payload tables. Production trial values are `56570237/100000000` (lower, outward sign −1) and `56574828/100000000` (upper, outward sign +1). The scale is `2^40` ticks per potential unit.

| Row | Exact internal support (also require `c≤D`) | Incoming physical legs |
|---|---|---|
| P5 | Three internal `D=2…9`, `c=2…6`; sum of degrees ≤18. The inherited adjacent-pair bound ≤21 is redundant for these headers. | One P4 gain transfer at each deleted end. |
| B5 | Hub `D≤8`, `c=3…6`; long-arm interior `D≤8`, `c=2…6`; sum ≤14. | Two P4 gain transfers at the hub (one per removed short tip), one S3 gain transfer at the arm interior (remove the far tip). |
| S4 | Hub `D≤8`, `c=4…6`. | Four S3 gain transfers at the hub, one per removed spoke. |

Each endpoint uses `G5`: exact through `D=10,c≤5`, arity overflow clipped at 6, a pure degree-11 tail knot at `D=11`, and weight `11/D` on that knot plus `1−11/D` on `INF` for `D>11`. Only internal support failures zero the entire motif. A gain uses the actual `D+1` lookup, and a promotion keeps `D` and uses actual `c+1` for the old vertex. The top row is

```
σ[-λh + Σ_x (D_x−c_x + Σ_{y~x} 1/D_y) Δ_gain,x h
       + Σ_x ((D_x−c_x)/D_x) Δ_promote,x h
       + Σ_(w,x) (1/D_w−r1(x)) Δ_gain,x h_deleted-w].
```

Every transferred term is subtracted once from its source. Equal-labelled tips are still distinct vertices. There is no outgoing higher transfer or newborn term at size five: a newly created five-core occurrence's weight remains in its transmitting four-core row and uses the promoter's post-event arity. Singleton and edge potential, including the separate FY edge swap, never enter the paid reward. The literal full-generator comparison sums raw gain/promotion and paired inter-level arrows; the omitted edge-to-edge FY swaps sum to zero and are not an oracle for individual FY edge rows.

The paid quantity is `max(0,−row)` for the exact physical allocation. The production center/radius projection and the nonnegative host cap are separate: the projected incoming bound dominates each exact incoming arrow and the original clipped cap still covers `row + max(0,−row)`. The receipt's family maximum caps are `26/27` ticks (P5), `31/32` (B5), and `47/45` (S4), lower/upper. The replay does not substitute these for row charges or rederive each production block cap.

## Explicit mask and lifetime comparison

Root at a nonfounder promotion, use actual births in each prescribed chronological order, and retain exactly the children in the current prefix. Additional children of any participant permanently exit that prefix; the designated birth starts exactly one next prefix. The first child must be born at the declared root degree or above. At payment, all five vertices remain inside the listed inclusive degree caps, the top row is supported, and no participant has an extra child. The root's own parent lies outside the motif. Root, order, and membership are determined by physical genealogy, rather than by state labels.

| Reward family | Parent indices for successive births (root index 0) | Degree caps in birth order | Root degree at first child birth |
|---|---|---|---|
| P5 | `0,1,2,3` | `20,5,5,5,20` | ≥10 |
| B5 | `0,0,0,3` | `8,10,10,5,14` | ≥7 |
| S4 | `0,0,0,0` | `8,9,9,9,9` | ≥7 |

Degrees never decrease. If a prefix exits a degree cap, has a disallowed early birth, or misses the first-birth threshold, no later time can earn this mask's reward on that registration. At completion all internal arities are fixed, and failure of internal support cannot reverse under later degree gains. Thus zero continuation for these exits is valid for the **ceiling** as well as the floor. This excludes extra children, alternate chronological orders, founders, and high-degree endpoints beyond the stated caps; it does not treat their physical rewards as zero in any broader account. Family masks are disjoint by motif shape, and all rewards are on size-five rows.

At every prefix the only exterior port is the live root-parent reciprocal `ρ=1/D_parent`. The transition rate at the root and the top-row reward use **the same** `ρ`, which can change at arbitrary external events. The replay verifies each Bellman inequality over all of `[0,1/2]` by checking endpoints and the reward's sign breakpoint when present. It uses the certified Gold1 endpoints `LO,HI` for the discount and directed rounding on a `2^20` subgrid of ticks. Since rates for external events never need to move any tracked coordinate, the inequalities apply directly by summing over physical registrations; no independent or frozen donor distribution is assumed. At each physical promotion a new value starts at `(2,1)`. Summing the Bellman inequalities and using stable-age balance gives `LO·V_floor(2,1) ≤ π[R] ≤ HI·V_ceiling(2,1)`, as in `../higher/README.md`. No reward from earlier prefix shapes is booked.

| Witness | Exact floor (ticks/core) | Exact ceiling (ticks/core) |
|---|---|---|
| Lower | `175229628157123409323189617/2618354562655846400000000` | `105731103061418050055580191/813887764221132800000000` |
| Upper | `130520562249956414848092417/2618354562655846400000000` | `8485689864607149695141173/81388776422113280000000` |

Sibling diagnostics (`../diagnostic/report.log`) measured `2.973/3.135e−6/core` across the entire five-only account on pooled trees. In particular, the large low-degree B5 exterior-child signal and high-degree B5/S4 tip bins are largely excluded here. Those sample results motivate a broader donor-branch registration, but supply no floor or ceiling.

## Executed evidence

Run from `.` with `R=round2/occupation-quintet/five-arms/replay.py`. All numerical jobs used CPU 13, nice 19 and wall deadlines.

- `python3 -Bu "$R"`: passed, CPU `5.83+0.21s`, wall `7.14s`. Per witness: 110 top rows (9/48/53), gain exits 4/24/7, arity exits 1/30/45, repeated-leg rows 7/29/47, 986 nonzero paired arrows and 784 nonzero newborn routes. Whole-potential jump identities were exact. Isolated internal degree/sum and endpoint reciprocal seams were also checked.
- `python3 -Bu "$R" --palm`: printed the lower floor and ceiling and upper floor above; the invocation was interrupted after the upper B5 ceiling, before the final S4 and total. Registration matched brute chronological subsets across physical birth snapshots (`P5=6,B5=1,S4=1`).
- `python3 -Bu "$R" --palm --side lower --direction ceiling` and the same command with `--side upper`: both exited zero and reproduced the exact ceilings above. CPU `19.30+0.10s / 20.25+0.07s`, wall `33.43s / 20.74s`. Terminal mask states per family: 13,376/8,424/8,192. The Bellman corner counts at the terminal layer were 27,046/19,836/20,257 on the lower witness and 27,050/19,700/20,293 on the upper witness.
