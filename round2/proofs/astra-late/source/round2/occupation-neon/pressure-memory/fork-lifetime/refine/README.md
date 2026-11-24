> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Same-mask Yt refinement

**Full-mask width at θ=1: `3.551665636e−6` per uniformly sampled core**, down
35.8% from `5.532089690e−6`. The `1e−6` target remains open. The mask is
unchanged: all 4,095 oriented Yt rows with every tagged degree ≤5 and every
tagged core-child count ≤3. The signed complement is unpriced.

| Method | θ=.9 width ≤ | θ=1 width ≤ |
|---|---:|---:|
| Original fork-lifetime certificate | `5.232571583e−6` | `5.532089690e−6` |
| Same rates, physical five-bin birth floors | `4.682815950e−6` | `4.950318726e−6` |
| Asymmetric fork/tip memory, same birth floors | `3.797441923e−6` | `4.014967954e−6` |
| Same conditional endpoints, sided physical birth law | `3.361340730e−6` | `3.551665636e−6` |

The final intervals are `[-1.129666819e−6,2.231673911e−6]` and
`[-1.361020379e−6,2.190645258e−6]`. They intersect two certifications of
the **same** conditional payoff: no birth-law or incidence charge is added
twice. `8d8f` was unavailable as a local Git object; the sided reprice uses
the delivered `birth-law` replay (`round2/occupation-neon/pressure-memory/birth-law/replay.py`, historical reference), with its
cutoff-24 sided physical constraints and exact reduced-cost repair.

## Where the `(3,2)` width lives

The 819-row submask additionally fixes the trunk at `(L,k,h)=(0,2,1)`.
Both its children are the tracked tips, so **neither a trunk-side receive
nor a trunk exterior parent enters its paid row**. The actual apex parent
is the largest isolated uncertainty. With five-bin physical birth floors
and the original fork-only memory, the θ=1 widths are:

| Submask diagnostic | Width ≤ |
|---|---:|
| Coupled total | `3.424114560e−6` |
| Intercept | `0.552839862e−6` |
| Apex-parent receive | `1.475912668e−6` |
| Apex-side-child receives | `0.775319888e−6` |
| Both tips' side-child receives | `1.299458608e−6` |

The apex-parent diagnostic's conditional width at the high observation
discount is `1.836008430e−6` for birth-parent degrees 3–4, versus
`1.286703608e−6` for degree 2 and `0.290452573e−6` for degree ≥20.
This identifies the broad *initial* receive range of the low-degree tail
as a source of uncertainty. On this submask alone, tracking both tips'
side donors to degree two narrows the coupled width to `2.531487690e−6`
under the same five-bin weighting. The isolated tip-source width drops
to `0.593173444e−6`, while the apex-parent contribution still has width
`1.427811765e−6`. These are separate projections, not additive gaps.
The parallel parent-time audit (`round2/occupation-neon/pressure-memory/birth-law/parent-time/README.md`, historical reference)
has a different conditional pressure refinement on the same submask;
its reported `2.600865763e−6` width and age/loss diagnostics are not stacked.

## Shared rates and the signed drift test

The replay keeps the trunk, both tips, and the apex live in one projected
state. Tip subtrees share the real trunk's degree and send rate, not two
independently adjustable copies. The tagged internal sends, the exterior
parent's send, and its own gain are different physical events; the latter
includes the apex's exact reciprocal. Every exterior receive coefficient
uses the **same jump** as its affected tagged degree and the same rate in
its successor and holding denominator. Side children attached to distinct
tagged roles are distinct donors on a tree. Symmetric completed tip
censuses can be merged after registration; prefix registrations are not
merged into independent physical promotions. Both tips' side donors keep
their actual birth and aging until the trunk exceeds degree three, then
their detailed pressures are relaxed (the full masked reward still
continues). Fork side donors keep the original cap three. The replay also
implements an independently adjustable apex-side cap and a monitored
parent `(D,k)` with the actual reciprocal and leaf stock. The completed
principal benchmark leaves apex-side donors opaque and monitors parent
degree 2; an attempted larger joint mother run did not complete and is
not used in the table.

For the signed apex-parent term on `(3,2)`, construct a killed potential
`Φ` with `Φ(D+1)−Φ(D)=−slope(D)` by summing the signed slopes backward
from apex degree 5 and zero beyond the mask. Solve the rewritten reward
`q + (A−μ)Φ`, including the entry impulse when the last edge registers and
all other drift/exit terms, using the common observation discount on
`−μΦ` and the holding denominator. The original exterior-parent slope
cancels exactly before rounding. The resulting interval equals the direct
one at displayed precision (`1.427811765e−6` width): a change of potential
alone does not constrain the allowed adaptive rates. It transfers the
uncertainty to the residual drift rather than removing it.

The next meaningful improvement needs a joint live apex-parent pressure
constraint, especially within birth degrees 3–4, together with retained
tip pressure. Birth mixing alone still leaves a broad within-bin payoff
range; another one-source antiderivative cannot shrink it. There is no
proved lower bound on what a stronger physical relaxation can achieve.

## Replay

From `.` (use `round2/.work/venv/bin/python -Bu` for `$PY`):

```sh
taskset -c 7 nice -n 10 "$PY" -Bu round2/occupation-neon/pressure-memory/fork-lifetime/refine/replay.py --caps 0 3 2 --mother 2 --tip32 --forecast
taskset -c 7 nice -n 10 "$PY" -Bu round2/occupation-neon/pressure-memory/fork-lifetime/refine/sided.py
taskset -c 7 nice -n 10 "$PY" -Bu round2/occupation-neon/pressure-memory/fork-lifetime/refine/replay.py --only32 --caps 0 0 2 --mother 2 --barrier
```

The full forecast authenticated both D18 payloads and the scanner, visited
all 4,095 rows, and checked 623,740 integer drift residuals on 155,935
retained states (65.32 s wall, 20.56 CPU s, 175,256 KiB peak RSS). The sided
replay prices its exact frozen tick vectors with 613 constraints and 5,554
columns (21.86 s wall, 7.05 CPU s, 111,860 KiB peak RSS). The barrier run
covered all 819 rows and checked 892,740 residuals before a later optional
submask pruning change. Its numerical result is reported only for that
completed run.
