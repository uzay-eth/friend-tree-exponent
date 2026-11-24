> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# Gold2 original upper: physical S4

**Correction:** The outputs in `finite.log` and `mixed.log` are withdrawn.
They divided the exterior reward slope by the common rate denominator a
second time. Do not use those floors in any Gold2 composition. The corrected
replay rounds `slope × GRID`, as the common denominator multiplies both the
reward and the gain rate. The exact-corner fixture detects the old formula.

Stationary comparison for the five-core **star with four arms (S4)**,
using the original upper payload only. For each nonfounder hub we include its
real mother and three unordered real children; the hub's mother must itself
be a nonfounder. The selected children are unordered as a set, but their
actual birth order determines registration. These are S4 rows, disjoint from
both H3 and B5 credits. The number `~1.45 ppm` is a finite-tree diagnostic;
even its mother-high subpool is not an occupation bound.

## What is paid and what stays in the host

At `sigma=+1`, original `s=28286843/50000000`, read the charged S4 row
from the production arrays using the exact physical tip projections. Four
incoming S3 transfers are evaluated as four actual arm deletions, including
when the arm states coincide, using `(1/D_deleted-center_hub)` times the
appropriate original S3 hub gain difference. A gain or promotion updates the
state of the real vertex. S3 newborn S4 weight stays in the S3 transmitting
row at actual post-promotion arity. The only increment is `(-H_S4)_+`.
At production hub knot `h` the verifier retains `max(0,K_S4(h))` per S4
vertex subset (`verify_quintet.cpp`, sections 10–11 and the host check).
`H_S4+(-H_S4)_+=max(0,H_S4)<=max(0,K_S4(h))`; every original local cap
and its multiplicity are retained. The reported maximum `50` just checks
the fixture rows; it is not substituted for the local caps in composition.

The original private staged pair has both payload checks `OK`, source
`b435b4c824e2ba5835f58b85e67a27de7d6a0d49`, original upper payload SHA
`f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`, and
receipt SHA `5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284`.
Replay uses the frozen Snapshot and independent physical wire reader with
shared inode/size/time checks; it does not scan or open the lower payload.

## Comparison

Start at the birth of the nonfounder mother in `(2,1)`. Each of her
promotions registers a hub; each hub promotion both continues its current
prefix and registers one successor prefix with the newborn `(2,1)`, up to
three marked children. Earlier children not selected remain in the exterior
stock. Summing over histories counts each crown once. Sorting marked child
state tuples quotients storage only: gain and promotion are summed over
each labelled arm. The hub can have one or two extra real children in a
paid row (`c_hub` up to 6, including its mother); extra births advance live
stock instead of terminating a prefix. All four marked arm reciprocals are
live. Remaining neighbors contribute a shared adverse port per marked vertex
in `[0,(c−marked_degree)/2]`; the same unknown sum is used in both the charged
row and gain hazard. This bounds all actual missing-degree histories without
replacing their reciprocal rates by a point estimate. It explicitly tracks
exterior **count**, not individual exterior degrees.

Each state solves an integer-rounded one-step discounted Bellman inequality.
For the clipped reward, the maximum of the zero and affine reward
subsolutions is valid with the same already-computed continuation. Reward
gates leave continuation running. Degree and arity caps, and the monotone
mixed-child degree cut, can be killed on exit. The outer lifetime
normalization is `LO` and discount `HI` from the original Gold2 enclosure;
divide ticks/core by `2^40` for normalized credit. See the frozen post-FY
proof for this rooted-birth occupation normalization. S4 shares the generic
affine solver but does not call the frozen H3 layer or use its numerical floor.

## Replays

From repo root, `r=round2/occupation-quintet-gold2/s4/replay.py`.

- `python3 -Bu "$r" --fixtures --mother 55 --mother-stock 9 --hub 8 --child 3`:
  [finite.log](finite.log) is a historical, invalid numeric replay. Its mask
  used mother `D=19…55,c≤9`, hub `D≤8,c≤6`, three children `D≤3,c≤2`.
- `python3 -Bu "$r" --fixtures --mother 125 --mother-stock 10 --hub 7 --child 9 --child-stock 1 --young 3`:
  The [corrected replay](mixed-corrected.log) exited 0, pricing 211,860 paid
  mask states. Its floor is `8942817855164742770193669/6201700363750000000000`
  ticks/core = **1,441.9945064481096**, or **1.3114863636002243e-9** normalized.
  [mixed.log](mixed.log) is the historical, invalid numeric replay. The mask pays
  mother `D=19…125,c≤10`, hub `D≤7,c≤6`, three children each `c=1`, two
  with `D≤3` and the third with `D≤9` (any birth rank).
- `python3 -Bu "$r" --fixtures-only`: 25 physical rows, including 10 with
  four identical nonzero deletion arrows and 5 with nonzero exterior-tip
  terms; 35 distinct crowns. The Bellman corner comparison checks physical
  rates at all 12 corners across three scales and rejects the old scaling.

The practical truncation to watch is the older included child. The preserved
tree split from 835c gives the larger mother-high old+young cells at mother
`D=180…1998` and second-highest arm `D=19…179`; small-child receipts cannot
be extrapolated to those cells. A next comparison should keep at least
one young child, allow one older included child across projection regions,
and bound the remaining state costs by exploiting sparsity and safe slices.
