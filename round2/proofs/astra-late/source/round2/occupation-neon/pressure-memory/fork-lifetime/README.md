> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# All-age fork registration for coupled D18 contrasts

The principal mask includes **every** `hC` and `Yt` row with an oriented apex
(`h=1`) and `D≤5, k≤3` at each tagged role. There are 455 and 4,095 canonical
physical labels, respectively. Bounds below are per uniformly sampled core
(the sum of incident crowns anchored there), not per selected crown.

| Family | θ | Stationary `lower − θ upper` | Width ≤ | Matched no-memory width ≤ |
|---|---:|---:|---:|---:|
| `hC` | .9 | `[4.828e−12, 7.707e−12]` | `2.878e−12` | `4.188e−12` |
| `hC` | 1 | `[−1.4e−14, 8.6e−14]` | `1.00e−13` | `1.43e−13` |
| `Yt` | .9 | `[−2.420385877e−6, 2.812185706e−6]` | `5.232571583e−6` | `6.669167750e−6` |
| `Yt` | 1 | `[−2.743394253e−6, 2.788695438e−6]` | `5.532089690e−6` | `7.054151729e−6` |

Widths are rounded outward independently of displayed endpoints. The `Yt`
reductions are 21.5% and 21.6% at `.9` and `1`. Crown incidence masses are
`[.016014405498677,.025523849989836]` for `hC` and
`[.003569191871076,.006394842454402]` for `Yt`. These are expected counts.
The signed complement (any excluded degree or arity, all other families,
and the founder if considering a finite population) is not priced here.
The limiting uniform-core law has zero founder mass. No signed global
contrast, scalar exponent sign, or intrinsic `V_s b` credit follows.

## Registration and payment

Start one apex at its physical core birth, in its live Palm environment.
For `hC` register the older tip at its promotion, and from there register
each *subsequent* tip promotion as a completed pair. For `Yt` first register
each trunk at its promotion, then the older and younger tips of that trunk
at their own promotions. A gain that turns a leaf into a core changes the
sender/parent's census at fixed degree and births a core at `(2,0)`; it is
this edge birth, not a later send along the core edge, that extends a
prefix. At every eligible promotion, one successor keeps the old prefix
with the newborn as a side arm; one successor extends the prefix with the
same newborn as a tag. Their continuation values add because these are two
rewards on the same physical jump, not two clocks or independent copies.
For every pair of physical tip IDs, the older tip is the unique prefix
choice and the younger tip is the unique completion. Canonical sorting is
used only to query the signed row; it contributes no factor of two even
when the census labels coincide. The replay checks this with interleaved
trunk/child births, and checks swapped-tip coefficients against the row
source. Registering at actual births gives the preexisting side children
their true elapsed lifetime at the last-edge birth, with no restart at
completion.

Only a complete prefix earns the entire masked `scan.Row.coefficients`
reward, including the intercept, all received exterior slopes, and all
recharges already in the scanner row. The main solve pairs coefficients
at each θ *before* any physical rate extremization. Source diagnostics
split the intercept, apex-parent and side-child terms; a separate census
ledger splits by exact fork `(D,k)=(3,2),(4,2),(4,3),(5,2),(5,3)`.
Those ledgers describe the same reward, so they are not added back to the
coupled interval. The code checks both partitions' summed intervals contain
the coupled certificate. On this mask the `hC` exterior slopes vanish;
for `Yt` at θ=1 the `(3,2)` fork subinterval alone has width at most
`3.857018488e−6`, larger than any other census. The parent and side-child
source projection widths are `2.382406220e−6` and `3.258366266e−6`.

## Live rates and irreversible boundaries

The recursion starts at apex core birth, not at a typical masked label.
No environment is regenerated when a role becomes tagged. Tagged roles
retain their actual censuses, tagged internal reciprocals, and actual
leaf-promotion clocks. With `--scope fork`, every side child of the fork
is recorded at its own promotion with its physical `(D,k)` and follows
both gain and promotion while within its donor cap. A send *from* a side
child changes its receiver, not the side child's reciprocal; its gain
does reduce its reciprocal. Once a side child crosses degree three its
reciprocal stays in `[0,1/4]` indefinitely by monotonicity. Side children
of other roles stay represented by their live counts with `[0,1/2]` per
child, including nonvirgin side children. This relaxes their rates without
restarting their dynamics. `--scope all` retains donor state at every
role for smaller masks. The apex parent starts with an arbitrary live
degree under the birth Palm law. The principal solve monitors degree two
and uses `[0,1/3]` above it; it does not assume a virgin parent or an
independent initial degree.

For each column and every retained state the replay verifies the integer
sign of both discounted drift residuals, using the **same** exterior rate
in the holding denominator, gain continuation and row slope. Interval
controllers can depend on the full unseen past; no assumption of conditional
independence of side branches is used. For the unknown common observation
discount, use `[.565673,.565775]` and multiply the discounted birth-Palm
integral by the corresponding stationary-age normalizer, rounding outward.
The two payloads' own discounts inside their charge rows are left intact.

The zero continuation at a tagged role's first cap crossing is exact for
the masked reward because each tagged `D` and `k` never falls. A prefix
whose required parent's arity is already at the cap also cannot complete
inside the mask. These zero boundaries make no claim about signed rows in
the complement. By contrast, when an **exterior** donor crosses its cap,
its future influence is paid by its perpetual tail interval. All other
unobserved jumps are integrated into adapted rate bounds. The counted
incidence has bounded multiplicity under the mask; these facts give the
usual discounted drift comparison for its full remaining lifetime.

## Authenticated replay

`863c` was not an available Git object in this checkout. The replay pins
the available `all-contrast/scan.py` by SHA-256
`6f2a074e9a695b4969cb58d39376e2a1f8c17ab41b04a3b2ec775a18f05ca02d`.
Its frozen-query checker passed all 41 independently regenerated physical
queries, checking both signed affine rows and both θ pairs. Each replay
hashes both original D18 payloads against `extract.PAYLOAD_SHA` before
using the rows. It also requires exact equality between the paid row set
and the canonical physical mask, verifies drift residuals, and checks
that the matched memory certificate nests in the same mask without memory.

From `.`, run the command below once for each family, replacing
`Yt` with `hC` for the other run:

```sh
taskset -c 7 nice -n 10 timeout 85s /usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M KiB status=%x' env OPENBLAS_NUM_THREADS=1 round2/.work/venv/bin/python -Bu round2/occupation-neon/pressure-memory/fork-lifetime/replay.py --family Yt --degree 5 --arity 3 --donor 3 --mother 2 --scope fork --ledger --matched
```

The completed `hC` replay checked 221,312 memory and 5,682 matched-control
residuals (`13.01` s wall, `5.64+0.83` s CPU, `70,604` KiB RSS). The `Yt`
replay checked 1,992,644 memory and 51,204 matched-control residuals
(`78.51` s wall, `28.06+1.54` s CPU, `329,172` KiB RSS). Both exited zero.
A smaller full-role donor-age check at `D≤4,k≤2` with donor/mother caps 2
covered 72 `hC` and 360 `Yt` rows and reduced the `Yt` θ=1 width from
`3.684576509e−6` to `1.899283317e−6`; the larger mask is a different
quantity, so that smaller width is not booked as a further improvement.
