> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../../../README.md#claim-21) for exact scope and citations.

# Frozen EDGE plus its unspent first-body and edge surplus

**Candidate outward bracket: `[0.565678056632067, 0.565770469400619]`.**
The full producer replay and 7b08's independent calculation agree on both
exact fractions. This is a calculation for review, not a board update.
The baseline is the exact [02d0bd8 EDGE freeze](../../../../edge/FREEZE.md).

## Physical identity and support

Use the unchanged original GMP allocation at each witness's `λ`. Let
`Q_x=q(x)+t(k−h)` be the frozen conservative all-eight potential charge.
Then `Φ≤ΣQ_x`, `Q_x≤M`, and the recharge sums to zero. Write `δ=ε(s−λ)`.
The original signed drift at `s` differs from that at `λ` by `−δΦ`.
Use the exact same correction from 02d0bd8, with census value `z`, edge
value `V`, newborn census offset `−z(J)`, and founder constant `−δM/s`.
Pointwise, the new signed drift decomposes as

```
ε(ℒ−s)Φ_EDGE = m(C−1) + Σ_nonfounder R_f + Σ_edges R_g
            + old reserved g and higher charges + unbooked nonnegative terms

R_f = F + δ(M−Q_x) + (A_ρ−s)z + b V(post-parent,J),
R_g = G − G_reserved + (A_edge−s)V,
m   = s z(J)−δM.
```

The unbooked terms include `δ(ΣQ_x−Φ)`, the unused founder shift
`δ(M−Q_founder)`, and unused original charged rows. The edge sum here is
over existing edges with nonfounder parent. Their adjoints move the old
surviving edge; they do not include the edge created at promotion. The
new edge's value is paid **once** in `R_f`. The conserved shift constant
and the founder have vanishing normalized contributions. The original
signed margin `m` and stationary reserve are carried forward once, even
though `m<0` at both trials. This does not price original first-body `F`
again or spend a higher-row credit from the reserve.

The default objectives pay `R_f` only at `D≤18` (170 rows) and `R_g`
only when both endpoint degrees are at most four (54 edge classes).
Every priced edge is both-inner, and its full uniform inner reserve is
subtracted. No high-parent mixed `g` or higher-family observation enters
the new price. The four edge corners and both first-body receive corners
are checked exactly. Outside the masks the objective is zero; the
unbounded arity, reciprocal and tail balances remain present. The
two LPs have 9,674 identical columns after removing just the 20 memory
rows. The shared receive rate enters the same census balance and `f`
coefficient; the edge uses the same two endpoint receive moments.

## Exact composition and scale

Use the frozen denominator upper `D` and reserve lower `Q` unmodified.
The checked new stationary reward lower `r` yields
`s+ε(m+Q+r)/D`, with positive `m+Q+r` on both sides. The script builds
all three quantities before checking that `s+ε(m+Q)/D` equals the
baseline's **exact** fraction (not its displayed outward endpoint).

| Witness | Flow-only reward lower | Memory reward lower | Exact inward gain, approx. | Reward needed for `1e−6` move | Multiple of this price |
|---|---:|---:|---:|---:|---:|
| Lower | `9.14435e−8` | `1.03096e−7` | `7.59460e−8` | `1.35749e−6` | 13.17× |
| Upper | `5.13399e−8` | `5.99909e−8` | `4.38793e−8` | `1.36718e−6` | 22.79× |

The present estimator does not demonstrate a `1e−6` move. The earlier
family ablation recovers only about `1.40e−9` lower / `1.08e−10` upper
by adding the small exact-edge mask to `f` with memory. Reaching the
target would require roughly another `1.25e−6` / `1.31e−6` of certified
unspent reward per core. Larger edge masks with conditional donor
constraints are a possible source; neither these prices nor the table
establish a ceiling on what those masks could recover.

Exact candidate endpoints before rounding:

- Lower: `2586244278310822397274177587149/4571936719109803018500000000000`.
- Upper: `6574835699364880591816419305857/11621030179129551777500000000000`.

## Authenticated replay

The manifest pins the exact 02d0bd8 inputs and both projection interfaces,
including the single root `pressure-memory/marginal/replay.py`. The
script hashes both full original payloads, rebuilds each frozen EDGE
witness, checks each rational dual and recomputes the two exact endpoints.
Floating proposal failures trigger a different proposal, not an unchecked
price. No frozen source or board file is edited.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 15 /usr/bin/time -f 'wall=%e cpu=%U+%S exit=%x' round2/.work/venv/bin/python -Bu round2/occupation-neon/interpolate/compose/residual/edge-projection/combined/replay.py --lower b3f18_lo.bin --upper b3f18_hi.bin
```

[`full.log`](full.log): both signs passed, `218.46 s` elapsed,
`103.05+1.87` CPU seconds, exit `0`.
