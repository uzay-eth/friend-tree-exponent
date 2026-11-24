> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Exact reduced-cost pricing of the incidence LP

`replay.py` consumes `IncidenceModel`; it makes no changes to the model. For
the six-family D18 contrast (lower-witness reward minus upper-witness reward,
with `theta=1`), uniform cap four and order four give these certified
stationary enclosures per core. Decimal endpoints are rounded outward below;
the logs give exact fractions.

| Objective | No shared exterior donors | Shared through order four | Width reduction |
|---|---:|---:|---:|
| All six families | `[-9.259574551e−6, 9.436025665e−6]` | `[-8.081492576e−6, 8.268281238e−6]` | 12.55% |
| Same six families, model's `--unpaid` gate | `[-9.446033743e−7, 8.612153958e−7]` | `[-7.248356335e−7, 6.574359272e−7]` | 23.45% |

The matched models use the same objective and 4,899 balance rows. The
no-sharing model has 24,802 columns and 15,011 inequality rows; the shared
model has 70,342 columns and 30,014 inequalities. End-to-end generation on one
niced CPU took 40.18 CPU / 40.43 wall seconds for the full contrast and 37.60
CPU / 38.50 wall for unpaid. Those times include both signs, authentication,
model construction, no-share optimization, all pricing passes and output.
No-share IPM took about 5 s per sign. Shared polish at depths 1–3 took about
1–1.2 s of IPM per sign; the remaining time was model and exact-arithmetic
overhead. Repairs for the two shared full endpoints were `5.59e−15` and
`9.26e−15` before stationary conversion; unpaid `4.39e−15` and `4.53e−15`.
Both enclosures contain zero. These are bounds on contrasts, not new exponent
endpoints or optimality claims for the shared LP.

## Certificate and search

For minimization, quantize the equality and nonpositive inequality prices
to rational `y,z`. For every original column compute its exact residual
`r_i = c_i − (A_eq^T y)_i − (A_ub^T z)_i`. With the model's certified
nonnegative exact limits `B_g` on total nonnegative mass in each group,

```
b_eq⋅y + b_ub⋅z − Σ_g B_g max(0, max_{i∈g}(−r_i))
```

is a lower bound. Negate a checked lower bound on `−c` to get an upper
bound on `c`. The proof uses the original rows and exact group limits, not
the solver's rounded objective or reduced costs. The checker iterates over
all columns, including omitted, extension-only, zero-cost and off-block
columns. For stationary conversion it tries all discount/endpoint corners.

The successful route first solves the no-sharing LP by IPM. An unchanged
row is transported by column/row identity. For a previously unshared port
with old upper row `p − n m/2 ≤ 0`, let `t` be the new unseen reciprocal
ceiling and `q` a known donor reciprocal. Combining the new upper row and
`(1−2t)` times its paired lower row, then dividing by `2t`, recovers the
old row plus nonnegative coefficients on resolved masses, because
`t−q+(1−2t)q = t(1−2q)`. The transported certificate is always repriced
on the destination model; changed cutoffs never bypass that check.

`polish` holds deeper multipliers fixed and proposes corrections to the
equality and inequality prices through each requested depth. Group
penalties enter the sparse floating-point LP explicitly; untouched columns
provide floors on those penalties. The depth-three pass adjusts 3,206
rows touching 11,137 columns. Each accepted move is selected by its exact
all-column bound. `blocks` instead adjusts one topology at a time; the
tested depth-four sweep did not materially improve the contrast.

`cg` starts with the balance and objective columns, excludes extension-only
columns, and adds the worst negative candidates in batches. Artificial
variables on nonzero equality RHSs and negative inequality RHSs ensure
an initially feasible restricted master. On cap three, it took five/four
rounds lower/upper (2,055/2,096 original columns in the final masters out
of 2,671), with all-column pricing at every round. Cap-three full IPM and
column generation both match the original solver within `3.71e−14` in
unnormalized endpoints. On cap four, full shared IPM, restricted shared
IPM and explicit all-row group dual IPM all hit their budgets without a
SciPy multiplier result. Native restricted simplex did expose an interim
candidate, but its exact bound was worse than the fallback. We do not use
any of those failed attempts for the reported endpoints.

The zero certificate is always available and is exactly repaired (zero
multipliers do **not** imply a zero bound). A saved candidate is accepted
only against a fingerprint of the exact original model and signed objective,
and is repriced. For any search, the returned result is at least as good as
the checked fallback. The default `ipm` is useful for small/no-share models;
use `--seed 0 --method polish` for the uniform shared cap-four case.

## Replay

Benchmarks used SciPy 1.16.3 with crossover and single-threaded HiGHS.
`solve.py` at the measured cut had SHA-256
`206b93131081a128c549255070e98896fefab2076fd40f34d2c619536e0cfee0`;
uniform-model certificates are keyed to matrix content, not that file hash.
From the repo root, with the original authenticated payloads available:

```sh
P=round2/.work/venv/bin/python
R=round2/occupation-neon/embedding-flow/pricing
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
taskset -c 0 nice -n 19 timeout 90s "$P" -Bu "$R/replay.py" --lower b3f18_lo.bin --upper b3f18_hi.bin --cap 4 --seed 0 --share-through 4 --method polish --seconds 24
```

Add `--unpaid` for the other objective. For exact verification without
optimization, replace `--seed 0 --share-through 4 --method polish --seconds 24`
with `--share-through 0 4 --method verify --read-certificate "$R/cap4-full.json"`
(or `cap4-unpaid.json` with `--unpaid`). Both replayed exactly in `verify.log`,
at 16.24 and 15.36 CPU seconds. Cap-three parity uses
`--cap 3 --share-through 0 4 --parity --seconds 12`; for the separate column
generation check use `--cap 3 --share-through 4 --method cg --batch 250
--rounds 8 --seconds 12 --parity`. See `cap4-full.log`, `cap4-unpaid.log`,
`verify.log` and `cap3-cg.log` for successful receipts.
