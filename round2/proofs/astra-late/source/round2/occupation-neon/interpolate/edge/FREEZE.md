> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Two-sided signed EDGE freeze

**Fresh clean-archive replay, 2026-11-24 UTC:** the signed EDGE composition certifies the outward enclosure **`[0.565677980686114, 0.565770513279953]`**. The prior release is `80fe02660b9f0a0d4e7fdc8240b0815615ff5896`. The original [producer note](README.md) is retained byte-for-byte; its review status is superseded by the [7b08 cold report](../../../../../../../reviews/cold/astra/checks/neon-edge-sol.md), verdict **CORRECT in scope**, and the fresh replay below.

The trial-point margins are negative on both sides. The enclosure uses the positive **net stationary numerator** after adding the disjoint reserved rewards to each signed margin. It inherits the original D18 GMP B verifications, physical mean and stable-age identities, and input interval `[0.565673570322, 0.565775]`. It does not rerun exhaustive GMP, include the separate memory prices, or certify eight decimals. The companion independently reconstructs charge rows, checks physical generator controls and composes the endpoint; it shares the producer’s exact Bellman solver.

## Clean source and payload authentication

The clean tree was extracted with `git archive 80fe02660b9f0a0d4e7fdc8240b0815615ff5896`, then overlaid with exactly seven files from the [SHA-256 manifest](inputs.sha256): the three declared producer files, both cold-review files, and `fringe/check.py` and `all-contrast/check/row_oracle.py` at their canonical paths. The 17 other manifest files matched the archive tree and Git bytes. Both that manifest and the seven-entry parent composition manifest passed `sha256sum --strict -c` from the archive root. The new fringe helper hash identifies the bytes used in **this** replay; no historical runtime hash is asserted.

The producer authenticated both external original D18 payloads in full. Use those same read-only files for the companion and regression. No memory-pricing source is required.

## Commands from the clean root

```sh
sha256sum --strict -c round2/occupation-neon/interpolate/edge/inputs.sha256
sha256sum --strict -c round2/occupation-neon/interpolate/compose/inputs.txt
ulimit -v 10485760
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
P=./round2/.venv/bin/python # Python 3.12.3, NumPy 2.1.3
R=round2/occupation-neon/interpolate/edge
A=round2/checks/neon-edge-sol.py
nice -n 15 taskset -c 10 timeout 200s "$P" -EBu "$R/replay.py" b3f18_lo.bin --degree 40 --arity 8 --challenge --mask disjoint --s .565677930574
nice -n 15 taskset -c 11 timeout 200s "$P" -EBu "$R/replay.py" b3f18_hi.bin --degree 32 --arity 7 --challenge --mask disjoint --s .565770634792
nice -n 15 taskset -c 8 timeout 30s "$P" -EBu "$R/regression.py" b3f18_hi.bin
nice -n 15 taskset -c 10 timeout 175s "$P" -EBu "$A" --full b3f18_lo.bin
nice -n 15 taskset -c 11 timeout 175s "$P" -EBu "$A" --full b3f18_hi.bin
```

| Run | Fresh result | Wall | Exit |
|---|---|---:|---:|
| Producer lower | 91,732 Bellman states; 119 whole-tree checks and 391 nonzero edge births; outward `0.565677980686114` | 66.18 s | 0 |
| Producer upper | 44,492 Bellman states; 119 whole-tree checks and 391 nonzero edge births; outward `0.565770513279953` | 43.24 s | 0 |
| Producer regression | Exact withdrawn seven-core residual reproduced, approximately `−0.0008223761382507083` | 2.08 s | 0 |
| Companion lower, full | 2,079 first-body rows, five stock floors, 43,776 inner charges, 91 port checks, 184 mixed reserves; omitted-source control gives three negative rows; outward `0.565677980686114` | 61.91 s | 0 |
| Companion upper, full | Same census, floor, inner and port counts; 105 mixed reserves; omitted-source control gives 88 negative rows; outward `0.565770513279953` | 42.32 s | 0 |

Exact fractions before outward rounding:

- Lower: `191976090122981726641370633087/339373453939523193250000000000`.
- Upper: `24172191945911554486185483421/42724375658564528593750000000`.
