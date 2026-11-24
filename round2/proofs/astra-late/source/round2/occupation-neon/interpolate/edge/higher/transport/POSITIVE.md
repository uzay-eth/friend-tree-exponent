> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../../README.md#claim-21) for exact scope and citations.

# Both-witness positive P4 candidate against f7cc2e1

**Candidate numerator floors per core: `8.99765651647766e−7` lower and `1.66895371114237e−7` upper.** Use the matched D18–120 row only. No new exponent endpoint or release. Baseline `f7cc2e1c77192ecea8f231519f4afd2bf50aab9b`.

| Cut `Da ka Db kb Dc kc Dt kt Ddon kdon` | Witness | Exact root interval ticks (`10⁻¹⁸`) | Physical interval per core (outward decimal) |
|---|---|---|---|
| `80 6 6 3 5 2 8 2 10 2` | lower | `1234658350768 2894236255229` | `[6.98412893253e−7, 1.63748651731e−6]` |
| same | upper | `223502327670 581796577711` | `[1.26429232200e−7, 3.29165958755e−7]` |
| `120 7 6 3 5 2 12 3 10 2` | lower | `1590610921235 3685563273776` | `[8.99765651647e−7, 2.08519956123e−6]` |
| same | upper | `295038602009 760905234812` | `[1.66895371114e−7, 4.30501159226e−7]` |

Multiply each lower root tick by `.565673/10^18` and each upper root tick by `.565775/10^18`. The exact D120 candidate lower floors are `179953130329553231/200000000000000000000000` and `166895371114237057/1000000000000000000000000`. The D80 alternative is nested; do not add it. The lower floor misses `1e−6` by about `1.0024e−7`; its physical ceiling does not rule out `1e−6` on this mask. The upper physical ceiling rules it out on **both tested masks**, by at least `5.6949e−7` even on D120. The wider lower-only exploration `120 7 7 3 6 3 14 3 10 2` reached 100,675,760 rows before the tool timed out and supplied no bound.

## Exact support and budgets

The paid atom is an oriented top-owned nonfounder `P4` with an actual terminal core child donor. It uses apex `D≥18,k≥2`, the declared monotone path caps, and donor arity ≤2 for both witnesses. Track each actual tagged donor from its promotion; above donor degree 10 retain its arity and bound its reciprocal in `[0,1/11]`. No witness-specific support selection or signed balancing is used. For each witness, each terminal off-path child sends the nonnegative excess `max(0,−b)/2 + b/D_donor` above the charged row's independent-port minimum; `b` is that witness's exact signed terminal-child coefficient. The exact row extraction checks the full boxed affine minimum is nonnegative, so this atom is a pathwise portion of the original `P4` charge, including when there is more than one terminal donor. The live donor rate shifts the same receiver and the payoff. Promotion produces one genuine new tagged prefix and also advances the old prefix, with post-promotion values. See [transport recurrence](README.md) for the proof and integer rounding.

At payment there are four path cores. The apex has at least one other core child, and the terminal has at least one core child; these are distinct from the path and each other. Hence the complete planted apex fringe has at least **six cores**. Its size is at least `Σ D_i−3`: count the four tagged cores, all their immediate leaves and every distinct off-path core child once. Since the degrees are at least `(18,2,2,2)`, this is at least **21 ordinary vertices**, excluding the apex's parent. The replay explicitly checks the size/core exclusion on every paid descriptor. Thus it misses the actual size-12 full-fringe payment (`round2/occupation-neon/interpolate/higher/README.md`, historical reference) and size-13–20 suffix (`round2/occupation-neon/interpolate/renewal/next/README.md`, historical reference) even without using the suffix's `D≤16` condition. It misses the `e911a1c` original P4 band (apex ≤12), six-part complement (apex ≤9 in the frozen run), and the in-progress path-fast low-apex replacement (apex ≤16). The `f7cc2e1` additions themselves spend `hV`, not `P4`. Other unbooked high-apex P4 results may overlap; they are not being combined here.

## Frozen replay

Source SHA-256, relative to this directory:

```text
a7b00b02e8933029d5a5d103f6e45a68cc4950d27a6b07ea6a4a7b79fe716ead  sweep.cpp
a92a5e830ac6ed092eeb2be1c73709244f66aef9df1a227cf969785467a7e870  replay.py
9d819bd1fd40efb5bd8b0ffd37de64fed074319473545e7f8100696cb7bec07f  check.py
```

The Python entry pins the frozen EDGE extractor/dependencies and checks both original payload SHA-256s. From the repo root:

```sh
ulimit -v 10485760
export OPENBLAS_NUM_THREADS=1
P=./round2/.venv/bin/python
R=round2/occupation-neon/interpolate/edge/higher/transport
B=<scratch-path>
g++ -std=c++17 -O3 -Wall -Wextra -Werror "$R/sweep.cpp" -o "$B"
nice -n 15 taskset -c 4 "$P" -Bu "$R/check.py" "$B"
nice -n 15 taskset -c 4 "$P" -Bu "$R/replay.py" b3f18_lo.bin b3f18_hi.bin --native "$B" --config 80 6 6 3 5 2 8 2 10 2 --objective both
nice -n 15 taskset -c 4 "$P" -Bu "$R/replay.py" b3f18_lo.bin b3f18_hi.bin --native "$B" --config 120 7 6 3 5 2 12 3 10 2 --objective both
```

The exact-corner checker passed 6,610 P4 and 7,020 Yb reference states. The D80 lower completed with 15,375,661 rows before the paired job's host timeout on CPU11; the standalone upper then finished in 257.36 s on CPU4 with the same row count. The full D120 replay finished both sides on CPU4 in 1,866.04 s / 304,724 KiB peak RSS, with 62,259,407 rows, 2,125,004 nonzero birth sources and 58,935,040 bytes of native value slabs **per side**. CPU3/7 were not used for these replays. Independent cold review of production operands and support is pending.
