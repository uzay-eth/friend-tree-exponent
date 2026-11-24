CORRECT

# Cold review: frozen five-pass quartet donor replacement

Scope. Reviews [`split-audit.md`](../five-next/split-audit.md), `split-reference.tsv` (`five-next/split-reference.tsv`),
`split-inputs.sha256` (`five-next/split-inputs.sha256`), `replay-split.sh` (`five-next/replay-split.sh`),
`split_integers.py` (`five-next/split_integers.py`) and the pinned source (`reserve.cpp`, `quartet.hpp`, `flow.hpp`,
`projection.hpp`). Conditioned on the separately reviewed raw certificate and mask63 facts. The billion-row production
replay (710/755 s per side) was **not** run; this is a short one-core check on the frozen transcripts and the pinned
small-split solver. No author file was changed and nothing was committed.

## 1. Frozen hashes, 50 counts and arithmetic

I recomputed all 25 SHA-256 digests in the manifest independently with `hashlib`: every pinned snapshot, source, log
and reference file matches. The transcript arithmetic is self-consistent: lower pass sum 29,374,764+5,855,805+
890,500+146,067+22,982 = 36,290,118; upper 26,942,162+4,889,420+732,894+119,891+18,704 = 32,703,071; increments
36,290,118−7,402,100 = 28,888,018 and 32,703,071−7,709,975 = 24,993,096, with the archived values taken from the
pinned `occupation-five-lower.txt` / `occupation-five-6.txt`. The 50 rank counts are 5 passes × 5 ranks × 2 sides and
match the frozen table exactly. `split_integers.py` also rejects a re-added archived reserve, a one-tick pass change
and a removed rank-five row; my checker independently rejects re-adding the baseline and changing the count convention.

## 2. Deletion, incoming/youngest and donor eligibility

On 398 real small friend trees I enumerated 1,038 connected chronological five-sets and every induced leaf deletion.
Each deletion gives a connected older quartet with a unique contact; 71 deletions remove the oldest member, as the
construction requires. Deleting the youngest of any connected set of size 2..5 recovers the unique older prefix, so
the incoming/youngest accounting is a bijection on these trees. The matched quartet family is the supported rooted
shape (`quartet.hpp` `deletions`), and the outside donor is a non-founder core neighbour; founder donors are excluded.
Sources from distinct deletions are distinct (five-set, donor, contact) incidences, so the sum over deletions charges
each original donor once.

## 3. Maximum physical weight and successive passes

The max-robustness step is the elementary identity max_j U_j(x′) ≥ U_{j*}(x′) for a maximizer j* at x; I verified it on
500 random candidate vectors, and `Hierarchy::combine` in `reserve.cpp` implements exactly `max(full, by_donor)` with
`full` from whole-source candidates and `by_donor` as the disjoint sum. The later passes use only the outside-gain
slack of the immediately preceding weights: `Sources` in `quartet.hpp` calls `prev->correction(n-1, part, x)` with the
signed gain on the maximum physical weight, and `donor()` applies the J difference with the true or pooled reciprocal.
The pass structure is visible in `reserve.cpp` (`prior=std::move(next)`), and the frozen totals confirm the telescoping
numerators. I did not recompute the full production passes (reserved for your slot); the small-split mode was run and
verified end to end.

## 4. No overlap with mask63 or the original replaced Q4

The replacement is reported as `replacement − archived`, never `replacement + archived`; the archived Q4 total is
subtracted once per side. The split uses the original quartet outside-donor incidences through `quartet.hpp` slopes,
whereas mask63's five added families use singleton/pair/edge/triple APIs (`reserve.hpp`, `credit.hpp`) and its fifth
family is an uncovered triple donor on a four-set. These are different source functions and different physical
incidences, so the split does not double-claim any mask63 source. The manifest pins the mask63-era inputs unchanged.

## 5. Signed integer Bellman bounds and physical controls

I compiled the pinned `reserve.cpp` in `small-split` mode and queried the exact scores on the 398 real trees and all
their successors, then independently computed the physical quartet source J(a,p)=max(a,0)/2−p·a from the frozen
certificate and checked the aggregated inequality source ≥ u·T−G+v·B with u=1133/2000 and exact Fractions. It holds on
all 398 trees for both witnesses (minimum slack 0). The per-row residual test in `flow.hpp`
(`need(residual(z)>=0)`) is the integer form of the same inequality after multiplying by 2000·2·∏D_j; my aggregation
checks its summed consequence. Planted errors are rejected: flipping the sign of a changes the computed source on 113
trees per side; dropping the 1/2 factor changes it on 123/124 trees; omitting an incidence exposes an under-count on
190 trees. The small-split solver uses UBSan in the author's replay; my pilot used the pinned source unmodified.

## 6. Verdict and limitations

Command, from the repository root on one core:

```
taskset -c 0 python3 -B round2/checks/five_split_hy4.py
```

Output: `FIVE SPLIT COLD CHECK: PASS` in about eight seconds. The full production replay, the six-site extension and
the high-degree residual are separate work and were not evaluated. The frozen lower run's known projection caveat is
disclosed in the audit and is not an arithmetic inconsistency in the pinned transcript. Within this scope no gap or
error was found in the frozen five-pass replacement.
