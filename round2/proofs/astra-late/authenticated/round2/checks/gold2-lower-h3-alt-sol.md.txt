# Cold audit: alternative original Gold2 lower H3

**GAP in both advertised 8790 prices.** The same-cut payment has exactly the support of 2417’s frozen lower H3 and could only replace it. The added mother-80–100 seam is genuinely disjoint. Both have verified, smaller floors after correcting the reward sign:

| Actual nonfounder mother → middle → child mask | Verified exact ticks/core | Approx. |
|---|---:|---:|
| Mother `19..79,c≤8`; middle `D≤5,c≤4`; child `D≤18,c≤3` | `5481957648397036499283/13005908321239040000000` | `0.4214974850656785` |
| Mother `80..100,c≤8`; same middle; child `D≤5,c≤3` | `229762209993424478073686403/16257385401548800000000` | `14132.789764062283` |
| Disjoint sum | `919076249761939897477242027/65029541606195200000000` | `14133.211261547349` |

## Reward correction and scope

The physical outward row is `P=A+Σ b_i ρ_i`; the payment is `[−P]_+`. In pinned `higher/cohort.py:244`, `vg-(-GRID*b//1)` is `vg+ceil(GRID*b)`. The affine reward branch needs `vg+floor(−GRID*b)` because the same gain rate `ρ_i` multiplies both the continuation and the *negative* row slope. The original residual check therefore checks a different objective even though it checks its corners exactly. In particular a positive slope can erroneously increase the reward at the largest exterior rate. The independent code uses the previously audited compact row evaluator and newly builds and enumerates the rate corners directly; it does not use the alternative producer's sweep. For each fixed set of successors it takes the maximum of the zero and correctly rounded affine lower branches. Rounding the intercept and slopes downward separately is conservative.

The pay masks are on actual current marks, not on tail corner labels. `D_mother>D_child` is automatic on both masks; the explicit threshold exclusion is redundant in `80..100`. The latter crosses the G3 endpoint seam `87→88`; it does not move the middle outside M3. At `D≥88` the endpoint interpolation acts only on the original row weights. Physical clocks keep the actual mother's degree and arity. Literal gain arrows and promoted newborns retain the paired M3 edge transfers, subtract supported outgoing P4/S3 transfers at pre-gain marks, and keep newborn weights at the post-promotion marks. The original local production cap still covers `P+[−P]_+`; these floors may therefore compose with the original singleton reserves, D18 edge reserve, FY gaps, and pointwise higher host-cap savings that retain the original positive H3 allocation. No sibling/upstream decomposition or second overlapping H3 credit is included.

The independent recursion retains actual post-promotion marks in the older prefix and newly created prefix, including the correct third-slot newborn. Extra children remain distinct. Degree and arity caps are permanent exits, and `c=D` has zero promotion rate. Summing over physical registrations cancels the internal seed terms; at each true core birth the only new unallocated value is the singleton at `(2,1)`. The physical stable-age birth rate is `μ`. Nonnegative values permit discounting at original Gold2 `HI` and multiplying by original `LO` once. These floors concern only the original **lower** witness.

## Separate issue: previously frozen H3 replay

The previously accepted frozen replay has a **different dimensional gap**. In `post-fy/replay.py:154–166`, `width` is the physical exterior rate width multiplied by `d`, the rate denominator. The signed reward variation `w` is expressed in ticks across the **physical** width. The shift must be `round(w*d/width)` (equivalently `round(−slope)` with the branch direction). The replay instead uses `round(w/width)`. The checker plants an exact example with reward `2−ρ`, `ρ∈[0,1]`, discount 1, continuation 0 and rate denominator 2: the true adverse value is `1/2`; the unscaled shift returns `3/4` and violates the physical residual. 2417 independently confirmed the dimensional error. It applies to the shared lower/upper H3 replay path; both derived endpoint claims need corrected accounting before reuse. This finding does not invalidate the original standalone Gold2 production interval or the independently verified FY and D18 credits. The prior three-credit proof at `99/175` is not established by these corrected floors.

## Evidence

Original source `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`; Gold2 lower SHA-256 `41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413`, rehashed through the pinned host reader on the same open inode. Producer `higher/cohort.py` SHA-256 `b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c`; frozen `post-fy/replay.py` SHA-256 `d8a9a8a65d6b652151c41de18fa7a0f4e8aec4a6ed8e89395412e4ddd5bc3c23`.

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 round2/.work/venv/bin/python -IBu round2/checks/gold2-lower-h3-alt-sol.py
```

Exited zero, reproducing the two corrected exact fractions above. It checked 1,514 literal arrows, 12 exact H3 rows and 12 rejected donor-denominator mutations, and compared the two compact row implementations at 21 marks across the cutoff and interpolation seam. The independent sweep counted 76,329 seam terminal states and 271,350 same-cut terminal states, including phase-valid zero-reward states. The producer's advertised seam rerun exited zero and reproduced its *uncorrected* arithmetic; a direct full producer rerun timed out at 220 seconds after 42 fixtures. Neither producer outcome is used as evidence for the corrected floor.
