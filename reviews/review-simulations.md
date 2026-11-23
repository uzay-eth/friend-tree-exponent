VERDICT: CONSISTENT / APPROVE (2 non-blocking notes)

# Cold review of the simulation package — Claude Opus 5.5

Reviewer: Claude Opus 5.5 (`claude-opus-5-5`, Claude Code harness), 2026-11-23. The review covers [`simulations/`](../simulations/README.md): the fast continuous-time (CT) generator, the fixed-size discrete generator, the archived trajectories, the analysis and fit scripts, the separate 48-run C program, and the claims in the simulations README. I wrote the [independent discrete cross-check](../simulations/cross-check/README.md) for this review without reusing any packaged code. I also read the separate [Node.js check](../simulations/independent-check/README.md).

## Verdict

The generators implement the stated growth rule. Every packaged number I recomputed agrees with its archived value. Every fresh comparison with independent code agrees within Monte Carlo error. Every negative control is rejected. Neither note below changes a published number.

## What I ran

| Check | Command or input | Outcome |
|---|---|---|
| Checksums | `sha256sum -c SHA256SUMS` | all files OK |
| Demo | `bash demo.sh` on a copy of the package | exit 0 in 2.70 s at load average about 12. The archived summary, checkpoint table, fixed-size summary and both fit tables were reproduced. |
| Pooled estimate | `analyze_ct.py 21` on the three archived batches | $`0.5657230\pm0.0000032`$ and, with raw conversions, $`0.5657253\pm0.0000032`$, byte-identical to `outputs/sim_summary.txt` |
| Tables and fits | `local_slopes.py --table`, `fit_window.py --table` | byte-identical to `local_slopes.tsv`, `fixed-size/summary.tsv` and `finite_window.txt` |
| 48-run C program | `gcc -O3 -fopenmp`, full replay | stdout byte-identical to `route-exact/replay.txt`. It gives $`0.5657242\pm0.0000175`$ at $`t=27`$ in 128.72 s with a peak of 1,252,068 KiB. |
| Exact small trees | independent enumeration, $`n\le21`$ | all tree counts and rational $`\mathbb E[C_n]`$ equal the packaged table; 2,144,505 trees at $`n=21`$, 987 s |
| Direct chain vs exact | $`10^8`$ trajectories to $`n=21`$ | max $`\lvert z\rvert`$ 1.36 over the means of $`C_n`$ and $`S_n`$. The variances also agree. |
| Fast CT vs direct chain | $`10^6`$ fresh `friend_ct` runs to $`t=12`$ vs the direct chain at exact Yule size draws | max $`\lvert z\rvert`$ 0.94 for $`\mathbb E[C(t)]`$ and $`\mathbb E[S(t)]`$ at $`t=2,4,\ldots,12`$ |
| CT population law | the same $`10^6`$ runs, plus $`10^7`$ runs to $`t=1`$ | mean population within 0.6 standard errors of $`3e^t`$ for the $`10^6`$ runs and 1.02 for the $`10^7`$ runs. The archived batches are within 1.14. |
| Fixed-size generator vs direct chain | 9,600 fresh seeds of `fixed_size.cpp` vs 12,000 direct trajectories at the nine published sizes | max $`\lvert z\rvert`$ 0.72 |
| Negative controls | attach to $`V`$; degree-biased neighbour; deterministic size $`3e^t`$ | rejected with $`\lvert z\rvert`$ from 39 to 7,445 |

The cross-check README gives the commands, seeds, timings and raw aggregate outputs for every row from the exact enumeration down.

## Code reading

**`sim/friend_ct.cpp`.** The lazy leaf update is exact. A Yule process started from $`L`$ individuals has $`L+\mathrm{NB}(L,e^{-\tau})`$ individuals after time $`\tau`$. The code draws a sum of geometric variables when $`L\le12`$ and a gamma–Poisson mixture otherwise, which is the same law. The Marsaglia–Tsang gamma sampler and the PTRS Poisson sampler are the standard exact rejection methods. The documented rounded-Gaussian branch for intensities above $`10^9`$ is the only approximation. A core firing selects a uniform core vertex, brings its leaf count up to date, and picks a uniform neighbour among $`L_u+c_u`$. A leaf choice converts that leaf into a core vertex with one leaf and one core neighbour. A core choice adds a leaf to the chosen neighbour. The neighbour index walks parent first, then the children list, and covers each core neighbour exactly once. Checkpoints update every vertex before recording, and the memoryless holding time is reused correctly across checkpoints.

**`fixed-size/fixed_size.cpp`.** A leaf stores its only neighbour. A nonleaf stores the index of its adjacency list. Each step draws a uniform vertex and then a uniform actual neighbour. When the target was a leaf, its list is created with both neighbours. The recorded core is the number of lists, which is the number of nonleaves.

**Stated identities.** In CT, the population is a Yule process from three individuals, independent of the jump chain. Hence $`N(t)-3`$ is negative binomial with shape 3 and success probability $`e^{-t}`$, and $`\mathbb E[C(t)]=\sum_n P(N(t)=n)\,\mathbb E[C_n]`$. A leaf's only neighbour is already core, so only core firings onto leaves create core vertices, and

```math
\frac{d}{dt}\,\mathbb E[C(t)]=\mathbb E[S(t)].
```

The fast-vs-direct comparison tests both statements without assuming either. The direct side never uses CT and observes $`S`$ independently.

## Notes

1. **The archived fixed-size sample sits about one standard error low.** At all nine sizes the 2,400 archived trajectories lie below the 12,000 direct trajectories, with $`z`$ between $`-1.62`$ and $`-0.81`$. The same trajectories contribute to every size, so this is one fluctuation. Rerunning the unchanged generator on 9,600 fresh seeds gives $`\lvert z\rvert\le0.72`$. A common relative offset leaves log–log slopes unchanged, so the fitted exponents are unaffected. No change is needed. The archived means are simply about 1.2% below the population means.

2. **The replay table mixes measured and projected times.** The "Projected full batch" column and the sentence beginning "Allow roughly four hours" extrapolate from partial replays of 9 to 51 trajectories. Tree sizes vary strongly between trajectories, so these figures are estimates rather than measurements. I suggest keeping the measured partial times and trajectory counts, and stating that full regeneration was not timed. The packager has said these projections will be replaced by measured full-batch replay results. That replacement was not available when this review was written, so I have not checked it.

The historical Python mixture reference gives $`3.82579\pm0.02720`$ at $`t=2`$ and $`12.12375\pm0.09352`$ at $`t=4`$. These lie 1.34 and 1.50 of their standard errors below the direct-chain values of this review. That fits chance and its stated errors, and it needs no change.

## Scope

This review tests the simulations as Monte Carlo estimates at the stated sizes and horizons. It does not test whether any finite-horizon slope equals the limiting exponent; the proofs address the limit separately. I did not run the full regeneration of the three archived production batches. I relied on the byte-level replay of each batch's first trajectory reported in the package, the checksums, and the agreement of fresh fast-generator runs with independent code.

The package does not record which model wrote the original CT generator. If it was a Claude model, this review is not cross-family for that file. The separate Node.js check, written by GPT Codex, then supplies the cross-family comparison. I read it and reran its exact-value mode, which gave a maximum $`\lvert z\rvert`$ of 1.167 as stated.
