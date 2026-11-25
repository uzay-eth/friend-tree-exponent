# The non-leaf growth exponent of the random friend tree

**Can Uzay** · [Paper (PDF)](paper/friend-tree-exponent.pdf) · [Supplement: complete proofs (PDF)](paper/friend-tree-exponent-supplement.pdf) · [Statement map](docs/proof-map.md) · [Check it yourself](#check-it-yourself)

Start with an edge. At each step, pick a vertex uniformly, pick one of its neighbours uniformly, and attach a new leaf to that neighbour. The non-leaves form the **core**, which Krapivsky and Redner call the *nucleus*. $`C_n`$ counts the non-leaves when the tree has $`n`$ vertices.

A leaf becomes a non-leaf only when its one neighbour is chosen and passes it the new leaf. So two trees with the same number of non-leaves can grow at different rates. For example, two adjacent non-leaves with 1 and 3 leaves give a promotion probability of $`5/24`$ at the next step, while 2 and 2 leaves give $`2/9`$.

**Results.**
- $`\log C_n/\log n`$ and $`\log\mathbb E C_n/\log n`$ converge to the same exponent $`\mu`$.
- $`C_n/n^\mu`$ converges almost surely to a positive random variable with a density.
- Exact computer-assisted certificates give $`0.565714341933\le\mu\le0.565736986679`$, rounded outwards. This excludes $`17/30`$, the value Strogatz asked about, and every rational with denominator below 426, but not $`241/426`$.
- The exact value of $`\mu`$, and whether it is rational, are open.
- No human referee has checked the proofs, which were developed and cross-checked with substantial help from AI models.

The paper states each result with the argument behind it and points to a complete proof in the supplement. Every result has a persistent tag, which the paper, the supplement and this repository all use.

| Tag | Result | Proof |
|---|---|---|
| EXP | common exponent in mean and almost surely; $`113/200\le\mu\le1133/2000`$ | supplement §EXP, §MOT3 |
| ENC | $`0.565714341933\le\mu\le0.565736986679`$; any rational $`\mu`$ has denominator $`\ge426`$ | supplement §ENC |
| AMP, DEN | $`C_n/n^\mu\to X`$ a.s., $`0<X<\infty`$, with a density | supplement §CORE, §DEN |
| DEG | fixed-degree fractions of the core converge a.s. to $`p_k\sim\mu Ak^{-1-\mu}`$ (size limit first) | supplement §CORE |
| MASS | limiting degree shares sum to one a.s. | supplement §MASS |
| ORD | founders' shares agree in law; strict stochastic order from label two onward | supplement Part II (birth-order module), §ORD |
| PREC | $`\mathbb E C_n`$ is not P-recursive; says nothing about whether $`\mu`$ is rational | paper, Appendix PREC |

These answer Questions 1–4 of [Addario-Berry et al., *Random friend trees*, §7](https://arxiv.org/html/2403.20185v1#S7). The fixed-degree law takes the size limit before the degree limit; exact-degree laws when the degree grows with $`n`$ are proved only for means (supplement §NOTES) and are otherwise open.

## What is new here, and what came before

- **Krapivsky and Redner** estimated $`\mu\approx0.566`$ from fitted slopes between 0.5652 and 0.5673 ([2017, §5.1](https://arxiv.org/abs/1706.01514)).
- **Steven Strogatz** publicized the question in a [2022 post](https://x.com/stevenstrogatz/status/1572287647715934210) relaying Redner's account: “it’s an open problem to prove that 17/30 is really the exponent.”
- **Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou and Lugosi** named the random friend tree, proved bounds on core growth and asked the questions of their Section 7 ([Random friend trees](https://arxiv.org/abs/2403.20185)).
- **[Christopher D. Long](https://github.com/long-mathematics/friend-tree-growth)** gave an earlier certified exclusion of 17/30, $`\mu\le0.5664485422`$, with a [characterization of $`\mu`$ by discounted generation sums and a terminating two-sided approximation](https://github.com/long-mathematics/friend-tree-growth/tree/090bbbe5faec83ec4124bd427eb5fda32c79bef1). He also proved $`\mathbb E C_n\asymp n^\mu`$.
- **Mohamed Sabba** reported the interval $`[0.565615,0.565809]`$ and an estimate of 0.5657237 ([interval](https://x.com/ma_sabba/status/2124768206141206762), [estimate](https://x.com/ma_sabba/status/2124791407667716208)).
- **Hartle, Krapivsky, Redner and Zhang** reported a finite-size fit of 0.569 ([Fig. 6](https://arxiv.org/html/2604.01540v1#S3.F6)), outside the certified interval.

**New here:**
- the enclosure to width about $`2.3\times10^{-5}$ and the denominator bound;
- the almost-sure limit $`X`$ and its density;
- the fixed-degree law;
- unit hub mass;
- the birth-order module;
- the non-P-recursiveness of the mean.

The three-vertex bounds of Theorem EXP first appeared in the author's [earlier preprint](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/paper/friend-tree-exponent.pdf).

## Check it yourself

From the repository root on x86-64 Linux:

1. **Compact replay.** `bash verify/compact_replay.sh` needs Python 3.10+ and `g++` with C++17. It runs:
   - the exact checks behind PREC;
   - the receipt arithmetic behind the lower endpoint of ENC, with planted defects;
   - both stochastic-order certificates of ORD, with their audits and controls;
   - the check that the printed enclosure matches the receipts.

   It ends with `=== ALL COMPACT CHECKS PASSED ===`.
2. **Three-vertex bounds.** `bash verify/run_all.sh` checks both sides of the numerical bounds of EXP ([setup](verify/README.md)). It takes about 11 minutes and 2.19 GiB of RAM.
3. **Selected supporting results.** `bash round2/verify/run_all.sh` replays 42 selected results ([setup](round2/verify/README.md)). It takes about 6 minutes and 605 MiB.
4. **Five-vertex upper endpoint.** Rebuild and verify the pair of five-vertex tables as in the [rebuild guide](round2/verify/gold/README.md). Allow 10 GiB of RAM to build, 20 GiB to verify and 12 GiB of disk.
5. **Refined lower endpoint.** The 13 GB lower witness is not stored here. [verify/refined-lower](verify/refined-lower/README.md) explains how to regenerate it byte for byte and how to replay it with the undefined-behaviour sanitizer.
6. **Paper and supplement.** [Building the paper](paper/README.md) rebuilds both PDFs from the [source bundle](paper/friend-tree-exponent-source.tar.gz). The [editorial checks](paper/checks/README.md) compare printed numbers with receipts. They also compare the transcribed proofs with their sources, formula by formula, and check every link.

A passing program run checks finitely many inequalities. The reductions from all trees to those inequalities are proved in the supplement and are not checked by the programs.

## Repository

The [statement map](docs/proof-map.md) links each tag to its proof, sources, checkers and review scope. Earlier repository material is kept unchanged, including every proof, review and research folder cited by the paper. The [credit record](round2/who-did-what.md) and [review archive](round2/reviews/README.md) say which model wrote and which checked each part. Can Uzay directed the work and is responsible for the manuscript.
