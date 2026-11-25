# The non-leaf growth exponent of the random friend tree

**Can Uzay** · [Paper (PDF)](paper/friend-tree-exponent.pdf) · [LaTeX source](paper/friend-tree-exponent.tex) · [Statement map](docs/proof-map.md) · [Building the paper](paper/README.md)

Start with an edge. At each step, pick a vertex uniformly, pick one of its neighbours uniformly, and attach a new leaf to that neighbour. The non-leaves form the **core**, called the *nucleus* by Krapivsky and Redner. It starts empty, and $`C_n`$ counts its vertices when the tree has $`n`$ vertices.

The core grows like $`n^\mu`$ with

```math
\frac{3843726496355222669}{6794465353700000000}\le\mu\le\frac{127993535256830043}{226242120050000000},
```

or $`0.565714341933\le\mu\le0.565736986679`$ rounded outwards (Theorem 2). This excludes $`17/30`$ and every rational number with denominator below 426. It does not exclude $`241/426`$.

## Results

- **The exponent (Theorems 1–2).** The core grows with the same exponent in expectation and almost surely, and exact computer-assisted drift certificates give the enclosure above.
- **A random multiplier with a density (Theorems 4 and 8).** $`C_n/n^\mu`$ converges almost surely to a finite, strictly positive random variable whose law has a density.
- **Degree proportions (Theorems 4–5).** For each fixed $`k\ge2`$ the fraction of non-leaves of degree $`k`$ converges almost surely to a deterministic $`p_k`$, and $`p_k\sim\mu A k^{-1-\mu}`$ as $`k\to\infty`$ for some $`A>0`$.
- **Limiting degree shares (Theorem 9).** Each vertex $`v`$ has a limiting degree share $`Z_v=\lim_n D_v(n)/n`$. Almost surely $`\sum_v Z_v=1`$, and the next attachment's label law converges to these shares in total variation.
- **Birth order (Theorems 10–12).** The two founders' shares have the same law. From label two onward the shares decrease strictly in first-order stochastic order:

  ```math
  \mathbb P(Z_u > x) > \mathbb P(Z_v > x) \qquad \text{for all } 2 \le u < v \text{ and } 0 \le x < 1.
  ```

  Given any tree on $`n\ge2`$ vertices, before labels $`n+1`$ and $`n+2`$ are placed, the earlier label already has the strictly larger conditional tail at every such threshold. Hub probabilities and all positive moments also decrease strictly from label two onward.
- **Arithmetic of the mean (Theorem 15).** The sequence $`\mathbb E C_n`$ satisfies no linear recurrence with polynomial coefficients, so its generating function is not D-finite (in particular not algebraic). This says nothing about whether $`\mu`$ is rational.
- **Supplementary, in mean only.** When $`k\to\infty`$ and $`k/n\to0`$ jointly, $`\mathbb E N_{=k}(n)\sim\mu A\,\mathbb E X\,n^\mu k^{-1-\mu}`$ (Section 6.3). This is separate from the fixed-degree law above.

These settle the growth and limit in Question 1, the fixed-degree law in Question 2, the birth-order comparisons in Question 3 and the unit total share in Question 4 of [Addario-Berry et al., *Random friend trees*, §7](https://arxiv.org/html/2403.20185v1#S7).

## What is still open

- The exact value of $`\mu`$, and whether it is rational, algebraic or transcendental.
- Almost-sure or in-probability laws for exact-degree counts when the degree grows with the tree.
- Numerical values of the constant $`A`$ and of the tail constants of the random multiplier.
- The remaining parts of Questions 4–8 of Addario-Berry et al.

## Kinds of proof

Most results are proved by hand. The enclosure of $`\mu`$ and the birth-order results reduce every tree to finitely many local inequalities, which exact integer or rational programs check; a passing run checks the inequalities, not the reduction. Several limit theorems use the computer-checked bound $`\mu\ge113/200`$ as an input. Theorem 15 is analytic apart from five exact probabilities at size seven. AI models from three families developed and checked the proofs under Can Uzay's direction; no human referee has reviewed them ([Section 8 of the paper](paper/friend-tree-exponent.pdf), [credit record](round2/who-did-what.md)).

## Check it yourself

From the repository root on x86-64 Linux:

1. **Compact replay.** `bash verify/compact_replay.sh` needs Python 3.10+ and `g++` with C++17. It runs the exact checks behind Theorem 15, the receipt arithmetic behind the lower endpoint of Theorem 2 with planted defects, both exact stochastic-order certificates of Theorem 12 with their audits and controls, and the printed-enclosure check. It ends with `=== ALL COMPACT CHECKS PASSED ===`.
2. **Three-vertex bounds.** `bash verify/run_all.sh` checks both sides of the first enclosure of Theorem 1 ([setup](verify/README.md)). About 11 minutes and 2.19 GiB peak RAM on a Ryzen 7 5800X; it ends with `ALL 7 VERIFICATION SUITES PASSED`.
3. **Selected supporting results.** `bash round2/verify/run_all.sh` replays 42 selected results ([setup](round2/verify/README.md)). About 6 minutes and 605 MiB.
4. **Five-vertex upper endpoint.** Rebuild and verify the Gold2 pair as in the [rebuild guide](round2/verify/gold/README.md). Allow 10 GiB RAM to build, 20 GiB to verify and 12 GiB of disk.
5. **Refined lower endpoint.** The 13 GB lower witness is not stored here; [verify/refined-lower](verify/refined-lower/README.md) explains how to regenerate it byte for byte and how to replay it with the undefined-behaviour sanitizer.

To rebuild the PDF and the source archive, see [Building the paper](paper/README.md).

## Prior work

- **Saramäki and Kaski** studied the family of random-walk attachment models that includes this rule ([2004 paper](https://arxiv.org/abs/cond-mat/0404088)).
- **Krapivsky and Redner** studied the rule as isotropic redirection and estimated $`\mu\approx0.566`$ ([2017 paper, §5.1](https://arxiv.org/abs/1706.01514)). Their range of fitted slopes included $`17/30`$.
- **Steven Strogatz** publicized the open problem in a [2022 post](https://x.com/stevenstrogatz/status/1572287647715934210) relaying a conversation with Redner: “it’s an open problem to prove that 17/30 is really the exponent.”
- **Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou and Lugosi** named the random friend tree and proved rigorous bounds on core growth ([Random friend trees](https://arxiv.org/abs/2403.20185)).
- **Hartle, Krapivsky, Redner and Zhang** solved two modified redirection models and reported a finite-size fit of $`0.569`$ for isotropic redirection ([Fig. 6](https://arxiv.org/html/2604.01540v1#S3.F6)).
- **[Christopher D. Long](https://github.com/long-mathematics/friend-tree-growth)** gave an earlier [generation-sum characterization and terminating two-sided approximation theorem](https://github.com/long-mathematics/friend-tree-growth/tree/090bbbe5faec83ec4124bd427eb5fda32c79bef1) with the certified bound $`\mu\le0.5664485422<17/30`$, proved constant-factor mean bounds, and gave a [signed hypergeometric product identity for the original tree](https://github.com/long-mathematics/friend-tree-growth/blob/5bd25a7d506cb64a8444b86334ab531d295f6c78/README.md#signed-hypergeometric-products-for-the-original-tree) in his transcendence programme.
- **Mohamed Sabba** reported the interval $`[0.565615,0.565809]`$ and a numerical estimate of $`0.5657237`$ in replies to Strogatz’s thread. [Interval](https://x.com/ma_sabba/status/2124768206141206762) · [Estimate](https://x.com/ma_sabba/status/2124791407667716208)

Earlier repository material, including all proofs, reviews and research folders cited by the paper, is kept unchanged.
