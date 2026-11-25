# The random friend tree exponent is not 17/30

**Can Uzay** · [Read the Round 3 paper (PDF)](paper/friend-tree-exponent.pdf) · [LaTeX source](paper/friend-tree-exponent.tex) · [Proof and checker map](docs/proof-map.md) · [Round 3 checkpoint (scoped status by question)](round3/README.md)

The exponent $`\mu`$ lies in the exact interval $`\frac{35083855339718511}{62017003637500000}\le\mu\le\frac{127993535256830043}{226242120050000000}`$, or $`0.565713486333\le\mu\le0.565736986679`$ rounded outwards (Theorem 2).

Start with an edge. At each step, pick a vertex uniformly, pick one of its neighbours uniformly, and attach a new leaf to that neighbour. The non-leaves form the **core**, called the *nucleus* in the 2017 paper. It starts empty, and $`C_n`$ counts its vertices when the tree has size $`n`$.

## New in Round 3: older vertices win at every threshold

Let $`Z_v`$ be the long-run share of all vertices that end up attached to vertex $`v`$ (its degree divided by the tree size). Addario-Berry et al. asked whether these shares are [stochastically ordered by birth label](https://arxiv.org/html/2403.20185v1#S7). They are. Apart from the two founders, whose shares have the same law, an older vertex is strictly more likely than a younger one to keep a share above any threshold:

```math
\mathbb P(Z_u > x) > \mathbb P(Z_v > x) \qquad \text{for all } 2 \le u < v \text{ and } 0 \le x < 1.
```

This holds even before the next two vertices arrive: given any tree on $`n\ge2`$ vertices, vertex $`n+1`$ has the strictly larger conditional tail than vertex $`n+2`$ at every threshold (Theorem 12 and Appendix D of the paper). The proof reduces every tree to finitely many inequalities, which two C++ programs check in exact 128-bit integer arithmetic.

- **Proof:** [Part 1, §§1–3](round3/q3-birth-order/PROOF.md) · [Part 2, §§4–10](round3/q3-birth-order/PROOF-PART-2.md) (one frozen proof, split only for GitHub's formula renderer)
- **Certificates:** [`cert_all_m.cpp`](round3/q3-birth-order/cert_all_m.cpp) (all $`m\ge200`$) · [`cert_small_m.cpp`](round3/q3-birth-order/cert_small_m.cpp) ($`4\le m\le199`$ and the founder seam) · [primitive audit](round3/q3-birth-order/audit_primitives.py) · [planted-defect controls](round3/q3-birth-order/planted_defects.py)
- **Reviews:** [mathematics](round3/q3-birth-order/REVIEW-Q3.md) and [implementation](round3/q3-birth-order/REVIEW-IMPLEMENTATION.md), both `HOLDS` on the exact bytes listed in the [provenance manifest](round3/PROVENANCE.md), with four recorded local clarifications.

Gemini 4 Neon wrote the proof and both certificates. GPT-6.2 Astra audited the argument and the certificate source bytes and wrote its own checks. That is one model family besides the author's; no human has refereed it. Round 3 also adds a short supplementary result on the mean number of non-leaves of degree $`k`$ when $`k\to\infty`$ and $`k/n\to0`$, expectation only ([argument](round3/q2-growing-degrees/CONTINUATION_DISCRETE_MEAN.md), [erratum](round3/q2-growing-degrees/ERRATUM_DISCRETE_MEAN.md), [review](round3/q2-growing-degrees/REVIEW.md)). The [Round 3 checkpoint](round3/README.md) records scoped progress and remaining gaps on all eight questions of their §7; outside Questions 2 and 3, that progress is not part of the paper.

## Results

- **The exponent (Theorem 1).** The core grows with the same exponent in expectation and almost surely.
- **A random multiplier (Theorem 4).** In an individual tree, $`C_n/n^\mu`$ converges almost surely to a finite, strictly positive, nonconstant random variable.
- **A density (Theorem 8).** The law of that random multiplier has a density.
- **Stable degree proportions (Theorems 4–5).** As $`n\to\infty`$, the fraction of non-leaves of each fixed degree $`k`$ converges almost surely to a deterministic $`p_k`$. As $`k\to\infty`$ these limits decay as $`p_k\sim a k^{-1-\mu}`$ for some $`a>0`$.
- **Hub mass (Theorem 9).** The sum of the fixed-label limits $`Z_v=\lim_{n\to\infty}D_v(n)/n`$ is one almost surely, and the next attachment’s conditional label law converges to these shares in total variation.
- **Birth order (Theorems 10–12).** The founders’ limiting shares have the same distribution. From label 2 onward, hub probabilities (Theorem 10), all positive moments (Theorem 11) and the tails at every threshold (Theorem 12) decrease strictly.

Krapivsky and Redner predicted the scaling of the normalized mean counts. Addario-Berry et al. state the fixed-degree almost-sure version in [*Random friend trees*, §7(2)](https://arxiv.org/html/2403.20185v1#S7). Theorems 4–5 give this almost-sure ratio for each fixed degree and then its large-degree tail. The [proof map](docs/proof-map.md) also covers the operator, drift, and conditional forest results.

## What is still open

- What is the exact value of $`\mu`$, and does it have a closed form?
- Is $`\mu`$ irrational?
- Do exact-degree counts follow the predicted scaling on individual trees, not just in mean, when the degree grows with the tree?
- Remaining parts of Questions 4–8 of Addario-Berry et al.: see the scoped status and open gaps in the [Round 3 checkpoint](round3/README.md#2-status-across-all-eight-questions-of-7).

Earlier versions: [Round 2 paper (PDF)](https://github.com/uzay-eth/friend-tree-exponent/blob/c250a1f2fdd53ced91b7b485e94ac8a16106b9c7/paper/friend-tree-exponent.pdf) · [Round 1 paper (PDF)](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/paper/friend-tree-exponent.pdf)

## Check it yourself

From the repository root, run the suites on x86-64 Linux with GCC/G++, OpenMP, GMP and CPython 3.12 with `venv`. The Round 3 check needs only `g++` with C++17 and Python 3.10 or newer. The Gold2 build also needs GNU `time`. Suite timings use one worker on a Ryzen 7 5800X, with peak RAM sampled across each process tree. [Round-one setup](verify/README.md) · [Round-two setup](round2/verify/README.md)

1. **Round-one bounds.** `bash verify/run_all.sh` checks both sides of the first enclosure. About 11 minutes, 2.19 GiB peak RAM. The measured run ended with `ALL 7 VERIFICATION SUITES PASSED (618s)`. The [three-vertex upper-witness replay](docs/proof-map.md#the-shorter-route-to-1730) already gives $`\mu\le1133/2000<17/30`$.
2. **Round-two checks.** `bash round2/verify/run_all.sh` replays 42 selected results. About 6 minutes, 605 MiB peak RAM. The measured run ended with `ALL 42 SELECTED RESULTS PASSED (381s)`.
3. **Theorem 2, exhaustive Gold2 replay.** On the Ryzen 7 5800X the two builds took 17 minutes in total and the strict verification took 2 minutes. Allow 10 GiB RAM to build, 20 GiB to verify and at least 12 GiB of free disk. Choose an output path that does not yet exist. [Gold rebuild guide](round2/verify/gold/README.md) · [Resource measurements](round2/verify/gold/RESULTS.md)

   ```sh
   bash round2/verify/gold/rebuild.sh "$PWD/gold-rebuild" gold2
   pair="$PWD/gold-rebuild/gold2-pair"
   mkdir "$pair"
   ln -s "$PWD/gold-rebuild/gold2-lower/quintet_lower.bin" "$PWD/gold-rebuild/gold2-upper/quintet_upper.bin" "$pair/"
   NEON2_THREADS=4 bash round2/verify/neon2/run.sh gold2-check "$pair"
   ```

   The builder prints `gold2-lower: match` and `gold2-upper: match` with both SHA-256 digests. The verifier prints the JSON receipt, both exact endpoints and `historical receipt matched byte for byte`.
4. **Round 3 checks.** `bash round3/verify_checkpoint.sh` checks the two certificate digests, compiles and runs both Q3 certificates with their primitive audits and planted-defect controls, then runs the exported checkers for Questions 2 and 5–8. About 2 minutes (1 min 38 s on the same Ryzen 7 5800X with g++ 13.3 and Python 3.12); its largest process used 301 MiB. It ends with `=== ALL ROUND 3 PACKAGED CHECKS PASSED ===`; the two certificates print `Boxes PASSED: 52, FAILED: 0` and `Boxes PASSED: 84, FAILED: 0`. [Round 3 verification notes](round3/README.md#5-reproducibility--verification)
5. **Printed numbers.** `python3 -B paper/checks/check_literals.py --enclosure` checks that the numbers printed in the paper match the certificate receipt. With Python 3.10 or newer and the standard library, this took 0.05 seconds and 21 MiB peak RAM. It ends with `PASS: Gold2 receipt matches Theorem 2 and the abstract.`

To rebuild the PDF and the arXiv source archive, see [Building the paper](paper/README.md).

## How it was done

Under Can Uzay’s direction, AI models developed the proofs and code, including the computer-assisted bounds used for the growth and degree laws. Two other model families reviewed the enclosure and replayed its original verifier, independent code checked selected parts, and the proofs await human refereeing. In Round 3 the stochastic-order proof had one outside model-family review, and Claude Opus 5.5 prepared this revision of the paper and pages without changing any proof, certificate or constant. [Contributions and review scope](round2/who-did-what.md) · [Round 3 provenance](round3/PROVENANCE.md) · [Computational dependencies](docs/proof-map.md#where-computation-enters)

## Prior work

- **Saramäki and Kaski** studied the family of random-walk attachment models that includes the undirected rule here ([2004 paper](https://arxiv.org/abs/cond-mat/0404088)).
- **Krapivsky and Redner** studied the rule as isotropic redirection and estimated $`\mu\approx0.566`$ ([2017 paper, §5.1](https://arxiv.org/abs/1706.01514)). Their slope range included both $`17/30`$ and the new enclosure, which rounds to their estimate.
- **Steven Strogatz** publicized the open problem in a [2022 post](https://x.com/stevenstrogatz/status/1572287647715934210) relaying a conversation with Redner: “it’s an open problem to prove that 17/30 is really the exponent.”
- **Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou and Lugosi** named the random friend tree and proved rigorous bounds on core growth ([Random friend trees](https://arxiv.org/abs/2403.20185)).
- **Hartle, Krapivsky, Redner and Zhang** solved two modified redirection models and reported a finite-size fit of $`0.569`$ for isotropic redirection ([Fig. 6](https://arxiv.org/html/2604.01540v1#S3.F6)).
- **[Christopher D. Long](https://github.com/long-mathematics/friend-tree-growth)** gave an earlier [generation-sum characterization and terminating two-sided approximation theorem](https://github.com/long-mathematics/friend-tree-growth/tree/090bbbe5faec83ec4124bd427eb5fda32c79bef1), and separately refined Uzay’s round-one witness to certify $`\mu\le0.5664485422<17/30`$. His [24 November draft](https://github.com/long-mathematics/friend-tree-growth/tree/5bd25a7d506cb64a8444b86334ab531d295f6c78) proves constant-factor mean bounds and gives a [signed hypergeometric product identity for the original tree](https://github.com/long-mathematics/friend-tree-growth/blob/5bd25a7d506cb64a8444b86334ab531d295f6c78/README.md#signed-hypergeometric-products-for-the-original-tree) in his transcendence programme.
- **Mohamed Sabba** reported the interval $`[0.565615,0.565809]`$ and a numerical estimate of $`0.5657237`$ in replies to Strogatz’s thread. [Interval](https://x.com/ma_sabba/status/2124768206141206762) · [Estimate](https://x.com/ma_sabba/status/2124791407667716208)
