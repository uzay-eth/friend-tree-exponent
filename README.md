# The random friend tree exponent is not 17/30

**Can Uzay** · [Read the paper (PDF)](paper/friend-tree-exponent.pdf) · [LaTeX source](paper/friend-tree-exponent.tex) · [Proof and checker map](docs/proof-map.md)

The exponent $`\mu`$ lies in the exact interval $`\frac{35083855339718511}{62017003637500000}\le\mu\le\frac{127993535256830043}{226242120050000000}`$, or $`0.565713486333\le\mu\le0.565736986679`$ rounded outwards (Theorem 2).

Start with an edge. At each step, pick a vertex uniformly, pick one of its neighbours uniformly, and attach a new leaf to that neighbour. The non-leaves form the **core**, called the *nucleus* in the 2017 paper. It starts empty, and $`C_n`$ counts its vertices when the tree has size $`n`$.

## Results

- **The exponent (Theorem 1).** The core grows with the same exponent in expectation and almost surely.
- **A random multiplier (Theorem 4).** In an individual tree, $`C_n/n^\mu`$ converges almost surely to a finite, strictly positive, nonconstant random variable.
- **A density (Theorem 8).** The law of that random multiplier has a density.
- **Stable degree proportions (Theorems 4–5).** As $`n\to\infty`$, the fraction of non-leaves of each fixed degree $`k`$ converges almost surely to a deterministic $`p_k`$. As $`k\to\infty`$ these limits decay as $`p_k\sim a k^{-1-\mu}`$ for some $`a>0`$.
- **Hub mass (Theorem 9).** The sum of the fixed-label limits $`Z_v=\lim_{n\to\infty}D_v(n)/n`$ is one almost surely, and the next attachment’s conditional label law converges to these shares in total variation.
- **Birth order (Theorems 10–11).** The founders’ limiting shares have the same distribution. From label 2 onward, hub probabilities (Theorem 10) and all positive moments (Theorem 11) decrease strictly.

Krapivsky and Redner predicted the scaling of the normalized mean counts. Addario-Berry et al. state the fixed-degree almost-sure version in [*Random friend trees*, §7(2)](https://arxiv.org/html/2403.20185v1#S7). Theorems 4–5 give this almost-sure ratio for each fixed degree and then its large-degree tail. The [proof map](docs/proof-map.md) also covers the operator, drift, and conditional forest results.

## What is still open

- What is the exact value of $`\mu`$, and does it have a closed form?
- Is $`\mu`$ irrational?
- Are the limiting shares stochastically ordered by birth label at every threshold?

## Check it yourself

From the repository root, run the suites on x86-64 Linux with GCC/G++, OpenMP, GMP and CPython 3.12 with `venv`. The Gold2 build also needs GNU `time`. Suite timings use one worker on a Ryzen 7 5800X, with peak RAM sampled across each process tree. [Round-one setup](verify/README.md) · [Round-two setup](round2/verify/README.md)

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
4. **Printed numbers.** `python3 -B paper/checks/check_literals.py --enclosure` checks that the numbers printed in the paper match the certificate receipt. With Python 3.10 or newer and the standard library, this took 0.05 seconds and 21 MiB peak RAM. It ends with `PASS: Gold2 receipt matches Theorem 2 and the abstract.`

## How it was done

Under Can Uzay’s direction, AI models developed the proofs and code, including the computer-assisted bounds used for the growth and degree laws. Two other model families reviewed the enclosure and replayed its original verifier, independent code checked selected parts, and the proofs await human refereeing. [Contributions and review scope](round2/who-did-what.md) · [Computational dependencies](docs/proof-map.md#where-computation-enters)

## Prior work

- **Saramäki and Kaski** studied the family of random-walk attachment models that includes the undirected rule here ([2004 paper](https://arxiv.org/abs/cond-mat/0404088)).
- **Krapivsky and Redner** studied the rule as isotropic redirection and estimated $`\mu\approx0.566`$ ([2017 paper, §5.1](https://arxiv.org/abs/1706.01514)). Their slope range included both $`17/30`$ and the new enclosure, which rounds to their estimate.
- **Steven Strogatz** publicized the open problem in a [2022 post](https://x.com/stevenstrogatz/status/1572287647715934210) relaying a conversation with Redner: “it’s an open problem to prove that 17/30 is really the exponent.”
- **Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou and Lugosi** named the random friend tree and proved rigorous bounds on core growth ([Random friend trees](https://arxiv.org/abs/2403.20185)).
- **Hartle, Krapivsky, Redner and Zhang** solved two modified redirection models and reported a finite-size fit of $`0.569`$ for isotropic redirection ([Fig. 6](https://arxiv.org/html/2604.01540v1#S3.F6)).
- **[Christopher D. Long](https://github.com/long-mathematics/friend-tree-growth)** gave an earlier [generation-sum characterization and terminating two-sided approximation theorem](https://github.com/long-mathematics/friend-tree-growth/tree/090bbbe5faec83ec4124bd427eb5fda32c79bef1), and separately refined Uzay’s round-one witness to certify $`\mu\le0.5664485422<17/30`$. His [24 November draft](https://github.com/long-mathematics/friend-tree-growth/tree/5bd25a7d506cb64a8444b86334ab531d295f6c78) proves constant-factor mean bounds and introduces a signed hypergeometric product in his transcendence programme.
- **Mohamed Sabba** reported the interval $`[0.565615,0.565809]`$ and a numerical estimate of $`0.5657237`$ in replies to Strogatz’s thread. [Interval](https://x.com/ma_sabba/status/2124768206141206762) · [Estimate](https://x.com/ma_sabba/status/2124791407667716208)
