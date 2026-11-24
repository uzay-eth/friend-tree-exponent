# The random friend tree’s core grows with exponent below 17/30

**Can Uzay** · [Read the paper (PDF)](paper/friend-tree-exponent.pdf) · [LaTeX source](paper/friend-tree-exponent.tex) · [Proof and checker map](docs/proof-map.md)

Start with an edge. At each step, pick a vertex uniformly, pick one of its neighbours uniformly, and attach a new leaf to that neighbour. The non-leaves form the **core**, called the *nucleus* in the 2017 paper. It starts empty, and $`C_n`$ counts its vertices when the tree has size $`n`$.

## Results

- **The exponent (Theorems 1–2).** The mean and almost-sure growth exponent $`\mu`$ lies in the computer-certified interval below, which excludes $`17/30`$.

  ```math
  \frac{35083855339718511}{62017003637500000}
  \le \mu \le
  \frac{127993535256830043}{226242120050000000}
  ```

  Outward-rounded decimals: $`0.565713486333\le\mu\le0.565736986679`$.
- **A random multiplier (Theorem 4).** In an individual tree, $`C_n/n^\mu`$ converges almost surely to a finite, strictly positive, nonconstant random variable.
- **A density (Theorem 8).** The law of that random multiplier has a density.
- **Stable degree proportions (Theorems 4–5).** As $`n\to\infty`$, the fraction of non-leaves of each fixed degree $`k`$ converges almost surely to a deterministic $`p_k`$. As $`k\to\infty`$ these limits decay as $`p_k\sim a k^{-1-\mu}`$ for some $`a>0`$.
- **Hub mass (Theorem 9).** The sum of the fixed-label limits $`Z_v=\lim_{n\to\infty}D_v(n)/n`$ is one almost surely, and the next attachment’s conditional label law converges to these shares in total variation.
- **Birth order (Theorems 10–12).** The founders’ limiting shares have the same distribution. From label 2 onward, hub probabilities (Theorem 10) and all positive moments (Theorem 11) decrease strictly. Theorem 12 separately verifies hub-probability order through label 8.

Krapivsky and Redner predicted the scaling of the normalized mean counts. Addario-Berry et al. state the fixed-degree almost-sure version in [*Random friend trees*, §7(2)](https://arxiv.org/html/2403.20185v1#S7). Theorems 4–5 give this almost-sure ratio for each fixed degree and then its large-degree tail. The [proof map](docs/proof-map.md) also covers the operator, drift, and conditional forest results.

## What is still open

- What is the exact value of $`\mu`$, and does it have a closed form?
- Is $`\mu`$ irrational?
- Are the limiting shares stochastically ordered by birth label at every threshold?

## Check it yourself

The [three-vertex upper-witness replay](docs/proof-map.md#the-shorter-route-to-1730) already gives $`\mu\le1133/2000<17/30`$. For the sharper interval, the command below checks the printed numbers, not the drift inequalities.

Download or clone the revision you are reading. From its top-level directory, run this command with Python 3.10 or newer (standard library only).

```sh
python3 -B paper/checks/check_literals.py --enclosure
```

The command prints this final line.

```text
PASS: Gold2 receipt matches Theorem 2 and the abstract.
```

It authenticates the stored receipt and recomputes both exact endpoints and their outward rounding. The [exhaustive build and replay](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/round2/verify/gold/README.md) uses the `gold2` preset: about 17 minutes to build with 10 GiB RAM and 2 minutes to verify with 20 GiB RAM on the recorded Ryzen 7 5800X, plus at least 12 GiB of scratch disk. The replay uses links to the snapshots in its pair directory.

## How it was done

Under Can Uzay’s direction, AI models developed the proofs and code, including the computer-assisted bounds used for the growth and degree laws. Two other model families reviewed the enclosure and replayed its original verifier, independent code checked selected parts, and the proofs await human refereeing. [Contributions and review scope](round2/who-did-what.md) · [Computational dependencies](docs/proof-map.md#where-computation-enters)

## Prior work

- **Saramäki and Kaski** studied the family of random-walk attachment models that includes the undirected rule here ([2004 paper](https://arxiv.org/abs/cond-mat/0404088)).
- **Krapivsky and Redner** studied the rule as isotropic redirection and estimated $`\mu\approx0.566`$ ([2017 paper, §5.1](https://arxiv.org/abs/1706.01514)). Their slope range included both $`17/30`$ and the new enclosure, which rounds to their estimate.
- **Steven Strogatz** publicized the open problem in a [2022 post](https://x.com/stevenstrogatz/status/1572287647715934210) relaying a conversation with Redner: “it’s an open problem to prove that 17/30 is really the exponent.”
- **Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou and Lugosi** named the random friend tree and proved rigorous bounds on core growth ([Random friend trees](https://arxiv.org/abs/2403.20185)).
- **Hartle, Krapivsky, Redner and Zhang** solved two modified redirection models and reported a finite-size fit of $`0.569`$ for isotropic redirection ([Fig. 6](https://arxiv.org/html/2604.01540v1#S3.F6)).
- **Christopher D. Long** gave an earlier [generation-sum characterization and terminating two-sided approximation theorem](https://github.com/long-mathematics/friend-tree-growth/tree/090bbbe5faec83ec4124bd427eb5fda32c79bef1), and separately refined Uzay’s round-one witness to certify $`\mu\le0.5664485422<17/30`$. His [24 November draft](https://github.com/long-mathematics/friend-tree-growth/tree/5bd25a7d506cb64a8444b86334ab531d295f6c78) gives constant-factor mean bounds and develops his transcendence programme.
- **Mohamed Sabba** reported $`[0.565615,0.565809]`$ and a separate numerical estimate of $`0.5657237`$ (as recorded in the [attribution record](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/round2/who-did-what.md), which does not link a primary result).
