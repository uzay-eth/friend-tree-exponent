# The random friend tree’s core grows with exponent below 17/30

**Can Uzay** · [Read the paper (PDF)](paper/friend-tree-exponent.pdf) · [LaTeX source](paper/friend-tree-exponent.tex) · [Proof and checker map](docs/proof-map.md)

Start with an edge. At each step, pick a vertex uniformly, pick one of its neighbours uniformly, and attach a new leaf to that neighbour. Let $`C_n`$ count the non-leaves in the tree on $`n`$ vertices.

## Results

- **The exponent (Theorems 1–2).** The mean and almost-sure growth exponent $`\mu`$ lies in the exact interval below, which excludes $`17/30`$.

  ```math
  \frac{35083855339718511}{62017003637500000}
  \le \mu \le
  \frac{127993535256830043}{226242120050000000}
  ```

  Outward-rounded decimals: $`0.565713486333\le\mu\le0.565736986679`$.
- **A random multiplier (Theorem 4).** In an individual tree, $`C_n/n^\mu`$ converges almost surely to a finite, strictly positive, nonconstant random variable.
- **A density (Theorem 8).** The law of that random multiplier has a density.
- **Stable degree proportions (Theorems 4–5).** As $`n\to\infty`$, the fraction of non-leaves of each fixed degree $`k`$ converges almost surely to a deterministic $`p_k`$; as $`k\to\infty`$ these limits decay as $`p_k\sim a k^{-1-\mu}`$ for some $`a>0`$.
- **Hub mass (Theorem 9).** The sum of the fixed-label limits $`Z_v=\lim_{n\to\infty}D_v(n)/n`$ is one almost surely, and the next attachment’s conditional label law converges to these shares in total variation.
- **Birth order (Theorems 10–12).** The founders tie in law; from label 2 onward, the chance of a positive share and every positive moment strictly decrease with birth label.

Krapivsky and Redner predicted the degree law; the fixed-degree almost-sure statement answers the formulation in [*Random friend trees*, §7(2)](https://arxiv.org/html/2403.20185v1#S7), not a joint limit with degree growing with tree size. The [proof map](docs/proof-map.md) also covers the operator, drift, and conditional forest results.

## What is still open

- What is the exact value of $`\mu`$, and does it have a closed form?
- Is $`\mu`$ irrational?
- Are the limiting shares stochastically ordered by birth label at every threshold?

## Check it yourself

Download or clone the revision you are reading. From its top-level directory, run this command with Python 3.10 or newer (standard library only).

```sh
python3 -B paper/checks/check_literals.py --enclosure
```

Observed final line:

```text
PASS: Gold2 receipt matches Theorem 2 and the abstract.
```

This checks the receipt’s SHA-256 digest, recomputes the exact and outward-rounded endpoints, and compares them with the manuscript. It does **not** rerun the local inequalities. The [exhaustive build and replay](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/round2/verify/gold/README.md) uses the `gold2` preset: about 17 minutes to build with 10 GiB RAM and 2 minutes to verify with 20 GiB RAM on the recorded Ryzen 7 5800X, plus at least 12 GiB of scratch disk; link the snapshots into the pair directory rather than copying them. The [proof map](docs/proof-map.md#the-shorter-route-to-1730) also gives a smaller replay that already excludes $`17/30`$.

## How it was done

The proofs use continuous-time clocks, and the enclosure comes from bounding the expected change of weighted connected groups of up to five non-leaves. Can Uzay directed the work and set the evidence criteria; AI models wrote the proofs and code, and two other model families reviewed the enclosure deduction and replayed the same verifier, with independently written checks for selected parts but no human refereeing ([who did what](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/round2/who-did-what.md)).

## Prior work

- **Saramäki and Kaski** introduced the undirected attachment rule ([2004 paper](https://arxiv.org/abs/cond-mat/0404088)).
- **Krapivsky and Redner** studied it as isotropic redirection and estimated $`\mu\approx0.566`$ ([2017 paper, §5.1](https://arxiv.org/abs/1706.01514)); the new interval rounds to their estimate.
- **Steven Strogatz**, after hearing about that work from Redner, asked whether the exponent was exactly $`17/30`$ ([2022 post](https://x.com/stevenstrogatz/status/1572287647715934210)).
- **Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou and Lugosi** named the random friend tree and proved rigorous bounds on core growth ([Random friend trees](https://arxiv.org/abs/2403.20185)).
- **Hartle, Krapivsky, Redner and Zhang** solved two modified redirection models and reported a finite-size fit of $`0.569`$ for isotropic redirection ([Fig. 6](https://arxiv.org/html/2604.01540v1#S3.F6)).
- **Christopher D. Long** gave an [earlier generation-sum characterization, terminating two-sided approximation theorem, and certified upper bound of 0.5664485422](https://github.com/long-mathematics/friend-tree-growth/tree/090bbbe5faec83ec4124bd427eb5fda32c79bef1), [constant-factor mean bounds](https://github.com/long-mathematics/friend-tree-growth/tree/790491725d9cf445789e14c24c277532f5d2e5c7), and a [transcendence programme](https://github.com/long-mathematics/friend-tree-growth/tree/5bd25a7d506cb64a8444b86334ab531d295f6c78); his kernel retains the full birth environment of the inert-newborn process, unlike the marked-fringe kernel here.
- **Mohamed Sabba** reported $`[0.565615,0.565809]`$ and a separate numerical estimate of $`0.5657237`$ ([attribution record](https://github.com/uzay-eth/friend-tree-exponent/blob/a0972085ad72ad5685ae7b424988648d1965e73c/round2/who-did-what.md)); neither is an input to the enclosure above.
