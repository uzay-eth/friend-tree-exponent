# Submission metadata (draft)

**Status** — Draft metadata for the paper dated 25 November 2026 and its supplement. Nothing has been uploaded or submitted anywhere. This sheet records what a submission would contain. The source bundle has **not** been tested as an arXiv upload, and arXiv processes one main document per submission. Distributing the supplement there would need a different layout, for example a separate submission or ancillary file.

**Title** — The non-leaf growth exponent of the random friend tree

**Authors** — Can Uzay

**Affiliation** — University of Nevada, Reno

**Abstract** — The abstract of the paper, as one paragraph. It has 1518 ASCII characters, including spaces and TeX. [arXiv metadata instructions](https://info.arxiv.org/help/prep.html#abstract-required)

```text
The random friend tree starts from one edge. At each step a uniformly chosen vertex picks a uniformly chosen neighbour, and a new leaf is attached to that neighbour. A leaf becomes a non-leaf only when its one neighbour is chosen and passes it the new leaf, so the number $C_n$ of non-leaves at size $n$ depends on how the leaves are spread over the non-leaves, not only on how many non-leaves there are. We prove that $C_n$ grows with a common exponent $\mu$ in mean and almost surely, and that $C_n/n^\mu$ converges almost surely to a positive finite random variable whose law has a density. Exact computer-assisted certificates give $0.565714341933\le\mu\le0.565736986679$, rounded outwards. They exclude $17/30$ and every rational with denominator below $426$, but not $241/426$; the exact value of $\mu$ and whether it is rational remain open. For each fixed $k\ge2$ the fraction of non-leaves of degree $k$ converges almost surely to a deterministic limit with a power-law tail, and the limiting degree shares of the vertices sum to one. The mean $\E C_n$ satisfies no linear recurrence with polynomial coefficients. An appendix proves the last result; a companion supplement gives complete proofs of the others and a module on birth order. There the two founders' limiting shares agree in law, and from label two onward an earlier label has a strictly larger tail at every threshold $0\le x<1$. The proofs were developed and cross-checked with substantial help from AI models; no human referee has checked them.
```

| Field | Entry | Reason |
| --- | --- | --- |
| Primary | math.PR (Probability) | The main theorems identify the mean and almost-sure exponent of a random growing tree and the stochastic order of its limiting degree shares. |
| Cross-list | cond-mat.stat-mech (Statistical Mechanics) | The rigorous interval tests the growth exponent measured for isotropic redirection in the statistical-mechanics literature. |
| MSC 2020 | 60J80, 05C80 (primary); 05C05, 60C05, 05A15, 11B37 | Branching processes, random graphs and trees, combinatorial probability, exact enumeration and recurrences. |

**Comments**

```text
18 pages, with a 80-page supplement of complete proofs. Proofs, exact certificate verifiers and simulation data: https://github.com/uzay-eth/friend-tree-exponent
```

**Licence recommendation.** Choose Creative Commons Attribution 4.0 International (CC BY 4.0) to match the repository’s licence for its text and proofs. It permits commercial reuse and adaptations with attribution, and the grant is irrevocable, so check compatibility with the intended journal. [Repository licence](../LICENSE) · [arXiv licences](https://info.arxiv.org/help/license/index.html)

## AI-use statement

Can Uzay directed the work and is the sole listed author. AI models developed the mathematical proofs, wrote and ran the numerical checkers, drafted text, and reviewed claims across model families. No human referee has reviewed the proofs. The paper’s “Proofs, checks and provenance” section, supplement §PROV, the [contribution record](../round2/who-did-what.md) and the [stochastic-order provenance manifest](../round3/PROVENANCE.md) give the scope of these checks. Claude Opus 5.5 prepared the paper and supplement from the frozen proofs and reviews without changing their mathematics, certificates or constants. The proofs of Theorem PREC and of the refined lower endpoint of Theorem ENC are write-ups of frozen notes whose SHA-256 digests are listed in supplement §VER, and the supplement reproduces three further reviewed Markdown proofs verbatim.

**Source** — [Source bundle](friend-tree-exponent-source.tar.gz) of both TeX files and their build script; pdfLaTeX, bibliographies embedded. See [Building the paper](README.md).
