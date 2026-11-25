# arXiv submission sheet

**Status** — Draft metadata for the paper dated 25 November 2026. No arXiv upload has been performed, and nothing has been submitted anywhere. This sheet only records what a submission would contain.

**Title** — The non-leaf growth exponent of the random friend tree

**Authors** — Can Uzay

**Affiliation** — University of Nevada, Reno

**Abstract** — Paste the contents of the box as a single paragraph. It has 1916 ASCII characters, including spaces and TeX. [arXiv metadata instructions](https://info.arxiv.org/help/prep.html#abstract-required)

```text
The random friend tree starts from one edge. Each step attaches a leaf to a uniformly chosen neighbour of a uniformly chosen vertex. Let $C_n$ count the non-leaves (the core) at size $n$. We prove that $C_n$ grows with a common exponent $\mu$ in mean and almost surely, and that $C_n/n^\mu$ converges almost surely to a positive finite limit with a density. Exact computer-assisted certificates give $0.565714341933\le\mu\le0.565736986679$, rounded outwards, excluding $17/30$ and every rational with denominator below $426$, but not $241/426$. For fixed $k\ge2$ the fraction of non-leaves of degree $k$ converges a.s.\ to a deterministic $p_k$, and $k^{1+\mu}p_k$ tends to a positive limit as $k\to\infty$. Let $Z_v$ be the limiting degree share of vertex $v$. Then $\sum_v Z_v=1$ almost surely, the two founders' shares agree in law, and $\mathbb{P}(Z_u>x)>\mathbb{P}(Z_v>x)$ for all $2\le u<v$ and $0\le x<1$. Given any current tree, the earlier of the next two labels has the strictly larger conditional tail at those thresholds. A supplementary result gives a mean exact-degree equivalent as $k\to\infty$ and $k/n\to0$. Finally, $\mathbb E C_n$ satisfies no linear recurrence with polynomial coefficients, so its generating function is not D-finite; the proof is $p$-adic, for primes $p\equiv2\pmod3$, and says nothing about whether $\mu$ is rational. This settles the growth and limit (Question~1), the fixed-degree law (Question~2), the birth-order comparisons (Question~3) and the unit total share (Question~4) of Addario-Berry et al. The exact value and arithmetic nature of $\mu$ and a.s.\ exact-degree laws across the full intermediate range remain open. The growth and degree results use a computer-checked lower bound on $\mu$, the birth-order results exact finite checks. AI models checked the proofs within documented scopes, each main result by a family other than its author's. No human referee has.
```

| Field | Entry | Reason |
| --- | --- | --- |
| Primary | math.PR (Probability) | The main theorems identify the mean and almost-sure exponent of a random growing tree and the stochastic order of its limiting degree shares. |
| Cross-list | cond-mat.stat-mech (Statistical Mechanics) | The rigorous interval tests the growth exponent measured for isotropic redirection in the statistical-mechanics literature. |
| MSC 2020 | 60J80, 05C80 (primary); 05C05, 60C05, 05A15, 11B37 | Branching processes, random graphs and trees, combinatorial probability, exact enumeration and recurrences. |

**Comments**

```text
55 pages. Proofs, exact certificate verifiers and simulation data: https://github.com/uzay-eth/friend-tree-exponent
```

**Licence recommendation.** Choose Creative Commons Attribution 4.0 International (CC BY 4.0) to match the repository’s licence for its text and proofs. It permits commercial reuse and adaptations with attribution, and the grant is irrevocable, so check compatibility with the intended journal. [Repository licence](../LICENSE) · [arXiv licences](https://info.arxiv.org/help/license/index.html)

## AI-use statement

Can Uzay directed the work and is the sole listed author. AI models developed the mathematical proofs, wrote and ran the numerical checkers, drafted text, and reviewed claims across model families. No human referee has reviewed the proofs. The paper’s “Provenance and checks” section, the [contribution record](../round2/who-did-what.md) and the [stochastic-order provenance manifest](../round3/PROVENANCE.md) give the scope of these checks. Claude Opus 5.5 prepared this manuscript from the frozen proofs and reviews without changing their mathematics, certificates or constants. The proofs of Theorem 15 and of the refined lower endpoint in Theorem 2 are write-ups of frozen notes whose SHA-256 digests are listed in Appendix G.

**Source upload** — [Single-file source archive](friend-tree-exponent-arxiv.tar.gz). Select the pdfLaTeX engine. The bibliography is embedded in the source.
