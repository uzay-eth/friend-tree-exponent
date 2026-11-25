# arXiv submission sheet

**Status** — Draft metadata for the Round 3 version of 25 November 2026. No arXiv upload has been performed, and nothing has been submitted anywhere. This sheet only records what a submission would contain.

**Title** — The non-leaf growth exponent of the random friend tree

**Authors** — Can Uzay

**Affiliation** — University of Nevada, Reno

**Abstract** — Paste the contents of the box as a single paragraph. It has 1667 ASCII characters, including spaces and TeX. [arXiv metadata instructions](https://info.arxiv.org/help/prep.html#abstract-required)

```text
In the random friend tree the number $C_n$ of non-leaves has a growth exponent in the interval $0.565713486333\le\mu\le0.565736986679$. The interval rounds to Krapivsky and Redner's estimate of $0.566$ and excludes $17/30$. The finite drift inequalities for this enclosure were replayed with one exhaustive verifier. From one edge, $C_n/n^\mu$ converges almost surely to a finite, strictly positive random variable with a density. For each fixed $k\ge2$, the fraction of non-leaves of degree $k$ converges almost surely to a deterministic $p_k$. As $k\to\infty$, $p_k\sim a k^{-(1+\mu)}$ for some $a>0$. This proves the almost-sure form of Krapivsky and Redner's degree prediction recorded by Addario-Berry et al. in Section 7(2). For each birth label $v$, let $Z_v$ be the limit of its degree divided by tree size. Almost surely $\sum_v Z_v=1$. The limiting shares are strictly ordered by birth label in the first-order stochastic sense: $Z_1$ and $Z_2$ agree in law, and $\mathbb{P}(Z_u>x)>\mathbb{P}(Z_v>x)$ for all $2\le u<v$ and $0\le x<1$. Indeed, given any tree on $n\ge2$ vertices, before labels $n+1$ and $n+2$ are placed, the earlier label has the strictly larger conditional tail at every such threshold. These results answer Questions 1 to 4 in that section, Question 2 for fixed degrees. A supplementary result gives the mean number of non-leaves of degree $k$ when $k\to\infty$ and $k/n\to0$ jointly. The growth and fixed-degree limits use a computer-checked lower bound on $\mu$. The birth-order results also use exact finite checks; the threshold order reduces to two certificates in exact integer arithmetic. No human referee has reviewed the proofs.
```

| Field | Entry | Reason |
| --- | --- | --- |
| Primary | math.PR (Probability) | The main theorems identify the mean and almost-sure exponent of a random growing tree and the stochastic order of its limiting degree shares. |
| Cross-list | cond-mat.stat-mech (Statistical Mechanics) | The rigorous interval tests the growth exponent measured for isotropic redirection in the statistical-mechanics literature. |
| MSC 2020 | 60C05 (primary), 05C80, 05C05 | Combinatorial probability, random graphs and trees. |

**Comments**

```text
44 pages. Round 3 version. Proofs, exact certificate verifiers and simulation data: https://github.com/uzay-eth/friend-tree-exponent
```

**Licence recommendation.** Choose Creative Commons Attribution 4.0 International (CC BY 4.0) to match the repository’s licence for its text and proofs. It permits commercial reuse and adaptations with attribution, and the grant for this version is irrevocable, so check compatibility with the intended journal. [Repository licence](../LICENSE) · [arXiv licences](https://info.arxiv.org/help/license/index.html)

## AI-use statement

Can Uzay directed the work and is the sole listed author. AI models developed the mathematical proofs, wrote and ran the numerical checkers, drafted text, and reviewed claims across model families. No human referee has reviewed the proofs. The paper’s “How this work was done” section, the [contribution record](../round2/who-did-what.md) and the [Round 3 provenance manifest](../round3/PROVENANCE.md) give the scope of these checks. Claude Opus 5.5 prepared the Round 3 revision of the manuscript from the reviewed Round 3 documents without changing any proof, certificate or constant.

**Source upload** — [Single-file source archive](friend-tree-exponent-arxiv.tar.gz). Select the pdfLaTeX engine. The bibliography is embedded in the source.
