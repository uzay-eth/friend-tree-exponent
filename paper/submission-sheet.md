# arXiv submission sheet

**Title** — The non-leaf growth exponent of the random friend tree

**Authors** — Can Uzay

**Affiliation** — University of Nevada, Reno

**Abstract** — Paste the contents of the box as a single paragraph. It has 743 ASCII characters, including spaces and TeX. [arXiv metadata instructions](https://info.arxiv.org/help/prep.html#abstract-required)

```text
In the random friend tree, each new vertex chooses an existing vertex uniformly and attaches to a uniformly chosen neighbour of that vertex. Let $C_n$ be the number of non-leaves at size $n$. We prove that $\log \mathbb{E}[C_n]/\log n$ and $\log C_n/\log n$ converge to the same deterministic exponent $\mu$, the latter almost surely. The proof plants finite seeds in a continuous-time embedding and controls the interaction across the edge joining each seed to the rest of the tree. Computer-assisted proofs with exact rational drift certificates give finite-size product bounds and imply $113/200 \le \mu \le 1133/2000 < 17/30$. An independent certificate based on degree bands gives the upper bound $5657/9983$, also strictly below $17/30$.
```

| Field | Entry | Reason |
| --- | --- | --- |
| Primary | math.PR (Probability) | The main theorem identifies the mean and almost-sure exponent of a random growing tree. |
| Cross-list | cond-mat.stat-mech (Statistical Mechanics) | The rigorous interval tests the growth exponent measured for isotropic redirection in the statistical-mechanics literature. |
| MSC 2020 | 60C05 (primary), 05C80, 05C05 | Combinatorial probability, random graphs and trees. |

**Comments**

```text
27 pages, 1 figure. Proofs, exact certificate verifiers and simulation data: https://github.com/uzay-eth/friend-tree-exponent
```

**Licence recommendation.** Choose Creative Commons Attribution 4.0 International (CC BY 4.0) to match the repository’s licence for its text and proofs. It permits commercial reuse and adaptations with attribution, and the grant for this version is irrevocable, so check compatibility with the intended journal. [Repository licence](../LICENSE) · [arXiv licences](https://info.arxiv.org/help/license/index.html)

## Authorship and endorsement

arXiv’s policy says “generative AI language tools should not be listed as an author; instead authors should refer to (1).” Item 1 requires disclosure of significant tool use, and item 2 requires signatories to “take full responsibility for all its contents”. Keep the disclosure covering proof development, checks, writing and model review. Can should submit only when he can take responsibility for the arguments and citations. [arXiv moderation policy, generative AI](https://info.arxiv.org/help/moderation/index.html#policy-for-authors-use-of-generative-ai-language-tools)

“arXiv requires that users be endorsed before submitting their first paper to arXiv or a new category.” For a first math.PR submission, automatic endorsement can come from a qualifying institutional account email **and** claimed authorship on an accepted arXiv Mathematics paper. The January policy update says “If you are a new submitter who does not meet BOTH of the requirements, you must follow the second path and seek personal endorsement.” With no prior arXiv Mathematics paper, the author needs personal endorsement from an eligible endorser. The submission system supplies the request link. This concerns account eligibility and leaves the public manuscript contact as supplied. [arXiv endorsement](https://info.arxiv.org/help/endorsement.html) · [Mathematics policy](https://blog.arxiv.org/2025/12/10/updated-endorsement-policy-for-arxiv-mathematics/) · [January 2026 update](https://blog.arxiv.org/2026/01/21/attention-authors-updated-endorsement-policy/)

**Source upload** — [Two-file source archive](friend-tree-exponent-arxiv.tar.gz). Select the pdfLaTeX engine. The bibliography is embedded in the source.
