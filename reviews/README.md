# Cold Reviews and Adversarial Attacks

Every mathematical claim included in the repository was reviewed cold by at least one model from a different model family than its author, and attacked by rival teams using independently written code, brute-force tree tests, and planted errors.

Two `GAP` verdicts occurred in the original reviews:
1. **Degree tail (`degree-tail.md`):** the first GPT-6.2 Sol cold review found a `GAP` in the summation-by-parts equality of part (c) at the endpoint $`k = 2`$ ($`m_0 = n`$) plus two overstated prose sentences. The proof was repaired and re-checked `CORRECT` (both the initial report and the re-check are preserved in [`review-degree-tail-gpt.md`](review-degree-tail-gpt.md)).
2. **Excluded tagged-neighbour route:** Claude Fable 5.2 returned `GAP` on an alternate upper-bound construction (`tag.cpp`) whose strict $`< 17/30`$ bound remained unverified ([`review-excluded-tagged-fable.md`](review-excluded-tagged-fable.md)). That construction is excluded from all results.

During review of the expanded existence proof, an exact counterexample settled a question about the mean within a planted branch. Its first derivative agrees with that of the free mean at the exhibited tree, but its second derivative is larger by $`113/2640`$. The proof uses a lower bound with an acceptance weight. [Counterexample, derivation and executable test](review-existence-counterexample-gpt.md) · [Weighted comparison](../proofs/existence.md)

*(The review files below are preserved as written during the cold audits, with only internal branch names, job-runner names, and local filesystem paths redacted where noted at the top of a file, and dead relative markdown links to the original working tree converted to plain code spans. Relative paths mentioned inside those historical transcripts refer to the original working tree; the packaged proofs and checkers live in [`proofs/`](../proofs/README.md) and [`verify/`](../verify/README.md).)*

## Index

| Target | Author (Family) | Reviewer / Attacker (Family) | Verdict | Report |
|---|---|---|---|---|
| [existence.md](../proofs/existence.md) | GPT-6.2 Astra (GPT) | Gemini 4 Neon (Gemini) | `CORRECT` | [review-existence-gemini.md](review-existence-gemini.md) |
| [almost-sure-exponent.md](../proofs/almost-sure-exponent.md) | GPT-6.2 Astra (GPT) | Gemini 4 Neon (Gemini) | `CORRECT` | [review-almost-sure-exponent-gemini.md](review-almost-sure-exponent-gemini.md) |
| [existence.md](../proofs/existence.md) & [almost-sure-exponent.md](../proofs/almost-sure-exponent.md) | GPT-6.2 Astra (GPT) | Claude Opus 5.5 (Claude) | `HOLDS` / `HOLDS` | [attack-existence-and-sample-claude.md](attack-existence-and-sample-claude.md) |
| [upper-bound-banded.md](../proofs/upper-bound-banded.md) | GPT-6.2 Astra (GPT) | Gemini 4 Neon (Gemini) | `CORRECT` + indep. `PASS` | [review-upper-bound-banded-gemini.md](review-upper-bound-banded-gemini.md) |
| [upper-bound-banded.md](../proofs/upper-bound-banded.md) | GPT-6.2 Astra (GPT) | Claude Opus 5.5 (Claude) | `HOLDS` | [attack-upper-bound-banded-claude.md](attack-upper-bound-banded-claude.md) |
| [bounds-transfer.md](../proofs/bounds-transfer.md) | GPT-6.2 Astra (GPT) | Gemini 4 Neon (Gemini) | `CORRECT` | [review-bounds-transfer-gemini.md](review-bounds-transfer-gemini.md) |
| [bounds-hub.md](../proofs/bounds-hub.md) | GPT-6.2 Astra (GPT) | Claude Fable 5.2 (Claude) | `CORRECT` | [review-bounds-hub-fable.md](review-bounds-hub-fable.md) |
| Excluded tagged-neighbour route | GPT-6.2 Astra (GPT) | Claude Fable 5.2 (Claude) | `GAP` (excluded) | [review-excluded-tagged-fable.md](review-excluded-tagged-fable.md) |
| [bracket-motifs.md](../proofs/bracket-motifs.md) | Gemini 4 Neon (Gemini) | GPT-6.2 Sol (GPT) | `CORRECT` | [review-bracket-motifs-gpt.md](review-bracket-motifs-gpt.md) |
| [bracket-motifs.md](../proofs/bracket-motifs.md) | Gemini 4 Neon (Gemini) | Claude Opus 5.5 (Claude) | `CORRECT` | [review-bracket-motifs-claude.md](review-bracket-motifs-claude.md) |
| [bracket-motifs.md](../proofs/bracket-motifs.md) | Gemini 4 Neon (Gemini) | Claude Opus 5.5 (Claude) | `HOLDS` / `HOLDS` | [attack-bracket-motifs-claude.md](attack-bracket-motifs-claude.md) |
| [almost-sure-lower.md](../proofs/almost-sure-lower.md) | Claude Fable 5.2 (Claude) | GPT-6.2 Sol (GPT) | `CORRECT` | [review-almost-sure-lower-gpt.md](review-almost-sure-lower-gpt.md) |
| [almost-sure-lower.md](../proofs/almost-sure-lower.md) | Claude Fable 5.2 (Claude) | Gemini 4 Neon (Gemini) | `CORRECT` | [review-almost-sure-lower-gemini.md](review-almost-sure-lower-gemini.md) |
| [cross-charging.md](../proofs/cross-charging.md) | Claude Fable 5.2 (Claude) | GPT-6.2 Sol (GPT) | `CORRECT` | [review-cross-charging-gpt.md](review-cross-charging-gpt.md) |
| [cross-charging.md](../proofs/cross-charging.md) | Claude Fable 5.2 (Claude) | Gemini 4 Neon (Gemini) | `CORRECT` | [review-cross-charging-gemini.md](review-cross-charging-gemini.md) |
| [degree-tail.md](../proofs/degree-tail.md) | Claude Fable 5.2 (Claude) | GPT-6.2 Sol (GPT) | v1 `GAP` -> repaired -> `CORRECT` | [review-degree-tail-gpt.md](review-degree-tail-gpt.md) |
| [degree-tail.md](../proofs/degree-tail.md) | Claude Fable 5.2 (Claude) | Gemini 4 Neon (Gemini) | `CORRECT` (noted $`k=2`$, $`\varepsilon<\mu`$) | [review-degree-tail-gemini.md](review-degree-tail-gemini.md) |
| [polynomial-scale-lower.md](../proofs/polynomial-scale-lower.md) | GPT-6.2 Sol / Carbon (GPT) | Claude Opus 5.5 (Claude) | `CORRECT` | [review-polynomial-scale-lower-claude.md](review-polynomial-scale-lower-claude.md) |
| Earlier baseline `[0.543, 0.586]` | Baseline | Gemini (Gemini) | `CORRECT` | [review-baseline-gemini.md](review-baseline-gemini.md) |
| Earlier baseline `[0.543, 0.586]` | Baseline | GLM (GLM) | `CORRECT` | [review-baseline-glm.md](review-baseline-glm.md) |
| Combined technical report ([REPORT.md](../REPORT.md)) | Gemini 4 Neon (Gemini) | Claude Opus 5.5 (Claude) | `CONSISTENT` / `APPROVE` | [review-combined-claude.md](review-combined-claude.md) |
| Combined technical report ([REPORT.md](../REPORT.md)) | Gemini 4 Neon (Gemini) | GPT-6.2 Sol (GPT) | `CONSISTENT` / `APPROVE` | [review-combined-gpt.md](review-combined-gpt.md) |
| Public repository packaging ([README.md](../README.md), [REPORT.md](../REPORT.md), [verify/run_all.sh](../verify/run_all.sh)) | Gemini 4 Neon (Gemini) | Claude Opus 5.5 (Claude) | `APPROVE` (after 2 revision rounds) | [review-public-packaging-claude.md](review-public-packaging-claude.md) |
| [Simulation sources, data and fits](../simulations/README.md) | Original continuous-time runs and new fixed-size experiment | Claude Opus 5.5 (Claude) | `CONSISTENT / APPROVE` | [review-simulations.md](review-simulations.md) |

## Reviews of the expanded proofs

| Text or test | Cold reviewer | Verdict and scope |
|---|---|---|
| [Motif drift and GMP correspondence](review-motif-writeup-gpt.md) | GPT-6.2 Astra | CORRECT for both shipped certificates, including the boundary at the first tail degree. |
| [Existence and sample paths](review-existence-and-sample-writeups-gpt.md) | GPT-6.2 Astra | CORRECT after revision, including the weighted marginal comparison, the preparation flags and the random first positive grid index. |
| [Unweighted branch comparison](review-existence-counterexample-gpt.md) | GPT-6.2 Astra | COUNTEREXAMPLE to a discarded comparison. The weighted lemma admits the exhibited tree. |
| [Banded upper bound](review-expanded-banded.md) | Claude Opus 5.5 | CORRECT for the history-dependent potential and its integer checker. |
| [Supporting results and exhaustive tree tests](review-supporting-results-audit.md) | Claude Opus 5.5 | CORRECT for the six supporting deductions and the direct tree and history tests reviewed. |
| [Degree tail and polynomial-scale typography](review-degree-typesetting.md) | Claude Opus 5.5 | CORRECT for both mathematical texts, including the empty sum and the fourteen numbered equations. |

## Paper

Claude Opus 5.5 checked the 24-page preprint against the proofs, certificate inputs, numerical data and cited literature, and rebuilt the PDF from source. [Cold review and verdict](review-preprint.md)

## Pre-submission referee round

The four referees read the version at `526f350` independently. Their assignments overlap. The first three each checked all three headline arguments, including the full motif drift lemma. The fourth also read the main proofs in support of the introduction and the verification section.

| Referee | Assigned close reading (original sections) | Report |
| --- | --- | --- |
| R1 · Claude Opus 5.5, max | §§1–5 | [Opus referee report](review-submission-opus.md) |
| R2 · GPT-6.2 Astra, ultra | §§1–5 | [Astra referee report](review-submission-astra.md) |
| R3 · Gemini 4 Neon, xhigh | §§2–4, 6–8, production and references | [Neon referee report](review-submission-neon.md) |
| R4 · Claude Fable 5.2, max | Abstract, §§1, 5–8, production and references, with §§2–4 and the supporting proofs | [Fable referee report](review-submission-fable.md) |

The README attributes the existence and degree-band proofs to GPT, the motif proof to Gemini, and several complementary bounds to Claude. The assignments therefore give the existence and degree-band arguments cold checks from Claude and Gemini, and the motif argument cold checks from Claude and GPT. All four reports also check the bibliography.

- [Independent editorial assessment](review-submission-editor-gpt.md)
- [Findings and rulings on the revision](review-submission-findings.md)
- [Fresh cross-family confirmation of the revised manuscript](review-submission-confirmation-claude.md)
- [Revised manuscript and source bundle](../paper/README.md)
