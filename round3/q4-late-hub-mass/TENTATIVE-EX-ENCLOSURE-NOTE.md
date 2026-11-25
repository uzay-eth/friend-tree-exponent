# Tentative Status Note — Q4 Numerical $`\mathbb{E}[X]`$ Enclosure (Full Draft Withheld)

> **Status: TENTATIVE / WITHHELD (neither proved in this checkpoint nor unexamined).** This note replaces the full draft `late-mass-neon/EX-ENCLOSURE.md` (original SHA-256 `f1ca28334b7a70578bea5ed0c1c4ab85dc4bb619953684bf13d2b4674ea069fb`), which is **withheld** from this public checkpoint because it contains a known lower-comparison defect at $`T_3`$ and the full cross-family repair/audit files referenced on the research board are not present on disk.

## What Was Proposed

In local research commit `08a7ee9314774eb4e1b03ab764bae4ce086662bb`, the `late-mass-neon` lane proposed a zero-floating-point (`__int128` + Python `fractions.Fraction`) renewal-resolvent certificate enclosing the limiting core-growth mean $`\mathbb{E}[X] = c_2`$ (under Round 1/2 Tier A and Tier B inputs) and induced bounds for $`\kappa = \mathbb{E}[X]\,c_{\mathrm{mass}}`$ (inheriting Tier C late-mass inputs):

- Proposed endpoints in frozen `EX-ENCLOSURE.md`: $`c_* \in [0.591238, 2.661580]`$, $`\mathbb{E}[X] = c_2 \in [0.330680, 1.488701]`$.
- The author’s end-to-end checker passed on the original worktree (including a concrete-tree `Fraction` oracle and planted-defect checks).

## What the Cross-Family Board Summary Reported — and What Is Missing

The rival lane’s final board summary (`late-mass-astra/RESULT.md`) reported a scoped **`HOLDS` with a repair**, specifically recording:

1. **Known $`T_3`$ defect in Proposition 4.5 of frozen `EX-ENCLOSURE.md`:** At the free three-vertex seed $`T_3`$, the printed lower residual in Proposition 4.5 exceeded its asserted right-hand side by $`38732423062131 / 10^{15}`$.
2. **Reported repair:** Substituting the actual direct-birth lower bound (using the certified virgin Perron coordinate) was reported to restore the signed inequality on every free tree with slack $`13804937869 / 10^{15}`$.
3. **Arithmetic documentation note:** The summary also noted that the frozen prose overstated universal overflow guards, while a signed-overflow build and explicit narrowing check reproduced identical output, and mentioned an author-checked tighter endgame ($`\mathbb{E}[X] \in [0.330713, 1.476087]`$).
4. **Missing archived audit/replay artifacts:** The standalone audit and replay files cited by `late-mass-astra/RESULT.md` (`REVIEW-EX-08A7EE9.md` and `EX-AUDIT-EVIDENCE.md`) are **absent** from the frozen board and worktree artifacts on disk (research seats were paused before those files were written or archived).

## Why the Full Candidate Is Withheld Here

Because the frozen `EX-ENCLOSURE.md` text carries an acknowledged Proposition 4.5 defect at $`T_3`$ and the complete repaired proof, standalone audit (`REVIEW-EX-08A7EE9.md`), and replay log (`EX-AUDIT-EVIDENCE.md`) cannot be backed by pinned artifacts on disk, publishing the defective full draft as usable would be misleading. Conversely, claiming the enclosure was never examined would ignore the recorded cross-family summary. It is therefore held back as a **tentative** item awaiting a complete, archived, independently checkable repaired proof and audit.
