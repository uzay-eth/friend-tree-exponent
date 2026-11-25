> **Public checkpoint export note:** Edited export of `walk-neon/ERRATUM-3.md` (frozen original SHA-256 `e77c8a93001aa9560a0edc8139f45fa750b9ced93305ef2bfa94ee92dd8dd0da`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.** Read alongside [`REVIEW-K9-NEON.md`](REVIEW-K9-NEON.md), which also records the separate `2s(2,1) = 221/250` printed-factor typo and unrestricted first-visit root-return proof still unaddressed in this erratum.

# Erratum 3 — Printed General Polynomial Display in Lemma 1(2) of Frozen `NINE-STEP.md`

Frozen `NINE-STEP.md` (SHA-256 `8d444e33dd93392276c316da805f2ef319b8b29c89b1ce076a694ef861594974`) is preserved verbatim. This erratum corrects a transcription error in the displayed general closed-form polynomial $`P_{2m+1}(a; j, s)`$ in Lemma 1(2) (line 61 of `NINE-STEP.md`), caught during cross-family review (`walk-astra/RESULT-K9-REVIEW.md`).

## 1. What was printed vs. the exact polynomial

In Lemma 1(2) of `NINE-STEP.md` (line 61), the sum $`\Psi_y = \sum_{p=0}^{m-1} \tilde{u}_p + j s \tilde{u}_m`$ of the exact decoupled recursion

```math
\tilde{u}_0 = a, \qquad \tilde{u}_p = (t - a) + a \tilde{u}_{p-1} = t - \sum_{r=1}^{p-1} \frac{1}{j} a^r - a^p + a^{p+1} \quad (p \ge 1,\; t = 1 - \tfrac{1}{j})
```

was mistakenly printed as degree $`m + 2`$ with leading term $`(m + j s) t`$:

```math
P_{2m+1}^{\mathrm{printed}}(a; j, s) = (m + j s) t - \sum_{r=1}^{m-1} \frac{m - r + j s}{j} a^r - s a^m + (1 - j s) a^{m+1} + j s a^{m+2}.
```

That display counted $`p = 0`$ as contributing $`t`$ instead of $`\tilde{u}_0 = a`$, shifting the powers of $`a`$ up by $`1`$.

The **exact** closed form of $`\sum_{p=0}^{m-1} \tilde{u}_p + j s \tilde{u}_m`$ has degree $`m + 1`$ and reads:

```math
P_{2m+1}(a; j, s) = (m - 1 + j s) t - \sum_{r=1}^{m-1} \frac{m - 1 - r + j s}{j} a^r + (1 - j s) a^m + j s a^{m+1}. \qquad (\text{P-exact})
```

## 2. Impact on Lemma 1(2), Proposition 1, and Theorems 1–3 of `NINE-STEP.md`

- **Descartes quasiconvexity (Lemma 1(2)) is completely unchanged.** Differentiating (P-exact) gives

```math
P_{2m+1}'(a; j, s) = -\sum_{r=1}^{m-1} r\,\frac{m - 1 - r + j s}{j} a^{r-1} + m(1 - j s) a^{m-1} + (m + 1) j s\, a^m.
```

Whenever $`0 \le j s \le 1`$, the coefficients of $`a^0, \dots, a^{m-2}`$ are strictly negative (or nonpositive when $`r = m - 1`$ and $`s = 0`$), while the top two coefficients (of $`a^{m-1}`$ and $`a^m`$) are nonnegative. By Descartes' rule of signs, $`P_{2m+1}'(a)`$ changes sign at most once on $`(0, \infty)`$ (from negative to positive), so $`P_{2m+1}(\cdot; j, s)`$ is quasiconvex on every $`[0, b_0] \subseteq [0, t]`$ and attains its maximum at an endpoint $`a \in \{0, b_0\}`$.
- **Endpoint formulas (line 62 of `NINE-STEP.md`) were already exact:** evaluating (P-exact) gives $`P_{2m+1}(0; j, s) = (m - 1) t + (j - 1) s`$ and $`P_{2m+1}(t; j, s) = \sum_{r=1}^m t^r + j t^{m+1} s`$, matching line 62 verbatim.
- **The pinned verifier [`verify_k9.py`](verify_k9.py) (commit `f2e72034`) was never affected:** it computes $`\Psi_y`$ directly from the exact recursion $`\tilde{u}_0, \dots, \tilde{u}_4`$, so all rational tables, Bernstein certificates, and finite-tree audits in Proposition 1 and Theorems 1–3 of `NINE-STEP.md` remain unchanged.

Exact symbolic equivalence between (P-exact) and the recursion $`\sum_{p=0}^{m-1} \tilde{u}_p + j s \tilde{u}_m`$ is verified in `verify_even.py` (author worktree script, not exported) (Part 3) for all $`m = 1, \dots, 8`$, $`j \in \{2, 3, 4, 7, 17\}`$, $`s \in \{0, \frac{1}{2j}, \frac{1}{j}\}`$.
