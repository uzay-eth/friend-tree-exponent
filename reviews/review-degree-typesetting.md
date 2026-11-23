# Cold review of the typeset degree-tail and polynomial-scale proofs

Reviewer: Claude Opus 5.5 (Claude family). [`degree-tail.md`](../proofs/degree-tail.md) was written by Claude Fable 5.2. [`polynomial-scale-lower.md`](../proofs/polynomial-scale-lower.md) was written by GPT-6.2 Sol under Gemini 4 Carbon. The typesetting was done by a GPT-family model.

Method: I compared both files word by word with the versions on `origin/main`. Every changed formula was checked against the old text, and every constant and endpoint was recomputed. I did not edit either proof.

Verdicts: **SAME** means the mathematics is unchanged and still correct. **CORRECT** means the step was rechecked from scratch. **NOTE** marks something that the typesetting did not cause.

## 1. `degree-tail.md`

| Step | Check | Verdict |
|---|---|---|
| Model and $`N_k(n)`$ | same definitions; $`N_2(n) = C_n`$ | SAME |
| Introduction | "no two-sided comparison up to constants" is now "on the logarithmic scale only", which is the same claim | SAME |
| Theorem D(a) | $`a_{m_0} + \frac{2n}{k}\sum_{m=m_0}^{n-1}(a_{m+1}-a_m)/(m+1)`$, for $`3 \le m_0 \le n`$ | SAME |
| D(b) | $`n \ge 8m'k`$, constant $`1/44160`$ | SAME |
| D(c) | $`n \ge 24k`$; bounds $`\frac{a}{44160}(n/(12k))^\delta`$ and $`A(2 + 1/(1-\lambda))(3n/k)^\lambda`$ | SAME |
| D(d) | limit as $`n/k \to \infty`$, uniform in $`k \ge 2`$ | SAME |
| Corollary | Gautschi gives $`\Gamma(m+0.55)/\Gamma(m) > m^{0.55}(m/(m+1))^{0.45} \ge m^{0.55}(3/4)^{0.45}`$ for $`m \ge 3`$, and $`0.5022 \cdot 0.75^{0.45} \approx 0.4412 \ge 0.44`$. The upper bound is $`3.974 \cdot (2 + 1/0.42) \cdot 3^{0.58} \approx 32.93 < 33`$. | CORRECT |
| Corollary, last paragraph | The old text said a density exponent of $`1 + \mu`$ is "not something proved here". The new sentence is conditional ("A pointwise power law … would give") but drops that explicit disclaimer. | NOTE (wording) |
| Conditional law of $`W_j`$ | $`\frac1j\sum_{u \sim w} 1/D_u(j)`$ | SAME |
| Lemma 1 | $`\mathbb E[D(j+1)] \le D(j)(j+1)/j`$; $`Q = D(D+1)`$ rises by $`2(D+1)`$; $`\prod_{j=m}^{n-1}(j+2)/j = n(n+1)/(m(m+1))`$ | CORRECT |
| Proof of (a) | Markov's inequality with $`D_v(m+1) = 2`$ gives $`2n/((m+1)k)`$ | CORRECT |
| Lemma 2 | drift $`\ge (L-1)/j`$; second moment at most $`(\ell+2)(\ell+3)(n/m)^2`$; Paley–Zygmund at half the mean; $`(\ell-1)^2/(4(\ell+2)(\ell+3))`$ increases in $`\ell`$ and equals $`4/120 = 1/30`$ at $`\ell = 3`$ | CORRECT |
| Eligibility and Lemma 3 | $`2I \le C-1`$ and $`n_3 \le n_1 - 2`$, so $`2E \ge (C+1)/2 + 2`$. The internal wording is replaced by a link to Lemma 2.1 of `existence.md`, which exists, and the anchor `#2-eligible-vertices` resolves. | CORRECT |
| Lemma 4, $`G`$ | $`\prod_{j=m'}^{8m'-1}(1 - 1/j) = (m'-1)/(8m'-1) \ge 2/23`$ | CORRECT |
| Lemma 4, Case A | $`1 - (m'-2)/(3m'-2) = 2m'/(3m'-2) \ge 2/3`$; $`1 - (3m'-2)/(8m'-2) = 5m'/(8m'-2) \ge 5/8`$; product $`(2/23)(2/3)(5/8) = 5/138`$ | CORRECT |
| Lemma 4, Case B | $`\frac{m'-2}{2m'-2}\big(1 - \frac{m'-3}{2m'-3}\big) = \frac{m'(m'-2)}{(2m'-2)(2m'-3)} \ge \frac14 \iff m' \ge 3`$; the two later blocks are each at least $`1/2`$; product $`(2/23)(1/4)(1/2)(1/4) = 1/368`$ | CORRECT |
| Proof of (b) | $`(\ell-1)n/(16m') \ge n/(8m') \ge k`$; $`4 \cdot 368 \cdot 30 = 44160`$ | CORRECT |
| Proof of (c), upper | $`m_0 = \lceil 2n/k \rceil`$ with $`48 \le 2n/k \le n`$; the sum is empty exactly when $`k = 2`$; summation by parts is correct; tail sum at most $`m_0^{\lambda-1}/(1-\lambda)`$; $`m_0 \le 3n/k`$ | CORRECT |
| Proof of (c), lower | $`m' = \lfloor n/(8k) \rfloor \ge 3`$, and $`n/(8k) - 1 \ge \frac23 \cdot \frac{n}{8k}`$ because $`n/(8k) \ge 3`$ | CORRECT |
| Proof of (d) | $`a_m \ge 1`$ gives constants for all $`m \ge 3`$ | SAME |
| Remarks 1–3 | The internal wording is replaced by a link to `almost-sure-exponent.md#2-one-block`, which resolves. The numbers are unchanged. | SAME |

## 2. `polynomial-scale-lower.md`

The only change is in 14 places: the equation tag on equations (1)–(14) became a plain label `\qquad (n)`. Each label is placed after the closing delimiter or at the end of its aligned row. The prose references to "(1)" through "(14)" still match their labels. No symbol, constant or inequality changed. **Verdict: SAME.**

## 3. Conventions

Neither file has an equation tag, an operator-name macro, `$$` or `\eqref`. Inline math uses dollar-backtick and display math uses fenced `math` blocks. `\begin{split}` appears in both files; it was already in `polynomial-scale-lower.md` on `origin/main`. Neither file has internal role wording or internal paths. All relative links resolve.

## 4. Points that the typesetting did not cause

These were already on `origin/main` and are not regressions.

1. `degree-tail.md` cites "Claim 2 of `RESULT.md`" twice, but no file of that name is in the release. Claim 2 is the cross-charging bound $`0.5022\,\Gamma(n+0.55)/\Gamma(n) \le \mathbb E[C_n] \le 3.974\,\Gamma(n+0.58)/\Gamma(n)`$, stated in [`REPORT.md`](../REPORT.md) and proved in [`cross-charging.md`](../proofs/cross-charging.md). The citation should point to one of those.
2. Remark 3 cites `outputs/hist_t27.txt`, which is not in the release. It is labelled as a numerical illustration only.
3. Optional: restore the old clause in the Corollary that a pointwise density exponent is not proved.

## Addendum: recheck of the notes

I reread `degree-tail.md` after the follow-up edits. Only these three passages changed, and each change addresses a note in §4.

1. The Corollary now cites the [cross-charging bounds](../proofs/cross-charging.md) instead of `RESULT.md`. The reference to `RESULT.md` in the last paragraph is also gone.
2. Remark 3 now links the [archived histogram](../simulations/outputs/hist_t27.txt), which is shipped. That file has the run's degree-histogram rows. I did not recompute the quoted fractions from it.
3. The Corollary now says that the theorem concerns cumulative counts and that the pointwise degree law remains open.

The final revision also changes two remarks.

4. Remark 1 no longer says that a block argument would give an almost-sure version. It now cites Theorem 4.8 of Addario-Berry et al. for fixed $`k \ge 2`$. By the summary in [`REPORT.md`](../REPORT.md), that theorem gives $`c_k < \liminf_n N_{k+1}(n)/N_k(n)`$ almost surely, and $`N_{k+1} \le N_k`$. Iterating from $`N_2 = C_n`$ gives $`N_k(n) \asymp C_n`$ almost surely. Together with the sample-path exponent of [`almost-sure-exponent.md`](../proofs/almost-sure-exponent.md), this gives $`\log N_k(n)/\log n \to \mu`$ almost surely for each fixed $`k`$. The remark leaves the uniform result for growing $`k`$ open, which agrees with `REPORT.md`. The deduction is correct. I checked Theorem 4.8 only against its statement in `REPORT.md`, not against the cited paper.
5. Remark 3 no longer reads the two finite-size slopes as approaching $`-\mu`$. It now calls them finite-size slopes from a single tree, which claims less and is accurate.

No other passage changed. No mathematical statement of Theorem D, its lemmas or its proofs changed. All links resolve.

## Verdict

**`degree-tail.md`: CORRECT.** The mathematics matches `origin/main`, every constant and endpoint checks out, the three notes in §4 are addressed, and the final remarks claim no more than is proved.

**`polynomial-scale-lower.md`: CORRECT.** Only the equation labels changed.

VERDICT: CORRECT
