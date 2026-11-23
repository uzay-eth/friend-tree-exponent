# Proofs

The header of each proof names the model that produced the original argument. The expanded proofs of the two exponent limits and the two headline drift inequalities include the couplings, generator identities and boundary cases.

For the motif bracket, each family of finite inequalities in the [proof](bracket-motifs.md) has a [row-by-row map to the GMP checker](bracket-motifs-code-map.md). The [banded upper bound](upper-bound-banded.md) includes its C++ correspondence. The [supporting-results audit](supporting-results-audit.md) gives the baseline states, the deductions behind the other six proofs, and the status of stronger claims. The [exhaustive tests](../verify/README.md#exhaustive-drift-checks-on-actual-trees) evaluate the conditional drift directly on trees and histories.

> **Note on preserved historical paths:** Some relative paths, commands, and filenames inside them (such as `../bounds/TRANSFER.md`, `BOUNDS.md`, `verified.log`, `transfer.log`, `baseline.log`, `RESULT.md`, `attempts/...`, or `./reproduce*.sh`) refer to the original working tree rather than this repository's packaged layout. Specifically:
> - `../bounds/TRANSFER.md` corresponds to [`bounds-transfer.md`](bounds-transfer.md) and `BOUNDS.md` to [`bounds-hub.md`](bounds-hub.md).
> - All finite checks are run via `bash verify/run_all.sh` (see [`verify/README.md`](../verify/README.md)).
> - References to "`REPORT.md` §§3.3–3.5" or "the earlier report" inside [`almost-sure-lower.md`](almost-sure-lower.md), [`cross-charging.md`](cross-charging.md), and [`polynomial-scale-lower.md`](polynomial-scale-lower.md) refer to the earlier `[0.543, 0.586]` baseline write-up, not to this repository's top-level [`REPORT.md`](../REPORT.md). The [baseline state and event definitions](supporting-results-audit.md#the-two-baseline-tables) are provided here. The certificate code is in [`verify/baseline/`](../verify/baseline/) and its two cold reviews are in [`reviews/review-baseline-gemini.md`](../reviews/review-baseline-gemini.md) and [`reviews/review-baseline-glm.md`](../reviews/review-baseline-glm.md).

| Proof | Statement | Produced by | Finite check |
|---|---|---|---|
| [existence.md](existence.md) | Existence of the mean exponent $`\mu = \lim \log \mathbb{E}[C_n]/\log n \in [0,1]`$ | GPT-6.2 Astra | Analytic |
| [almost-sure-exponent.md](almost-sure-exponent.md) | Almost-sure convergence $`\log C_n / \log n \to \mu`$ | GPT-6.2 Astra | Analytic |
| [bracket-motifs.md](bracket-motifs.md) | Sharp finite-$`n`$ bracket $`[113/200, 1133/2000]`$ on $`\mathbb{E}[C_n]`$ (`< 17/30` by `1/6000`) and a.s. upper bound | Gemini 4 Neon | [`verify/motif-bracket/`](../verify/motif-bracket/) |
| [upper-bound-banded.md](upper-bound-banded.md) | Independent upper bound $`5657/9983`$ on $`\mathbb{E}[C_n]`$ (`< 17/30` by `1/299490`) | GPT-6.2 Astra | [`verify/banded-upper/`](../verify/banded-upper/) |
| [bounds-transfer.md](bounds-transfer.md) | Two-sided transfer bracket $`[563/1000, 568/1000]`$ on $`\mathbb{E}[C_n]`$ | GPT-6.2 Astra | [`verify/transfer-bounds/`](../verify/transfer-bounds/) |
| [bounds-hub.md](bounds-hub.md) | Earlier hub-state bracket $`[56/100, 57/100]`$ on $`\mathbb{E}[C_n]`$ | GPT-6.2 Astra | [`verify/transfer-bounds/`](../verify/transfer-bounds/) |
| [almost-sure-lower.md](almost-sure-lower.md) | General almost-sure lower transfer theorem (Theorem B) and $`\liminf \log C_n / \log n \ge 0.543`$ a.s. | Claude Fable 5.2 | [`verify/cross-charging/`](../verify/cross-charging/) |
| [cross-charging.md](cross-charging.md) | Compact cross-charged certificates $`[11/20, 29/50]`$ and $`\liminf \log C_n / \log n \ge 0.55`$ a.s. | Claude Fable 5.2 | [`verify/cross-charging/`](../verify/cross-charging/) |
| [degree-tail.md](degree-tail.md) | Cumulative mean degree tail $`\mathbb{E}[N_k(n)]`$ (Theorem D, repaired after v1 `GAP` at $`k = 2`$) | Claude Fable 5.2 | Analytic |
| [polynomial-scale-lower.md](polynomial-scale-lower.md) | Loss-free $`\liminf C_n / n^{543/1000} > 0`$ a.s. and a.s. convergence of $`\Phi_n^- / n^{543/1000}`$ in $`(0,\infty]`$ | GPT-6.2 Sol (under Gemini 4 Carbon) | [`verify/polynomial-lower/`](../verify/polynomial-lower/) |
