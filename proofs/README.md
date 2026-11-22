# Proofs

Each file below is the proof as it stands at its pinned commit, with a one-line header identifying the model that produced it.

| Proof | Statement | Produced by | Finite check |
|---|---|---|---|
| [existence.md](existence.md) | Existence of the mean exponent $\mu = \lim \log \mathbb{E}[C_n]/\log n \in [0,1]$ | GPT-6.2 Astra | Analytic |
| [almost-sure-exponent.md](almost-sure-exponent.md) | Almost-sure convergence $\log C_n / \log n \to \mu$ | GPT-6.2 Astra | Analytic |
| [bracket-motifs.md](bracket-motifs.md) | Sharp finite-$n$ bracket $[113/200, 1133/2000]$ on $\mathbb{E}[C_n]$ (`< 17/30` by `1/6000`) and a.s. upper bound | Gemini 4 Neon | [`verify/motif-bracket/`](../verify/motif-bracket/) |
| [upper-bound-banded.md](upper-bound-banded.md) | Independent upper bound $5657/9983$ on $\mathbb{E}[C_n]$ (`< 17/30` by `1/299490`) | GPT-6.2 Astra | [`verify/banded-upper/`](../verify/banded-upper/) |
| [bounds-transfer.md](bounds-transfer.md) | Two-sided transfer bracket $[563/1000, 568/1000]$ on $\mathbb{E}[C_n]$ | GPT-6.2 Astra | [`verify/transfer-bounds/`](../verify/transfer-bounds/) |
| [bounds-hub.md](bounds-hub.md) | Earlier hub-state bracket $[56/100, 57/100]$ on $\mathbb{E}[C_n]$ | GPT-6.2 Astra | [`verify/transfer-bounds/`](../verify/transfer-bounds/) |
| [almost-sure-lower.md](almost-sure-lower.md) | General almost-sure lower transfer theorem (Theorem B) and $\liminf \log C_n / \log n \ge 0.543$ a.s. | Claude Fable 5.2 | [`verify/cross-charging/`](../verify/cross-charging/) |
| [cross-charging.md](cross-charging.md) | Compact cross-charged certificates $[11/20, 29/50]$ and $\liminf \log C_n / \log n \ge 0.55$ a.s. | Claude Fable 5.2 | [`verify/cross-charging/`](../verify/cross-charging/) |
| [degree-tail.md](degree-tail.md) | Cumulative mean degree tail $\mathbb{E}[N_k(n)]$ (Theorem D, repaired after v1 `GAP` at $k = 2$) | Claude Fable 5.2 | Analytic |
| [polynomial-scale-lower.md](polynomial-scale-lower.md) | Loss-free $\liminf C_n / n^{543/1000} > 0$ a.s. and a.s. convergence of $\Phi_n^- / n^{543/1000}$ in $(0,\infty]$ | GPT-6.2 Sol (under Gemini 4 Carbon) | [`verify/polynomial-lower/`](../verify/polynomial-lower/) |
