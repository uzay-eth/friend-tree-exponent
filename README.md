# The random friend tree exponent is not 17/30

```math
\frac{113}{200} = 0.5650 \;\le\; \mu \;\le\; \frac{1133}{2000} = 0.5665 \;<\; \frac{17}{30} = 0.5666\overline{6}.
```

Here $`C_n`$ counts the non-leaves in the tree on $`n`$ vertices. The same exponent $`\mu`$ governs $`\mathbb{E}[C_n] = n^{\mu + o(1)}`$ and $`\log C_n / \log n \to \mu`$ almost surely.

**Can Uzay**

Start from a single edge $`T_2`$; at each step pick a uniform vertex $`V_n`$, pick a uniform neighbour $`W_n`$ of $`V_n`$, and attach a new leaf to $`W_n`$. Addario-Berry et al. call this model the random friend tree. [Random friend trees](https://arxiv.org/abs/2403.20185)

In Krapivsky and Redner’s paper and in Strogatz’s thread, $`N_k`$ denotes the number of vertices of degree **exactly** $`k`$. Here $`N_k(n)`$ counts vertices of degree **at least** $`k`$, so $`N_2(n)=C_n`$. The exact-degree count is $`N_k(n)-N_{k+1}(n)\le C_n`$ for every $`k\ge2`$. Thus both finite-size mean upper bounds below apply to the expectation of each exact-degree count. The cumulative lower bounds give no matching lower bound for that difference, and its proposed pointwise scaling remains open. [Krapivsky and Redner, §5.1](https://arxiv.org/pdf/1706.01514#page=12) · [Strogatz’s thread](https://x.com/stevenstrogatz/status/1572287647715934210)

The certified gap below $`17/30`$ is at least $`1/6000`$. For every $`n \ge 3`$,

```math
\frac{7304352308}{8763936965}\prod_{j=3}^{n-1}\!\left(1 + \frac{113}{200j}\right) \;\le\; \mathbb{E}[C_n] \;\le\; \frac{14606753031}{2989767809}\prod_{j=3}^{n-1}\!\left(1 + \frac{1133}{2000j}\right).
```

An independent upper-bound construction also separates $`\mu`$ from $`17/30`$, giving $`\mathbb{E}[C_n] \le \frac{1539738114735}{3908512529}\prod_{j=3}^{n-1}\!\left(1 + \frac{5657}{9983j}\right)`$ for all $`n \ge 3`$ ($`\mu \le 5657/9983 = 17/30 - 1/299490`$). The exact value of $`\mu`$ remains open.

Krapivsky and Redner estimated $`\mu\approx0.566`$ for isotropic redirection. After speaking with Redner, Steven Strogatz asked in 2022 whether the exact exponent could be $`17/30`$. [Krapivsky and Redner, §5.1](https://arxiv.org/pdf/1706.01514#page=12) · [Strogatz’s post](https://x.com/stevenstrogatz/status/1572287647715934210)

Hartle, Krapivsky, Redner and Zhang obtained $`0.569`$ by nonlinear least squares on mean core sizes from $`N=10^2`$ to $`10^6`$. Their fitted slope exceeds the asymptotic upper bound $`0.5665`$, which constrains the limiting logarithmic growth rate. [Hartle et al., Fig. 6](https://arxiv.org/html/2604.01540v1#S3.F6)

The [technical report](REPORT.md), [proofs](proofs/README.md) and [review reports](reviews/README.md) give the arguments and checks.

## Check it yourself

Every finite check runs in exact integer or rational arithmetic (GMP, 128-bit integers with UBSan, or Python `fractions.Fraction`). On Ubuntu/Debian you need `build-essential`, `libgmp-dev`, `python3-venv`, `python3-pip`, and `gzip`.

```bash
git clone https://github.com/uzay-eth/friend-tree-exponent.git
cd friend-tree-exponent
bash verify/run_all.sh
```

The script creates a local virtual environment at `verify/.venv` from [`verify/requirements.txt`](verify/requirements.txt) (`numpy==2.1.3`, `mpmath==1.3.0`, `scipy==1.14.1`), builds and runs all seven verification suites (including negative controls), and ends with one `PASS` line per result:

```text
PASS  motif-bracket (Neon)          E[C_n] product bounds, exponents 113/200 (lower) and 1133/2000 (upper, < 17/30 by 1/6000); 17,329,108 GMP checks/side, 6/6 neg controls (21s)
PASS  banded-upper (Astra)          E[C_n] <= 1539738114735/3908512529 prod(1+5657/(9983j)) (< 17/30 by 1/299490); 226,545,617 rows (69s)
PASS  banded-upper-audit (indep.)   independent verifier 5657/9983; 226,545,617 rows, 436 histories, 50 adversarial trees, 14/14 neg controls (158s)
PASS  transfer-bounds (Astra)       [563/1000, 568/1000] (9,639,742 states) and hub [56/100, 57/100] (89s)
PASS  almost-sure-lower (Fable)     liminf log C_n/log n >= 543/1000 a.s.; 1,205,767 exact checks, (H1)-(H5), 1/1 neg control (101s)
PASS  cross-charging (Fable)        [11/20, 29/50] + 11/20 a.s.; 5,853/5,851 exact checks, 200 trees, 2/2 neg controls (2s)
PASS  polynomial-lower (Carbon/Sol) liminf C_n/n^(543/1000) > 0 a.s.; 319,086 upper checks, 102 trees, 84 branches, 15/15 neg controls (148s)
ALL 7 VERIFICATION SUITES PASSED (594s)
```

All seven suites passed in 594 s (9 min 54 s) on an AMD Ryzen 7 5800X and in 1,342 s (22 min 22 s) from a fresh clone on a separate 8-vCPU AMD EPYC 7543P VM. The Ryzen system had 8 cores, 16 threads and 62 GiB RAM and ran Ubuntu Linux (`6.8.0-149-generic`), GCC/G++ `13.3.0`, GMP `2:6.3.0`, Python `3.12.3` and GNU gzip `1.12`. See [per-suite details](verify/README.md).

## Who did what

I chose the question, set up the system that worked on it, and decided what counted as an acceptable answer: exact arithmetic only, every claim attacked by someone other than its author, and nothing stated more strongly than its proof. The mathematics was produced with AI models. Four models each led an independent attempt in a separate copy of the repository, on the same night, and each had to try to break the others' claims.

**Existence and almost-sure exponent (GPT-6.2 Astra).** Astra proved that the limit $`\mu = \lim \log \mathbb{E}[C_n]/\log n`$ exists ([`proofs/existence.md`](proofs/existence.md)) and that $`\log C_n/\log n \to \mu`$ almost surely ([`proofs/almost-sure-exponent.md`](proofs/almost-sure-exponent.md)). The argument embeds the chain in continuous time and shows that a positive fraction of core vertices can plant any prescribed finite tree with uniformly positive probability. Without this result the drift bounds below would be statements about the mean (and almost-sure $`\limsup`$) only.

**First upper bound below 17/30 (GPT-6.2 Astra).** Astra built a positive weight on every local state of a core vertex and verified, row by row in 128-bit integers, that $`n`$ times the expected one-step drift of the total weight is at most $`\lambda = 5657/9983`$ times the weight ([`proofs/upper-bound-banded.md`](proofs/upper-bound-banded.md); building on [`proofs/bounds-transfer.md`](proofs/bounds-transfer.md) and [`proofs/bounds-hub.md`](proofs/bounds-hub.md)). Two exact cancellations along core edges make this possible. The certificate has 226,545,617 finite rows plus two integer inequalities covering all degrees $`D \ge 2400`$. Combined with existence it gives $`\mu \le 5657/9983 = 17/30 - 1/299490`$. A second verifier ([`verify/banded-upper/independent/verify_banded.cpp`](verify/banded-upper/independent/verify_banded.cpp)), written from the specification before its author saw the first, reproduces it.

**The sharp bracket (Gemini 4 Neon).** Neon proved, for every $`n \ge 3`$, explicit two-sided product bounds on $`\mathbb{E}[C_n]`$ with exponents $`113/200 = 0.5650`$ below and $`1133/2000 = 0.5665`$ above, using directed 2- and 3-vertex local motifs and an analytic treatment of the $`1/D`$ tail ([`proofs/bracket-motifs.md`](proofs/bracket-motifs.md)). Each side is 17,329,108 exact GMP checks. Combined with Astra's existence theorem, this separates $`\mu`$ from $`17/30`$ by at least $`1/6000`$.

**Almost-sure lower bounds, cross-charging, and the degree tail (Claude Fable 5.2).** Fable proved a general transfer theorem showing that lower drift certificates of this kind bound $`C_n`$ from below almost surely and not only in mean ([`proofs/almost-sure-lower.md`](proofs/almost-sure-lower.md)). That gives $`\liminf \log C_n/\log n \ge 0.543`$ on the earlier baseline table and $`\ge 0.55`$ on a much smaller cross-charged certificate (5,666 table values; 5,851 lower and 5,853 upper exact inequalities; [`proofs/cross-charging.md`](proofs/cross-charging.md)), independently of Astra and Neon. Fable also proved the cumulative mean degree-tail law: combined with existence of $`\mu \in (0,1)`$, $`\lim_{n/k \to \infty} \log \mathbb{E}[N_k(n)]/\log(n/k) = \mu`$ uniformly over $`k \ge 2`$ ([`proofs/degree-tail.md`](proofs/degree-tail.md)).

**A polynomial-scale lower bound (Gemini 4 Carbon, proof by GPT-6.2 Sol).** Using both earlier baseline tables (`543/1000` lower and `293/500` upper), Carbon proved $`\liminf C_n/n^{543/1000} > 0`$ almost surely and almost-sure convergence of the normalized lower potential in $`(0,\infty]`$ with no $`-\varepsilon`$ loss in the exponent, via an $`L^2`$-bounded centered martingale on saturated branches ([`proofs/polynomial-scale-lower.md`](proofs/polynomial-scale-lower.md)). As a standalone bound on $`C_n`$, $`C_n/n^{0.543} \to \infty`$ a.s. already follows from Fable's $`0.55`$ bound; Carbon's distinct contribution is eliminating the exponent loss at a certificate's own exponent. (Both this result and Fable's $`0.543`$ application rest on the earlier `[0.543, 0.586]` baseline tables, whose exact verification code and cold reviews are in [`verify/baseline/`](verify/baseline/) and [`reviews/`](reviews/README.md); the earlier baseline write-up itself is not packaged here.)

**How the results were combined.** Neon's Theorem 1 bounds $`\mathbb{E}[C_n]`$ (and gives $`\limsup \log C_n/\log n \le 1133/2000`$ a.s.); Astra's existence theorem turns those bounds into bounds on the single exponent $`\mu`$, and Astra's sample-path theorem carries them to individual trees. Before combining, the definitions were matched line by line across all five works (same starting edge $`T_2`$, same indexing $`n`$, same core count $`C_n`$, mean versus almost-sure statements). The combined report ([`REPORT.md`](REPORT.md)) was written and all certificates were re-run from clean copies by a separate Gemini 4 Neon session, and the combined statements were checked cold by Claude Opus 5.5 and GPT-6.2 Sol.

**How the claims were tested.** Every proof was reviewed cold by at least one model from a different family than its author. Three adversarial attack suites were then run by Claude Opus 5.5 (working under the Fable lead): certificates were recomputed with independently written code, drift inequalities were tested by brute force (Neon's on 34,776 trees including all 32,972 rooted trees with at most 14 vertices; Astra's across 2,336 random histories and about 50 corner histories, covering 133,864 checked times and 1,532,992 vertex rows), and deliberate errors were planted to confirm the checks catch them. One gap was found in a claimed proof (in the first version of Fable's degree-tail theorem, at the endpoint $`k = 2`$ in part (c), plus two overstated sentences), then repaired and re-verified (`CORRECT`). An alternate tagged-neighbour upper-bound route by Astra also drew a `GAP` verdict and was excluded. Every review and attack report, with its verdict, is in [`reviews/`](reviews/README.md).

**Earlier work this builds on.** Krapivsky and Redner ([arXiv:1706.01514](https://arxiv.org/abs/1706.01514)) introduced isotropic redirection and measured $`\mu \approx 0.566`$. Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou, and Lugosi ([arXiv:2403.20185](https://arxiv.org/abs/2403.20185)) named the random friend tree, proved $`(\sqrt{13}-3)/6 \le \liminf \log C_n/\log n \le \limsup \log C_n/\log n \le \sqrt{2} - 1/2`$ almost surely (approximately $`[0.1009, 0.9142]`$), and left existence and the value of $`\mu`$ open. Hartle, Krapivsky, Redner, and Zhang ([arXiv:2604.01540](https://arxiv.org/abs/2604.01540)) solved two modified local redirection models analytically while noting that isotropic redirection remains open.

## What is still open

- **Exact value or closed form of $`\mu`$.** High-precision continuous-time simulations (heuristic, not proof) sit near $`0.565723 \pm 0.000003`$ inside our certified interval $`[0.5650, 0.5665]`$. Whether $`\mu`$ has an exact formula, operator characterization, or closed form is open.
- **Sharp asymptotics.** Neither $`\mathbb{E}[C_n] \asymp n^\mu`$ nor convergence of $`C_n / \mathbb{E}[C_n]`$ is proved; the existence proof allows subpolynomial factors $`n^{o(1)}`$.
- **Pointwise degree law and uniform almost-sure degree tail.** The cumulative mean limit for $`\mathbb{E}[N_k(n)]`$ is proved, and fixed-$`k`$ almost-sure convergence $`\log N_k(n)/\log(n/k) \to \mu`$ follows from Astra's sample-path theorem together with Theorem 4.8 of Addario-Berry et al. Both an almost-sure limit uniform over growing $`k`$ and the Krapivsky–Redner pointwise scaling $`\mathbb{E}[\#\{v : D_v(n) = k\}] \sim c\,n^\mu k^{-(1+\mu)}`$ remain open.

Preprint-level work, not peer reviewed. If you find an error, please open an issue.
