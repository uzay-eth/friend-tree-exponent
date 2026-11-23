# Editorial report on the random friend tree manuscript

**Recommendation. Revise before external review.** The paper establishes a deterministic logarithmic growth exponent for the non-leaf population, with matching limits in mean and almost surely. That is the main contribution for a probability journal. The exact drift certificates then locate the exponent in a narrow interval and settle the comparison with $`17/30`$. The paper should put the probabilistic argument on the shortest route from the model to the theorem and keep the full computer-assisted proofs in the same submission.

## Length and emphasis

The supplied PDF is 24 pages. The common-exponent proof starts on page 2, the motif proof on page 7, the degree-band proof on page 14, and the verification section on page 20. Aim for **12–14 pages of main text and 7–9 pages of proof appendices**, plus references. The main gain is a shorter first read. The technical appendices should remain in the PDF so a referee can check how each finite certificate implies an inequality for an arbitrary tree or growth history.

Put the model, the common-exponent theorem, and the sharp interval on page 1. State the two motif product bounds together, then state the separate product bound for every $`n\ge3`$ as a proposition. The degree-band estimate has exponent $`5657/9983\approx0.5666633`$ and prefactor about $`393.94`$. The motif upper estimate has exponent $`1133/2000=0.5665`$ and prefactor about $`4.89`$. The larger exponent and prefactor make the banded bound weaker at every size. Its value is a second construction that independently puts the upper exponent below $`17/30`$.

Reduce the literature review to a short paragraph after the result. Identify the previously known almost-sure bounds, the numerical estimate near $`0.566`$, and the question about $`17/30`$. Distinguish the historical discussion of exact degree counts from the paper’s theorems about the core and cumulative counts. Use “logarithmic growth exponent” early so readers know which asymptotic the theorem establishes.

## Proposed order and proof allocation

1. **Model and results, about two pages.** Give the attachment rule and initial edge, define $`C_n`$, and state the common exponent and the two-sided product bounds with their exact constants. Display the sharp bracket and its separation from $`17/30`$. State the degree-band proposition and give a brief map of the probabilistic and certificate arguments.

2. **One continuous-time proof of the common exponent, about four to five pages.** Keep the clock embedding and the independence of population waiting times from the jump chain. Use the disjoint-site count of at least $`C/6`$ at the outset and state a single conditional block lemma. For fixed error parameters and a sufficiently long comparison time $`s`$, the block endpoint dominates a sum of independent copies of a contribution $`Y_s`$ from the selected sites, with

```math
\mathbb E Y_s \ge c\,e^{-\eta s}\bigl(F(s)-k-O(e^{-s})\bigr),
\qquad F(s)=\mathbb E C(s).
```

Here $`F`$ is the mean from a single edge, $`k`$ is the fixed seed size, $`c>0`$ is independent of $`s`$, and $`\eta`$ can be made arbitrarily small. Prove the lemma in three steps. Pack the sites and plant a fixed seed with positive probability. Couple each planted branch to a free reference across its cut edge. Restrict to good reference histories and select a root that carries a positive share of the good expected core count. Keep the conditional product of the root-choice acceptance factors inside the reference expectation. Spell out how fresh clocks and thinning give independent site inputs even though the branches share an exterior. A small diagram of the cut edge could replace some of the description.

   Take expectations of the block lemma to obtain the mean iteration, which can replace the separate $`C/4`$ count and marginal-site argument. For the pathwise lower bound, use the same lemma with conditional Chebyshev. Each $`Y_s`$ is bounded by its reference population and has finite second moment for fixed $`s`$, so the truncation and its constants can be dropped. Keep the proof that the core tends to infinity and the martingale argument that controls the fraction of unsuccessful blocks. Finish with Markov and Borel–Cantelli for the pathwise upper bound. Keep both transfers to fixed size in one subsection. The mean uses independence and Yule tails, and the pathwise limit uses the asymptotics of the hitting times. State the order of limits in the proof. Fix the seed and preparation parameters before taking $`s`$ large, then send the loss $`\eta`$ to zero.

3. **Common drift argument and motif bounds, about two to three pages.** Define the scaled generator and the two possible updates once. State the elementary product lemma that iterates a drift bound and a one-sided comparison with $`C_n`$. Specify the comparison in each direction and the initial value at the three-vertex path. Describe the motif potential as a sum over vertices, edges, directed two-edge paths, and sibling pairs. Show where the newborn terms enter, and explain how the redistribution of rates and the corner bounds give the two global drift inequalities. State the certified values used in the product lemma. Put the full generator decomposition, projection rules, all row formulas, the signed corner bounds, and the per-core comparison in Appendix A. Keep the special quadratic bound for a tail inverse degree multiplying a jump that depends on that same endpoint.

4. **Degree-band proof, about one to two pages.** Explain why the weight depends on history, what the recorded labels guarantee, and how the two transfers cancel in the conditional drift. Describe how trial rewards remove the newborn dependency during reverse-order construction, and how the final comparison with the constructed newborn weights closes the argument. Give the finite-row inequality, the two affine tail endpoint tests, and the resulting product constant. Put the full history invariant and its proof, the band table, the row specification, and the recursion and verification details in Appendix B. Use one short dyadic argument to record the direct almost-sure upper bound from either product.

5. **Consequences and discussion, about one page.** Keep the mean exponent for the cumulative degree count as a named corollary, uniformly as $`n/k\to\infty`$ with logarithmic denominator $`\log(n/k)`$. Give its proof in the supplement, with a short explanation in the text that the upper estimate comes from the times of promotion and the lower estimate from growth of eligible old vertices. State the almost-sure consequence for each fixed $`k`$. Move the preliminary pathwise exponents and the polynomial-scale bound to the supplement. Keep the open questions on the exact exponent, normalization, and exact degree counts. Put the simulation protocol, fit table, comparisons, and figure in the supplement. A brief main-text sentence about finite-size slopes can retain the link to the numerical literature. Describe sampling errors at the measured sizes and disclose the rounded Gaussian draws in the simulation methods.

Give the two certificate lemmas a common reproducibility paragraph. Link the immutable data and checker inputs, identify the precise commands and expected summaries, and report the arithmetic conventions and resource requirements. In the appendices, present the universal reduction, the exhaustive arithmetic checks, and their reported outcomes next to the respective lemma. The direct small-tree tests and deliberate perturbations are useful checks of the encoding and implementation. The finite-row proofs establish the universal claims through the reduction. Fix the undefined endpoint-measure symbol in the motif table description by using $`\zeta_s`$, the measure defined there.

## Title and arXiv abstract

**The non-leaf growth exponent of the random friend tree**

**Original abstract.** It accurately describes the attachment rule, the common deterministic exponent in mean and almost surely, the exact-arithmetic bounds, and the finite-size products. It omits the comparison with $`17/30`$, which both upper certificates settle, and leaves the interaction across the cut edge implicit in its description of the proof. The closing reference to simulations uses space that could explain either point.

The original TeX body has **690 characters**. Its only manuscript-defined macro is `\E`. For the arXiv field, replace it with standard `\mathbb{E}` and keep the inline TeX. The resulting one-paragraph ASCII text has **698 characters**, including spaces and TeX source and excluding the abstract environment and trailing newline. It leaves **1,222 characters** under the 1,920-character limit.

**Proposed abstract.**

Use the abstract for the common exponent, the continuous-time idea, and the two certificate constructions. The exact product prefactors, degree-tail corollary, and numerical fits can stay in the paper. The following is one paste-ready paragraph in ASCII, with inline TeX. Its length is **702 characters**, counting spaces, punctuation, and TeX source and excluding the surrounding fence and trailing newline. It leaves **1,218 characters** under the 1,920-character limit.

```text
In the random friend tree, each new vertex chooses a vertex uniformly and attaches to a uniformly chosen neighbor of that vertex. Let $C_n$ be the number of non-leaves at size $n$. We prove that $\log \mathbb{E}[C_n]/\log n$ and $\log C_n/\log n$ converge to the same deterministic exponent $\mu$, the latter almost surely. The proof plants finite seeds in a continuous-time embedding and controls the interaction across the edge joining each seed to the rest of the tree. Exact rational drift certificates give finite-size product bounds and imply $113/200 \le \mu \le 1133/2000 < 17/30$. An independent certificate based on degree bands gives the upper bound $5657/9983$, also strictly below $17/30$.
```

## How this work was produced

Keep a visible disclosure before the references. The current section is specific about AI contributions, which is worth preserving. Shorten the main text to a paragraph covering the author’s role, the use of AI in the mathematics and writing, and the checks that were performed. Move the model-by-model contribution record to a versioned supplement. Attribute human verification, interpretation, and responsibility only as the author can personally attest. Keep the exact-arithmetic verification distinct from the AI-assisted review in the description.

> The author selected the problem, coordinated independent attempts, and set the criteria for proofs and exact-arithmetic checks. Generative AI systems developed proofs and certificates, wrote verification and simulation code, drafted the manuscript, and reviewed the arguments in separate runs. The certificate inequalities were checked with exact arithmetic. Independently implemented checks and direct computations on small trees tested parts of the verification. The supplement provides the certificate inputs, programs, and instructions for reproducing the checks.
