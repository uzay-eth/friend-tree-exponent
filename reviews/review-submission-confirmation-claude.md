# Confirming referee report on the revised manuscript

**Claude Opus 5.5 · maximum reasoning.** Fresh cross-family read of the complete revision.

**Verdict: release. No false or unjustifiable proof step remains. All 23 findings are confirmed fixed at the release head `d85b546`, and the supplement pin `2919520` is public and complete. The optional items in the last section do not block release.**

## Scope

- **Manuscript.** `paper/friend-tree-exponent.tex` at the release head `d85b546` of the submission branch. I read every section, including Appendices A and B, and the complete diff against `526f350`. That covers each intermediate state I reviewed (`5d181e8`, `2919520`) and the final re-pin.
- **Findings.** Every ID in `reviews/review-submission-findings.md` at the head: U1–U4, S1–S8, C1–C5, E1–E6, 23 items. I checked each against the revised text and, where useful, against the code, the supplement and the cited sources.
- **Also checked:**
  - Proposition 6.1 along both its pathwise and its mean route;
  - the three headline results;
  - the arXiv metadata and two-file source bundle;
  - the supplement pin, through the GitHub API;
  - the four original referee reports and the editor's report, to confirm that each point they raised maps to a finding;
  - possible regressions.
- **Materials.** The clone of the revised repository, the public GitHub record of the pinned commit, and the cited public sources for the citation items. I used no correspondence and no other working directories.

## Headline results

| Result | Status | Basis of this confirmation |
| --- | --- | --- |
| Theorem 1.1, common deterministic exponent of $`\mathbb E C_n`$ and of $`C_n`$ almost surely | **Stands** | Line-by-line reading of Section 2, including the new Coupling principle and Marginal site bound (U1–U3) |
| Theorem 1.2, motif products and $`113/200\le\mu\le1133/2000=17/30-1/6000`$ | **Stands** | Appendix A logic, archived arrays, my own rebuild, and the GMP checker rerun |
| Theorem 1.2 and Proposition 4.1, banded product with exponent $`5657/9983=17/30-1/299490`$ | **Stands** | Appendix B logic, primary and independent checkers rerun |
| Proposition 6.1, fixed exact degree $`d\ge2`$, pathwise and in mean | **Stands** | Line-by-line reading (section below) |

The prefactors reduce as printed, each in lowest terms:

```math
\frac{14608704616}{17527873930}=\frac{7304352308}{8763936965},\qquad
\frac{14606753031}{2989767809},\qquad
\frac{1539738114735}{3908512529}.
```

The banded tail margins are as printed:

```math
5657A-9983(B+a_*)=6235,\qquad 5657(A+B)-9983s=4755 .
```

Exclusion of $`17/30`$ has two independent supports, the motif and banded certificates. As §5 now states, the numerical interval $`[0.5650,0.5665]`$ rests on one exhaustive motif implementation, together with the independent direct and boundary tests.

## Independent reproduction

These runs used the verifier and data files of the release head. Between the run commit and `d85b546` these files are identical, except that the head adds the two archives, which I checked separately. Platform: Linux x86_64 with GCC 13.3.0, the reference compiler of Appendix A.

| Check | What I ran | Outcome |
| --- | --- | --- |
| Motif arrays rebuilt | `cert_build.c` with `-O3 -Wall -Wextra -fopenmp -lm`, arguments `1 1133 2000` and `2 113 200` | SHA-256 `3a4aef888a5b567a8ada05f9546f98c979d361da9f9b33ba6e794eb91fb30abe` (upper) and `9e966fba9aee6eae3da659b693d3bb9c481c1cea66635b5218d0887f182a60d9` (lower), equal to the Appendix A digests |
| Archived arrays | `xz -t`, then `xz -dc … \| sha256sum` on `verify/motif-bracket/reference/{upper,lower}.bin.xz` | Both archives pass the integrity test. Both decompress to exactly the Appendix A digests and are byte-identical (`cmp`) to my rebuild. |
| Motif GMP checker | `cert_verify_gmp.c` on both arrays, built plain and with `-fsanitize=undefined -fno-sanitize-recover=all` | `17329108 exact GMP inequality checks, 0 failures` per direction in both builds. Counts $`484138/160545/13519800/3164625`$. The checker prints `Phi(T3)=14606753031/10000000000, min W(u)=2989767809/10000000000` (upper) and `Phi(T3)=14608704616/10000000000, max W(u)=17527873930/10000000000` (lower), and the header `Domains: f(800,150) D0=24(KT=12,NG=601)`. |
| Banded primary | `check.cpp` with UBSan on `geometry.txt` and the decompressed recipe, `5657 9983` | `verified rows 226545617`, `minimum integer drift slack 3`, `minimum newborn reserve 4340525`, `weight range [3908512529, 1764716280715], centre 1539738114735`, `PASS`. Wall time 73 s. |
| Banded independent | `verify_banded.cpp` with UBSan, same inputs | `verified 226545617 finite states + analytic D>=2400 tail`, `min finite slack = 0`, newborn reserve `4341454`, prefactor `1539738113358/3908512428`, `PASS`. All seven band counts equal the Appendix B table. Wall time 160 s. |
| Direct tests | `bash verify/tree-drift/run.sh` | `[PY-EXACT 1] … ALL 19 sizes n=3..21: PASS`; `[PY-EXACT 2] … ALL 32972 non-isomorphic rooted trees … failures=0 for BOTH UPPER & LOWER: PASS`; `PASS banded: 23116 attachment histories through n=9, … failures=0`. Wall time 2 min 4 s. |
| Recipe digest | `sha256sum verify/banded-upper/recipe.txt.gz` | `584bb93c2797138ae3de54cb85c7d8e487bb7d7cc5320c021b5a64ec0cf0f60a`, as printed. The file has 12 697 rewards, 63 309 signed entries and 18 nonnegative entries. |

I recomputed the printed counts as well.

- **Index counts.** Inner states 572. Endpoint, parent and child indices 601, 555 and 290. Active centres 84 and 75.
- **Vertex rows.**

```math
484138=572+2\cdot241024+302+4\cdot302+4\cdot2 .
```

- **Edge pairs.** The 160 545 pairs equal $`555\cdot290-27\cdot15`$.
- **Motif test trees.** $`32972=A(14)-A(2)`$, where $`A`$ is the rooted-tree count A000081 with $`A(14)=32973`$.
- **Banded histories.** $`23116=\sum_{n=3}^{9}(n-1)!/2`$.

## Proposition 6.1 (fixed exact degree)

**Confirmed. Both routes are complete. I found no false or unjustified step.**

- **Seed and parameters.** The star on $`d`$ vertices rooted at its centre is a $`k=d`$ seed with $`M=2`$ and $`h=d+1`$. For $`d=2`$ the last stage is empty. The text notes, correctly, that preparation needs only $`M,k\ge2`$. The condition $`k\ge6M/\epsilon`$ belongs to the free-tree comparison.
- **Protection.** After preparation the root's neighbours are the base and $`d-1`$ leaves. A vertex gains a neighbour only through a ring of one of its neighbours. Silence of the base and the $`d`$ branch clocks for one unit therefore keeps the root at degree $`d`$ on $`[l+h,l+h+1]`$.
- **Independence.** Each flag is a function of one site's reserved clocks, its post-selection clocks, its marks and its coins. The protection unit uses increments on a new interval. The flags are therefore conditionally independent given $`\mathcal F_l`$, each with probability exactly $`q=p(2,d)e^{-(d+1)}`$.
- **Pathwise route.** With $`r\ge C(l)/6`$ sites, the Chernoff bound with $`\delta=1/2`$ gives

```math
\Pr\Big(\inf_{[l+h,l+h+1]}X_d<\tfrac{q}{12}C(l)\,\Big|\,\mathcal F_l\Big)\le e^{-rq/8}\le e^{-qC(l)/48}.
```

- The promotion intensity is at most $`C(t)`$ and each promotion adds one core vertex. The Yule domination is therefore a valid thinning coupling.
  - Chebyshev gives $`\Pr(C(l+H)>2ce^H\mid\mathcal F_l)<1/c`$.
  - Both failure events are path functionals. Separate couplings with overlapping windows are therefore legitimate for the first Borel–Cantelli lemma, after restriction to $`\{C(l)\ge e^{\mu l/2}\}`$.
  - Since $`h`$ is an integer, the unit windows tile all large times.
- **Mean route.** From size $`m=\lfloor ne^{-h-1/2}\rfloor`$:
  - $`\Pr(\tau_n\notin[h,H]\mid T_m)=O(1/m)`$ uniformly in the tree;
  - the flag count is at most $`C_m\le m`$;
  - $`\mathbb E C_n\le(n/m)\mathbb E C_m`$.

  Hence

```math
\liminf_{n\to\infty}\frac{\mathbb E N_{=d}(n)}{\mathbb E C_n}\ge\frac{q\,e^{-h-1/2}}{6}>0 .
```

- **Question 2 of Addario-Berry et al.** I checked §7 in the arXiv version. It asks whether, for fixed $`k`$, the number of degree-$`k`$ vertices is of the same order as the number of degree-$`(k+1)`$ vertices, and separately asks about large $`k`$. The comparability claim is correctly restricted to degrees at least two (S7), and the large-$`k`$ part is correctly left open.

## Findings by ID

### U1. Simultaneous preparations and independent contributions

**Confirmed fixed.** Lemma 2.3 makes each stage an explicit sufficient event on reserved inputs. I rechecked each part of it:

- **Common probability.** Thinning gives the same $`p`$ at both site types.
- **Selected leaf.** Its fresh clock is attached after selection. A leaf of $`u`$ lies in no other site's reserved set.
- **References and tests.** The free references and exterior tests are generated at every site, independently of preparation.
- **Failed tests.** These follow the ordinary law.
- **Site count.** It gives $`C\le6A_0+3J`$.
- **Stage probabilities.** Stage 1 is $`e^{-2}(1-e^{-1})/2`$ and stage 3 is $`e^{-3}M/(M+2)`$. Each final unit is $`e^{-(j+2)}/(k+1)`$, thinned by $`e^{-(k-1-j)}`$.

The ruling's three stated differences from the supplement's §2 each still give the stated lower bound:

- the smaller per-unit factor $`1/(k+1)`$ against $`1/d_j`$;
- the looser count of $`2A_0-1`$ bare paths against $`2A-3`$;
- the additional good-history restriction on $`Y`$.

### U2. Cut-edge comparison and iteration

**Confirmed fixed.**

- **Coupling principle.** It is now a named statement with a proof (strong Markov residuals, fresh marks, induction on the jump index, non-explosion).
- **Root acceptance.** The root accepts with probability $`d/(d+1)`$.
- **Base marks.** Since $`D_u\ge M`$, a base mark outside $`[0,1/M]`$ cannot choose the root.
- **First discrepancy.** The actual choice and birth are made, and then the drivers switch.
- **Acceptance product.** Given the free history, the acceptance marks are independent, so (1) keeps $`\prod_jd_j/(d_j+1)`$ inside the free expectation.
- **Good-history estimates.** I rechecked $`b\log(b/\epsilon)-b+\epsilon\ge6`$, $`8\log8-7>5`$, the Poisson lower tail with mean $`k\epsilon s/M\ge6s`$, and $`\eta=b\log2+(A+1)/M`$.
- **Marginal site bound.** It is a named statement. Its summands are path functionals on disjoint outward components.

### U3. Moment and time limits

**Confirmed fixed.**

- **Intensity.** It is bounded by the integrable $`N(t)`$, so $`F(t)\ge\mathbb E(t-\tau)^+\to\infty`$.
- **Order of limits.** It is stated: $`M,\epsilon,k`$ first, then $`s\to\infty`$, then $`\eta\downarrow0`$.
- **Holding times.** Non-explosion and independence of the holding times are stated at the start of Section 2.
- **Grid index.** It is deterministic in the mean iteration and random in the pathwise iteration.
- **Second moment.** $`\mathbb E Y^2\le k(k+1)e^{2s}`$, and conditional Chebyshev gives the following without truncation:

```math
q_j\le\frac{24\,\mathbb E Y^2}{C(t_{j-1})(\mathbb E Y)^2}.
```

- **Rest of the argument.** The martingale and Kronecker step, Markov with Borel–Cantelli, and the hitting-time transfer are correct.

### U4. Baseline drift at $`0.543`$

**Confirmed fixed.** The audit's baseline section now treats the two unbounded regimes explicitly. Both rational verifiers implement exactly that argument.

- **Leaf tail.** `certs9.py` and `certs6.py` evaluate the row at $`L=L_{\max}+1`$.
  - The focal gain increment there is zero, by the clamp.
  - Each remaining term is bounded at the endpoints $`L/D\in\{61/(61+c),1\}`$ and $`1/D\in\{0,1/(61+c)\}`$.
  - The verifiers first assert $`L_{\max}\ge K_0+1`$ and $`L_{\max}+1\ge D_{\rm cap}`$ (respectively $`P_{\max}`$), which makes the newborn type constant.

  Each term is affine in a single parameter, so this bound covers every $`L>60`$.
- **Core-degree overflow.** For $`c>10`$ the focal weight is constant forever.
  - Lower table: the weight is $`0`$, so the row is the nonnegative newborn charge.
  - Upper table: the verifier asserts $`\max(\text{newborn})\le\lambda G`$.

The (H2) deduction and the polynomial-scale hypothesis (1) now link to this argument. No table changed, and the headline results do not use these tables.

### S1. Motif potential at endpoints

**Confirmed fixed.**

- **Measure.** Appendix A defines $`\zeta_s`$ with $`\sum_a\zeta_s(a)i(a)=1/D`$ and averages independently at each endpoint.
- **Implementation.** The brute-force potential in `verify_python_cross.py` uses exactly this measure: a point mass for inner states, and $`(1-t)\delta_\infty+t\delta_{q_{\kappa h}}`$ with $`t=D_0/D`$ for tail states, product-averaged for $`g`$, $`h_V`$ and $`h_C`$.
- **Corner rules.** I re-derived the four rules of Lemma A.2.
- **New half-weight example.** It is correct. At $`D=48`$, $`t=1/2`$ and the credit is $`1/48`$, so the product is $`(d_0+d_1)/96`$. The two corners average to $`\max(d_0,d_1)/48`$, which bounds it for either sign (and the minimum does so for the lower bound).

### S2. Motif certificate identity

**Confirmed fixed.**

- **Digests.** Appendix A prints both digests and the reference build. My rebuild reproduces both digests.
- **Archives.** The supplement now archives both arrays losslessly. Both decompress to the printed digests.
- **Floating-point proposal.** Exact acceptance is separated from the floating-point proposal, and any other digest requires a fresh exact check.
- **Banded constants.** Their dependence on the weight-selection rule (R1) is stated in Appendix B and §5, with the independent integer-rule constants given.

### S3. Roles of the tree and row checks

**Confirmed fixed.** `check_tree` tests only the drift inequality and never compares with the core count, so the revised docstring is accurate. Section 5 separates three things from the universal reduction:

- the drift test on 32 972 shapes;
- the exact means for sizes 3–21;
- the identity and comparison test on selected trees.

### S4. Exact and cumulative degrees

**Confirmed fixed.** Proposition 6.1 is correct; see above. The paper no longer suggests that fixed exact degrees are out of reach.

### S5. Motif header and integer arithmetic

**Confirmed fixed for the shipped proof.**

- **Where the condition is used.** $`D_0\le K_F+1`$ is used in vertex case 5.
- **Where it is checked.** §5 checks it ($`24\le151`$), and the checker prints the header.
- **Sanitizer run.** My UBSan run of the GMP checker found no runtime errors.

### S6. Numerical arithmetic

**Confirmed fixed.** The two generators differ as stated:

- The separate C program uses $`\lfloor xm/2^{64}\rfloor`$ without rejection. Each index then gets $`\lfloor2^{64}/m\rfloor`$ or $`\lceil2^{64}/m\rceil`$ words, a deviation below $`2^{-64}`$.
- The primary C++ generator uses Lemire's rejection method.

The rounded-Gaussian thresholds $`10^9`$ and $`10^{14}`$ match the code.

### S7. Lower endpoint in the Question 2 sentence

**Confirmed fixed.** For $`n\ge3`$, $`N_{=1}(n)=n-C_n`$ while $`N_{=2}(n)\le C_n=n^{\mu+o(1)}`$ with $`\mu<1`$, so degrees one and two are not comparable. The paper, the degree-tail supplement and the pinned README now restrict the claim to degrees at least two.

### S8. Finite-size fits and the figure caption

**Confirmed fixed.**

- **Reproduction check.** On the archived fixed-size sample of 2 400 runs at 9 sizes, my code first reproduces both Table 1 fixed-size entries exactly: $`0.5660272\pm0.0002417`$ and $`0.5675481\pm0.0004706`$.
- **Weighted fit.** With weights $`1/\mathrm{se}_i^2`$ it gives $`0.5675819`$. The delete-one jackknife error is $`0.0004739`$ with the weights held fixed, or $`0.0004775`$ if they are recomputed in each deleted sample. Either way this matches the printed $`0.56758\pm0.00047`$ to its precision.
- **Exact slope.** From the archived exact means, $`20(a_{21}-a_{20})/a_{20}=0.583346`$, matching $`\approx0.5833`$.
- **Caption.** The caption now says "region around" the published range.

### C1. Addario-Berry and coauthors

**Confirmed fixed.** In the arXiv version (only v1 is posted):

- Theorem 4.8 gives $`c_k<\liminf_nX_n^{\ge k+1}/X_n^{\ge k}\le1`$ for $`k\ge2`$.
- Proposition 6.20 gives $`(\sqrt{13}-3)/6`$.
- Proposition 6.22 gives $`\sqrt2-1/2`$.
- The abstract states $`n^{9/10}`$.

The revision cites the journal record and both locators. I did not see the journal text, so the journal numbering, Theorem 4.6, rests on the R2 and R4 reports.

### C2. Exact degrees and the historical conjecture

**Confirmed fixed.** The syndicated post, dated 20 September 2022 at 18:12 UTC, contains the quoted sentence verbatim. The degree notation and scaling law are attributed to Krapivsky and Redner, Eq. (1) and §5.1.

### C3. Neighbouring growth models and numerical fits

**Confirmed fixed.** The text now has [SK, ES], [HIH], [LSS] with its triangle seed, Cannings–Jordan Theorem 2.1, and the Hartle et al. 0.569 as a finite-size fit. I confirmed the three added references bibliographically:

- Krapivsky–Redner, *Phys. Rev. E* 63 (2001), 066123 (ancestor redirection);
- Gabel–Krapivsky–Redner, *Phys. Rev. E* 88 (2013), 050802(R) (enhanced redirection);
- the survey, *Indian J. Phys.* 98 (2024), 3777–3787, arXiv:2305.10628.

### C4. Supplement pin and long-term archiving

**Confirmed fixed.** The bibliography pins `2919520d6673198777ab13a3ae80b6a2c650d37e`.

- **Public record.** The GitHub API returns this commit, and its recursive tree has 161 blobs. All are identical to the local commit, with none missing, none extra and no mismatches.
- **Contents.** It contains the two archived arrays, all four original referee reports, the supporting audit and the recipe. It also contains the offline build notes in `paper/README.md` and the README qualified to $`k\ge2`$.
- **Difference from the release head.** Only `paper/friend-tree-exponent.tex`, the PDF and the bundle differ, so every supplement path named in the paper resolves at the pin.

During this read-back the previous pin `bb536bc` lacked those build notes, which were added later in `e4d4081`. I reported this, and the re-pin closes it. The DOI deposit is an owner decision, recorded as deferred.

### C5. Degree-tail corollary and historical coefficients

**Confirmed fixed.**

- **Reference.** "Claim 2" is replaced by a reference to the cross-charging products.
- **Prefactors.** The finite-size table in `REPORT.md` gives $`0.50220838\ge0.5022`$ and $`3.97372002\le3.974`$.
- **Coefficients.** Gautschi's inequality gives $`0.5022\,(3/4)^{0.45}=0.4412`$ and $`3.974\,(2+1/0.42)\,3^{0.58}=32.93`$, as printed.
- **Higher-precision value.** The article uses the replayed $`0.5657242\pm0.0000175`$.

### E1. Title, abstract and opening

**Confirmed fixed.**

- **Title.** It names the observable, and the PDF metadata title matches.
- **Order.** The theorems precede the history.
- **Abstract.** It now spells "neighbour", consistent with the body. It states the cut-edge comparison and the double exclusion of $`17/30`$. The sheet copy is identical, ASCII only and 743 characters, under arXiv's 1 920-character limit.

### E2. Length and certificate structure

**Confirmed fixed.**

- **Layout.** The main text runs to page 14. Appendices A and B (pages 15–26) keep the complete reductions.
- **Worked three-vertex transition.** It is correct: rates $`2`$ and $`1`$, and $`C_4=2`$ with probability $`1/3`$.

### E3. Marginal versus simultaneous preparations

**Confirmed fixed.** The marginal argument is kept as the named Marginal site bound. The block lemma is stated separately with the second-moment route.

### E4. Later sections and production note

**Confirmed fixed.**

- **Earlier certificates.** The account is one paragraph.
- **Simulations.** The figure, protocol, horizons and limitations are retained.
- **Production note.** It is one paragraph and states that the sole listed author is responsible.

### E5. Notation and navigation

**Confirmed fixed.**

- **Local cutoffs.** $`H`$ and $`T`$ are identified as local to Appendix B.
- **Scaled drift.** The two-event expression is named the scaled drift.
- **Degree-count notation.** $`N_k`$ is identified with $`N_{\ge k}`$.
- **References.** These lead directly to the lemmas.
- **Section 6.** Its opening now separates the in-paper proof of Proposition 6.1 from the supplementary bounds.

### E6. Worked examples, direct tests and independent verification

**Confirmed fixed.**

- **Named statements.** Both are present.
- **Strict inequality.** The planting proof states $`L_u/D_u\ge(M+1)/(M+3)>M/(M+2)\ge U`$.
- **Final-stage bound.** The degree bound is correct: at most $`k-1`$, including the exterior edge.
- **Appendix B example.** The two-step degree-band history is correct. After the promotion both vertices are $`(2,1)`$ and record each other. The gain makes the centre $`(3,1)`$ and updates its recorded degree to three.
- **Test coverage.** §5 states that the small exhaustive tests contain no tail state and no band crossing (a tail state needs 25 vertices, a crossing needs degree 16). It cites R1's tests, which the R1 report records: 135 trees per direction up to $`n=3500`$, and 16 histories up to size 220 crossing 16, 32, 80 and 200.
- **Build notes.** The offline replay commands match the ones I ran.

## arXiv metadata and source bundle

- **Bundle contents.** The bundle holds `friend-tree-exponent.tex` and `local_slopes.png` at top level, both byte-identical to the repository files at `d85b546`.
- **Clean `pdflatex` build.** In an empty directory, three passes each exit 0. The log reports `Output written on friend-tree-exponent.pdf (27 pages, …)` with no overfull or underfull boxes, no undefined references and no warnings (TeX Live 2023, pdfTeX 1.40.25). Numbering: Theorems 1.1–1.2, Lemmas 2.1–2.3, Propositions 4.1 and 6.1, Lemmas A.1–A.3, Lemmas B.1–B.2.
- **Tectonic rebuild.** Tectonic 0.17.0 with `SOURCE_DATE_EPOCH=1795392000` gives a PDF byte-identical to the committed `paper/friend-tree-exponent.pdf` (SHA-256 `5b44197b148b5e967a18a2d12683b12d6654625f6e078c1d6f0693e0fa057911`, 27 pages).
- **Submission sheet.**
  - The title, author and affiliation match the PDF.
  - The comment line's 27 pages and 1 figure match the build.
  - The categories (math.PR, cross-listed to cond-mat.stat-mech) and MSC codes 60C05, 05C80 and 05C05 are appropriate.
  - The quoted arXiv policy sentences on generative AI tools and on endorsement match the current policy pages.

## Regressions

None found. The revision changes no constant, certificate input or verifier. Since `526f350` the verifier tree gained only the two archives and one docstring. Every new or altered mathematical sentence was rechecked above:

- Section 2: the named statements, the strict inequality and the degree bound in Lemma 2.2;
- the two appendix examples;
- the U4 argument;
- the S8 numbers.

The bundle builds cleanly, and the committed PDF matches its source.

## Owner items (not paper defects)

- Submission requires the author's personal responsibility for the arguments and citations, which arXiv policy and R4 both stress. The production note and submission sheet state this.
- The DOI or journal deposit of the supplement is deferred (C4).
- arXiv endorsement for a first math.PR submission is still needed (submission sheet).

## Optional improvements (non-blocking)

1. Three supplement landing texts still say the 0.566/0.569 difference "remains unexplained": `README.md` (results paragraph), `REPORT.md` (simulation estimates) and `simulations/README.md` (slope paragraph). §7 now attributes it to the fitting weights.
2. `REPORT.md` (first paragraph) still says the cumulative lower bounds give no lower bound for exact-degree counts. It never mentions Proposition 6.1, although `README.md` does.
3. The degree-tail link to "the two prefactors in the finite-size bounds table" opens the top of `REPORT.md` rather than the table.
4. [GKR] is a Rapid Communication, usually cited as 050802(R).
5. $`h`$, $`b`$ and $`A`$ still carry different local meanings in Section 2 and the appendices. Each is defined where it is used.

## Read-back of the integrated difference

This report reads back the complete integrated revision from all four referees, the editor and the independent degree reader:

- the diff `5d181e8 → cfae347 → 2919520 → d85b546`;
- the final findings file with its 23 IDs.

The one gap found during the read-back was the missing build notes at the earlier pin. It was reported during the read-back and is closed by the re-pin to `2919520`. According to the integration notes, only the review index and this confirmation text are still to change. Both lie outside the paper and the pinned supplement paths it cites.
