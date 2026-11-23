# Pre-submission referee findings

The four independent reports concern the manuscript at `526f350`. R1 is Claude Opus 5.5, R2 GPT-6.2 Astra, R3 Gemini 4 Neon and R4 Claude Fable 5.2. A separate GPT-6.2 Sol editor considered the structure. The reports use original section and line numbers. This list records the rulings on the revision. Related instances of the same issue are grouped.

### U1. Original §2, simultaneous preparations and independent contributions

**Unjustified step · fixed.** Raised by R1, R2, R3, editor. A selected leaf can depend on exterior rings. Allocating one clock per starting site alone does not establish independent success. In the independent-contributions lemma each stage has an explicit sufficient event on reserved clocks, uniform marks and thinning coins. The identity of a leaf is decided before attaching its fresh clock. The proof now treats failed sites and conditions only on the history revealed so far. The expanded sample-path proof supplies the same construction.

### U2. Original §2, cut-edge comparison and iteration

**Unjustified step · fixed.** Raised by R1, R2, R3. The lower bound needs the root-acceptance product inside the free-tree expectation. The revision defines the free history and records which marks remain independent. It performs the actual birth at the first discrepancy, then switches drivers. The marginal summands are actual path observables on disjoint cut sides and may be evaluated by separate couplings.

### U3. Original §2, moment and time limits

**Unjustified step · fixed.** Raised by R1, R2, editor. The revision supplies the integrable promotion intensity, the quantifier order in the block iteration, non-explosion and independence of the size clock, and the distinction between the first positive random grid index and a deterministic index. The independent-block lemma gives a finite second moment for each fixed block directly from the Yule law. Conditional Chebyshev then proves the needed vanishing failure probability without truncation.

### S1. Original §3, motif potential at endpoints

**Wrong or imprecise statement · fixed.** Raised by R2, R3, editor. The old endpoint measure was never defined. Appendix A uses the previously defined endpoint measure at each end of a path and each child of a cherry. The checker uses precisely these interpolations.

### S2. Original §§3, 5, motif certificate identity

**Wrong or imprecise statement · fixed.** Raised by R1, R2, R4. The public proof already listed the reference array digests. The article did not. R1 reproduced identical arrays across thread counts; R4 found a three-byte difference in the upper array with native compiler optimizations. These observations are compatible. Appendix A gives both digests and the reference build, separates the floating-point proposal from exact acceptance and requires a new exact check for different arrays. The matching supplement now makes that qualification.

### S3. Original §5, roles of the tree and row checks

**Wrong or imprecise statement · fixed.** Raised by R2. The Python shape loop tests drift on 32,972 rooted shapes and reads archived means. It does not test core comparison on every such tree or prove the infinite tail from these examples. The separate identity test checks terms and comparisons on selected trees. Section 5 and the Python module description now distinguish those checks from the universal reduction, its finite arithmetic and its tail inequalities. Executable tests and acceptance rules are unchanged.

### S4. Original §§1, 6, 8, exact and cumulative degrees

**Wrong or imprecise statement · fixed.** Raised by R2. The old upper estimate is correct, but leaves the impression that even fixed exact-degree growth is inaccessible. Consecutive cumulative lower bounds cannot be subtracted. The new fixed-degree proposition protects roots of disjoint planted stars over overlapping time windows and bounds core growth during each window by a Yule process. This proves eventual pathwise comparability and a positive lower ratio of means for each fixed degree. The matching proof is in the degree-tail supplement. Limits for proportions and their tail in the degree remain open.

### C1. Original §§1, 6, 8, Addario-Berry and coauthors

**Citation · fixed.** Raised by R1, R2, R3, R4. The 2024 preprint is genuine. The cumulative comparison is Theorem 4.8 there and Theorem 4.6 in the published article. Cite both locators and the journal record, and attribute the normalization and degree questions to §7. The quoted upper exponent 0.9142 is what Proposition 6.22 yields; the source abstract states 0.9. The revision identifies both places without inferring an unprinted improvement.

### C2. Original §1, exact degrees and the historical conjecture

**Citation · fixed.** Raised by R2, R4. The September 2022 post itself gives the fraction. Attribute the degree notation to Krapivsky and Redner and the scaling prediction to their equation (1) and §5.1. Give a short exact quotation from the post. The accessible opening post and attached figure suffice for this attribution; the rest of the social-media thread is unnecessary.

### C3. Original §§1, 7, neighbouring growth models and numerical fits

**Citation · fixed.** Raised by R2, R3, R4. Identify the one-step random-walk literature, the mediation-driven rule and the friend-of-a-friend specialization, noting the triangle seed in the last paper. Add a precise Cannings–Jordan locator. Separate the finite-size slope of Hartle and coauthors from the asymptotic exponent. Primary papers and edition-specific details appear in the reports. The redirection survey covers the related model family. The models with different redirection rules and the general branching-process references suggested by R3 provide background; the article uses no theorem from them.

### E1. Title, abstract and opening

**Exposition · fixed.** Raised by Editor, R2. Put the common mean and pathwise exponent and the product bounds before the historical discussion. The title names the observable. The short ASCII abstract describes the cut-edge comparison and explicitly states the exclusion of 17/30 by both constructions. The original abstract was already within the arXiv length limit.

### E2. Original §§3–5, length and certificate structure

**Exposition · fixed.** Raised by Editor, R1, R2. The main text explains the drift, sign cancellations, product iteration and independent upper bound. Appendices A and B retain the complete motif drift and corner bounds, routing history and integer row construction. A worked three-vertex transition, range references, digests and precise checker scopes connect the code to the mathematics. Both proofs remain in the submission.

### E3. Original §2, marginal versus simultaneous preparations

**Exposition · fixed.** Raised by Editor, R2. The editor suggested deriving both limits from a single block lemma. We kept the short marginal expectation argument because it does not require simultaneous independence. The pathwise proof now states its stronger independent-block lemma separately and uses the second-moment simplification.

### E4. Original §§6–8 and production note

**Exposition · fixed.** Raised by Editor, R4. Shorten the account of earlier certificates, retain the empirical figure, protocol and limitations in the submission, and reduce the production note to a single paragraph naming the human and model contributions. The different size and time horizons and their sampling uncertainties remain explicit.

### S5. Original §§3, 5, motif header and integer arithmetic

**Wrong or imprecise statement · fixed.** Raised by R3. An overflow vertex must be a tail vertex. Section 5 explicitly checks the header inequality 24 ≤ 151 and points to the per-header and intermediate-range audit. The independently instrumented C replays also passed. This resolves the shipped proof; a request to validate arbitrary future binary inputs would change the scope of the verifier and is outside the mathematical claim.

### S6. Original §7, numerical arithmetic

**Wrong or imprecise statement · fixed.** Raised by R3. The independent simulator selects an index by truncating a 64-bit product. For m indices each individual probability differs from 1/m by less than 2 to the power −64; the revised methods disclose this alongside the rounded Gaussian draws. These are numerical measurements. The certificate proofs do not use them.

### C4. Supplement pin and long-term archiving

**Citation · fixed.** Raised by R3, R4. The cited original commit exists and supplies the data used in the paper. The revised link identifies the immutable commit containing the updated proofs and checker documentation. A permanent journal or DOI deposit would help long-term preservation; no deposit was made in this round. That publication decision can be made with the eventual venue.

### E5. Original §§3, 4, 6, notation and navigation

**Exposition · fixed.** Raised by R3. In the degree-band section, identify H and T as local cutoffs and describe the two-event expression as the scaled drift. In the degree discussion, identify the supplement’s cumulative-count notation. Main-text references now lead directly to the decomposition and the certificate lemmas, with the complete reductions in the appendices.

## Disposition

| Severity | Fixed | Open | Requires owner ruling |
| --- | ---: | ---: | ---: |
| Error in a proof | 0 | 0 | 0 |
| Unjustified step | 3 | 0 | 0 |
| Wrong or imprecise statement | 6 | 0 | 0 |
| Citation | 4 | 0 | 0 |
| Exposition | 5 | 0 | 0 |
| **Total** | **18** | **0** | **0** |

No false step was identified in the existence proof, the two-sided motif bound or the independent degree-band upper bound. Outstanding research questions in §8 concern sharper theorems. This is a model review. A human probability referee has yet to read the arguments.
