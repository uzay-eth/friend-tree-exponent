VERDICT: ACCEPT

# Cold review of the preprint "Growth of the non-leaf population in the random friend tree"

Reviewer: Claude Opus 5.5 (`claude-opus-5-5`), 23 November 2026. I read the paper cold as an independent reviewer and wrote none of its text, proofs or code.

**Reviewed PDF.** [`paper/friend-tree-exponent.pdf`](https://github.com/uzay-eth/friend-tree-exponent/blob/1db5264b672cf4942b2278cc7078d4a96285258f/paper/friend-tree-exponent.pdf), 24 pages, SHA-256 `e95e17a2f6aa2c55ac1f73e6e360870083d2fc7253aa08e537a8fef597fad905`.
**Its source.** [`paper/friend-tree-exponent.tex`](https://github.com/uzay-eth/friend-tree-exponent/blob/1db5264b672cf4942b2278cc7078d4a96285258f/paper/friend-tree-exponent.tex), SHA-256 `f5140fab890ff03b8dd59f50ec803d6c22216790e6e9b93cfd2077dae2cf096e`.
**Cited supplement.** Revision `c4653036b993bf81799c51b60cdd18a500142147` of this repository. The proofs, simulations and verifiers match that revision byte for byte.

**Note added 23 November 2026.** The byte-for-byte match above held when this review was written. Since then, under `proofs/`, `verify/` and `simulations/`, the only file that differs from the cited revision is `simulations/README.md`. Its replay paragraph and timing table now report full replays of the three batches and their measured times. The proofs and verifiers still match the cited revision byte for byte. The PDF and source digests above are unchanged, and so is the verdict.

## Steps

1. **Build: PASS.** Rebuilding from the source and figure with Tectonic 0.17.0 and `SOURCE_DATE_EPOCH=1795392000` gives a byte-identical PDF (the SHA-256 above), 24 pages, with no TeX warnings.
2. **Cross-references and citations: PASS.** No `??` appears in the PDF text. Every `\ref`/`\eqref` has a label. Every `\cite` has a bibliography entry and every entry is cited.
3. **Figure: PASS.** The included `simulations/local_slopes.png` passes `sha256sum -c SHA256SUMS`. Its title now reads "Growth of the mean core size at fixed time". The data, panels, sample sizes (912 and 91 trajectories) and the reference line agree with the caption and with `local_slopes.tsv`.
4. **Existence of the mean exponent: PASS.** Section 2 agrees with [`proofs/existence.md`](../proofs/existence.md): eligible sites at least $`C/4`$; planting probability $`p(M, k)`$ with stage factors $`e^{-2}(1-e^{-1})/2`$, $`e^{-1}\Pr(\mathrm{Pois}(1)\ge M)`$, $`e^{-3}M/(M+2)`$ and $`(e^{-(k+1)}/(k+1))^{k-2}`$; the cut-edge comparison with acceptance weight $`d/(d+1)`$ and exterior rate $`1/M`$; $`b\log(b/\epsilon)-b+\epsilon\ge6`$; $`\eta=b\log2+(A+1)/M`$; the order of limits; and de-Poissonisation.
5. **Path-defined $`X_u`$ (first-moment step): PASS.**
   - **Definition.** $`X_u`$ is now a functional of the actual tree path alone. It counts endpoint cores on every side $`w`$ of an edge $`uw`$ whose $`w`$-component at time $`t+h`$ contains no vertex that was core at time $`t`$.
   - **Disjointness.** Each such side is the component of $`w`$ after removing the connected old-core subtree, and later leaf additions keep those components apart. Hence $`\sum_u X_u\le C(t+h+s)`$ on every path.
   - **Per-site bound.** A successful planted branch is one of the components counted by $`X_u`$. So $`\mathbb E[X_u\mid\mathcal F_t]`$ may be bounded with that site's own law-preserving coupling.
   - **Consequence.** Summing the bounds no longer needs simultaneous success or any comparison of success events across sites.
6. **Sample-path limit: PASS.** The section is unchanged since my previous reading and agrees with [`proofs/almost-sure-exponent.md`](../proofs/almost-sure-exponent.md).
   - **Independent copies.** The block argument uses at least $`C/6`$ disjoint sites. Thinned flags give each site common success probability $`p`$, and given the past at the block start the copies of $`Y`$ are independent.
   - **Failure probability.** Chebyshev gives $`16R^2/(m(\mathbb EY)^2)\le96R^2/(C(\mathbb EY)^2)`$.
   - **Limits.** The martingale/Kronecker count of failed blocks gives the lower limit, and Markov plus Borel–Cantelli gives the upper limit.
7. **Motif certificate: PASS.** Section 3 agrees with [`proofs/bracket-motifs.md`](../proofs/bracket-motifs.md).
   - **Proof structure.** States and rates, tail interpolation, tables, allocation and FY splits, newborn bookkeeping, corner lemma, the five vertex cases, and the edge, path and cherry rules all match.
   - **Row counts, recomputed:**
     - paths: $`555\cdot84\cdot290=13{,}519{,}800`$;
     - cherries: $`75\cdot290\cdot291/2=3{,}164{,}625`$;
     - edges: $`555\cdot290-27\cdot15=160{,}545`$;
     - vertex entries: $`2\cdot801\cdot151+2=241{,}904`$.
   - **Checker scaling.** The integer multiplier stated for path and cherry rows matches `cert_verify_gmp.c`.
   - **Constants.** The comparison constants and both prefactors match the proof.
8. **Banded certificate: PASS.** Section 4 agrees with [`proofs/upper-bound-banded.md`](../proofs/upper-bound-banded.md).
   - **Inputs.** The band table is identical to `geometry.txt`, and the band counts sum to $`226{,}545{,}617`$. The `recipe.txt.gz` SHA-256 matches. The recipe has 12,697 rewards, 63,309 signed entries and 18 arrow entries.
   - **Tail margins.** $`5657A-9983(B+a_*)=6235`$ and $`5657(A+B)-9983s=4755`$.
   - **Proof steps.** The history invariant, the two exact cancellations, the rate intervals and the integer acceptance test all match the proof.
   - **Checker output.** Weights, minimum and root agree with the checker output.
9. **Direct tree tests: PASS.**
   - **Motif test.** $`32{,}972`$ shapes is the number of rooted trees on 3 to 14 vertices with non-leaf root.
   - **Banded test.** $`23{,}116=1+3+12+60+360+2520+20160`$ histories is every attachment sequence from $`T_3`$ through size 9.
   - Both counts agree with the test output, and `verify/tree-drift/run.sh` is present at the cited revision.
10. **Related results and baseline tables: PASS.**
    - **Statements.** The Section 6 statements match Theorem B with (H1)–(H5), the cross-charging certificates, the polynomial-scale theorem and Theorem D(a)–(d).
    - **Baseline references.** The baseline tables are now pointed to `verify/baseline/certs9.py` (lower), and to `certs6.py` and `find_upper_lin.py` (upper). Their states and rates are in [`proofs/supporting-results-audit.md`](../proofs/supporting-results-audit.md).
    - **Baseline section against the code.** That section matches the verifiers:
      - creator label $`P=\max(3, \min(D_{\rm creator},20))`$;
      - pendant gain $`i+1/D`$ and pendant conversion $`i/(i+1)`$;
      - class-B gain in $`[1/D,1/D+1/2]`$;
      - degree growth in $`[0, d]`$;
      - caps 60 and 10.
    - **Comparison constants.** These are now stated: weight at least one (upper), and $`\sup f`$ (lower).
11. **Prior work: PASS.** I checked each citation against the public source.
    - **Krapivsky–Redner, §5.1:** $`\mu\approx0.566`$, fitted slopes 0.5652 to 0.5673, and $`N_k`$ counting exact degree.
    - **Addario-Berry et al.:**
      - Proposition 6.20 gives $`(\sqrt{13}-3)/6`$.
      - Proposition 6.22 gives $`\sqrt2-1/2`$; both bounds are almost sure.
      - Theorem 4.8 gives comparability of successive cumulative counts.
      - Saramäki–Kaski are credited with introducing the model, and Cannings–Jordan with $`n-o(n)`$ leaves.
    - **Hartle et al., Fig. 6:** 0.569 from nonlinear least squares, $`N=10^2`$ to $`10^6`$, 2400 realizations per point. The citation Phys. Rev. E 114, 034314 (2026) is correct.
    - **Strogatz post:** dated 20 September 2022.
12. **Simulations: PASS.** Every number I checked agrees with the archived outputs:
    - the fit table and adjacent slopes;
    - the final checkpoint $`0.5657203\pm0.0000038`$ at mean size $`8.117\times10^{13}`$;
    - the pooled estimates $`0.5657230`$ and $`0.5657253\pm0.0000032`$ over $`22{,}223{,}928{,}834`$ firings;
    - the 48-run value $`0.5657242\pm0.0000175`$;
    - the fit window $`t=4, \ldots,12`$.

    The text correctly says that these estimates give no bound on $`\mu`$.
13. **Supplement availability: PASS. This closes required change R1 of my previous review.**
    - **Pinned citation.** The bibliography now cites the immutable revision `c4653036b993bf81799c51b60cdd18a500142147`, which is publicly reachable.
    - **Cited files present.** At that revision every file the paper names exists: `verify/run_all.sh`, `verify/tree-drift/run.sh`, the four proof files under `proofs/`, `proofs/supporting-results-audit.md`, `proofs/bracket-motifs-code-map.md`, the three baseline scripts and `simulations/`.
14. **AI disclosure: PASS. This closes required change R2 of my previous review.**
    - **What it now states.** GPT-6.2 Astra drafted the manuscript and its LaTeX. Claude Opus 5.5 expanded the motif, existence and sample-path proofs. The archived continuous-time simulations came from earlier model-led attempts, and Astra recovered their seeds and analysed the fitting window. GPT models packaged the direct tree tests. Claude Opus 5.5 and GPT-6.2 Sol tested the simulations against independent implementations.
    - **Unchanged.** The statement that no human referee has read the proofs remains.
15. **Standalone nature: PASS.** Theorems 1.1 and 1.2 are proved in the paper. The only parts left to named checkers are the finite arithmetic and enumeration, and the paper gives their inputs, counts and digests. The Section 6 results are presented as recorded results whose proofs are in the supplement.
16. **Contact address: PASS.** The owner explicitly selected the contact address in the title block.

## Limits of this review

- **Certificate suites not rerun.** I did not rerun `bash verify/run_all.sh` or `bash verify/tree-drift/run.sh` on the cited revision. For the certificate and direct-test figures I relied on recorded outputs from runs of the same checkers on the same inputs, together with my own recomputation of the counts, hashes and integer margins listed above.
- **No line-by-line code audit.** I read checker code only where the paper states a specific implementation detail.
- **Attributions not recorded in the supplement.** The proof files do not record that Claude Opus 5.5 expanded the motif write-up. The owner confirmed it from the drafting record, which I did not see myself. Which attempts produced the archived simulations is taken on the owner's word.
- **Scope.** The review covers the mathematics and presentation of this PDF. It is not a substitute for human refereeing.

## Non-blocking notes

- **Overloaded notation.** Several symbols carry different meanings in different sections (for example $`M`$, $`H`$, $`J`$, $`h`$, $`A`$, $`B`$, $`T`$, $`m`$). Each use is defined locally, so readers can follow it, but a notation note would help.
- **Unused labels.** Seven labels are never referenced. This is harmless.
