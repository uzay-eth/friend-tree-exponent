> Packaging note: internal branch names on lines 6 and 42 were redacted for public release; all mathematical content and commit hashes are unchanged.

# Cold cross-family review of the combined report — Claude family

* **Referee:** Claude Opus 5.5 (`claude-opus-5-5`, Claude Code harness), cold, no prior involvement with any lead.
* **Target:** combined report at commit `a8781f2dfd71b551364c145f4f6f357c1d0054e0`; `FINAL.md`, `OZET-TR.md`.
* **Scope:** cross-work compatibility of the stated theorems (I did not re-review the underlying proofs from zero), and a PR review of the two report files. I edited no existing file and made no commit.

## Part 1: Cold check of the combined statement

- **P1 (model equivalence):** All five works use the same model: start from the single edge T_2, n=|V(T_n)|, the same C_n (C_2=0, C_3=1). Neon's rooting at the centre of the deterministic T_3 adds no randomness. The Astra/Baseline/Carbon history marks only need the tower property on F_n. Astra's rate-1 clocks give total rate |V| and a uniform source independent of the holding time, so the jump chain has the exact law of (T_n) from two vertices (N(t) is a sum of two Geom(e^{-t})). §1.1 probabilities are correct for n≥3. **Consistent.**
- **P2 (Astra 1 + Neon 1 ⇒ μ∈[113/200,1133/2000]):** Valid. Astra Claim 1 makes the limit exist, so it equals both liminf and limsup. Neon's product bounds then give 0.5650 ≤ μ ≤ 0.5665. The margin is 17/30−1133/2000 = 1/6000 exactly (Fraction). The fallback 563/1000 ≤ μ ≤ 5657/9983 (Astra only) is also valid, with margin 1/299490 exactly. **Consistent.**
- **P3 (Astra 1–2 a.s. with Neon 1 and Astra 3–4):** The a.s. limit equals the same deterministic μ, so log C_n/log n → μ ∈ [0.5650,0.5665] a.s. Neon Cor. 2 (limsup ≤ 0.5665) and the a.s. lower bounds 0.543/0.55 are implied by this, and none contradicts it. Every certified bracket is nested: [0.5650,0.5665] ⊂ [0.563,0.56666] ⊂ [0.563,0.568] ⊂ [0.55,0.58] ⊂ [0.543,0.586]. I computed exact E[C_n] for n=3..15 by full rational enumeration of unlabelled trees. All 9 certified finite-n bounds bracket it (Neon, Astra 3/4, Fable 2, Baseline; `ALL OK`). **Consistent.**
- **P4 (product ↔ Γ ↔ Gautschi arithmetic, §1.3 and §2.1):** I recomputed everything independently (Python `fractions`, `math.gamma`):
  - The identity P_a(3,n) = 2Γ(n+a)/(Γ(3+a)Γ(n)) checks out numerically for n=3…199999.
  - The Gautschi sandwich (3/4)^{1−a}·2/Γ(3+a)·n^a < P_a(3,n) < 2/Γ(3+a)·n^a holds for n≥3.
  - Neon's constants come out at 0.41166674 (>0.41) and 2.73020494 (<2.74).
  - All 9 c_Γ entries match to ≤1e-8. The last-digit differences, e.g. Fable upper 3.97372001 vs 3.97372002, come from 8-digit inputs.
  - All 9 margins are exact (1/600, 1/6000, 1/299490, 11/3000, −1/750, 1/60, −1/75, 71/3000, −29/1500).
  - In 3(b), 2.74·(2+1/0.4335)·3^0.5665 = 21.98847 and 2.73020·(…) = 21.90986, both <22.
  - Fable's own corollary gives 0.4412 and 32.925.
  - The only defect is the labelling of the 21.90986 parenthetical (see N3). **Consistent.**
- **P5 (Carbon vs Fable/Astra/Neon):** The non-upgrade claim is accurate as a statement about hypotheses. Carbon needs radius-1 per-vertex tables with degree-H saturation, f_+≥1 and λ<2δ. Fable 0.55 uses cross-charging transfers, Astra 0.560/0.563 use a distance-2 witness, pooled history marks and (0.563) signed transfers, and Neon 0.5650 uses motif terms with a non-saturating 1/D tail. The stated implication (Astra 1–2 + Neon 1 ⇒ C_n/n^α → ∞ for α<0.565) is correct. **But the stated relationship is incomplete:** Fable Claim 2 (liminf log C_n/log n ≥ 0.55 a.s.) does not use Astra or Neon, and on its own it already gives C_n/n^{0.543} → ∞ a.s. (because 0.55 > 0.543). So §1.3.4 overstates Carbon's "independent value" when it says Carbon proves this "without Astra Claims 1–2 or Neon Theorem 1", and §2.1 4(c) is a consequence of 4(b). Carbon's genuinely new content is the method (no loss at the certificate's own exponent) and the a.s. convergence of Φ_n^−/n^{0.543} in (0,∞]. **No mathematical inconsistency; a wording fix is needed (B1).**
- **P6 (dependencies, domains, overclaims in §2.1/§2.2):**
  - Domains match. Neon, Astra 3/4 and Fable 2 hold for n≥3. Fable D(c) only needs a_m bounds for m≥3, since m′=⌊n/8k⌋≥3 and m0=⌈2n/k⌉≥3 for n≥24k, k≥2. Carbon roots at T_3. Astra converts to discrete size via hitting times. No n≥2 vs n≥3 mismatch.
  - D(d) correctly needs μ∈(0,1) from Astra 1 plus a positive / sub-1 bound.
  - §2.2 is appropriately cautious: no regular variation, no pointwise tail, simulations labelled heuristic. 0.000944/0.000003 ≈ 315σ, so "~300 standard errors" is fine.
  - The only overclaim is the Carbon point in P5. **Consistent.**

VERDICT: CONSISTENT

## Part 2: PR review of `FINAL.md` and `OZET-TR.md`

### Evidence (commands I ran, outputs abridged)

```text
# arithmetic (Python, exact Fractions + math.gamma): all c_Gamma, margins, Gautschi, D(c) constants -> as in P4
# exact E[C_n], n=3..15, rational enumeration of unlabelled trees; all 9 certificate bounds checked -> ALL OK
#   e.g. E[C_15] = 3.179932316303627; Neon lower/E = 0.672, upper/E = 3.949
git diff --stat <pinned> HEAD -- $(git ls-tree -r --name-only <pinned>)    # for 45bd14f ec15a72 9d18e22 9ad27a2 1cbd617
  -> only .gitignore differs (expected: merge keeps Fable superset); git diff 3ac8946 HEAD adds only FINAL.md, OZET-TR.md
pinned commits: neon=9d18e22, astra=ec15a72, fable=9ad27a2, carbon=1cbd617, baseline=45bd14f
markdown link check over FINAL.md + OZET-TR.md: all 33 relative links resolve (0 MISSING)
fresh clone in /tmp at a8781f2: bash attempts/strogatz-exact/cert/run_all.sh
  -> 17329108 checks, 0 failures (both sides); PY-EXACT 1/2 PASS; 6/6 negative controls rejected; ELAPSED=21.59 EXIT=0
  -> SHA-256 upper 3a4aef88...0abe, lower 9e966fba...60d9  (byte-identical to FINAL §4.4)
fresh clone: bash attempts/strogatz-astra/banded/reproduce.sh
  -> verified rows 226545617, signed transfers 63309, minimum integer drift slack 3; PASS upper 5657/9983; ELAPSED=69.81 EXIT=0
grep: x-strogatz/run_all.sh:3 "mpmath not required"; find_upper_lin.py:19 and certs9.py:220 import mpmath  (Finding #3 confirmed)
grep: verify_as_lower.py calls certs9.lower_table / certs6.upper_lin_table, never *.verify  (Finding #2 confirmed)
verdict lines in all 20 cited checks/*.md match the §5.1 table (CORRECT / HOLDS / GAP as stated)
```

Not re-run by me: the Fable `run_all.sh`, the Astra audit, transfer and baseline scripts, the Carbon scripts, and the Baseline `certs9`/`find_upper_lin`.

### Findings

- **Executive Summary:** All numbers are exact and match the sources: the fractions, 17,329,108, 226,545,617, the 1/6000 and 1/299490 margins, and the simulation 0.565723±0.000003. The caveats (AI-only review, μ open) are present. The one overclaim is calling the Carbon result an "independent" bound without noting that it follows from the 0.55 bound in the same sentence (B1).
- **§4 clean-clone claims:** Both runs I repeated (Neon, Astra banded) match the stated outputs, hashes and row counts. Findings #2 and #3 are confirmed in the code.
- **§5 review history:** Every cited file exists and its verdict line matches. "Board attacks" (`attacks/*.md`, `fable-on-carbon-1.md`, `carbon.md`) are not in this repository. They are corroborated only by the leads' RESULT.md files (N5).
- **Other leads' files:** These are not modified. Only `.gitignore` differs, as disclosed.
- **OZET-TR:** It is mostly faithful and plain-language, with the key caveats intact (human review missing, μ open, simulation ≠ proof). One dependency caveat is dropped (B3), and there are small precision slips (N6–N8).

### Required fixes (blocking)

- **B1 — Carbon's value is overstated (FINAL Exec Summary, §1.3 item 4, §2.1 4(c); OZET item 3).** Fable Claim 2 alone, which is independent of Astra and Neon, gives liminf log C_n/log n ≥ 0.55 a.s. That already implies C_n/n^{0.543} → ∞ a.s.
  - Fix: say explicitly that Carbon's stated conclusion is implied by §2.1 4(b).
  - Present Carbon's contribution as the loss-free method and the a.s. convergence of Φ_n^−/n^{0.543} in (0,∞].
  - Replace "Carbon's independent value is proving … without Astra Claims 1–2 or Neon Theorem 1" accordingly.
- **B2 — FINAL §3, Carbon paragraph, misstates the proof.** It says "the normalized submartingale Φ_n^−/P_δ(3,n) is L²-bounded on the planted branch and converges a.s. to a positive limit." Carbon proves L²-boundedness only of the martingale part Z_n in M = X_m + A_n + Z_n (eq. 9–10). Carbon explicitly warns that "the normalized lower sum can have unbounded second moments" (eq. 11 and the text after it). The limit lies in (0,∞], with positive probability on the planted branch. The second-moment bound is K·N_n/n with N_n ≤ Y_n, not "C_1 Φ_n^+/n" in Φ form.
  - Fix: reword to match eqs. (6), (9)–(12).
- **B3 — OZET-TR item 3 drops a dependency caveat.** It says Fable "kanıtlar" (proves) that E[N_k(n)] grows with the same μ. In FINAL §2.1 3(a) this needs Astra Claim 1 (existence of μ) plus a positive and a sub-1 mean bound; Fable D(d) alone is conditional.
  - Fix: add "(Astra'nın varlık teoremiyle birlikte)" ("together with Astra's existence theorem") or equivalent.

### Minor nits (non-blocking)

- **N1:** §3 Astra says that Claim 1 "plants independent fresh edge copies" and that it concludes "by Fekete's lemma". `existence.md` explicitly needs no independence between plantings (that is Claim 2's C/6 construction). Existence also comes from an approximate-supermultiplicativity iteration with η→0, not from Fekete's lemma.
- **N2:** §1.3 item 4's "cannot be upgraded directly" should read "the existing sharper certificates do not satisfy Carbon's hypotheses". This is a check of hypotheses, not an impossibility result. OZET's "0,565'e doğrudan yükseltilemez" ("cannot be raised directly to 0.565") says the same thing more strongly.
- **N3:** In §2.1 3(b) the parenthetical says "= 21.90986… < 21.98847 < 22 with A = 2.74". In fact 21.98847 is the value at A=2.74, and 21.90986 is the value at the exact A=2.73020.
- **N4:** §4.4 is titled "Verbatim" but contains elisions ("...") and an editorial line in the Baseline block. Rename it to "abridged".
- **N5:** §5.1 cites "board attacks" whose files are not in this repository, while the header says every link resolves. Either state that these are external and corroborated by the leads' RESULT.md files, or include them.
- **N6:** OZET says the simulation sits "tam ortasında" ("exactly in the middle") of [0.5650,0.5665]. The midpoint is 0.56575 and the estimate is 0.565723. Use "içinde, ortasına yakın" ("inside, near the middle").
- **N7:** OZET says "İki küçük paketleme bulgusu" ("two small packaging findings") but then lists three items (numpy, the board-summary misstatement, mpmath). FINAL §4.3 has four.
- **N8:** The Exec Summary and the OZET call 17/30 Strogatz's "conjecture". §6 more accurately says he "publicized the question". Consider "proposed value" or "question".

PR REVIEW VERDICT: CHANGES_REQUESTED

## Re-check after revision (58dca869325aed21c611569c0faac8b3413c712e)

This was a fast, focused check of `git diff a8781f2..58dca86 -- FINAL.md OZET-TR.md`, read against the revised files at `58dca86`. Scope and isolation:
- `git diff --name-status 3ac8946 58dca86` shows only additions under `attempts/strogatz-final/`: FINAL.md, OZET-TR.md and the two review files.
- Each pinned lead/base commit (`45bd14f ec15a72 9d18e22 9ad27a2 1cbd617`) differs from `58dca86` only in `.gitignore` (disclosed). **No other leads' files changed.**

- **B1 — resolved.** The changes are in the Exec Summary, §1.3 item 4, §2.1 4(c) and OZET item 3. All of them now say that $\liminf C_n/n^{0.543}>0$ a.s. already follows from Fable's $0.55$ (and from Astra 1–2 + Neon 1). Carbon's distinct content is now described as the loss-free method and the a.s. convergence of $\Phi_n^-/n^{0.543}$ in $(0,\infty]$, with its dependence on both Baseline §3.5 tables stated.
- **B2 — resolved.** §3 now says that only the centred martingale part $Z_n$ in $M_n=X_m+A_n+Z_n$ is $L^2$-bounded, and that $M_n$ can have unbounded second moments with a limit in $(0,\infty]$. The second-moment bound $(H+3)B_F^2N_j/j$ via $N_j\le Y_j$ matches Carbon eqs. (6), (9)–(12).
- **B3 — resolved.** OZET item 3 now says Fable's degree-tail limit holds "Astra'nın varlık teoremi ve $(0,1)$ içindeki üs sınırlarıyla birleştirildiğinde" ("when combined with Astra's existence theorem and exponent bounds in $(0,1)$"), in mean and on the log scale, and $N_k(n)$ is defined.
- **Nits N1–N8:** All eight are addressed:
  - N1: §3 Astra no longer cites Fekete's lemma or independence between plantings.
  - N2: §1.3 now speaks of hypotheses not being satisfied, not an impossibility.
  - N3: the constant is relabelled (21.988467 at A=2.74, 21.909862 at A=2.7302).
  - N4: §4.4 is renamed "Excerpts".
  - N5: the board logs are marked as external.
  - N6: OZET now says "aralığın içinde (orta noktaya yakın)" ("inside the interval, near the midpoint").
  - N7: OZET now says "Üç ... bulgu" ("three findings"), with three items listed.
  - N8: "conjecture" no longer appears in either file; the wording is now "proposed value" / "gündeme getirdiği".
- **New content from the GPT review (checked for consistency only):**
  - The §2.2/§8 fixed-$k$ a.s. statement is logically sound *given* the quoted Addario-Berry et al. Theorem 4.8 ($c_k<\liminf N_{k+1}/N_k$ a.s.). That gives $N_k\asymp C_n$ a.s. for fixed $k$, and with Astra Claim 2 and $\log(n/k)/\log n\to1$ it gives the limit $\mu$. I did not re-read the paper to confirm the theorem's exact wording.
  - The literature endpoints are confirmed: $(\sqrt{13}-3)/6=0.1009252$ and $\sqrt2-1/2=0.9142136$ (Python).
- **No new overclaim found** in the diff. The theorem statements and numbers in §2.1 are unchanged from what I verified in Part 1.

VERDICT: CONSISTENT
PR REVIEW VERDICT: APPROVE
