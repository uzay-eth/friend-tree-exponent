# Cold review of the supporting-results audit and the tree-drift tests

Reviewer: Claude Opus 5.5 (Claude family). The audit was written by GPT-6.2 Astra (GPT family).

Objects:

- [`proofs/supporting-results-audit.md`](../proofs/supporting-results-audit.md), including the new section "The two baseline tables".
- [`verify/tree-drift/`](../verify/tree-drift/) and its section in [`verify/README.md`](../verify/README.md).

Both were read against their definitions: the linked proofs and the verifiers [`certs6.py`](../verify/baseline/certs6.py), [`certs9.py`](../verify/baseline/certs9.py), [`verify_python_cross.py`](../verify/motif-bracket/verify_python_cross.py), [`test_banded_audit.py`](../verify/banded-upper/independent/test_banded_audit.py) and [`check.cpp`](../verify/banded-upper/check.cpp).

**Verdict on the audit: CORRECT, with three required changes.** I found no wrong deduction. The changes concern unshipped artefacts and wording.

**Verdict on the tree-drift tests: CORRECT.** Both tests compute the true drift from the attachment rule. The counts in `verify/README.md` match a fresh run.

## 1. The two baseline tables (new section)

| Claim | Source | Verdict |
|---|---|---|
| Classes: pendant $`i`$ is $`(i, 1)`$; class B is $`(0, 2)`$; the remaining neighbours have degree at least 3 and never return to a recorded class | `certs6.py` docstring and `Model`; a degree-2 core vertex is $`(1,1)`$ or $`(0,2)`$ | CORRECT |
| $`B_0 = b/2 + \sum_i m_i/(i+1)`$ | `rho_bounds`, `m[l]/(l+2)` with 0-based `l` | CORRECT |
| Upper table: $`K_0 = 1`$, $`P = \min(D_{\mathrm{creator}}, 20)`$ if the creator is in the remaining class, otherwise $`P = 3`$ | `nb` uses `max(3, min(D, Pmax))`; the runner passes `K = 1`, `Pmax = 20`. The two expressions agree because a creator in the remaining class has $`D \ge 3`$. | CORRECT |
| Remaining credit in $`[0, o/3]`$, improved to $`(o-1)/3 + 1/P`$ when $`P > 3`$ | `rho_bounds`; the creator stays a remaining neighbour and its degree never decreases | CORRECT |
| Lower table: $`K_0 = 2`$, exact $`n_3, n_4, n_5`$, other neighbours in $`[0, 1/6]`$ | `certs9.py` with `K0 = 2`, `Dcap = 6` (runner arguments `60 10 2 6`) | CORRECT |
| Event rates: pendant gain $`i + 1/D`$, conversion $`i/(i+1)`$, class-B gain $`1/D + 1/D_x`$, remaining gain in $`[0, d]`$ | `certs6.py` E3–E5, `certs9.py` E7 | CORRECT |
| Leaf count frozen at 60; for $`c > 10`$ the lower table returns zero and the upper table a positive constant $`G`$ | `certs9.py` 119–121, `certs6.py` 111–113, 197–199, 276–278 | CORRECT |
| Tail handled by the endpoints of the intervals for $`L/D`$ and $`1/D`$ | `certs9.py` 193–196 | CORRECT |

The section replaces the earlier link to an unpublished source. It now states every event the two tables use, and each statement is checked by a shipped verifier.

## 2. The other sections

| Section | What I checked | Verdict |
|---|---|---|
| Hub 1–6 | rates and ownership against [`bounds-hub.md`](../proofs/bounds-hub.md). Tail: the upper drift of $`A(1+\lambda q)`$ is at most $`A\lambda(1-q) + qs \le \lambda A(1 + \lambda q)`$, which holds exactly when $`s \le A\lambda(1+\lambda)`$. The lower drift of $`Bq`$ is at least $`q(s - B/D) \ge \delta Bq`$ when $`s \ge B(\delta + 1/H)`$. Prefactors $`15493968692/17611447440`$ and $`15349885606/810439425`$ match the proof and the log. | CORRECT |
| Transfer 1–5 | witness rules and the cancellation against [`bounds-transfer.md`](../proofs/bounds-transfer.md); 6,445 pairs; $`H = 208`$; $`30 \cdot 8723588643 - 17 \cdot 15361771189 = 557549077`$ (recomputed) | CORRECT |
| Almost-sure lower 1–5 | Lemma 1: $`n \ge \max(3, M^2+M)`$ and probability at least $`1/(2n)`$. Lemma 2: the bounded-increment dichotomy. Lemma 3: $`\kappa = 2B\delta + 4B^2R/f_{\min}`$ follows from (1) and the union bound. Theorem B: $`K_0 \ge 2B`$ keeps the stopped denominator positive. Theorem A: $`R = 60 + 5 + 1 + 60 = 126`$ and 7 newborn classes. All against [`almost-sure-lower.md`](../proofs/almost-sure-lower.md). | CORRECT |
| Cross-charging 1–6 | summary of [`cross-charging.md`](../proofs/cross-charging.md); the counts 5,853 and 5,851 and the root-row failures at $`579/1000`$ and $`552/1000 = 69/125`$ match the runner log | CORRECT |
| Degree tail 1–8 | Lemma 2 minimum $`4/120 = 1/30`$ at $`l = 3`$. Lemma 3: $`2E \ge (C+1)/2 + 2`$, so $`E \ge (C+5)/4`$. Lemma 4: $`(m-1)/(8m-1) \ge 2/23`$ and $`m(m-2)/((2m-2)(2m-3)) \ge 1/4 \iff m \ge 3`$. D(b): $`4 \cdot 368 \cdot 30 = 44160`$. D(d): $`0.5022 \cdot 0.75^{0.45} \approx 0.4412 > 0.44`$ and $`3.974 \cdot (2 + 1/0.42) \cdot 3^{0.58} \approx 32.93 < 33`$, which is tight but true. | CORRECT |
| Polynomial-scale 1–10 | spot checks only: $`\sum_u \rho_u = \sum_w c_w/D_w \le C_n`$, the bound $`N_n + 1 \le 2N_n`$ in (5), and the hypothesis on $`H`$. I did not re-derive items 4–9. They agree with the statements of [`polynomial-scale-lower.md`](../proofs/polynomial-scale-lower.md). | CONSISTENT (not re-derived) |
| Fresh verification | each quoted output line occurs in the runner log, and the hub state counts are 3,550,712 and 5,960,583 | CORRECT |

All relative links in the audit resolve to files in the release tree.

## 3. Tree-drift tests

- **Motif.** `check_tree` computes $`\sum_W \big(\sum_{V \sim W} 1/\deg V\big)(\Phi(T + W) - \Phi(T))`$, which is $`n`$ times the true expected increment. It evaluates $`\Phi`$ from the true rooted states and does not use the decomposition of the proof. The generator lists every non-isomorphic rooted tree whose root has degree at least 2. Since $`\Phi`$ depends only on the rooted shape, this covers every pair $`(T, r)`$ of that size. The run covers 32,972 trees for $`n = 3, \dots, 14`$.
- **Banded.** The test starts from $`T_3`$ and branches over every attachment target. That gives $`(n-1)!/2`$ histories of size $`n`$, and $`1 + 3 + 12 + 60 + 360 + 2520 + 20160 = 23{,}116`$ through $`n = 9`$, as reported. The drift is $`\sum_V (1/\deg V) \sum_{W \sim V} (\Phi(\text{step}(W)) - \Phi)`$, and the test requires $`9983 \cdot \text{drift} \le 5657\,\Phi`$. The label updates in `BandedHistoryTree.step` follow the definitions: marks are frozen once pooled, all but the largest are reset, witnesses are kept, renewed or dropped by the stated rules, and a newborn's witness is the minimum label. Weights come from the production checker through a compile-time observer, after every finite row has passed.
- **Scope.** Histories through $`n = 9`$ reach only small degrees, so they test witnesses and exact labels but not the pool marks above $`h_c`$, the boundary or the tail. `verify/README.md` describes the test without claiming more.

## Required changes

1. **Unshipped artefacts in "Fresh verification".** The quoted commands use `build/support-audit/x/upper.pkl` and `build/support-audit/x/lower.pkl`. These paths are local and ignored by `.gitignore`. Use the runner's own invocations, `xverify.py "$TMP/xc_upper.pkl" 29/50` and `xverify.py "$TMP/xc_lower.pkl" 11/20` (`verify/run_all.sh` 151–152), or say that the files are local. The closing sentence ("I also independently enumerated the 6,445 eligible endpoint pairs…") cites scripts that are not in the release. Either ship them or say they are not part of the release.
2. **First person.** Change "I also independently enumerated…" to "The audit also enumerated…", in line with the other public files.
3. **Optional precision.** In the upper-table item, write $`P = \max(3, \min(D_{\mathrm{creator}}, 20))`$ as in `certs6.py`, or add that $`D_{\mathrm{creator}} \ge 3`$ in that case.

## Evidence

```text
$ bash verify/tree-drift/run.sh
[PY-EXACT 1] Theorem 1 exact rational bounds bracket exact E[C_n] for ALL 19 sizes n=3..21: PASS
[PY-EXACT 2] Pure-Python Fraction brute-force drift passed on ALL 32972 non-isomorphic rooted trees with core root (n=3..14), failures=0 for BOTH UPPER & LOWER: PASS
verified rows 226545617, signed transfers 63309, minimum integer drift slack 3
PASS upper 5657/9983
PASS banded: 23116 attachment histories through n=9, 304 distinct weight queries, failures=0
exit status 0
```
