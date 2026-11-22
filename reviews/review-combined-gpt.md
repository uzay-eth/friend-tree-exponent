> Packaging note: internal branch name on line 5 was redacted for public release; all mathematical content and commit hashes are unchanged.

# Cold GPT-family referee report

Reviewed the combined report at commit `a8781f2dfd71b551364c145f4f6f357c1d0054e0`. This review checks the interfaces between the pinned results and the report's claims. The underlying drift certificates were taken as inputs.

## Part 1 — Cross-work consistency

P1 — All five use a single edge $`T_2`$, size index $`n`$, and the same core count $`C_2=0`$, $`C_3=1`$. Mean product/Gamma bounds start at $`n=3`$. Fable has $`N_2(n)=C_n`$. In §1.2 add $`3\le m_0\le n`$ in D(a) and $`m'\ge3`$ in D(b), and label Carbon's conclusion asymptotic.

P2 — D(d) composes with Astra Claim 1 and either Neon Theorem 1 or Fable Claim 2. Both exponent brackets, $`[0.5650,0.5665]`$ and $`[0.55,0.58]`$, lie inside $`(0,1)`$. Say “positive lower exponent and upper exponent below one.”

P3 — Both inequalities hold on their entire stated domains, including $`k=2`$ and $`n=24k`$. Gautschi gives $`0.411666741459\ldots`$ and $`2.730204939524\ldots`$. D(c) gives $`21.988467893197\ldots<22`$ with $`A=2.74`$, and $`21.909862647653\ldots`$ with the sharper Gamma coefficient. Correct §2.1(3b).

P4 — Baseline `certs9` supplies Fable Theorem A via Theorem B at $`543/1000`$. Fable's cross-charged lower table supplies $`11/20`$ a.s. via Theorem B, and $`29/50`$ bounds the mean exponent above. Both a.s. lower bounds pair with Neon Corollary 2 independently of Astra Claims 1–2.

P5 — The product identity, every §1.3 prefactor at its printed precision, and all nine signed margins check out. The two upper margins are $`1/6000`$ and $`1/299490`$, and the width is $`3/2000`$.

P6 — §2.1(2c) misattributes the Neon interval to Astra-only inputs, and §2.1(1a) overstates equivalence and finite-size sharpness. §2.2 overlooks the fixed-$`k`$ almost-sure limit from Astra Claim 2 plus Addario-Berry et al., Theorem 4.8, with $`k=2`$ immediate.

VERDICT: INCONSISTENT

## Part 2 — PR review

### Required fixes

1. **Executive Summary and combined statement.** The summary has exactly five sentences, and the rational constants and primary row counts agree with the cited certificates. The three middle sentences are dense. Attribute existence and sample-path convergence to Astra, the tighter exponent bracket to Astra Claim 1 plus Neon, the uniform tail to Fable D(d) with the existence and interior-bracket inputs, and Carbon's polynomial-scale result to Carbon plus both Baseline tables. Say that two independently derived upper *exponents* fall below $`17/30`$, and that $`1/6000`$ is the certified minimum gap. Give Astra's range $`n\ge3`$ explicitly. Limit the reproducibility claim to the certificate commands and checks recorded in §4. Analytic proofs were reviewed by reading, and Carbon's local replay and the full Baseline table replays are separate checks. Report the simulation as a finite-size effective-slope estimate with a sampling standard error.

2. **§§1.2–2.2, scope and dependencies.** Apply P1–P3 and P6. §2.1(1a) should say the exact products *imply* the power bounds, and “tightest” should describe the exponent bracket. At $`n=3`$ Fable's lower product coefficient $`0.882192\ldots`$ already exceeds Neon's $`0.833455\ldots`$. Put Astra's independent, weaker interval on its own line in §2.1(2c). In §2.2 and §8(4), specify the open sample-path question as a logarithmic limit uniform over growing $`k`$, or as convergence of normalized counts. For each fixed $`k\ge2`$, Theorem 4.8 gives $`N_k(n)\asymp C_n`$ almost surely, and Astra Claim 2 then gives $`\log N_k(n)/\log(n/k)\to\mu`$ almost surely. A countable intersection gives the fixed-threshold statements simultaneously. [Addario-Berry et al., Theorem 4.8](https://arxiv.org/html/2403.20185v1#S4)

3. **§3, certificate methods.** Astra's banded weight above the cutoff is $`A+B L/D`$, and its labels store observed degrees, possibly capped core degrees, pooled marks and a witness. The signed transfers add and subtract the same gain-rate term across each edge. The second cancellation is the transfer of $`a(D_u,c_u)/D_w`$ on a core edge with $`D_w\ge T`$. The trial newborn reward and its reserve are a separate check. For Neon, use $`f(L,k,h)`$, exact endpoint states below $`D_0=24`$, and multilinear $`D_0/D`$ interpolation for $`g,h_V,h_C`$ at tail endpoints. The centering interpolants have **nine** knots in the pinned header, while $`K_T=12`$ is the tail child-count parameter. The counts 84 and 75 enumerate active middle and apex states. A step affects motifs within a bounded graph radius, with the number of affected motifs potentially growing with degree. Fable D(a) uses a conditional first-moment recursion and Markov's inequality, followed by summation by parts in D(c). D(b) activates eligible vertices already in the core at $`m'`$ and then uses conditional moments and Paley–Zygmund.

4. **§3, almost-sure methods.** Astra Claim 1 plants rooted $`k`$-vertex seeds at eligible sites, uses marginal first moments, and iterates the resulting inequality. Claim 2 uses conditionally independent lower trials within a block and a martingale-difference law across blocks. Fable's transfer needs all (H1)–(H5), and its stopped reciprocal is $`1/(\Phi_{n\wedge\tau}p_{n\wedge\tau})`$ with $`p_m=\prod_{j=\sigma}^{m-1}(1-(\delta-\varepsilon)/j)`$. Carbon controls the $`L^2`$ norm of the centered martingale part. The normalized lower sum has a nondecreasing predictable part and may have an infinite limit. Carbon gets probability one from a uniform lower bound on the conditional success probability at every starting size, using approximation in the full history sigma-field. The trials can share history.

5. **§5, audit provenance and limits.** The Fable sequence is recorded correctly. The v1 `GAP` concerned $`k=2`$ and two prose overstatements. Commit `b432127` repairs them, the GPT referee records `CORRECT` after re-reading, and `b1e74f0` incorporates Gemini's additional $`\varepsilon<\mu`$ edit. Astra's `tag.cpp` `GAP` concerns the unverified *strict* upper bound. The partial referee report says a positive reserve had not been produced and the referee did not replay the full $`17/30`$ table. Replace “flawed alternate construction” with “alternate construction whose strict bound remained unverified.” Tie each purported independent full pass to a source and output or narrow the count. The Opus Astra term-by-term replay extends `check.cpp`. The linked independent full reconstruction is Gemini's `verify_banded.cpp`. For Neon, GPT's two arithmetic versions share row formulas, and the cold Claude and Fable-delegated Opus referees each report a separate full replay. The extra Carbon full-replay claims in the two counts need linked evidence. Keep the human-referee, formalization, high-degree tree coverage, and Carbon/Baseline scope limits. Describe the cross-family coverage at the level of the underlying claims, since the combined statements have their own review. `Fable GPT re-check` · `Astra tagged audit` · `Astra Opus audit` · `Neon GPT audit`

6. **§6, literature.** Correct the titles of arXiv:2305.10628 (*The Magic of Networks Grown by Redirection*) and arXiv:2604.01540 (*Anomalous scaling in redirection networks*). Hartle et al. cite the earlier $`0.566`$ estimate and obtain about $`0.569`$ from their own finite-size IR fit in Fig. 6. Their analytical exponents concern modified models. Addario-Berry et al., Propositions 6.20 and 6.22, give almost-sure statements for lower powers strictly below $`(\sqrt{13}-3)/6`$ and upper powers strictly above $`\sqrt{2}-1/2`$, respectively. State the resulting log-liminf and log-limsup bounds with exact endpoints. The printed $`0.9142`$ is below the exact endpoint, so the cited proposition does not give $`C_n=O(n^{0.9142})`$ a.s. Explain the discrepancy with the arXiv theorem's displayed $`0.9`$ ceiling. Phrase the literature gate as a dated search result and cite the Strogatz post through the primary link and the Baseline's preserved record. This review could not retrieve the post directly. [Redirection review](https://arxiv.org/abs/2305.10628) · [Hartle et al., Fig. 6](https://arxiv.org/html/2604.01540v1) · [Random friend trees, §6.5](https://arxiv.org/html/2403.20185v1#S6.SS5) · [Journal record](https://doi.org/10.1002/rsa.70082)

7. **§8, open problems.** Apply the fixed-$`k`$ correction above. Present a CMJ, multitype, or operator representation as a research question. The attachment rate is determined locally from leaf count and neighbour degrees, while descendant dynamics can remain coupled across edges. Drop the blanket claim that finite-motif certificates cannot identify an exponent, since a matching upper and lower certificate would do so.

8. **Turkish summary.** The interval, margins, Fable repair, and exclusion of the tagged route are correctly noted. In line 11 specify “eşit olasılıkla” for both choices. Use “üst sınırın $`5657/9983`$ olduğunu” or write $`\mu\le5657/9983`$ explicitly in line 18. In line 22 write “Fable’ın sonucu Astra’nın varlık teoremi ve üs sınırlarıyla birleştirildiğinde, $`n/k\to\infty`$ iken bu limit $`k\ge2`$ için tekdüze olarak $`\mu`$’ye gider.” Define $`N_k(n)`$ and keep the mean and logarithmic qualifications. In line 23 say “$`n^{0{,}543}`$ ölçeğinde pozitif bir alt sınır” and retain the Baseline dependency. In line 30 replace “tam ortasında” with “aralığın içinde” (midpoint $`0.56575`$), identify Fable's simulation error as sampling error, and give Neon's $`1.6\times10^{12}`$ as a mean size at fixed time. In line 31 carry over the open normalization and growing-$`k`$ caveats. In lines 37–40, carry over the limits of the replays, absence of formal proof verification, the scope of the tagged `GAP`, and the qualified literature result.

### Minor nits

- $`T_3`$ has deterministic shape and a unique centre. Its label can be random, and rooting requires no extra random choice. The Baseline's $`0.50`$ remains a valid weakening of $`0.502`$.
- Label the §1.3 prefactor decimals as approximations, and identify the table as covering the selected certificates. §4.4 includes abridged logs, so “excerpts” would be more accurate than “verbatim outputs.”
- §8(1) can say “an order of magnitude narrower” without promising four or five determined decimal digits. In Turkish, “derecesi en az $`k`$ olan düğümler” is clearer here than “derece-kuyruğu.”

### Numerical evidence

`python3` with `Fraction` and `math.gamma` recomputed the table. The Neon and Astra inputs were the cited exact ratios, the Fable inputs came from the hexadecimal certificate weights, and the Baseline inputs came from the full-precision values printed in its `outputs/run_all.log`. The columns below are the computed $`c_{\rm prod}`$, $`2c_{\rm prod}/\Gamma(3+a)`$, and $`17/30-a`$.

```text
Neon lower       0.833455596175    0.466545636804     1/600
Neon upper       4.885581076574    2.730204939524     1/6000
Astra 4 author 393.944781629994  220.107368872299     1/299490
Astra 4 audit  393.944791457626  220.107374363257     1/299490
Astra 3 lower    0.876152855091    0.491550237250    11/3000
Astra 3 upper   23.995087667543   13.386542041760    -1/750
Fable lower      0.882192339679    0.502208383199     1/60
Fable upper      7.219817464634    3.973720017481    -1/75
Baseline lower   0.875518155482    0.502326215141    71/3000
Baseline upper   3.954623733153    2.161877845316   -29/1500
```

For the universal bounds, apply Gautschi at $`n=512`$ to enclose the Gamma coefficients using an exact finite product. For arbitrary $`m\ge3`$ the lower Gamma ratio contributes at least $`(3/4)^{1-\delta}m^\delta`$. D(c) uses $`m_0=\lceil2n/k\rceil\in[3,n]`$ and $`m'=\lfloor n/(8k)\rfloor\ge3`$, with $`m'\ge n/(8k)-1\ge n/(12k)`$ for $`n\ge24k`$. At $`k=2`$ the upper increment sum is empty. The lower coefficient is $`41/4416000`$. I ran the following exact comparisons.

```sh
python3 - <<'PY'
from fractions import Fraction as Q
from math import prod

def enclosed_coefficient(c, a, n, scale=10**12):
    exact_product = c * prod(1 + a/j for j in range(3,n))
    scaled = exact_product * scale
    floor = scaled.numerator // scaled.denominator
    return Q(floor,scale), Q(floor+1,scale)

n = 512
delta = Q(113,200)
lam = Q(1133,2000)
lo, _ = enclosed_coefficient(Q(7304352308,8763936965),delta,n)
_, hi = enclosed_coefficient(Q(14606753031,2989767809),lam,n)
lower_test = (lo / Q('0.41'))**200 * Q(3,4)**87 > n**113
upper_test = (hi / Q('2.74'))**2000 * Q(n+1,n)**867 < n**1133
tail_test = (Q('2.74') * (2 + 1/(1-lam)) / 22)**2000 * 3**1133 < 1
assert lower_test and upper_test and tail_test
print(f'Gautschi rational enclosure anchor n={n}; all m>=3: lower={lower_test}, upper={upper_test}')
print(f'D(c) coefficient <22, exact Fraction test={tail_test}')
PY
```

```text
Gautschi rational enclosure anchor n=512; all m>=3: lower=True, upper=True
D(c) coefficient <22, exact Fraction test=True
```

PR REVIEW VERDICT: CHANGES_REQUESTED

## Re-check after revision (58dca869325aed21c611569c0faac8b3413c712e)

- **P1–P6 resolved in the combined statement.** P1 records the full Fable domains and Carbon's asymptotic scope. P2 states the interior-exponent requirement. P3 attributes each upper coefficient to the correct input. P4 and P5 still match the pinned sources. P6 separates the Astra-only bracket from the tighter Neon bracket, makes the product-to-power implication one-way, and recovers the fixed-$`k`$ almost-sure limit in §2.2.

VERDICT: CONSISTENT

- **Required fixes 1–4.** The five-sentence summary now gives the dependencies, exact margins, replay limits, and sampling caveat. Several sentences are still dense. The combined theorem and the Astra, Neon, Fable, and Carbon method corrections are in place.
- **Required fix 5.** The two GAP histories and replay scopes are now distinguished correctly. Qualify the new claim in §§5.1–5.2 that every external attack is corroborated by the leads' `RESULT.md` files. For example, the pinned reports do not corroborate the new Carbon full-replay entries. Attribute those entries to the external board and give a reference readers can check, or omit them from the corroborated count.
- **Required fixes 6–7 still need two mathematical edits.** In `§6.1`, Theorem 4.8 gives a positive lower bound on the liminf of $`N_{k+1}/N_k`$. The upper bound is $`\le1`$ by nesting. Replace the asserted strict $`\limsup<1`$. The fixed-$`k`$ consequence uses only the lower bound. In `§8(3)`, the unqualified sum in $`L_u+\sum_{v\sim u}1/D_v`$ counts leaves twice. Use the scaled rate $`L_u+\sum_{v\sim u,\,D_v\ge2}1/D_v`$ from §1.1. The titles, published exponent endpoints, and open-problem framing are otherwise corrected. [Addario-Berry et al., Theorem 4.8](https://arxiv.org/pdf/2403.20185#page=8)
- **Required fix 8.** The Turkish dependency, limit, simulation, and audit caveats are substantially corrected. Remaining wording edits are to state both choices are uniform (“eşit olasılıkla”), explicitly apply “neredeyse kesin” to Carbon's liminf, and specify that the tagged GAP concerned a bound strictly below $`17/30`$.

`git diff --check a8781f2..58dca869325aed21c611569c0faac8b3413c712e -- attempts/strogatz-final/FINAL.md attempts/strogatz-final/OZET-TR.md` exited 0 with no output. The Python `Fraction` check above again returned `lower=True, upper=True` for Gautschi and `True` for the D(c) coefficient.

PR REVIEW VERDICT: CHANGES_REQUESTED

## Final re-check (4d5e7995505279d6cdc5f638afb010cfa0f80aa5)

The focused diff resolves all four outstanding items. §6.1 now states the positive liminf bound from Theorem 4.8 and the non-strict upper bound from nesting. §8(3) uses the scaled rate with a sum over core neighbours. §§5.1–5.2 distinguish attacks corroborated in the pinned reports from later Carbon entries reported on the external board. The Turkish summary now specifies both uniform choices, Carbon's almost-sure liminf, and the scope of the tagged GAP.

`git diff --check 58dca86..4d5e7995505279d6cdc5f638afb010cfa0f80aa5 -- attempts/strogatz-final/FINAL.md attempts/strogatz-final/OZET-TR.md` exited 0 with no output. This re-check covered the revision diff. The external board logs were outside its scope.

VERDICT: CONSISTENT

PR REVIEW VERDICT: APPROVE
