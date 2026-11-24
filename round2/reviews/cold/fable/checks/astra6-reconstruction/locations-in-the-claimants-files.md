# Statement locations (comparison only)

All paths in the table are relative to `ASTRA_REPO`, except the board at `BOARD/astra.md`. Brackets give section, equation number where present, and source lines. The five principal board entries are claims 2, 3, 5, 6, 7, lines 14–22.

| Statement | Source location |
|---|---|
| Model, notation and exponent inputs | `proofs/bracket-motifs.md` [§1, Theorem 1 and remark, lines 5–36]; `proofs/existence.md` [Theorem 1 and §7–8]; `proofs/almost-sure-exponent.md` [Theorem 2, lines 9–13]; `round2/operator.md` [§2, $`T_3`$ time exponent, line 41]. |
| S2.1 | `round2/fringe/renewal.md` [§1, lines 9–15; §2, (1), lines 17–23]. |
| S2.2 | `round2/fringe/renewal.md` [§3, weight and (4)–(5), lines 35–55]. |
| S2.3 | `round2/fringe/renewal.md` [§2, (2)–(3), lines 25–33; §3, (6) and score class, lines 55–59]. |
| S2.4 | `round2/fringe/renewal.md` [Proposed result, lines 3–7; §4, (7), lines 61–81; §5, lines 83–103]; board claim 2. |
| S5.1 | `round2/fringe/finite-kernel.md` [Construction, (1), lines 3–9; §1, lines 11–17]; board claim 5. |
| S5.2 | `round2/fringe/finite-kernel.md` [§2, (2), lines 19–36; §3, (3)–(4), lines 38–55; §4, (5)–(7) and last unnumbered slope estimate, lines 57–94]. |
| S3.1 | `round2/fringe/asymptotics.md` [Statement, (1)–(2), lines 3–15; §1, definitions, lines 17–25; §2, (4), lines 35–50; §3, (5) and fixed-size constant, lines 52–76; §4, line 80]; board claim 3. |
| S3.2 | `round2/operator.md` [Statement, (1), lines 3–11; §1, (2)–(3), lines 13–37; §2, tilts, (4)–(5), lines 39–55; §4, line 75; §5, (8), lines 79–88]; `round2/fringe/asymptotics.md` [§4, line 80]. |
| S3.3 | `round2/degree-law.md` [Statement, (1), lines 3–9; Birth times, (2), lines 11–21; Tail constant, lines 23–29; Point probabilities, crossing identity, lines 31–37]. |
| S6 (theorem) | `round2/fringe/almost-sure.md` [Statement, (1), lines 3–11; §4, lines 115–117; §5, lines 119–134]; board claim 6. |
| S6.1 | `round2/fringe/almost-sure.md` [§1, kernel and $`RB_0g`$ bound, integral and decay, pole and residue, (2), lines 13–30]. |
| S6.2 | `round2/fringe/almost-sure.md` [§1, newborn renewal expansion and (3), lines 32–41]. |
| S6.3 | `round2/fringe/almost-sure.md` [§2, initial conditions, weight/immigration bounds, (4) and degree-count assertion, lines 43–57]; `round2/fringe/renewal.md` [§2, (2)–(3); §3, (6)]. |
| S6.4 | `round2/fringe/almost-sure.md` [§2, two-sided comparisons and (5), lines 59–64]. |
| S6.5 | `round2/fringe/almost-sure.md` [§3, (6) and expected combined dirty count, lines 68–73]. |
| S6.6 | `round2/fringe/almost-sure.md` [§3, (7), one-unit estimates and (8), lines 75–87]. |
| S6.7 | `round2/fringe/almost-sure.md` [§3, integer-time recurrence, (9), conditional probability lower bound and positivity, lines 89–106]. |
| S6.8 | `round2/fringe/almost-sure.md` [§4, conditional means and centered error, all real times, deterministic ratio, lines 108–117]. |
| S6.9 | `round2/fringe/almost-sure.md` [§5, discrete mean, holding-time transform, stopped second moment and conditional mean of $`X`$, lines 119–134]. |
| S7 (theorem) | `round2/hub-mass.md` [Statement, (1), lines 3–10; localized variant, lines 30–41]; board claim 7. |
| S7.1 | `round2/hub-mass.md` [A bound in expectation, (2), degree sum and attachment formula, lines 12–28]; fixed-label limits attributed there to Addario-Berry et al., Proposition 3.1. |
| S7.2 | `round2/hub-mass.md` [Using only an almost-sure power upper bound, stopped sum and Fatou estimate, lines 30–41]. |

The specified S3 sources (`asymptotics.md`, `degree-law.md`, `operator.md`) and S7 source contain no separately labelled lemma/proposition blocks. S3.1–3.3 record the pertinent claims in those files; the S7.1–2 labels and S6.x labels are extraction labels. The citations in the proof maps name any external results used as inputs.

## Ambiguities and notation collisions

- **Admissible types (S2.1/S6.2):** renewal §1 describes physically assigned marks and finite auxiliary types, including initial mark zero, without separately enumerating which compatible but unreachable internal mark patterns belong to the score domain. The description uses the stated reciprocal inequality to describe admissibility; no assertion of reachability is made. Almost-sure §2 explicitly quantifies over any compatible capped marks on finite starting physical trees.
- **Uniform constants (S6):** $`C, C_0, \ldots, C_{12}`$, the kernel's $`O`$-constant, and the choice of $`\sigma\in(\rho, \mu)`$ are not given numeric values. The absolute $`C`$ in the generator bound and generic $`C`$ in kernel/semigroup bounds need not denote the same number. The statements retain each source quantifier: e.g. kernel decay may depend on $`\varepsilon`$; (4) is uniform in $`T`$ and $`f`$. Constants absorbed by the founder's contribution and in the conditional mean estimate are not named in source. Simultaneity in $`k`$ is a countable almost-sure assertion; the uniform score bound applies to all threshold scores as well.
- **Names:** exact-mark $`d_f`$ includes founder immigration, while capped $`d(f)`$ does not. Both sources use $`g, \ell`$ for (different) eigenvectors; neither fixes their separate scales, but the paired coefficients and limits are invariant. $`D^{\rm ker}`$, $`W^{\rm deg}`$, $`H_*`$, $`e_*`$, $`Z^{\rm tilt}`$ and $`B^{\rm sp}`$ distinguish the sources' $`D, W, H, e, Z, B`$ from other objects named the same way. The source time argument and discrete size both sometimes use subscripts on $`T`$; attached qualifiers distinguish them.
- **Core versus all neighbours (S3.3):** `degree-law.md` uses $`L+\sum_{w\sim o}1/D_w`$ and $`\rho_o\le c_o/2`$ in a leaf-marked core. Accordingly the sum is over core neighbours; physical leaves already contribute $`L`$. The un-subscripted $`\rho`$ in its line 27 is that core-neighbour gain contribution, not the residual growth parameter.
- **The law $`\pi`$ (S3.2):** the operator file constructs a weak subsequential limit, not a uniquely specified law on all spatial scores. S3 identifies the degree marginal of the constructed limit. The existence of moments for the limiting law is read from the bounds (5) and weak convergence, not a separate labelled proposition.
- **Truncation upper matrices (S5.2):** §3 claims both nested matrix bounds; its supersolution explanation tracks a scalar reward bound per label. The defined upper matrix includes *all* old exit rewards from a state, even when a larger cutoff can accumulate direct reward there before exiting again. Comparison requires the supersolution against direct reward plus remaining boundary reward, not merely bounding the remaining mass. See exact formula in (3).
- **Power upper bound in S7:** the preprint exponent and big-$`O`$ constant are not specified in `hub-mass.md`; the localized assertion only needs some deterministic $`\lambda<1`$ and a pathwise finite constant. Under that hypothesis one necessarily has $`\lambda\ge0`$ since $`C_n\ge1`$. This is the relevant range in the source bound $`C_{r\wedge\tau}\le K r^\lambda+1`$; for a formally negative $`\lambda`$ the hypothesis is empty. Localization counts promotions after $`m`$ regardless of vertex birth label, an upper bound for the late-label mass.
- **One-ring terminology:** a core promotion changes the existing leaf into a core and adds a different new leaf to it. Thus a core's promotion time and its original birth label are not the same; for the initial edge, the first ring creates the founder and the ordinary promotion description applies from size three.
