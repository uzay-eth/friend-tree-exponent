CORRECT

# Cold review of `round2/hub-mass.md`

Reviewer: Gemini 4 Carbon  
Target: `../hub-mass.md` (`hub-mass.md`) (referencing `../../proofs/bracket-motifs.md` (`proofs/bracket-motifs.md`) and Addario-Berry et al.)  
Checker: `hub_mass_carbon.py` (`python3 -B round2/checks/hub_mass_carbon.py`)

## Verdict summary

The deduction in `../hub-mass.md` (`hub-mass.md`) is **CORRECT**. Both routes—the unconditional expectation route using $`a_n = \mathbb E C_n \le K n^\lambda`$ ($`\lambda < 1`$) and the stopped-process localization using only $`C_n = O(n^\lambda)`$ almost surely ($`\lambda < 1`$)—rigorously prove $`\sum_{v \ge 1} Z_v = 1`$ almost surely, almost-sure total-variation convergence of $`\mathbb P(W_n = \cdot \mid T_n)`$ to $`(Z_v)_{v \ge 1}`$, and almost-sure weak convergence of $`\mathcal L(D_{W_n}(n)/n \mid T_n)`$ to $`\sum_{v \ge 1} Z_v \delta_{Z_v}`$.

---

## Audit of each step

### 1. Preprint and journal hypotheses and statements

- **Model and fixed-vertex limits.** In both the preprint (`arXiv:2403.20185`) and the published journal article (*Random Structures & Algorithms* 69(1), e70082, 2026), $`T_2`$ is an edge on $`\{1, 2\}`$ and step $`n \to n + 1`$ attaches vertex $`n + 1`$ to a uniform neighbour $`W_n`$ of a uniform vertex $`V_n \in \{1, \dots, n\}`$. Almost-sure convergence $`D_v(n)/n \to Z_v`$ and $`L_v(n)/n \to Z_v`$ for every fixed $`v \ge 1`$ is Theorem 3.1 in the preprint (Proposition 3.1 in the journal version), and $`\sum_{v \ge 1} Z_v \le 1`$ a.s. is observed immediately below it in Section 3 (from $`\sum_{v \le k} L_v(n) \le n - 1`$).
- **Exact match to the two Section 3 conjectures.** Section 3 of Addario-Berry et al. conjectures (i) $`\sum_{i \ge 1} Z_i = 1`$ almost surely (also Section 7, Question 4), and (ii) almost-sure Prokhorov (weak) convergence of the conditional law of $`D_n(W_n)/n`$ given $`T_n`$ to $`\sum_{i \ge 1} Z_i \delta_{Z_i}`$ on $`[0, 1]`$. Both conjectures—together with total-variation convergence of $`W_n`$ given $`T_n`$—are established here.
- **Upper-bound inputs.** Section 1 uses $`a_n \le K n^\lambda`$ for $`n \ge 3`$ with $`\lambda < 1`$, supplied in the repository by `../../proofs/bracket-motifs.md` (`proofs/bracket-motifs.md`) ($`\lambda = 1133/2000`$). Section 2 uses only $`C_n = O(n^\lambda)`$ a.s. for deterministic $`\lambda < 1`$, which is already proved in Addario-Berry et al. itself (Proposition 6.22 in both versions gives $`n^{-\lambda} C_n \to 0`$ a.s. for every $`\lambda \in (\sqrt{2} - 1/2, 1)`$; Theorem 4.7 in the preprint / Theorem 4.5 in the journal states the resulting polynomial upper bound). Note that $`C_n \ge 1`$ for $`n \ge 3`$ forces $`\lambda > 0`$ (indeed $`\lambda > 0.1`$).

### 2. Supermartingale from promotion onward and Equation (2)

- **Supermartingale property.** At size $`j \ge 3`$, $`\mathbb P(W_j = v \mid \mathcal F_j) = j^{-1}\sum_{u \sim v} D_u(j)^{-1} \le D_v(j)/j`$ since every neighbour has degree $`\ge 1`$, giving $`\mathbb E[D_v(j+1)/(j+1) \mid \mathcal F_j] \le D_v(j)/j`$.
- **Promotion indexing and Equation (2).** For $`n \ge m \ge 3`$, any vertex with label $`v > m`$ arrives as a leaf at size $`v \ge m + 1`$ ($`D_v(v) = 1`$). If $`D_v(n) \ge 2`$, it first enters the core at a unique size $`r \in \{v + 1, \dots, n\} \subseteq \{m + 1, \dots, n\}`$ via step $`T_{r-1} \to T_r`$ ($`r \ge 5`$). At size $`r`$, $`T_{r-1}`$ had size $`r - 1 \ge 3`$, so attaching leaf $`r`$ to leaf $`v = W_{r-1}`$ promotes **only** $`v`$, giving $`C_r - C_{r-1} = 1`$ and $`D_v(r) = 2`$. Conditioning on $`\mathcal F_r`$ bounds $`\mathbb E[D_v(n)/n \mid \mathcal F_r] \le 2/r`$; summing over $`r = m + 1, \dots, n`$ (including promotions of earlier-born leaves $`u \le m`$, which only increases the bound) gives the first inequality in (2). Setting $`r = j + 1`$ ($`j \ge m \ge 3`$), $`\mathbb P(C_{j+1} - C_j = 1 \mid T_j) = j^{-1}\sum_{u \in K(T_j)} L_u(j)/D_u(j) \le C_j/j`$ yields $`a_{j+1} - a_j \le a_j / j \le K j^{\lambda - 1}`$, giving the second inequality. Bounding the $`j = m`$ term by $`m^{\lambda - 2}`$ and $`\sum_{j \ge m+1} j^{\lambda - 2}`$ by $`\int_m^\infty x^{\lambda - 2}\, dx = m^{\lambda - 1}/(1 - \lambda)`$ ($`\lambda < 1`$) gives the third.
- **Deduction of (1).** In every $`n`$-vertex tree ($`n \ge 2`$), $`\sum_{v: D_v(n) \ge 2} D_v(n) = 2(n - 1) - (n - C_n) = n + C_n - 2`$. For each fixed $`v`$, $`D_v(n)\mathbf 1_{D_v(n) \ge 2}/n \to Z_v`$ a.s. (both on eternal-leaf paths where $`Z_v = 0`$ and on promoted paths where $`D_v(n) \ge 2`$ eventually) and lies in $`[0, 1]`$. Bounded convergence as $`n \to \infty`$ gives $`\mathbb E\sum_{v \le m} Z_v \ge 1 - 2K(m^{\lambda-2} + m^{\lambda-1}/(1-\lambda))`$; sending $`m \to \infty`$ gives $`\mathbb E\sum_{v \ge 1} Z_v \ge 1`$, which combined with $`\sum_{v \ge 1} Z_v \le 1`$ a.s. forces $`\sum_{v \ge 1} Z_v = 1`$ a.s.

### 3. Total-variation and empirical-degree weak convergence

- For $`n \ge 3`$ and every $`v \in \{1, \dots, n\}`$ (both cores and leaves, since a leaf's sole neighbour in $`T_n`$ has degree $`\ge 2`$), $`p_n(v) := \mathbb P(W_n = v \mid T_n) = n^{-1}(L_v(n) + \sum_{w \sim v, \, D_w(n) \ge 2} D_w(n)^{-1})`$. The core-neighbour sum is bounded by $`C_n/(2n) \le C_n/n \to 0`$ a.s., while $`L_v(n)/n \to Z_v`$ a.s. On the single probability-one event $`\Omega_0`$ where $`p_n(v) \to Z_v`$, $`D_v(n)/n \to Z_v`$ for all $`v \ge 1`$, and $`\sum_{v \ge 1} Z_v = 1 = \sum_{v \ge 1} p_n(v)`$, truncating at $`m`$ gives $`\frac{1}{2}\sum_{v \ge 1}|p_n(v) - Z_v| \le \frac{1}{2}\sum_{v \le m}|p_n(v) - Z_v| + (1 - \sum_{v \le m} p_n(v))/2 + (1 - \sum_{v \le m} Z_v)/2 \to 0`$ (first $`n \to \infty`$, then $`m \to \infty`$). On the same $`\Omega_0`$, for every $`f \in C_b([0, 1])`$, $`|\sum_v p_n(v)f(D_v(n)/n) - \sum_v Z_v f(Z_v)| \to 0`$, proving almost-sure weak convergence.

### 4. Stopped localization under $`C_n = O(n^\lambda)`$ almost surely

- **Predictability of stopping at $`\tau`$.** With $`\tau = \inf\{n \ge 3 : C_n > K n^\lambda\}`$, defining $`X_n`$ on $`\{\tau > m\}`$ to include promotions **up to and including step $`\tau`$** and freezing at $`\tau`$ makes the active indicator for step $`T_{r-1} \to T_r`$ ($`r \ge m + 1`$) equal to $`\mathbf 1_{\{\tau \ge r\}} = 1 - \mathbf 1_{\{\tau \le r-1\}} \in \mathcal F_{r-1}`$. Thus existing summands satisfy $`\mathbf 1_{\{\tau \ge r\}}\mathbb E[\Delta(D_u(r)/r) \mid \mathcal F_{r-1}] \le 0`$, while a promotion at $`r \le \tau`$ adds $`(2/r)(C_{r \wedge \tau} - C_{(r-1) \wedge \tau})`$, proving $`\mathbb E[X_n] \le 2\,\mathbb E\sum_{r=m+1}^n (C_{r \wedge \tau} - C_{(r-1) \wedge \tau})/r`$.
- **Unit-jump bound and partial summation.** For $`r \ge 3`$: if $`\tau = 3`$, $`C_{r \wedge \tau} = C_3 = 1 \le K r^\lambda + 1`$; if $`r < \tau`$, $`C_{r \wedge \tau} = C_r \le K r^\lambda`$; if $`r \ge \tau \ge 4`$, $`C_{r \wedge \tau} = C_\tau \le C_{\tau - 1} + 1 \le K(\tau - 1)^\lambda + 1 \le K r^\lambda + 1`$ ($`\lambda > 0`$). Summation by parts bounds $`2\sum_{r=m+1}^n (B_r - B_{r-1})/r`$ ($`B_r = C_{r \wedge \tau} \le (K + 1)r^\lambda`$) by $`2(K + 1)(1 + (1 - \lambda)^{-1})m^{\lambda - 1}`$ uniformly in $`n \ge m`$.
- **Order of limits.** On $`\{\tau = \infty\}`$, $`X_n \ge \sum_{v > m, \, D_v(n) \ge 2} D_v(n)/n`$ and $`C_n/n \le K n^{\lambda - 1} \to 0`$. Fatou's lemma as $`n \to \infty`$ gives $`\mathbb E[\mathbf 1_{\{\tau = \infty\}}(1 - \sum_{v \le m} Z_v)] \le O_\lambda((K + 1)m^{\lambda - 1})`$. Sending $`m \to \infty`$ at fixed $`K`$ forces $`\sum_{v \ge 1} Z_v = 1`$ a.s. on $`\{\tau = \infty\}`$, and exhausting $`K \to \infty`$ covers the probability-one event $`\{C_n = O(n^\lambda)\}`$.

---

## Finite-size enumeration and negative controls

`hub_mass_carbon.py` (`OMP_NUM_THREADS=1 taskset -c 0 python3 -B round2/checks/hub_mass_carbon.py`, `1.26 s`) enumerates the exact rational law of all labelled histories $`(T_2, \dots, T_7)`$, verifying the core degree identity, the $`W_n`$ formula, the one-step supermartingale and promotion bounds ($`a_3, \dots, a_7 = 1, 4/3, 19/12, 647/360, 17189/8640`$), Equation (2), and the stopped bound across thresholds triggering $`\tau = 3`$, $`3 < \tau < 7`$, and $`\tau = \infty`$, while rejecting five planted errors (wrong promotion prefactor, unshifted denominator, omitted core-neighbour credit, omitted $`+1`$ overshoot at $`\tau`$, and dropping the stopping step $`r = \tau`$).
