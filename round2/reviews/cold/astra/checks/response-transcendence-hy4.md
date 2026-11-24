CORRECT

# Cold audit: prime-denominator obstruction and algebraic-series argument

Scope. Audit of §2 of `isolated-transfer.md` only: the prime-denominator calculation for the newborn derivatives, the inference to a transcendental Laplace function, and the formal-series/polynomial argument. The verdict concerns the scalar response as a function of s; it makes no statement about the arithmetic of the critical root. Starting point is the reviewed fringe/transfer algebra; no excursion number is used here and no conclusion about μ is drawn. No author file was changed and nothing was committed. Author: hy4.

## 1. Prime-denominator step

For the constant mark A=1/2, `round2/checks/response_transcendence_hy4.py` reimplements from scratch the tree representation, the intrinsic transitions (leaf gains, promotions, child copy gains), the residual generator and the newborn reward b=L/D (D the **root degree**, not total size). It obtains independently:

- intrinsic derivatives through k=5: 1/2, −1/12, 11/72, −121/1080, 8789/129600, −184243/4536000;
- constant-gauge derivatives: 1/2, −1/3, 19/36, −127/160, 321353/259200, −26334359/13608000.

These reproduce the shipped values. For p=3,5,7 and k=p−2 the script enumerates **every** forward path of length k from the newborn and groups by final (size, degree, number of core children). At k=p−2 every reachable final state has size p; the only one with degree p is the pure star with p−1 leaves and no core child, and its total contribution is exactly (p−2)!(p−1)/p. All other final states have degree <p, so neither their reward nor their rates carry a factor p. The full derivative includes the diagonal (killing and total-rate) terms, whose denominators are likewise prime to p before the last action. Consequently d_{p−2}=(p−2)!(p−1)/p + (rational with denominator prime to p), so the reduced denominator of d_{p−2} is divisible by p for each of p=3,5,7, and the reduced denominator of the difference is not. This is exactly the stated claim, and it gives infinitely many distinct denominator primes because there are infinitely many odd primes.

The degree-versus-size distinction is essential: replacing b=L/D by L/(total size) makes non-star final states at p=5 carry denominator 5, and the uniqueness argument fails. That is included as a negative control and is rejected by the correct definition.

## 2. Algebraic-series step

The second part of the argument is audited as a formal statement about the response function. A rational Laplace transform forces a constant-coefficient recurrence on the derivatives, possibly after finitely many terms; rational elimination gives a rational solution of a maximal independent finite subsystem, hence only finitely many denominator primes. The stronger claim uses a possible polynomial relation P(s,F(s))=0 over C(s). The complete large-s expansion is s⁻¹H(1/s) with H(z)=Σd_kz^k; a polynomial relation then gives a nonzero formal relation for H, and rational coefficients reduce linear dependence over C to dependence over Q by rational elimination. For irreducible P over Q in characteristic zero, P_y(z,H) has finite order r with leading coefficient c≠0; for n>r the coefficient of z^{n+r} gives c·d_n=−[z^{n+r}]P(z,H_{<n}), because powers (H−H_{<n})^j with j≥2 start at degree ≥2n>n+r. By induction the denominator primes of every d_n lie in the finite set supplied by c and the initial block. This contradicts §1.

The mechanism was exercised on a Catalan toy, P(z,H)=H−z−zH², whose coefficients are reproduced by the stated truncation rule (1,1,2,5,14,42,132), and a planted new prime 7 in the eighth coefficient is correctly excluded by the recurrence. No error was found in the formal manipulation.

## 3. Extensions stated in the text

- For L⁰, P_A and Q_A the other transitions either evaluate the reward at a newborn or decrease subtree size, so the unique minimal path to a degree-p star is still the p−2 intrinsic gains. The finite-mark residual row works after excluding the finitely many primes dividing the mark denominators. Both are structurally correct.
- For the active/dormant return, the same conclusion follows from step counting, not from the one sentence in the text: starting from a newborn of size 2, every step that can appear on a minimal length-(p−2) path ending at a degree-p star must increase the rewarded particle's size by one. Dormant emission creates a fresh size-2 particle and activation preserves size, so neither can occur on such a minimal path; a path using either needs more than p−2 steps to reach degree p. I verified this counting but did **not** run a full block-system derivative computation for \widetilde K, so this extension is reasoning-level here rather than replay-level.

## 4. Checker output

```
p 3 unique terminal star 2/3 and d_{p-2} has p in reduced denominator
p 5 unique terminal star 24/5 and d_{p-2} has p in reduced denominator
p 7 unique terminal star 720/7 and d_{p-2} has p in reduced denominator
negative control: degree->size reward breaks uniqueness at p=5
formal recurrence c*d_n = -[z^(n+r)]P(z,H_<n) reproduces Catalan; planted new prime 7 is excluded
SCOPE: response as a function of s only; nothing claimed about the critical root
RESPONSE TRANSCENDENCE COLD CHECK: PASS
```

Commands, one core:

```
taskset -c 0 python3 -B round2/checks/response_transcendence_hy4.py
```

## 5. Scope boundary and limitations

The proof concerns the scalar response as a function of s. It says nothing about whether the critical root μ is algebraic or transcendental; no step here evaluates the response at μ or connects the denominator primes to μ. It also does not decide the existence of a classical integral or an ODE with variable coefficients for the complete response, exactly as the source states. The asymptotic-to-formal step assumes the complete descending-power expansion justified by Taylor expansion on a short interval plus the exponential tail bound; that justification is a stated prerequisite and was not reproved here.

## 6. Verdict

The prime-denominator computation and the algebraic-series argument are correct for the response function, with the small code-level limitation on the dormant extension noted in §3. No claim is made about the arithmetic of μ.
