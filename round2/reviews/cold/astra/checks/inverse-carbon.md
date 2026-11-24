CORRECT

# Cold review of `round2/inverse/search.py` and `round2/outputs/inverse.txt`

Reviewer: Gemini 4 Carbon  
Targets: `../inverse/search.py` (`inverse/search.py`) and `../outputs/inverse.txt` (`outputs/inverse.txt`)  
Checker: `inverse_carbon.py` (`python3 -B round2/checks/inverse_carbon.py`)

## Verdict summary

The exhaustive low-height search script `../inverse/search.py` (`inverse/search.py`) and its recorded output `../outputs/inverse.txt` (`outputs/inverse.txt`) are **CORRECT**. All six primitive relation counts, the exact rational enclosure of the Hartle–Krapivsky–Redner–Zhang (HKRZ) integral, the outward integer pruning bounds, and the Sturm / bisection root tests are mathematically sound and complete (zero false negatives on $`[a, b]`$, including multiple and endpoint roots), and `../outputs/inverse.txt` (`outputs/inverse.txt`) properly qualifies the search window as unproved.

---

## Audit of the four requested components

### 1. Count of primitive polynomials and HKRZ relations

- **Möbius inversion.** In `mobius(n)`, trial division removes smaller prime factors sequentially so any composite remainder satisfies $`n \ge p^2`$; square factors return $`0`$ and a prime remainder $`> 1`$ flips sign once. By Möbius inversion, `primitive_count(m, H)` $`= \sum_{k=1}^H \mu(k)((2\lfloor H/k\rfloor + 1)^m - 1)`$ equals $`P(m, H) := \#\{v \in [-H, H]^m \setminus \{0\} : \gcd(v) = 1\}`$.
- **Sign-normalized degree-$`d`$ polynomials.** Since $`\gcd(c_0, \dots, c_{d-1}, 0) = \gcd(c_0, \dots, c_{d-1})`$, primitive $`(d+1)`$-tuples with $`c_d \ne 0`$ number $`P(d+1, H) - P(d, H)`$, and negation bijects $`c_d > 0`$ with $`c_d < 0`$. Hence `(primitive_count(d+1, H) - primitive_count(d, H)) // 2` exactly counts primitive integer polynomials of exact degree $`d`$, height $`\le H`$, and $`c_d \ge 1`$. *(Note: "primitive" means $`\gcd(c_0, \dots, c_d) = 1`$ in $`\mathbb Z[x]`$, so reducible polynomials are included in both the count and the scan; all five reported polynomial matches in `../outputs/inverse.txt` (`outputs/inverse.txt`) are in fact $`\mathbb Q`$-irreducible.)*
- **HKRZ relations.** Excluding pure quadratics $`(f_0, f_1) = (0, 0)`$ ($`P(3, H)`$ primitive tuples) from $`P(5, H)`$ and halving by sign normalization ($`f_1 > 0`$, or $`f_1 = 0`$ and $`f_0 > 0`$) gives `(primitive_count(5, H) - primitive_count(3, H)) // 2` $`= 4{,}661{,}184`$ at $`H = 12`$. All six counts in `../outputs/inverse.txt` (`outputs/inverse.txt`) (`30401831`, `3355121`, `6261487`, `13418113`, `11116167`, `4661184`) match an independent Möbius-free GCD recursion.

### 2. HKRZ interval arithmetic and constant-term pruning

- **Exact alternating enclosure of $`F(x) = \int_0^1 t^x e^{t-1}\, dt`$.** Repeated integration by parts for $`x > -1`$ gives
  $`\displaystyle F(x) = \sum_{k=0}^{j-1}\frac{(-1)^k}{(x+1)(x+2)\cdots(x+1+k)} + R_{j-1}(x), \qquad R_{j-1}(x) = \frac{(-1)^j}{(x+1)\cdots(x+j)}\int_0^1 t^{x+j}e^{t-1}\,dt.`$
  For $`x \in (0, 1)`$, $`0 < \int_0^1 t^{x+j}e^{t-1}\,dt < \frac{1}{x+1+j}`$, so $`R_{j-1}(x)`$ has the exact sign of the $`j`$-th term $`a_j`$ and $`|R_{j-1}(x)| < |a_j|`$. Thus `min(total, new_total) < F(x) < max(total, new_total)` in exact `Fraction` arithmetic, terminating at $`j = 32 < 80`$ with width $`< 10^{-35}`$.
- **Monotonicity and positivity on $`[a, b] \subset (0, 1)`$.** Since $`\log t < 0`$ on $`(0, 1)`$, $`x \mapsto F(x)`$ is strictly decreasing, so `lower_f = hkrz(b)[0] < F(x) < hkrz(a)[1] = upper_f`. Because $`0 < a \le x \le b`$ and $`\mathrm{lower\_f} > 1/3 > 0`$, `[a * lower_f, b * upper_f]` validly encloses $`x F(x)`$ on $`[a, b]`$.
- **Integer pruning bounds.** In `search_polynomials`, `den = lcm(a.denominator, b.denominator)` makes `a**j * scale` and `b**j * scale` exact integers. In `search_hkrz`, all four positive lower endpoints are floored via `int(x * UNIT)` and upper endpoints ceilinged via `ceil(x * UNIT)`, widening `[lo, hi]` conservatively. For every signed integer pair `(lo, hi)` and positive `scale` / `UNIT`, Python's floor division satisfies `-(hi // scale) == ceil(-hi / scale)` and `(-lo) // scale == floor(-lo / scale)`, so no valid constant term $`c_0`$ is ever pruned.

### 3. Root completeness

- **Polynomials (`poly_has_root`).** Line 66 (`eval_poly(coeff, a) * eval_poly(coeff, b) <= 0`) catches all endpoint roots ($`P(a) = 0`$ or $`P(b) = 0`$) and odd-parity sign changes in exact `Fraction` arithmetic. When $`P(a)P(b) > 0`$, the exact rational Euclidean Sturm sequence $`[P_0, P_1, \dots, P_m]`$ counts **distinct** real roots in $`(a, b)`$ even if $`P`$ is not squarefree: since $`P_m = c\gcd(P, P')`$ is nonvanishing at $`a`$ and $`b`$, dividing the entire tuple by $`P_m(x)`$ preserves all sign variations while $`\tilde P_0/\tilde P_1 = P/P' \sim (x - r)/\mu_r`$ loses one variation across every root $`r \in (a, b)`$ of any multiplicity $`\mu_r \ge 1`$. Thus multiple (even-order tangential) roots and even-count root clusters inside $`(a, b)`$ cannot be missed.
- **HKRZ relations (`hkrz_root`).** A call returns `"outside"` only if `lo > 0` or `hi < 0` on the whole interval (or on both bisected halves). Any zero $`x^* \in [a, b]`$—including a tangential root or an endpoint zero—forces $`\mathrm{lo} \le 0 \le \mathrm{hi}`$ on every subinterval containing $`x^*`$, so it can never be classified `"outside"` (a non-sign-changing zero bisects to `depth == 20` and reports `"unresolved"`). Conversely, `"root"` requires strict opposite point-enclosure signs `left * right < 0`, certifying a root by the Intermediate Value Theorem. Because `../outputs/inverse.txt` (`outputs/inverse.txt`) reports `0 roots, 0 unresolved`, none of the $`4{,}661{,}184`$ primitive HKRZ relations of height $`\le 12`$ vanishes anywhere on $`[11314469/20000000, 11314477/20000000]`$.

### 4. Meaning of the statistical window

- **Heuristic status vs rigorous bounds.** The default window $`I = [0.56572345, 0.56572385] = [11314469/20000000, 11314477/20000000]`$ (width $`|I| = 4 \times 10^{-7}`$) is a provisional search band around the finite-time ($`t \in [28, 30]`$) motif-potential logarithmic-derivative estimates in `../numerics/RESULTS.md` (`numerics/RESULTS.md`) ($`\approx 0.56572363\text{--}0.56572368`$ with trajectory jackknife SEs $`\sim 3\text{--}7 \times 10^{-8}`$; see `../numerics/README.md` (`numerics/README.md`)). While it spans several Monte Carlo standard errors, deterministic pre-asymptotic relaxation bias between finite $`t \le 30`$ and $`t \to \infty`$ is unproved, so rigorous membership remains only the certified interval $`[113/200, 5657/9983]`$ (or $`1133/2000`$). Line 1 of `../outputs/inverse.txt` (`outputs/inverse.txt`) states this explicitly (`Membership of mu in this window is unproved.`).
- **Interpretation of the recorded hits and exclusions.** Conditionally on $`\mu \in I`$, `../outputs/inverse.txt` (`outputs/inverse.txt`) rigorously excludes all quadratics of height $`\le 100`$ and all HKRZ relations of height $`\le 12`$. Meanwhile, in degrees $`1, 3, 4, 5`$ where the candidate counts are $`3.04\times 10^7`$, $`6.26\times 10^6`$, $`1.34\times 10^7`$, and $`1.11\times 10^7`$, the product $`(\text{candidate count}) \times |I|`$ is of order $`1`$ (for example, $`\frac{3}{\pi^2}\cdot 5000^2 \cdot 4\times 10^{-7} \approx 3.04`$ expected reduced fractions of denominator $`\le 5000`$ in a generic interval of width $`4\times 10^{-7}`$), so finding $`2, 1, 1, 1`$ roots is fully consistent with ordinary root-density background rather than evidence for an algebraic closed form.

---

## Independent verification

`inverse_carbon.py` (`OMP_NUM_THREADS=1 taskset -c 0 python3 -B round2/checks/inverse_carbon.py`, `3.39 s`) verifies all six counts without Möbius inversion, cross-checks small heights against unpruned brute force, checks all five reported polynomials in `../outputs/inverse.txt` (`outputs/inverse.txt`), validates `hkrz(x)` against an independent 55-digit positive Taylor series $`e^{-1}\sum_{m\ge0}\frac{1}{m!(x+1+m)}`$, confirms `search_hkrz` recovers the solvable DAN (`(-1, 0, 1, 1, 0)`) and PAN (`(-1, 1, 0, 1, 0)`) equations near $`0.771192`$ and $`0.549735`$, verifies detection of double roots $`(3x-2)^2`$, paired close roots, and endpoint roots, and catches four planted errors.
