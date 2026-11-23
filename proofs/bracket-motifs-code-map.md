# Code map for the motif-bracket proof

This file links [`bracket-motifs.md`](bracket-motifs.md) to the checker that proves its finite inequalities. Line numbers refer to two files:

- `V` is [`verify/motif-bracket/cert_verify_gmp.c`](../verify/motif-bracket/cert_verify_gmp.c).
- `F` is [`verify/motif-bracket/cert_format.h`](../verify/motif-bracket/cert_format.h).

The builder `cert_build.c` is untrusted. The proof depends only on what `V` checks and on the header values of the two shipped certificates.

## 1. What the checker establishes

The checker reads one binary certificate and prints `VERIFIED OK` with exit status `0` only if every check below passes. The other statuses are as follows.

- Status `1`: the file cannot be opened, a header check fails (`V` 278–288), the payload is shorter than the layout requires (`V` 323), or some tail–tail entry of $`g`$ is nonzero (`V` 341–344). The tail–tail test is one of the mathematical checks, but it exits at once with `1`, not `2`. The loader does nothing else to validate a file. In particular it ignores trailing bytes.
- Status `2`: a finite inequality fails, the comparison fails, or a brute-force tree fails (`V` 1100–1111).
- Built-in controls: `--override-pq 1132 2000` on the UPPER file, `--override-pq 114 200` on the LOWER file, and `--perturb f|g|hV|hC` each return `2`, since the change makes some finite inequality fail. `verify/run_all.sh` requires status `2` for six such runs.

Exit status `0` establishes the following:

1. Every finite inequality of §5.2 holds in exact rational arithmetic for the tables in the file, for the direction (`mode`) and the exponent `P/Q` given in the header.
2. $`g[a, b] = 0`$ for every pair of tail states.
3. $`c_0 > 0`$ (UPPER) or $`C_0 > 0`$ (LOWER), and $`\Phi(T_3) > 0`$.

Together with the header conditions of §3 below, these are exactly the hypotheses of §5 and §6 of the proof. The reader must still confirm, from the printed lines, that `Mode` and `Target` are `UPPER 1133/2000` or `LOWER 113/200`, and that the `Domains` line shows the parameters of §3.1.

## 2. Line-by-line correspondence

### 2.1 Layout and indexing (proof §3.2, §9)

| Proof item | Code |
|---|---|
| Inner states $`\mathcal{I}`$, enumerated `h = 1, 0`; `D = 2..D0-1`; `k = 0..D-h` | `F` 92–97, 107–112 |
| $`\mathrm{INF}`$, $`T_1(\kappa, h)`$, $`\kappa \in \{0, \dots, K_T+1\}`$ | `F` 98–100, 113–118 |
| Clamp $`\kappa = \min(k, K_T+1)`$ | `F` 56–59 |
| $`\mathcal{P}`$ ($`\mathrm{INF}`$, or $`k \ge 1`$), $`\mathcal{W}`$ ($`\mathrm{INF}`$, or $`h = 1`$) | `F` 125–130 |
| $`H_u`$: $`h = 1`$, $`1 \le k \le K_{hu}`$, $`D < D_{hu}`$ | `F` 63–67, 136–140 |
| $`C_u`$: $`2 \le k \le K_{cu}`$, $`D < D_{cu}`$ | `F` 68–72, 147–152 |
| $`\pi_f`$: $`\mathrm{VSTAR}_h`$ if $`k > K_F`$, else $`(\min(L, L_F), k, h)`$ | `F` 51–55 |
| $`h_V`$ offset over $`\mathcal{P} \times H_u \times \mathcal{W}`$ | `F` 74–77 |
| $`h_C`$ offset over unordered pairs (symmetry holds by storage) | `F` 78–84 |
| $`\mathrm{inv}`$: $`0`$ at $`\mathrm{INF}`$, $`1/D`$ for inner states, $`1/D_0`$ at $`T_1`$ (`gD = D0`) | `F` 113, 117; used at `V` 361, 368, 376, 544–545, 589–590, 717, 725, 749, 753 |
| Reading the payload in the order of §9, and rejecting a short file | `V` 296–323 |

### 2.2 Header checks (proof §3.1)

| Requirement | Code |
|---|---|
| magic | `V` 278 |
| `mode` is 1 or 2 | `V` 282 |
| $`P, Q, S, S_r, Q_K > 0`$ | `V` 283 |
| $`D_{hu} \le D_0`$, $`D_{cu} \le D_0`$, $`K_T + 1 \le K_F`$, $`D_0 \le L_F`$ | `V` 284 |
| 9 knots, $`\xi_0 = 0`$, $`\xi_1 = 1/D_0`$, $`\xi_8 = 1/2`$ | `V` 285–287 |
| knots strictly increasing | `V` 288 |

§3 below lists the conditions that are not checked.

### 2.3 Tools (proof §5.1)

| Proof item | Code |
|---|---|
| $`A`$ at $`x = \mathrm{inv}(v)`$, exact linear interpolation | `V` 25–60 |
| Lemma 5.1(1): $`A^*_{\mathop{\mathrm{ext}}}(y) = \mathop{\mathrm{ext}}_j(\xi_j y - a_j)`$, with numerator $`\mathrm{knot}_j\,y - Q_K a_j`$ over $`Q_K S`$ | `V` 628–673 |
| Lemma 5.1(3): $`A(t/D_0)`$ is affine in $`t`$ | needs $`\xi_1 = 1/D_0`$ (`V` 285) |
| Exact $`g`$ at an actual tail state, $`((D - D_0)\,g[\mathrm{INF}, \cdot] + D_0\,g[T_1, \cdot])/D`$ | `V` 63–91 |
| Choosing between two rates by the sign of the jump (fact (b), Lemma 5.2) | `V` 255–258; inline at `V` 567–568, 575–576, 609–610, 616–617 |
| Lemma 5.2: $`\gamma_{\min} = (D_0 - \kappa - h)/(D_0+1)`$ and $`b_{\min} = (D_0 - \kappa - h)/D_0`$ for $`\kappa \le K_T`$, $`\gamma_{\min} = 0`$ for $`\kappa = K_T + 1`$ | `V` 566, 608, 721, 729, 862, 870 |
| Lemma 5.3: $`\frac{1}{D_0}\mathop{\mathrm{ext}}(d_{\mathrm{INF}}, d_{T_1})`$ | `V` 580–584, 621–624, 776–779, 801–804, 906–909, 928–931 |

### 2.4 Received FY pieces $`A_{\mathrm{tot}}[a, b]`$ (proof §4.3, last two lines of $`E`$)

| Piece | Condition | Code |
|---|---|---|
| $`k_b A_{pc}[b](\mathrm{inv}(a))`$ | $`b`$ inner, $`h_b = 1`$, $`k_b \ge 1`$, $`b \notin H_u`$; every $`a \in \mathcal{P}`$ | `V` 358–365 |
| $`A_{cp}[a](\mathrm{inv}(b))`$ | $`a`$ inner, $`h_a = 1`$, $`k_a \ge 1`$, $`a \notin H_u`$; every $`b \in \mathcal{W}`$ | `V` 366–371 |
| $`(k_a - 1)A_{cc}[a](\mathrm{inv}(b))`$ | $`a`$ inner, $`k_a \ge 2`$, $`a \notin C_u`$; every $`b \in \mathcal{W}`$ | `V` 373–380 |

### 2.5 Vertex family (proof §5.2, cases (A)–(E))

Common values: $`\lambda`$ at `V` 393, $`f(\mathrm{nb})`$ at `V` 386, 394, and the right side $`\lambda f(\pi_f(s))`$ at `V` 399–400 and 426–427. The loop over all $`(L, k, h)`$ with $`L \le L_F`$, $`k \le K_F`$ and $`D \ge 2`$ is `V` 420–495.

| Case | Inequality | Code |
|---|---|---|
| (E) $`\mathrm{VSTAR}_h`$ | corners $`b \in \{0, 1\}`$ and $`t \in \{0, D_0/(K_F+1+h)\}`$ of $`b\,(f(\mathrm{nb}) + g(s^{\mathrm{bud}} \to \mathrm{nb}))`$; $`t = 0`$ is $`g[\mathrm{INF}, \mathrm{nb}]`$, and the other $`t`$ is the exact tail value at $`D = K_F + 1 + h`$ | `V` 396–417 |
| (C) $`L = L_F`$ | $`\frac{L_F}{D}\big(f(\pi_f(L_F - 1, k+1, h)) - f + f(\mathrm{nb}) + g(\cdot \to \mathrm{nb})\big)`$, with $`g`$ exact at $`t = D_0/D`$ | `V` 429–439 |
| (D) $`L \ge L_F + 1`$ | the four corners $`b \in \{\frac{L_F+1}{L_F+1+k+h}, 1\}`$, $`t \in \{0, \frac{D_0}{L_F+1+k+h}\}`$; $`f(\pi_f(s^{\mathrm{bud}})) = f(L_F, k+1, h)`$ or $`\mathrm{VSTAR}_h`$ | `V` 441–457 |
| budding part $`B(s)`$ for (A) and (B) | $`[L \ge 1]\frac{L}{D}(f(\pi_f(s^{\mathrm{bud}})) - f + f(\mathrm{nb}) + g(s^{\mathrm{bud}} \to \mathrm{nb}))`$ | `V` 461–469 |
| (A) inner | rate $`L + h\,r_{\mathrm{par}} + k\,r_{\mathrm{ch}}`$ times $`\Delta^L f`$, plus $`B(s)`$ | `V` 471–479 |
| (B) tail, $`L < L_F`$ | rates $`L`$ and $`L + (h+k)/2`$ | `V` 480–490 |

Count: the 302 states with $`L = L_F`$ give 5 checks each, the 572 inner states 1 each, the 241,024 other tail states 2 each, and $`\mathrm{VSTAR}`$ gives 8. The total is $`1{,}510 + 572 + 482{,}048 + 8 = 484{,}138`$, as printed.

### 2.6 Edge family (proof §5.2)

The tail–tail test (§5.2, "Both ends tail") is `V` 341–344, over all tail pairs; such pairs are skipped at `V` 527. The loop over $`(a, b) \in \mathcal{P} \times \mathcal{W}`$ is `V` 515–686.

| Summand of $`\mathrm{LHS}_g(a, b)`$ | Code |
|---|---|
| $`a^+`$, $`a^{\mathrm{bud}}`$ (inner, $`T_1`$ at $`D_a = D_0 - 1`$, $`T_1(\kappa+1)`$, or $`\mathrm{INF}`$ for leaf gain at a tail end) | `V` 517–518; for $`b`$: 528–529 |
| membership flags $`a \in H_u`$, $`a \in C_u`$, $`b \in H_u`$; the FY conditions; $`a^{\mathrm{bud}} \in C_u`$, $`b^{\mathrm{bud}} \in H_u`$ | `V` 519–523, 530–532 |
| $`\Delta_p g`$, $`\Delta_w g`$, budding jumps | `V` 534–537 |
| received pieces $`A_{\mathrm{tot}}[a, b]`$ | `V` 539 |
| transfer $`(\mathrm{inv}(b) - r_{\mathrm{ch}}[a])\Delta^L f(a)`$ | `V` 543–548 (with $`\Delta^L f`$ at `V` 502–504) |
| kept $`p`$-rate $`L_a + [b \in \mathcal{I}]\mathrm{inv}(b) + [a \in H_u]x_{\mathrm{par}} + [a \in C_u](k_a - 1)x_{cc}`$ times $`\Delta_p g`$ | `V` 550–555 |
| budding at inner $`a`$ with the newborn cherry $`h_C[a^{\mathrm{bud}}, b, \mathrm{nb}]`$ | `V` 557–563 |
| tail $`a = T_1`$: $`\gamma`$ term, $`\beta`$ term if $`\kappa \le K_T`$, quadratic term for $`\iota_p\,\Delta^L_w g`$ | `V` 564–585 |
| transfer $`(\mathrm{inv}(a) - r_{\mathrm{par}}[b])\Delta^L f(b)`$ | `V` 588–593 |
| kept $`w`$-rate $`L_b + [a \in \mathcal{I}]\mathrm{inv}(a) + [b \in H_u]k_b x_{\mathrm{ch}}`$ times $`\Delta_w g`$ | `V` 595–599 |
| budding at inner $`b`$ with the newborn path $`h_V[a, b^{\mathrm{bud}}, \mathrm{nb}]`$ | `V` 601–606 |
| tail $`b = T_1`$: $`\gamma`$, $`\beta`$, quadratic term for $`\iota_w\,\Delta^L_p g`$ | `V` 607–625 |
| $`A^*_{pc}[a](\Delta_p g)`$ | `V` 628–643 |
| $`(k_a - 1)A^*_{cc}[a](\Delta_p g)`$ | `V` 644–658 |
| $`k_b A^*_{cp}[b](\Delta_w g)`$ | `V` 659–673 |
| comparison with $`\lambda g[a, b]`$ | `V` 675–678 |

At $`a = \mathrm{INF}`$ or $`b = \mathrm{INF}`$ neither tail block runs, matching "At $`a = \mathrm{INF}`$ all these summands are $`0`$". The transfer at the other end then uses $`\mathrm{inv}(\mathrm{INF}) = 0`$ (`V` 544, 589).

### 2.7 Path family (proof §5.2)

Loop over $`(a, s, b) \in \mathcal{P} \times H_u \times \mathcal{W}`$: `V` 706–830. Every term is multiplied by $`D_{\mathrm{all}} = 2D_aD_sD_b(D_0+1)S_r`$ (with $`D_{\mathrm{INF}} := 1`$) and compared as integers (`V` 742–743, 819–823).

| Proof item | Code |
|---|---|
| $`s^+`$ (inner or $`T_1(k_s, 1)`$); $`s^+, s^{\mathrm{bud}}`$ outside $`H_u`$ give $`0`$ | `V` 710–712, 733, 736 |
| $`\Delta^L_s g(a \to s)`$, $`\Delta^L_s g(s \to b)`$ | `V` 720, 728 |
| sources $`(\mathrm{inv}(a) - x_{\mathrm{par}}[s])\Delta^L_s g(s \to b) + (\mathrm{inv}(b) - x_{\mathrm{ch}}[s])\Delta^L_s g(a \to s)`$ | `V` 748–754 |
| top inner: rates $`L_a + \mathrm{inv}(s)`$ and $`+\,(h_a + k_a - 1)/2`$; exact budding | `V` 757–766 |
| top tail: $`\gamma`$, $`\beta`$, quadratic term for $`\iota_p\,\Delta^L_s h_V`$ | `V` 767–780 |
| bottom inner: rates $`L_b + \mathrm{inv}(s)`$ and $`+\,k_b/2`$; exact budding | `V` 783–791 |
| bottom tail: $`\gamma`$, $`\beta`$, quadratic term for $`\iota_w\,\Delta^L_s h_V`$ | `V` 792–805 |
| middle: rates $`r_0 = L_s + [a \in \mathcal{I}]\mathrm{inv}(a) + [b \in \mathcal{I}]\mathrm{inv}(b)`$ and $`r_0 + (k_s - 1)/2`$; exact budding | `V` 808–817 |

When both ends are tail, both quadratic terms are added at the corner, as §5.2 states (`V` 776–779 and 801–804).

### 2.8 Cherry family (proof §5.2)

Loop over $`s \in C_u`$ and unordered $`(b_1, b_2)`$, with $`b_1 = b_2`$ allowed: `V` 844–956.

| Proof item | Code |
|---|---|
| apex state, $`s^+`$, $`s^{\mathrm{bud}}`$ (exits give $`0`$) | `V` 845–854, 873–874 |
| sources $`(\mathrm{inv}(b_1) - x_{cc})\Delta^L_s g(s \to b_2) + (\mathrm{inv}(b_2) - x_{cc})\Delta^L_s g(s \to b_1)`$ | `V` 861, 869, 884–888 |
| child $`w_1`$ inner, or tail with quadratic term | `V` 891–910 |
| child $`w_2`$ inner, or tail with quadratic term | `V` 913–932 |
| apex: rates $`r_0`$ and $`r_0 + (h_s + k_s - 2)/2`$; exact budding | `V` 935–944 |
| comparison with $`\lambda h_C`$ | `V` 946–950 |

### 2.9 Comparison and start (proof §6, §7)

| Proof item | Code |
|---|---|
| $`g_{\mathrm{ext}}`$, $`h_{V,\mathrm{ext}}`$ (start at $`0`$, extremum over $`\mathcal{P}`$, resp. $`\mathcal{P} \times H_u`$) | `V` 966–981 |
| $`h_{C,\mathrm{ext}}`$ | `V` 982–991 |
| $`\underline W`$ / $`\overline W`$ for $`L \le L_F`$, $`k \le K_F`$: inner bottom state | `V` 997 |
| tail bottom state: extremum over **both** candidate endpoints $`\mathrm{INF}`$ and $`T_1(\kappa, 1)`$, including the endpoint of weight $`0`$ at $`D = D_0`$ (proof §6, $`E(u)`$) | `V` 998–1002 |
| cherry part $`\binom{k}{2}h_{C,\mathrm{ext}}`$ | `V` 1004–1005 |
| $`\mathrm{VSTAR}_0`$ (root, no assigned edge or path) and $`\mathrm{VSTAR}_1`$ ($`\kappa = K_T + 1`$) | `V` 1009–1017 |
| $`\Phi(T_3) = f(2, 0, 0)`$, and positivity of it and of $`c_0`$ or $`C_0`$ | `V` 1018–1019 |
| printed prefactor $`\Phi(T_3)/c_0`$ or $`\Phi(T_3)/C_0`$ | `V` 1020–1025 |

### 2.10 Not used by the proof

The 37 brute-force trees (`V` 1027–1098) and the controls `--override-pq` and `--perturb` (`V` 279, 326–332) are sanity checks. The proof does not depend on them. The trees are evaluated directly by `V` 144–240, which does not use the decomposition of §4.

## 3. Header requirements

| Requirement (proof §3.1) | Checked by `V`? | Shipped value | What fails for another binary if it is false |
|---|---|---|---|
| $`D_0 \le L_F`$ | yes, `V` 284 | $`24 \le 800`$ | — |
| $`K_T + 1 \le K_F`$ | yes, `V` 284 | $`13 \le 150`$ | — |
| $`D_{hu} \le D_0`$, $`D_{cu} \le D_0`$ | yes, `V` 284 | $`16, 12 \le 24`$ | — |
| $`\xi_1 = 1/D_0`$, $`\xi_8 = 1/2`$ | yes, `V` 285 | $`10/240`$, $`120/240`$ | — |
| $`D_0 \ge 3`$ | implied: strictly increasing knots give $`1/D_0 = \xi_1 < \xi_8 = 1/2`$ | $`24`$ | — |
| $`D_0 \le K_F + 1`$ | **no** | $`24 \le 151`$ | f-case (E) and §6. An overflow state could be inner. Then $`\mathrm{VSTAR}_h`$ would stand for states with different $`g`$ entries, inner states with $`k > K_F`$ would have no vertex check, and $`\mathrm{VSTAR}`$ could lie in $`C_u`$. The checks would no longer cover every tree. |
| $`K_T \ge 0`$ | **no** | $`12`$ | The argument works for $`K_T = -1`$. For $`K_T \le -2`$ the index range of $`T_1`$ is empty or negative and the indexing in `F` 56–59 is invalid. |
| $`D_{hu} \ge 2`$, $`K_{hu} \ge 1`$, $`D_{cu} \ge 2`$, $`K_{cu} \ge 2`$ | **no** | $`16, 8, 12, 6`$ | Not needed by the mathematics: smaller values only make $`H_u`$ or $`C_u`$ empty. Negative values give negative allocation sizes in `F` 133, 143. |
| $`D_0 \le 64`$ (`CERT_MAX_D0`) | **no** | $`24`$ | Writes past the array `g_lut` (`F` 43, 96). This is a memory-safety bound, not a mathematical one. |
| int64 range of all intermediates | **no** | see §4 | Signed overflow is undefined behaviour and could hide a failed inequality. |

None of the unchecked requirements matters for the two shipped certificates. Their headers satisfy all of them, and the checker prints the values it uses on its `Domains` line:
```text
Domains: f(800,150) D0=24(KT=12,NG=601) hV(16,8;act=84) hC(12,6;act=75)
```
That line shows $`L_F, K_F, D_0, K_T, D_{hu}, K_{hu}, D_{cu}, K_{cu}`$, so a reader can confirm $`D_0 \le K_F + 1`$ and the other conditions from the checker's own output. For an arbitrary binary, exit status `0` proves the theorem only if the unchecked rows above are also confirmed.

## 4. Integer ranges for the shipped certificates

All inequalities are decided in GMP, but some intermediates are formed in `int64_t` or `long` first. The largest table entries are $`|f| \le 16{,}784{,}520{,}759`$, $`|g| \le 2{,}632{,}796{,}304`$, $`|h_V| \le 1{,}080{,}087{,}147`$, $`|h_C| \le 1{,}390{,}614{,}866`$, $`|r_{\mathrm{par}}|, |r_{\mathrm{ch}}|, |x_{\cdot}| \le 416{,}667`$, and knot values at most $`163{,}678{,}802`$ in absolute value (both files).

| Intermediate | Code | Bound |
|---|---|---|
| vertex rate numerator $`L S_r + h r_{\mathrm{par}} + k r_{\mathrm{ch}}`$ | `V` 474 | $`3.3 \times 10^{7}`$ |
| vertex-table difference $`\Delta^L f`$ stored as `dLf_q` | `V` 502–504 | $`2\max\lvert f\rvert \le 3.4 \times 10^{10}`$ |
| edge rate numerators, $`S_r - D r_{\mathrm{ch}}`$ | `V` 545, 551, 590, 595 | $`3.3 \times 10^{7}`$ |
| sums and differences of table entries | `V` 537, 561, 574, 581, 604, 615, 621, 720, 728, 762–815, 861–943 | $`6.7 \times 10^{9}`$ |
| denominators $`D_0 S`$, $`Q_K S`$ | `V` 583, 623, 637 | $`2.4 \times 10^{12}`$ |
| $`D_{\mathrm{all}}`$ and $`c_s = D_{\mathrm{all}}/D - D_{\mathrm{base}}x`$ | `V` 742–753, 880–887 | $`6.2 \times 10^{11}`$ |
| path and cherry rates such as $`L_a D_{\mathrm{all}} + D_{\mathrm{all}}/D_s + (h_a + k_a - 1)D_{\mathrm{all}}/2`$ | `V` 759–815, 892–941 | $`1.6 \times 10^{13}`$ |
| comparison sum $`W`$ | `V` 995–1016 | $`4.2 \times 10^{10}`$ |
| piecewise-linear evaluation $`x Q_K`$, $`\mathrm{knot}_j\,x_{\mathrm{den}}`$ | `V` 27–41 | $`2{,}880`$ |
| tree evaluator $`p_{\mathrm{den}} w_{\mathrm{den}}`$ and $`p_w w_w`$ | `V` 99–112, 124–137 | at most $`n_{\max}^2 \le 804^2 < 6.5 \times 10^{5}`$, where $`n_{\max} = 804`$ is the largest tree evaluated (the tree with $`L = L_F + 2`$ at the root, plus one attached leaf). These factors are actual degrees, not degrees capped at $`D_0`$. |

The largest bound, $`1.6 \times 10^{13}`$, is below $`2^{63} - 1 \approx 9.2 \times 10^{18}`$ by a factor of about $`6 \times 10^{5}`$. The checker itself does not test these ranges. A build with `-fsanitize=undefined -fno-sanitize-recover=all` ran both shipped files with no runtime error (proof §8). That build traps on signed overflow.

## 5. Corrections made to the expanded draft

These corrections were made while checking the draft against the checker. None of them changes a checked inequality.

1. §1, first line: restored the editorial note that the title and introduction were edited, and added that §§2–9 were reorganised.
2. §1, remark: replaced internal role wording with links to the existence and almost-sure proofs.
3. §3.3 and §5: replaced the two display labels with prose references to (3.1) and (5.0).
4. §2.3: $`C_u`$ does not fix $`h`$; only $`H_u`$ fixes $`h = 1`$.
5. §3.1: the condition $`D_0 \le K_F + 1`$ is also used in §6.
6. §4.2: the "kept" row for $`g`$ applies when the event vertex $`x`$ is tail. A tail $`p`$ does not stop the splits at an inner $`w`$, and the checker does not stop them either (`V` 530–531, 595, 659).
7. §5.2, cherries: "leaf end" changed to "child end"; these vertices are core children.
8. §6: the extremum runs over the candidate endpoints $`E(u)`$, which include $`\mathrm{INF}`$ even when its weight is $`0`$ at $`D = D_0`$. The draft wrote $`\mathrm{supp}\,\mu_u`$, which omits $`\mathrm{INF}`$ there. The checker uses both endpoints (`V` 998–1002). The new bound is never sharper and is still valid.
9. §8: the constructed checker tree has $`L = L_F + 2`$ leaves at the root, not $`L_F + 1`$ (`V` 1069–1076).
10. §4.3 and §8: the audit script that checks Identity (I) on trees is not part of the release. The text now says so and no longer uses the first person.
