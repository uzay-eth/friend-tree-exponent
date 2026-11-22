VERDICT: CORRECT

The summation-by-parts equality in part (c), line 129, fails at the admissible endpoint \(k=2\), where \(m_0=n\). Its left side is zero and its right side is \(a_n/[n(n+1)]>0\). The inequality used next still holds. Handling that endpoint directly with part (a) proves part (c) with the stated constants. Theorem D and the numerical Corollary follow with this repair.

## Written proof

I refer to the line numbers in `degree_tail.md`. Only the written proof and the supplied Gamma-ratio bounds are used to judge the mathematical claims.

### Model and Lemma 1, lines 43–54

- **VALID** — Given the tree history up to \(T_j\), one must first select a neighbor \(u\) of \(w\), then select \(w\) among the \(D_u(j)\) neighbors of \(u\). Summing the conditional probabilities gives the attachment formula. Each term in the sum before the factor \(1/j\) is at most one.
- **VALID** — With \(d=D_v(j)\) and \(p=\Pr(W_j=v\mid\mathcal F_j)\le d/j\), the conditional first moment at \(j+1\) is \(d+p\le d(j+1)/j\). The tower property gives the bound from time \(m\) with \(\prod_{j=m}^{n-1}(j+1)/j=n/m\).
- **VALID** — The increment of \(Q=D(D+1)\) on a hit is \((d+1)(d+2)-d(d+1)=2(d+1)\). Thus the next conditional moment is at most \(Q(j)+2(d+1)d/j=Q(j)(j+2)/j\). The tower property and \(\prod_{j=m}^{n-1}(j+2)/j=n(n+1)/[m(m+1)]\) give the result. Both products are one when \(n=m\).

### Part (a), lines 56–60

- **VALID** — Degrees never decrease and each vertex enters as a leaf, including both initial vertices at time two. Every vertex counted in \(N_k(n)\), for \(k\ge2\), is either core at \(m_0\) or has a unique promotion step \(m\in[m_0,n-1]\). This includes vertices born after \(m_0\). The first group contributes at most \(C_{m_0}\). At step \(m\), the attachment creates one new leaf and promotes at most one previously existing leaf, so \(C_{m+1}-C_m\) is the indicator of promotion.
- **VALID** — On the promotion event, the identity of the promoted vertex and its new degree two are \(\mathcal F_{m+1}\)-measurable. Lemma 1 applies to that random vertex by conditioning separately on each possible label. Markov gives \(2n/[(m+1)k]\) as an upper bound for its conditional probability of reaching degree \(k\). Multiplication by the promotion indicator and averaging yields the summand in (a). The sum is empty if \(m_0=n\), and (a) then gives \(E[N_k(n)]\le a_n\).

### Lemma 2, lines 62–77

- **VALID** — Attachment to \(u\) adds a leaf neighbor. Its conditional probability is \((L_u+\rho_u)/j\), where \(\rho_u\) is the sum of reciprocal degrees of the core neighbors. A leaf neighbor of \(u\) can be promoted only if \(V_j=u\) and \(W_j\) is that leaf. The total probability of losing one leaf neighbor is \(L_u/(jD_u)\le1/j\). These changes are disjoint and exhaust the changes to \(L_u\).
- **VALID** — The drift bound yields \(E[L_u(j+1)-1\mid\mathcal F_j]\ge (1+1/j)(L_u(j)-1)\). The tower property preserves this inequality, and iteration multiplies by positive deterministic factors. Even on histories with \(L_u(j)-1<0\), the iteration gives \(E[L_u(n)-1\mid\mathcal F_m]\ge(l-1)n/m\).
- **VALID** — Pointwise \(L_u(n)^2\le D_u(n)(D_u(n)+1)\). At time \(m\), \(D_u(m)=l+c_u(m)\le l+2\). Lemma 1 and \((n+1)/(m+1)\le n/m\) give the stated upper bound \((l+2)(l+3)(n/m)^2\).
- **VALID** — Conditional Paley–Zygmund at half the conditional mean gives the stated bound for that threshold. The threshold in the lemma is at most half the conditional mean, so the event in the lemma includes the Paley–Zygmund event. The lower bound is one quarter the product of \((l-1)/(l+2)\) and \((l-1)/(l+3)\). Both factors increase for \(l\ge3\), and the product with the factor \(1/4\) is \(1/30\) at \(l=3\).

### Lemma 3, lines 79–90

- **VALID** — For \(C=1\) the single core vertex has no core neighbors and at least two leaf neighbors. It is eligible. For \(C\ge2\), the path between two core vertices stays inside the core because interior vertices of a path have degree at least two. The induced core is therefore a tree.
- **VALID** — A core-tree leaf has one core neighbor and at least one leaf neighbor in the original tree. A core vertex that is ineligible either has core degree at least three or is an isolated \((L,c)=(0,2)\) vertex among those of its kind. The latter \(I\) vertices share no edge with each other and each uses two distinct core-tree edges, hence \(2I\le C-1\).
- **VALID** — Summing \(c-2\) in the core tree gives \(-2\), so \(n_1=2+\sum_{c\ge3}(c-2)\ge 2+n_3\). Adding \(E\ge C-n_3-(C-1)/2\) and \(E\ge n_1\ge n_3+2\) gives \(E\ge(C+5)/4\) for \(C\ge2\), which implies \(E\ge C/4\).

### Lemma 4, lines 92–119

Write \(r=m'\) for the block calculations below.

- **VALID** — Label vertices by birth order, with the first two labels assigned to the initial edge. At time \(j\) the set of labels is deterministic, so independent draws uniform on the first \(j\) labels realize the \(V_j\). Use fresh independent randomness at each step to select the neighbor of \(V_j\). Thus \(V_j\) is independent of the tree history at \(j\), and the future \(V_j\) are jointly independent given \(\mathcal F_r\). Given \(\mathcal F_r\), the focus \(u\) is fixed. Conditioning further on \(G\) separately removes label \(u\) from the uniform law at each time and preserves independence and the fresh neighbor-selection randomness.
- **VALID** — The product for \(G\) telescopes to \((r-1)/(8r-1)\ge2/23\), since \(23(r-1)-2(8r-1)=7(r-3)\). On \(G\) no leaf of \(u\) can be promoted, as only \(u\) could select such a leaf as a neighbor. Every new neighbor of \(u\) starts as a leaf and likewise cannot be promoted on \(G\). Consequently \(L_u\) is nondecreasing, each of its leaves persists, and \(c_u\) is constant throughout the block.
- **VALID** — In case A fix an existing leaf \(x\) given \(\mathcal F_r\). Each draw of that leaf on \(G\) attaches to \(u\). The no-hit product for a fixed non-\(u\) label on \([s,t)\) given \(G\) is \(\prod_{j=s}^{t-1}(j-2)/(j-1)=(s-2)/(t-2)\). Thus the hit probabilities on \([r,3r)\) and \([3r,8r)\) are \(2r/(3r-2)\ge2/3\) and \(5r/(8r-2)\ge5/8\). The draws on the two blocks are independent given the initial history and \(G\), and the two attachments bring \(L_u\) to at least three. The product of the lower bounds is \(5/138\).
- **VALID** — In case B, \(z\) and \(w\) are distinct from \(u\) and from each other. Conditional on the initial history and \(G\), \(H_1\) is the event of no \(z\) minus the event of neither \(z\) nor \(w\) on \([r,2r)\). The two no-hit products are \((r-2)/(2r-2)\) and \((r-2)(r-3)/[(2r-2)(2r-3)]\), respectively. Their difference is \(r(r-2)/[(2r-2)(2r-3)]\). Comparing four times the numerator to the denominator leaves \(2(r-3)\ge0\), giving the bound \(1/4\).
- **VALID** — Before the first hit \(\tau\) of \(w\), neither \(u\) nor \(z\) was drawn as \(V_j\). Any first attachment to \(w\) would have required drawing one of these two original neighbors, hence its degree is still two at \(\tau\). Conditional on the chooser sequence specifying \(G,H_1,\tau\), and the prior neighbor draws, the neighbor choice at \(\tau\) is fresh and sends the attachment to \(u\) with probability \(1/2\).
- **VALID** — On that success, \(x\) is the label assigned at birth on the step from \(T_\tau\) to \(T_{\tau+1}\). It is present by \(T_{2r}\) and stays a leaf of \(u\) on \(G\). Conditional on \(G\) and the full choices before \(2r\) that produced \(x\), its label is fixed and distinct from \(u\). Future draws still have the independent uniform laws on labels other than \(u\). The hit probabilities on \([2r,4r)\) and \([4r,8r)\) are \(2r/(4r-2)\) and \(4r/(8r-2)\), both at least \(1/2\). These block events are conditionally independent, and their probabilities are the same for every allowed early history. This justifies multiplying by the earlier event probability even though the label \(x\) is random. The two hits add two leaves after the first one, and \(c_u=2\) throughout on \(G\). The product bound is \((2/23)(1/4)(1/2)(1/2)^2=1/368\).

### Part (b), lines 121–124

- **VALID** — For each eligible vertex in \(T_{m'}\), Lemma 4 supplies an event measurable at time \(m=8m'\) with conditional probability at least \(1/368\). On that event \(l\ge3\) and \(c\le2\), so Lemma 2 applies conditional on \(\mathcal F_{8m'}\). Since \(n\ge8m'k\), its threshold \((l-1)n/(2m)\) is at least \(n/(8m')\ge k\). Reaching at least \(k\) leaf neighbors implies reaching degree at least \(k\).
- **VALID** — The tower property gives the product \(1/(368\cdot30)\) for each eligible vertex. Count those distinct vertices by linearity of conditional expectation and apply Lemma 3. This introduces the factor \(1/4\) and gives \(4\cdot368\cdot30=44160\). Joint independence across eligible vertices is unnecessary.

### Part (c), lines 126–134

- **VALID** — Under \(n\ge24k\), \(2n/k\ge48\) and \(k\ge2\) gives \(2n/k\le n\). Thus \(3\le m_0=\lceil 2n/k\rceil\le n\), and (a) applies.
- **FALSE** — The equality on line 129, over its stated range, fails for \(m_0=n\). There are no increments on the left, and the right equals \(a_n(b_n-b_{n+1})>0\). This happens for every \(k=2\) in the range of (c), including \(n=48\). For \(m_0<n\) the equality follows by expanding the differences and collecting the coefficients of \(a_j\).
- **NEEDS REPAIR** — Replace the equality on line 129 by \(\le\). For \(m_0<n\) the right side is still equal to the left, and for \(m_0=n\) it is positive. An exact formula for the whole range uses \(a_n b_{n+1}\) as the first boundary term and sums the interior terms through \(j=n\). Alternatively, for \(k=2\) apply (a) directly to get \(E[N_2(n)]=a_n\le A n^\lambda\le A(3n/k)^\lambda\).
- **VALID** — With that repair, drop \(-a_{m_0}b_{m_0+1}\le0\), apply \(a_j\le A j^\lambda\), and use \(b_j-b_{j+1}=1/[j(j+1)]\le j^{-2}\). The sum through \(n-1\) can be extended to all \(j>m_0\) since the terms are nonnegative. The function \(x^{\lambda-2}\) is positive and decreasing, so that sum is at most \(\int_{m_0}^\infty x^{\lambda-2}\,dx=m_0^{\lambda-1}/(1-\lambda)\). Also \(n^{\lambda-1}\le m_0^{\lambda-1}\). This proves line 130.
- **VALID** — Inserting the bound in (a), the first term contributes at most \(A m_0^\lambda\). The remaining terms contribute at most \(A m_0^\lambda(1+1/(1-\lambda))\) since \(2n/k\le m_0\). The upper constant is therefore \(A(2+1/(1-\lambda))\) on \(m_0^\lambda\). The ceiling gives \(m_0\le2n/k+1\le3n/k\).
- **VALID** — With \(y=n/(8k)\ge3\), \(m'=\lfloor y\rfloor\) satisfies \(m'\ge3\), \(8m'k\le n\), and \(m'\ge y-1\ge(2/3)y=n/(12k)\). Applying (b) with \(a_{m'}\ge a(m')^\delta\) gives the lower constant \(a/44160\) on \((n/(12k))^\delta\).

### Part (d), lines 136–138

- **VALID** — Take \(0<\varepsilon<\min(\mu,1-\mu)\) when applying (c). The limit of \(\log a_m/\log m\) gives the required power bounds for sufficiently large \(m\). Enlarging the upper constant and shrinking the positive lower constant covers the finitely many smaller \(m\), since \(a_m\ge1\) for \(m\ge3\).
- **VALID** — With \(r=n/k\ge24\), (c) bounds \(E[N_k(n)]\) between positive constants times \(r^{\mu-\varepsilon}\) and \(r^{\mu+\varepsilon}\), and the constants only depend on \(\varepsilon\) and the sequence \(a_m\). Taking logs and dividing by \(\log r\) makes their contributions \(O_\varepsilon(1/\log r)\), uniformly in \(k\). Let \(r\to\infty\), then \(\varepsilon\downarrow0\). This establishes (d) once the endpoint in (c) is handled.

### Numerical Corollary, lines 31–35

- **VALID** — For \(x>0\), \(0<s<1\), strict concavity of \(t^s\) and the mean \(x\) of a gamma random variable of shape \(x\) and rate one give \(\Gamma(x+s)/\Gamma(x)<x^s\). Applying the upper bound at \(x=m+s\) with exponent \(1-s\) gives \(\Gamma(m+1)/\Gamma(m+s)<(m+s)^{1-s}\). Thus \(\Gamma(m+s)/\Gamma(m)>m(m+s)^{s-1}>m(m+1)^{s-1}\). This verifies the version of Gautschi used here.
- **VALID** — For \(m\ge3\), \(m/(m+1)\ge3/4\), so the assumed lower bound is at least \(0.5022(3/4)^{0.45}m^{0.55}\). The coefficient is \(0.44121907205027\ldots>0.44\). The given upper bound and Gautschi give \(a_m\le3.974\,m^{0.58}\).
- **VALID** — Substitution into (c) gives the written lower prefactor \(0.44/44160\). The coefficient of \((n/k)^{0.58}\) in the upper bound is \(3.974(2+1/0.42)3^{0.58}=32.92507122519528\ldots<33\). Hence the numerical inequalities follow with the repair to (c).

### Claims around the theorem

- **NEEDS REPAIR** — If “of the order of” on lines 10–11 is intended as two-sided constant-factor comparison with \(E[C_{\lfloor n/k\rfloor}]\), qualify it as a statement about the logarithmic exponent under (d). The displayed log-scale conclusion does not establish that comparison.
- **NEEDS REPAIR** — On lines 37–39 the conditional claim about the cumulative exponent is an application of (d) if the stated external exponent inputs hold. A cumulative logarithmic exponent alone does not determine a pointwise degree-density exponent. For example, a distribution putting mass proportional to \(2^{-j\mu}\) at degrees \(2^j\) has a tail of logarithmic exponent \(-\mu\) and zero mass at the intervening degrees. Describe \(1+\mu\) as the density exponent suggested by the power-law heuristic unless a separate density estimate is supplied. This does not affect the numerical inequalities.

## Numerical checks

I wrote `referee.py` and `boundary.py`. The Python test commands were pinned to logical CPU 0. Exact enumeration and product checks use rational arithmetic. Monte Carlo intervals are pointwise and serve as checks on the sampled configurations.

1. **Lemma 3.** Enumerate all unlabeled trees by adding a leaf to each vertex of every tree at the previous order, then deduplicating by canonical rooted encodings minimized over all roots. Removing a leaf proves enumeration coverage by induction. The two trees of orders one and two are core-free. Every tree of order at least two can also arise from the process, since each required attachment has positive probability through some neighbor. In all 2,286 trees with a core through order 13 the inequality holds. The smallest ratio is \(3/7\), attained by 14 trees at orders 10–13, with edge lists in the output. None attains \(E=C/4\). The smallest additive slack is \(3/4\), attained by the stars with \(C=1\), since the stronger bound for \(C\ge2\) gives slack at least \(5/4\). Omitting the adjacent-pair eligibility clause caused seven failures.
2. **Lemma 1.** Generate 224 simulated trees at seven orders up to 300. Enumerate each ordered \((V,W)\) pair with probability \(1/(jD_V(j))\), then evaluate both moments exactly for every vertex. Both inequalities hold on all 14,688 vertex checks. Using \(j+1\) in place of \(j\) as the growth denominator caused 109 failures in each moment.
3. **Lemma 2.** Use 768 paths per starting configuration at ratios \(2,8,50\), with path-shaped or fan-shaped branches when both are available. The case \(l=3,c=0\) forces the four-vertex star, since all neighbors of \(u\) are leaves and the tree is connected. I included it at \(m=4\) and also used \(m=5,l=4,c=0\) to cover the requested lower end of the starting-size range. With \(l=3\), the minimum orders for \(c=1\) and \(c=2\) are six and eight, respectively, because each core branch needs at least two vertices. The smallest pointwise Wilson lower endpoint was \(0.790215>1/30\). Switching from attachment to \(W\) to attachment to \(V\) yielded zero successes in 512 paths for the test at \(m=8,n/m=50\).
4. **Lemma 4.** Check every exact product and stated floor for \(3\le m'\le200\), and simulate 12,000 paths per activation configuration. Case A covers \(m'=3,4,10,50\). Case B covers \(m'=6,10,50\). Case B is structurally impossible at orders three and four, and also at order five. The adjacent pair of core vertices of full degree two requires a distinct other core neighbor on each side, and each of those two outer core neighbors needs at least one additional vertex. A six-vertex path attains that minimum. The exact-product checks still include algebraic values at \(m'=3,4,5\). The smallest empirical activation rate was \(0.1426\), with pointwise Wilson lower endpoint \(0.1364\), above \(1/368\). Every sampled occurrence of the full proof event produced activation. Ending the first case A block at \(2m'\) caused 196 product floor failures. Omitting the requirement \(W_\tau=u\) in the case B proof event certified 16 paths that did not activate.
5. **Parts (a) and (b).** Use 192 paths through \(n=200{,}000\), with four reported orders and each listed pair satisfying \(n\ge24k\). The simulation initializes the three-vertex path, which is the unique unlabeled outcome of the first step from one edge. For each pair use \(m_0=\lceil 2n/k\rceil\) and \(m'=\lfloor n/(8k)\rfloor\). The upper estimator averages \(C_{m_0}+(2n/k)\sum_{m=m_0}^{n-1}(C_{m+1}-C_m)/(m+1)\) path by path. The lower estimator averages \(C_{m'}/44160\). The output reports \(E[N_k(n)]\), \(E[C_{m_0}]\), \(E[C_{m'}]\), both estimated bounds, and paired gaps with normal-approximation 95% half-widths. The \(k=2\) upper comparison is an exact pathwise equality. Omitting the initial core term in (a) and multiplying by 44160 in (b) are the planted errors.
6. **Part (c) endpoint.** The additional exact-arithmetic check uses an illustrative positive sequence to isolate the algebra. The discrepancy at \(m_0=n\) is \(a_n/[n(n+1)]\) for every positive \(a_n\), including the model's sequence.

## Commands and captured output

The first command's saved output reaches the last completed Lemma 4 check and the parts (a) and (b) headers. The service restart interrupted that invocation before it printed the comparisons. The continuation reran only parts (a) and (b), then ran the constants subcommand. The boundary check had already finished.

### Original invocation (saved partial output)

```sh
taskset -c 0 env PYTHONDONTWRITEBYTECODE=1 attempts/strogatz-fable/.venv/bin/python -u checks/fable-degree-tail-gpt/referee.py all > checks/fable-degree-tail-gpt/run.txt 2>&1
```

```text
Lemma 3: exhaustive unlabeled trees by adding one leaf at every vertex
n trees core-free min(E/C) minimizers
1 1 1 - 0
2 1 1 - 0
3 1 0 1 1
4 2 0 1 2
5 3 0 2/3 1
6 6 0 2/3 1
7 11 0 2/3 2
8 23 0 3/5 2
9 47 0 1/2 1
10 106 0 3/7 1
11 235 0 3/7 2
12 551 0 3/7 4
13 1301 0 3/7 7
Global minimum 3/7; exact C/4 equalities: 0
  minimizer n=10, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8), (6, 9)]
  minimizer n=11, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 6), (3, 7), (5, 8), (6, 9), (7, 10)]
  minimizer n=11, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (4, 8), (5, 9), (6, 10)]
  minimizer n=12, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 6), (2, 7), (3, 8), (6, 9), (7, 10), (8, 11)]
  minimizer n=12, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 6), (3, 7), (5, 8), (5, 9), (6, 10), (7, 11)]
  minimizer n=12, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (5, 6), (5, 7), (6, 8), (7, 9), (8, 10), (9, 11)]
  minimizer n=12, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (4, 8), (5, 9), (5, 10), (6, 11)]
  minimizer n=13, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (1, 7), (2, 8), (3, 9), (7, 10), (8, 11), (9, 12)]
  minimizer n=13, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 6), (2, 7), (3, 8), (6, 9), (6, 10), (7, 11), (8, 12)]
  minimizer n=13, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 6), (6, 7), (6, 8), (7, 9), (8, 10), (9, 11), (10, 12)]
  minimizer n=13, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 6), (3, 7), (5, 8), (5, 9), (5, 10), (6, 11), (7, 12)]
  minimizer n=13, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 6), (3, 7), (5, 8), (5, 9), (6, 10), (6, 11), (7, 12)]
  minimizer n=13, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (5, 6), (5, 7), (6, 8), (7, 9), (8, 10), (8, 11), (9, 12)]
  minimizer n=13, E=3, C=7, edges=[(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (4, 8), (5, 9), (5, 10), (6, 11), (6, 12)]
Orders where case B is feasible in this enumeration: [6, 7, 8, 9, 10, 11, 12, 13]
Mutation caught: omit adjacent-pair clause, 7 failures; first (11, '((((((((((()))))))))))', 2, 9)
Lemma 1: 224 trees, 14688 vertex checks per moment, 28928 ordered outcomes, 83 equalities in both bounds
Mutation caught: use j+1 as growth denominator; 109 first- and 109 second-moment failures; first (j,v,d,E[D'])=(3, 0, 2, '8/3')
Lemma 2: seed=413, 768 independent paths per row, shared across ratios
m l c shape   hits/trials (estimate; pointwise Wilson 95% lower), ratios 2, 8, 50
  4 3 0 star  752/768 (0.9792; 0.9664) 706/768 (0.9193; 0.8979) 691/768 (0.8997; 0.8765)
  5 4 0 star  756/768 (0.9844; 0.9729) 728/768 (0.9479; 0.9299) 717/768 (0.9336; 0.9137)
  6 3 1 path  763/768 (0.9935; 0.9849) 720/768 (0.9375; 0.9181) 711/768 (0.9258; 0.9050)
  8 3 2 path  763/768 (0.9935; 0.9849) 732/768 (0.9531; 0.9358) 726/768 (0.9453; 0.9269)
 20 3 1 path  755/768 (0.9831; 0.9713) 708/768 (0.9219; 0.9007) 679/768 (0.8841; 0.8595)
 20 3 1 fan   753/768 (0.9805; 0.9680) 676/768 (0.8802; 0.8553) 654/768 (0.8516; 0.8247)
 50 3 2 path  764/768 (0.9948; 0.9867) 737/768 (0.9596; 0.9433) 721/768 (0.9388; 0.9196)
 50 3 2 fan   763/768 (0.9935; 0.9849) 691/768 (0.8997; 0.8765) 665/768 (0.8659; 0.8400)
200 3 1 path  755/768 (0.9831; 0.9713) 701/768 (0.9128; 0.8907) 684/768 (0.8906; 0.8666)
200 3 1 fan   751/768 (0.9779; 0.9648) 669/768 (0.8711; 0.8455) 629/768 (0.8190; 0.7902)
200 3 2 path  764/768 (0.9948; 0.9867) 733/768 (0.9544; 0.9373) 713/768 (0.9284; 0.9079)
200 3 2 fan   762/768 (0.9922; 0.9831) 679/768 (0.8841; 0.8595) 640/768 (0.8333; 0.8053)
Smallest pointwise Wilson lower bound: 0.790215; target=0.033333
Mutation caught: attach to V instead of its friend; m=8, n/m=50, 0/512
Lemma 4: exact rational products vs closed forms, all m'=3..200
  G             min=2/23         at m'=  3; claimed floor 2/23
  A1            min=200/299      at m'=200; claimed floor 2/3
  A2            min=500/799      at m'=200; claimed floor 5/8
  H1            min=1/4          at m'=  3; claimed floor 1/4
  B1            min=200/399      at m'=200; claimed floor 1/2
  B2            min=400/799      at m'=200; claimed floor 1/2
  A certificate min=90/1771      at m'=  3; claimed floor 5/138
  B certificate min=9/2530       at m'=  3; claimed floor 1/368
Mutation caught: A1 ends at 2m' instead of 3m'; 196 failures, first at m'=5 with probability 5/8
Lemma 4: seed=419, 12000 independent paths per row
case m' shape    activation/trials (estimate; pointwise Wilson 95% lower) certificate/trials
A  3 path   8377/12000 (0.6981; 0.6898)  573/12000
A  4 path   5792/12000 (0.4827; 0.4737)  628/12000
A  4 star   8426/12000 (0.7022; 0.6939)  654/12000
A 10 path   5865/12000 (0.4888; 0.4798)  597/12000
A 10 star   8015/12000 (0.6679; 0.6594)  665/12000
A 50 path   5722/12000 (0.4768; 0.4679)  615/12000
A 50 star   7890/12000 (0.6575; 0.6490)  573/12000
B  6 path   2694/12000 (0.2245; 0.2171)  51/12000
B 10 path   2669/12000 (0.2224; 0.2151)  51/12000
B 10 fan    2006/12000 (0.1672; 0.1606)  64/12000
B 50 path   2622/12000 (0.2185; 0.2112)  55/12000
B 50 fan    1711/12000 (0.1426; 0.1364)  53/12000
Mutation caught: omit W_tau=u from case B certificate; 16 unsuccessful activations among 105 claimed certificates (2684 actual activations in 12000 trials)
Parts (a),(b): seed=521, 192 independent paths, all from one edge
n k m0 m'   mean_N mean_C_m0 mean_C_m'   upper_a lower_b   paired_a_gap+-95%halfwidth paired_b_gap+-95%halfwidth
```

### Continuation (all 27 comparisons passed)

```sh
taskset -c 0 env PYTHONDONTWRITEBYTECODE=1 attempts/strogatz-fable/.venv/bin/python -u checks/fable-degree-tail-gpt/referee.py parts-ab > checks/fable-degree-tail-gpt/parts-ab.txt 2>&1
```

```text
Parts (a),(b): seed=521, 192 independent paths, all from one edge
n k m0 m'   mean_N mean_C_m0 mean_C_m'   upper_a lower_b   paired_a_gap+-95%halfwidth paired_b_gap+-95%halfwidth
   500   2    500    31    24.8906   24.8906    5.1302    24.8906 0.0001162  0.0000+-0.0000 24.8905+-1.8377
   500   3    334    20    15.6562   19.8177    4.0573    23.9511 0.0000919  8.2949+-0.6771 15.6562+-1.1854
   500   6    167    10     8.7188   13.2552    2.5938    19.9290 0.0000587  11.2103+-0.9072 8.7187+-0.6485
   500  12     84     5     5.3854    8.9948    1.6354    15.3065 0.0000370  9.9211+-0.8106 5.3854+-0.3754
  5000   2   5000   312    93.2708   93.2708   19.1302    93.2708 0.0004332  0.0000+-0.0000 93.2704+-6.9042
  5000   3   3334   208    57.6198   73.9792   14.9896    89.7072 0.0003394  32.0874+-2.4225 57.6195+-4.3231
  5000   6   1667   104    31.5208   49.6250   10.1302    74.6662 0.0002294  43.1454+-3.2907 31.5206+-2.3198
  5000  12    834    52    19.3958   33.5000    6.7604    57.2871 0.0001531  37.8913+-2.8448 19.3957+-1.4383
  5000  24    417    26    12.2448   22.5208    4.6823    42.0972 0.0001060  29.8524+-2.2611 12.2447+-0.9067
  5000  50    200    12     7.7917   14.5833    2.9635    29.5090 0.0000671  21.7173+-1.6580 7.7916+-0.5646
  5000 100    100     6     5.3333    9.9219    1.8542    20.6726 0.0000420  15.3392+-1.1948 5.3333+-0.3680
 50000   2  50000  3125   343.5052  343.5052   71.4062   343.5052 0.0016170  0.0000+-0.0000 343.5036+-26.0406
 50000   3  33334  2083   211.5938  273.1250   56.2917   330.5657 0.0012747  118.9719+-9.0626 211.5925+-16.0710
 50000   6  16667  1041   115.6719  183.4427   37.7292   275.3882 0.0008544  159.7164+-12.1359 115.6710+-8.7669
 50000  12   8334   520    70.8385  124.2604   25.4635   211.9689 0.0005766  141.1303+-10.7980 70.8380+-5.2836
 50000  24   4167   260    44.9427   84.0312   17.2552   156.3004 0.0003907  111.3577+-8.4435 44.9423+-3.3458
 50000  50   2000   125    28.5000   55.0625   11.2656   109.8152 0.0002551  81.3152+-6.1807 28.4997+-2.0873
 50000 100   1000    62    19.1354   36.7969    7.4635    77.0563 0.0001690  57.9209+-4.3634 19.1352+-1.4186
 50000 500    200    12     7.6042   14.5833    2.9635    32.5736 0.0000671  24.9695+-1.9066 7.6041+-0.5468
200000   2 200000 12500   751.4219  751.4219  156.0885   751.4219 0.0035346  0.0000+-0.0000 751.4183+-56.4211
200000   3 133334  8333   463.9948  598.4792  124.2552   723.2763 0.0028138  259.2815+-19.4478 463.9920+-34.9527
200000   6  66667  4166   254.5833  404.1615   84.0312   603.8177 0.0019029  349.2344+-26.1706 254.5814+-19.3492
200000  12  33334  2083   155.2396  273.1250   56.2917   465.4152 0.0012747  310.1756+-23.4557 155.2383+-11.7002
200000  24  16667  1041    98.6927  183.4427   37.7292   342.8130 0.0008544  244.1203+-18.6150 98.6919+-7.3150
200000  50   8000   500    62.8229  121.4948   24.8906   240.7695 0.0005636  177.9466+-13.4710 62.8224+-4.7461
200000 100   4000   250    41.4479   82.0573   16.8542   169.5525 0.0003817  128.1046+-9.6798 41.4475+-3.0875
200000 500    800    50    16.5417   32.6615    6.6302    71.8128 0.0001501  55.2711+-4.1440 16.5415+-1.2432
Mutation caught: omit a_m0 in (a), 12/27 paired comparisons fail
Mutation caught: multiply by 44160 in (b), 27/27 paired comparisons fail
```

### Continuation (coefficient checks passed)

```sh
taskset -c 0 env PYTHONDONTWRITEBYTECODE=1 attempts/strogatz-fable/.venv/bin/python -u checks/fable-degree-tail-gpt/referee.py constants > checks/fable-degree-tail-gpt/constants.txt 2>&1
```

```text
Lower corollary prefactor before rounding: 0.44121907205027
Upper corollary prefactor before rounding: 32.92507122519528
4*368*30 = 44160
```

### Boundary arithmetic (completed before the restart)

```sh
taskset -c 0 env PYTHONDONTWRITEBYTECODE=1 attempts/strogatz-fable/.venv/bin/python -u checks/fable-degree-tail-gpt/boundary.py > checks/fable-degree-tail-gpt/boundary.txt 2>&1
```

```text
Part (c): exact fractions on illustrative a_m=(m+1)/4
n=48, k=2, m0=48: written RHS minus LHS=1/192, corrected RHS minus LHS=0
n=72, k=3, m0=48: written RHS minus LHS=0, corrected RHS minus LHS=0
n=96, k=4, m0=48: written RHS minus LHS=0, corrected RHS minus LHS=0
Corrected identity passes all 4465 pairs 3 <= m0 <= n <= 96
```

The continuation completed with exit status zero. Omitting the initial term in (a) failed 12 of 27 paired comparisons. Multiplying by 44160 in (b) failed all 27.

## Re-check after repair

I re-read the current proof. The output of `git diff -- attempts/strogatz-fable/proofs/degree_tail.md` showed three changed passages.

- **Introduction, lines 9–13 — VALID.** The introduction describes the bounds in (a) and (b) in their stated ranges, with \(C_{n/(8k)}\) read as \(C_{\lfloor n/(8k)\rfloor}\). Its logarithmic comparison refers to (d) and uses the limit assumption there. The prose now explicitly excludes a two-sided comparison up to constants.
- **Paragraph after the Corollary, lines 39–43 — VALID.** The sharper exponent bounds and existence of \(\mu\) are identified as inputs. Under those inputs, (d) gives the cumulative logarithmic exponent. The paragraph distinguishes the pointwise density prediction under a power law from the cumulative result actually proved.
- **Proof of (c), lines 130–145 — VALID.** In the range \(n\ge24k\), \(m_0=\lceil2n/k\rceil=n\) for \(k=2\). For \(k\ge3\), \(2n/k\le2n/3\le n-1\), so \(m_0<n\). At \(m_0=n\) the increment sum is zero and its new upper bound in line 138 is \(a_n b_n>0\). For \(m_0<n\) the summation-by-parts equality applies, and dropping \(-a_{m_0}b_{m_0+1}\) gives the same inequality. The estimate of \(b_j-b_{j+1}\), extension of the positive sum, integral bound, and final ceiling and floor estimates yield both original constants. Both endpoint cases are covered.

The four lemmas and all parts of Theorem D are valid in the current text. The numerical Corollary follows from the supplied Gamma-ratio bounds. The earlier verdict **GAP** referred to the earlier text, which used the summation-by-parts equality at \(m_0=n\). The earlier numerical runs are recorded above.
