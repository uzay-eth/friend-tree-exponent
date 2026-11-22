VERDICT: CORRECT

Theorem B transfers the drift bound to an almost-sure lower exponent. The
tracked-degree certificate meets its hypotheses, and exact enumeration on the
rebuilt table found no drift violation. The reciprocal and random start need
only the local justifications below.

## 1. Theorem B

Line numbers refer to `attempts/strogatz-fable/proofs/as_lower.md`.

### Lemma 1, lines 59–69

**Valid as written with the standing $n\ge3$ convention.** Degrees never
decrease, and $C_n$ increases by one exactly when a leaf is the attachment
target. For $n\ge\max(3,M^2+M)$ with $C_n\le M$, each of the at least $M^2$
leaves has a core neighbour, so some core $u$ has $L_u\ge M$ and $c_u\le M-1$.
Its conversion probability is $L_u/[n(L_u+c_u)]\ge M/[n(2M-1)]\ge 1/(2n)$.
The step-$n$ conversion is in $\mathcal F_{n+1}$ and its probability is
conditional on $\mathcal F_n$, the indexing needed for conditional
Borel–Cantelli. On the event that $C_n$ stays at most $M$, the harmonic lower
bound forces infinitely many conversions, contradicting that bound. A
countable union over $M$ and monotonicity give $C_n\to\infty$ almost surely.
The maximum in the threshold above makes $n\ge3$ explicit for $M=1$.
[Durrett, Theorem 4.3.4](https://sites.math.duke.edu/~rtd/PTE/PTE5_011119.pdf#page=233)

### Lemma 2, lines 71–86

- **Doob decomposition, lines 73–76. Valid as written.** The potential is
  integrable since $0\le\Phi_n\le nF$. Each increment of $A$ is
  $\mathbf E_n\Delta_n\ge\delta\Phi_n/n$. The increment
  $\Delta_n-\mathbf E_n\Delta_n$ of $M$ has conditional mean zero and absolute
  value at most $2B$.
- **Dichotomy, lines 78–82. Valid as written.** The bound $2B$ is deterministic,
  as the bounded-increment dichotomy requires. On
  $E=\{\sup_n\Phi_n<\infty\}$, nonnegativity of $\Phi_3$ and $A_n$ gives
  $M_n\le\sup_n\Phi_n$. The dichotomy for $M$ holds on the whole probability
  space, and its finite-limit alternative holds almost surely on $E$.
  Then $A_n=\Phi_n-\Phi_3-M_n$ is bounded above and increasing there, so
  $A_n$, and then $\Phi_n$, converge finitely.
  [Durrett, Theorems 4.3.1–4.3.2](https://sites.math.duke.edu/~rtd/PTE/PTE5_011119.pdf#page=232)
- **$\Phi_\infty=0$ and the contradiction, lines 82–86. Valid as written.**
  The finite $A_\infty$ bounds $\delta\sum_k\Phi_k/k$. Any positive limit of
  $\Phi_k$ would give a positive lower bound for its tail and make the series
  diverge. Nonnegativity leaves zero as the only possible limit. On the
  probability-one event of infinitely many conversions, (H4) gives a
  subsequence with $\Phi_{n_i+1}\ge\eta$, which rules out convergence to zero.
  Thus $E$ has probability zero.

### Lemma 3, lines 88–118

- **Reciprocal identity and expansion, lines 94–102. Valid on the stated
  event.** With $P=\Phi_n\ge2B>0$ and $\Delta=\Delta_n$, we have
  $P+\Delta\ge P/2$. Multiplying the identity by $P^2(P+\Delta)$ gives $P^2$
  on both sides. Also $\Delta^2\le B|\Delta|=B(\Delta+2\Delta^-)$.
  Substitution gives coefficient $-(1-2B/P)/P^2$ on $\mathbf E_n\Delta$
  and $4B/P^3$ on $\mathbf E_n\Delta^-$, as in (1). To define the
  conditional expectation globally, set $1/0=+\infty$ and use nonnegative
  conditional expectation. On the stated event the reciprocal is finite
  for every one-step outcome.
- **Drift substitution, lines 104–105. Valid as written.** Since
  $1-2B/P\ge0$, replacing $\mathbf E_n\Delta$ by its lower bound
  $\delta P/n$ increases the right-hand side.
- **Union bound, lines 107–111. Valid as written.** A negative total change
  requires a strictly decreasing old contribution, which was positive
  before the step. There are at most $P/f_{\min}$ such candidates from the
  old tree. The transition law given $\mathcal F_n$ is the one given $T_n$.
  Hypothesis (H5) gives $\mathbb P_n(\Delta<0)\le RP/(nf_{\min})$ by the
  conditional union bound. Multiplication by $B$ bounds
  $\mathbf E_n\Delta^-$.
- **Constants, lines 89 and 113–118. Valid as written.** The two corrections
  are $2B\delta/(nP^2)$ and $4B^2R/(nf_{\min}P^2)$. This gives
  $\kappa=2B\delta+4B^2R/f_{\min}$. Taking
  $K_0=\max(2B,\kappa/\varepsilon)$ preserves the expansion bound and gives
  $\kappa/P\le\varepsilon$. Here $K_0$ is the stopping threshold.

### Stopping argument, lines 120–145

- **Start and stop, lines 120–129. Valid as written.** Both times are stopping
  times since $\Phi$ is adapted, and Lemma 2 ensures $\sigma<\infty$ almost
  surely. Since $K_1>K_0$, a finite $\tau$ is strictly greater than $\sigma$.
  The preceding value is at least $K_0$, so $\Phi_\tau\ge K_0-B>0$. This
  implicit bound keeps the denominator positive at the crossing step.
  With $a=\delta-\varepsilon$, all product factors are in $(0,1)$ and
  $p_\sigma=1$.
- **Supermartingale, lines 130–133. Valid, with the random-time step
  compressed.** On $\{\sigma\le n<\tau\}$, $\min(n+1,\tau)=n+1$, including
  a crossing on the next step. The factor $p_{n+1}=p_n(1-a/n)$ is already
  $\mathcal F_n$-measurable there. Lemma 3 gives the written conditional
  inequality. The value stays constant on $\{\tau\le n\}$. The shifted
  process is adapted, and
  $0<G_{\sigma+k}\le (K_0-B)^{-1}(1-a/3)^{-k}$ since its denominator
  product has at most $k$ factors. This supplies integrability. For
  $A\in\mathcal F_{\sigma+k}$, split into
  $A\cap\{\sigma+k=n\}\in\mathcal F_n$, apply the deterministic-time
  inequality to each piece and sum. This proves the claimed
  supermartingale inequality in the shifted filtration.
- **Maximal inequality, lines 135–137. Valid as written.** Its right-hand
  side is $K_0\mathbf E[G_\sigma]\le K_0/K_1$. On $\{\tau<\infty\}$,
  $G_\tau>1/K_0$ since $0<\Phi_\tau<K_0$ and $p_\tau\le1$.
- **Growth and limits, lines 138–145. Valid as written.** Nonnegative
  supermartingale convergence gives a finite limit and hence finite
  pathwise supremum $G^*$. For each realized finite $\sigma$ and
  $n\ge\sigma$, $\log(1-a/j)\le-a/j$ and
  $\sum_{j=\sigma}^{n-1}1/j\ge\int_\sigma^n dt/t=\log(n/\sigma)$, giving
  $p_n\le(\sigma/n)^a$. On $\{\tau=\infty\}$, the resulting bound
  $\Phi_n\ge(n/\sigma)^a/G^*$ gives the stated liminf. With the convention
  $\log 0=-\infty$, the liminf event is defined on every path and is
  independent of $K_1$. Its probability is at least $1-K_0/K_1$ for
  arbitrarily large $K_1$, so it has probability one. Finally
  $C_n\ge\Phi_n/F$ and a countable sequence $\varepsilon\downarrow0$
  give the claimed liminf endpoint.

These details require no change to Theorem B or its hypotheses.

## 2. Tracked-degree hypotheses

The parameters are $L_{\max}=60$, $c_{\max}=10$, pendant cutoff 2 and degree
cap 6. The [independent checker](fable-as-lower-gpt/check.py) derives the
eight state coordinates from adjacency lists. It recomputes the post-step
states by adding an edge to each possible target $W$, then enumerates all
$2(n-1)$ ordered pairs $(V,W)$ with exact probability $1/(nD_V)$ apiece
and checks that their probabilities sum to 1. It reuses the post-step
state of a target across different sources. I used 24 simulated histories
from a fixed PRNG seed, with checkpoints through $n=150$. Targeted stars
and subdivided stars exercise the leaf cap, the $c=10/11$ boundary and a
degree-5 target surrounded by five core vertices.

- **(H3), $b_0=6$. Holds.** On a leaf gain the target's $L$ changes and
  only its core neighbours can see its class change. A top target stays
  top. Every remaining target has at most five core neighbours, for at
  most six changed old states. On a conversion the sponsor's total degree
  is fixed. If its old $c\ge2$, it had a leaf and was other or top, and its
  new $c\ge3$ leaves its class unchanged. If its old $c=1$, at most its
  sole old core neighbour also sees a change. If $c=0$, there are no old
  core neighbours. The last case is implicit in the write-up. Only the
  selected old leaf becomes core, so the newborn count is at most one.
  All tested moves obeyed the stronger bound on raw-state changes. The
  check rejected $b_0=2$.
- **(H4). Holds.** Immediately after a conversion the selected old leaf
  has one new leaf and one core neighbour, its sponsor. The sponsor has
  one of the seven classes $P_1,P_2,B,O_3,O_4,O_5,\mathrm{top}$. The
  possible newborn states have one of the six counters equal to 1,
  or all six zero for a top sponsor, always with $L=c=1$. This list is
  exhaustive for all trees. Every state occurred among the tested moves
  and its exact rebuilt table value was positive. The transcript gives
  the minimum $\eta$.
- **(H5), $R=126$. Holds.** For a vertex with $c_v\le10$, a change of the
  kind counted can come from its own leaf gain below the cap, its own
  conversion, a leaf gain by an old core neighbour with degree at most
  5, or an old core neighbour's conversion. Their total rates in units
  of $1/n$ are bounded by $60+10/2$, $1$, $10\cdot5$ and $10\cdot1$,
  respectively, which sum to 126. The neighbour conversion bound counts
  all such events, including ones leaving that neighbour's class
  unchanged. The checker measured the exact mass of raw-state changes at
  every eligible vertex except pure increments of $L\ge60$, where it also
  checked equality of the table values. This bounds the event in (H5),
  since positive contributions have $c\le10$. The check rejected $R=1$.

## 3. (H2) and the rebuilt table

I built `certs9.Model(60,10,2,6)` and called `certs9.lower_table(M, 0.5433)`.
The checker uses those two entry points for the table. Its potential
calculation uses `Fraction` to represent table entries as exact binary
rationals. It calculates $\Phi(T)$ and each full post-step
$\Phi(T')-\Phi(T)$ and totals the exact probability-weighted drift. Every
tested tree satisfied $n\mathbf E[\Delta_n\mid T_n]\ge(543/1000)\Phi(T_n)$.
The full table also passed the sign and support checks for (H1).

For the universal implication of the earlier certificate, I traced the
drift identity and the verifier against the tree rules. Each old core
contribution is charged its exact own and neighbour events, and each birth
is charged once to its sponsor. The tracked degrees give exact credit
except for top neighbours, whose individual credit is in $[0,1/6]$. A
tracked other neighbour of degree $d$ gains a leaf at rate at most $d$,
and conversion preserves its degree and other class. The pendant and
class-B transitions also agree with the event list. The rate intervals
in the verifier cover the true rates. Bounding their linear terms
individually is valid even when the true rates are dependent. At $L>60$,
the own gain term vanishes, conversion uses the capped entry and the
newborn sees a top sponsor. The ranges for $L/D$ and $1/D$ cover the
whole tail. At $c>10$, $c$ cannot decrease, so the old contribution
stays zero and only the nonnegative newborn term remains in that local
summand. The supplied `outputs/as_lower_check.log` records 1,205,767
exact checks. I re-tested (H2) at the tree level and used the supplied
certificate verification for the universal quantifier over states.

## 4. Theorem A and literature

Theorem A follows by applying Theorem B to the certified drift at
$\delta=543/1000$ and the other verified hypotheses. Remark 1 should
state that future certificates also need (H3) and (H4), or arguments
that supply their roles. Both hypotheses enter the proof.

The comparison with Addario-Berry et al. is fair at the level of
liminf exponents. The arXiv and published versions of Proposition 6.20
both concern $X_n^{\ge2}=C_n$ in the same friend tree and give
$n^{-s}C_n\to\infty$ almost surely for every $0<s<(\sqrt{13}-3)/6$.
Taking a countable sequence of such $s$ increasing to the endpoint
gives the almost-sure liminf lower bound $(\sqrt{13}-3)/6\approx0.1009$.
The proposition's divergence claim is for strictly smaller $s$.
The new $0.543$ is a strict improvement in the same liminf formulation.
I found no stronger published almost-sure lower exponent in the sources
checked.
[arXiv v1](https://arxiv.org/html/2403.20185v1)
[Published article, Proposition 6.20](https://onlinelibrary.wiley.com/doi/full/10.1002/rsa.70082)

## Commands and output

I read the proof with `nl -ba attempts/strogatz-fable/proofs/as_lower.md`,
the background with `nl -ba attempts/x-strogatz/REPORT.md | sed -n '104,315p'`,
and the code and supplied log with `nl -ba attempts/x-strogatz/certs9.py`,
`nl -ba attempts/strogatz-fable/as_lower_check.py`, and
`cat attempts/strogatz-fable/outputs/as_lower_check.log`. I checked the
Durrett theorems and both versions of Proposition 6.20 at the sources
linked above.

Run from `.`, with the process pinned to CPU 0 and the
numerical libraries restricted to one thread.

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 taskset -c 0 attempts/strogatz-fable/.venv/bin/python -B -u checks/fable-as-lower-gpt/check.py
```

The first run exited 1 with `AssertionError` at `assert degree >= 3`
in the checker. I had tried to classify the sole core of the three-vertex
star as seen by a core neighbour, although it has no such neighbour.
The checker now skips class assignment for any vertex with $c=0$ and still
computes its own state. The same command then exited 0 with this output.

```text
random seed=314159; 24 histories; sample sizes=(3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 40, 60, 61, 62, 80, 100, 125, 150)
table: 1186328 entries, Perron value 1.000032724
H1: F=2878110131726845/4503599627370496, fmin=8326880206819883/2305843009213693952, zero entries=18
H4: f(1, 1, 0, 0, 0, 0, 0, 0) = 808493090670503/2251799813685248
H4: f(1, 1, 0, 0, 0, 0, 0, 1) = 818168068835477/2251799813685248
H4: f(1, 1, 0, 0, 0, 0, 1, 0) = 6648278195080043/18014398509481984
H4: f(1, 1, 0, 0, 0, 1, 0, 0) = 3399918893151513/9007199254740992
H4: f(1, 1, 0, 0, 1, 0, 0, 0) = 3601309680998547/9007199254740992
H4: f(1, 1, 0, 1, 0, 0, 0, 0) = 3409606259239179/9007199254740992
H4: f(1, 1, 1, 0, 0, 0, 0, 0) = 7144730011107121/18014398509481984
sample: 432 random and 18 targeted trees; 38692 directed moves
coverage: {'excluded_capped_growth_targets': 49, 'trees_with_c_above_10': 3}
H3: maximum changed=6 (subdivided-star-c5-L0, (V,W,changed vertices)=(1, 0, [0, 1, 2, 3, 4, 5])); at most one newborn for every move
H4: all 7 possible newborn states realized; eta=808493090670503/2251799813685248
H5: maximum n*p=80 (subdivided-star-c10-L59, v=0)
H2: all exact tree drifts passed; minimum slack=814400475490811/4503599627370496000 (random-00-n3)
negative control rejected: subdivided-star-c5-L0, (V,W,changed vertices)=(1, 0, [0, 1, 2, 3, 4, 5]): 6 > 2
negative control rejected: subdivided-star-c10-L59, v=0: 80 > 1
```
