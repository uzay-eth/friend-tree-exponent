> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; Claims 3 and 10; Gemini reviewed the base sequential law, with separate scopes for addenda. See the [release index](../../../../../README.md#claim-20) for exact scope and citations.

# A simultaneous low-forest limit at polylogarithmic cutoffs

The polylogarithmic diagonal was suggested in the Neon-2 plan. This proof uses the explicit cylinder rate from [age-effective-sol](../../../../../../../reviews/cold/astra/checks/age-effective-sol.md) and a census of bounded downward-fringe tests to certify one. All observables below are on the actual finite graph: delete using **original** degrees, retain ordinary component shapes, and count original cores before forgetting their kinds. The orientation toward the original founder supplies the top and its deleted parent, when present. No stable-age or future mark is an observation of the finite graph.

Let `M_{n,N}` be the measure obtained by counting each surviving core-containing component once, dividing by the **physical** `C_n`, and recording `(log(D_parent/N), F, K)`, where `F` is its entire ordinary top-rooted tree and `K` its original surviving core count. Use a special parent mark for the possible founder component. Define `M_N` by the same count at the top under `π`, with no founder. Use the mass norm `||λ||=sup_{|f|≤1}|λf|`. Put

```
eta = 2^(-348),  Nmax(n)=floor((log n)^(1/4)).
```

For **any deterministic** `N(n)→∞` with `2≤N(n)≤Nmax(n)`,

```math
\left\|M_{n,N(n)}-M_{N(n)}\right\|=o\big(N(n)^{-\mu}\log N(n)\big)\quad\text{a.s.} \qquad\text{(1)}
```

The conclusion is uniform over the available cutoffs after the explicit errors below and the stationary truncation error are separated. It does not follow from taking an arbitrary diagonal of the sequential result.

## 1. Uniform physical cylinder estimate at deterministic sizes

Consider a bounded score at a core reading only its downward core descendants through depth `k`, including core degrees and leaf counts at the boundary. Root at the physical founder. This is a bounded test on a single planted descendant fringe: the root slot restores its physical degree, and no ancestor reconstruction is needed. Capped birth marks are optional input and are not returned here. The marked physical immigration identity and the explicit renewal remainder of [age-effective-sol](../../../../../../../reviews/cold/astra/checks/age-effective-sol.md) therefore give for the sum over all cores, including the one founder,

```
|e^(-mu h) E_T sum f(T_h) - Phi(T) pi(f)|
    <= K 4^k ||f|| C(T) e^(-eta h).
```

Here `Phi` is the leading expected core-count coefficient; only a finite absolute `K` is needed. The oriented `pi` is the downward-fringe law of the same spatial limit. For downward tests the estimate directly follows from the weighted mean expansion for bounded fringe scores, with a bounded founder term; the factor `4^k` is harmless. The same estimate holds for `f−pi(f)`, whose leading term vanishes.

In a one-birth dirty-core coupling, discard dirty cores and their first `k` ancestors in each graph. All remaining reads agree. At elapsed `u` the expected combined number of dirty cores is at most `3e^{2u}` (a dirty core's rate-one initiations can create one dirty core or infect one pair). A promoted leaf can change both the old parent's mark and its own core status; the initial dirty set counts both. This gives a jump bound `J_f e^{2u}`, with `J_f=6(k+1)||f||` for the backward conditional expectation at lag `u`, and a direct bound `J_f` at lag zero. Boundary reads of the full degree are part of the boundary vertex's data. If a degree convention needs one additional level, increase `k` by one. This argument never enumerates all the vertices in a large ball.

Set `m=floor(n^(99/100))`, and restart the physical tree at size `m` for elapsed `h=log(n/m)` with an independent Yule clock. Let `S` be its terminal size. The clock is independent of the jump chain from `m`; hence to compare the scores use the same jump chain at sizes `S` and `n`. Uniformly in the initial graph, the clock satisfies

```
P(|S−n|>n^(11/20)) ≤ 2 exp(−c n^(9/100)).
```

Indeed `S` is a sum of `m` independent geometric variables of parameter `m/n`, with mean `n`; after scaling each by `m/n`, its mgf is bounded in a fixed interval around zero and the centered log-mgf is at most a constant times the squared argument there. Chernoff gives an exponent of order `m n^(11/10)/n²`. In particular `S≤2n` outside that error, and the direct jump bound pays at most `J_f n^(11/20)` in the clock transfer.

For the restarted continuous-time backward martingale stop at the size boundary `2m e^s` at elapsed time `s`. The chance of reaching it by `h` is at most `exp(−c m)`: `exp(θ e^{-s}|T_s|)` is a submartingale, and the same geometric mgf at time `h` gives the maximal bound for small fixed `θ>0`. Before stopping the bracket is at most `2m J_f² e^{4h}/3`, or `O(J_f² n^(103/100))`, with jump bound `O(J_f n^(2/100))`. The bounded-jump martingale exponential inequality at deviation `J_f n^(11/20)` gives `O(exp(−c n^(7/100)))`, uniformly in `T_m`. For the centered score use its own norm in `J_f` (at most twice the original norm). The inequality follows from the exponential compensated-jump supermartingale using `e^v−1−v≤v²/[2(1−|v|/3)]` for `|v|<3`. These estimates bound failure probabilities of the original jump chain; no common choice of auxiliary clocks across `n` is needed. Centering first at `pi(f)` and dividing by the physical `C_n`, the a.s. remainder, for any family of `exp(n^{o(1)})` such tests at size `n`, is bounded simultaneously by

```
O_path( ||f|| 4^k n^(-eta/100) + J_f n^(-(mu−11/20)) ).       (2)
```

Here `C_n/n^mu→X>0` and `e^{mu h} C_m/C_n→1` a.s. were used only in the final normalization. The probability union and Borel–Cantelli precede that division. This supplies a stable cylinder bound for tests depending on `n,N`; no convergence rate for `X` was used.

## 2. A finite marked census with a quantitative exceptional mass

Use `j=3 ceil( (132 Nmax^3+100 log log n)/3 )`, depth `k=j+3`, `U=(Nmax+1)^(j+3)`, and parent cap `B=U^32`. All have the required scale: `j=o(log n)` and `U,B=n^{o(1)}`. For any `2≤N≤Nmax`, admit a nonfounder component if its deleted parent has degree `≤B` and its surviving core height from the top is less than `j`. Its entire ordinary size is at most `U` by the low degree cutoff. The admitted graph census at a high parent `o` reads only these data:

1. Read the exact degree of `o`; reject it unless `N<D_o≤B`.
2. For each low core child explore downward, stopping at high cores and rejecting the component if a low vertex is reached at core depth `j` from that child. Restore the real leaf marks of every retained core, and record the original core kinds, the tree and the relative degree of `o`. No ordinary high-core descendant occurs as a retained leaf. A low child with no surviving neighbours is allowed.

Thus this is a downward cylinder of depth at most `k`; the ordinary shape has at most `U` vertices. Recording all retained original kinds and, optionally, all original degrees of retained cores costs at most `exp(O(U log(U+N)))` types (bound by labelled encodings), times at most `B` possibilities for the parent degree. The same bound with a factor `U` covers an additionally distinguished surviving core in the component. Counting at the parent counts every admitted component exactly once; a single parent's contribution to any coordinate is at most `B`, or at most `BU` if summing over every surviving core. These bounds apply under `π` as well by the stationary mass transport between the high parent and the component's top or its surviving cores. Under `π` the unique end determines the orientation. With the union over `N≤Nmax`, (2) applies to a finite list of cylinder coordinate tests, but **do not sum (2) over all types**. Apply it instead to all signed sums of the coordinate tests (coefficients `±1`, each test evaluated at a parent): their number is at most `exp(n^{o(1)})`, with sup norm still at most `BU`. Their supremum is the mass norm of the census difference. The right side of (2) is therefore `O_path(n^{-eta/200})`, for all sufficiently large `n`, for either count and uniformly in `N`. This also controls all real bounded functions of the admitted marks simultaneously, including functions varying with `n,N`.

For the discarded mass start with the **surviving-core** measure (unnormalized, per original core). A component of top height at least `j` has diameter at least `j`, so every surviving core in it has radius at least `j/2`. Apply the radius estimate at the multiple of three `j0=3 floor(j/6)≥64N³`. Then [low-forest §2](../../low-forest/README.md) bounds this mass by `5 exp(−.08 j0)=O((log n)^-4)`. For a short component with parent degree `>B`, reverse the paths from a surviving core to its deleted parent. The first step from the high endpoint has `c_o` choices and subsequent paths of length at most `j` through low cores have at most `U` choices altogether. By mass transport and the [audited core moments](../../../../../../../reviews/cold/astra/checks/age-effective-sol.md#effective-local-convergence-versus-a-limit-without-a-modulus), the mass is at most `U π(c_o 1_{D_o>B})≤K U B^-1/8=K U^-3`. This term is also `O((log n)^-4)` because `j≥100 log log n` and `Nmax+1≥3` eventually. The component-count mass of either discarded class is bounded by its surviving-core mass; no estimate pretends that the two normalizations are the same.

To bound the same discarded mass in the finite tree, use the **total surviving-core mass** and subtract the admitted surviving-core census. The total is the radius-zero empirical score `1{D≤N}` on original cores, covered by (2). There is no need to prove a finite-tree radius bound from scratch or assume a rate for an uncapped observable. The remainder is at most the stationary discarded survivor mass plus `O_path(n^-eta/200)`. This argument includes the finite founder component (which can have more than one surviving core); it is not treated as a one-count exception in the **survivor** measure. In the component-count measure its mass is at most the discarded survivor mass as well. In both measures the norm error is consequently

```
O_path(n^(-eta/200)) + O((log n)^(-4)), uniformly for 2≤N≤Nmax(n).    (3)
```

Since `d_N=M_N(1)∼A N^-mu log N`, (3) is `o(d_N)` on the stated diagonal. In fact the same proof permits `Nmax=floor((log n)^α)` for any fixed `0<α<1/3`, since `j log(Nmax+1)=o(log n)` and the discarded mass still beats `d_N`. Floors of the cutoff do not affect the conclusion. This proof uses an explicit radius modulus, not a diagonal extracted from fixed-`N` convergence. It does not prove cutoffs at positive powers of `n`; the census and tail budget would need a different argument there.

## 3. Finite-graph consequences and sampling conventions

* **Uniform component.** The actual number of core-containing components is asymptotic a.s. to `A C_n N^-mu log N`, equivalently to `AX n^mu N^-mu log N`. The conditional law given the actual tree is within `o(1)` in TV of `M_N/d_N`, jointly for the **exact lattice parent mark**, entire component and original `K`. The times and marks in the stationary construction are only proof devices. By [forest-age](../README.md), [its marked version](../marks/README.md), and [inert](../inert/README.md), this law is within `o(1)` of the product of a Pareto(`mu`) relative parent degree and the inert mixture at uniform age in `[0,log N]`, for bounded tests uniformly 1-Lipschitz in log relative degree and arbitrary in the entire component and `K`. In that mixture `K` counts cores with one restored root slot. The no-further-cut success in the stationary coupling preserves all original kinds, so adding `K` does not require pretending it is a function of the bare observed shape on bad histories. The exact parent lattice prevents TV to the continuous Pareto target; the shape and `K` projection alone has TV approximation. In particular `(D_parent/N, log K/log N, log |F|/log N)` tends weakly to `(R,mu V,V)`, where `R` is independent of the single uniform `V`. The founder mark has vanishing conditional mass. Bounds (1)–(3) quantify the finite-graph transfer, not the error of the separate stationary large-`N` approximation.
* **Uniform surviving core.** Count a component once for every surviving original core in it, optionally distinguishing that core as the observation root. The same (3) gives a TV comparison with the **stationary cutoff-dependent survivor law**, normalized by `π(D≤N)→1`. Retaining the orientation/top when reporting the parent mark is allowed. This law is not the Pareto–inert mixture above: core-size bias cannot be passed through the latter's TV error without weighted control. For local graph tests the survivor law tends to the full decorated `π`, since any fixed ball is eventually uncut under `π`; thus the finite diagonal has that same local limit. Entire components are compared to a moving stationary law, not to a finite limiting component type.
* **Uniform surviving ordinary vertex.** Most are ordinary isolated leaves, as the mass of all vertices in core-containing components is at most `(N+1)C_n= o(n)` along this schedule, while at least `n−C_n` ordinary leaves survive. This sampling law differs from both preceding laws.

The graph census supports bounded marks from the full observed component, including unrooted projections and original core counts. It does not introduce promotion ages into the empirical observations or assert a simultaneous asymptotic for uniformly chosen sibling tuples or largest components.

## Exact finite control

`python3 -B round2/open/forest-age/diagonal/check.py` compares a direct connected-component search with the parent census on growing labelled trees. It checks the exact rational parent ratio, the shape with original core flags, the survivor weighting, subtraction of the residual from the total surviving-core count, and the one-birth bound over the full signed census. Output:

```text
3924 exact parent/survivor censuses; 2879 nonempty residuals; 516 reversed-edge controls; 2516 nonzero one-birth changes
```

The reversed-edge controls are cases where an indiscriminate scan of all low neighbours of a high parent would inspect a low **ancestor**. The replay checks finite identities, not the stationary tail estimate or the rate in (2).
