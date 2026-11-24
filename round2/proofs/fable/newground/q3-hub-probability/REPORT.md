# Random friend tree: hub probability versus birth label

**Author:** Fable's GPT-6.2 Sol Ultrafast proof author. **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** Astra (GPT) and Neon (Gemini) audited the original range, and the Claude lead supplied a separate enumeration. The strongest scale conclusion uses the mean core-count theorem; the first certificate also uses the Addario-Berry et al. edge theorem. See [F3](../../README.md#f3).

Write `p_u = P(Z_u>0)`, and take `F_m` to be the history through size `m`.

## PROVED — finite-state bounds

Use the rate-one clocks of the continuous-time embedding, with one fresh clock at each birth. The next ringing vertex is uniform, and its mark picks a uniform neighbour. Total population is a Yule process; from any finite starting state `e^{-t} N(t)` converges a.s. to a finite strictly positive limit (the sum of independent exponential Yule-family limits). Time in the following lemmas is measured from the conditioning state. For a core vertex write `l=L_v`, `c=D_v-l`. The leaf-count birth rate is `l + sum_{core w~v} 1/D_w`; the loss rate is `l/(l+c)`. Only a loss increases `c`.

**Lower classifier.** From `l>0`, compare the leaf count below by a birth-death chain `Y` with rates `k`, `k/(k+c)` at `k>0`, with `c` frozen at its starting value and zero absorbing. Couple at equality by putting all lower-chain births among true births and all true losses among lower-chain deaths; off equality unit jumps cannot cross. Nonexplosion follows from the pure-birth upper bound for `Y`, and from the Yule population for the true process. With `a_0=1`, `a_k=1/[(c+1)...(c+k)]`, nonabsorption has probability

```
B(l,c) = (sum_{k=0}^{l-1} a_k) / (sum_{k=0}^infty a_k).
```

This is the usual finite gambler's-ruin recursion on `[0,H]` followed by `H→infty`; staying in any finite positive band forever has probability zero, by the positive probability of a consecutive sequence of downward jumps. On nonabsorption `Y(t)` is unbounded. At the first time when `Y>=K`, the true leaf count is at least `K`. Alternatively, to see directly that nonabsorption implies `Z_v>0`, use the following expectation argument on the **true** leaf count: at any size `m'` with `L_v>=K`, the Cauchy–Schwarz bound below alone is insufficient if `D_v` is much larger than `L_v`. Instead use the simpler comparison `Y_1` with births `k` and deaths 1 on positive integers, zero absorbing. The scaled `e^{-t}Y_1(t)` is a nonnegative supermartingale with bounded second moment (dominate by pure Yule of initial size `k`). Its limit has mean at least `k-1` by the drift bound for `e^{-t}(Y_1-1)` (also at zero), and second moment at most `k(k+1)`. Thus the probability of strictly positive limit from `k` is at least `(k-1)^2/[k(k+1)]`. From a true leaf count at least `K`, couple a fresh `Y_1` from `K`; a positive scaled limit gives `Z_v>0` since `e^{-t}N(t)` has a finite limit. Hence on the event that the original chain reaches `K` the unconditional probability mass with `Z_v=0` is at most `1-(K-1)^2/[K(K+1)]`; let `K→infty`. This proves `p_v(F_m)≥B(l,c)` and also justifies the zero-bound convention `B(0,c)=0`. For rational computation truncate at index 24 and upper-bound the omitted tail by `a_24/(c+24)`.

**Upper classifiers.** A hub that is currently a leaf must be promoted. Fix a marked leaf at its parent `w`, whose current `k=L_w` includes that leaf. If `c_w≥1`, the probability of any future promotion of the marked leaf is at most `1/k`. For a proof use `f(k)=1/k`, with terminal payoff 1 upon marked promotion. At `k≥2` the birth contributes at most `k(f(k+1)-f(k))=-1/(k+1)` and the combined promotion contribution is exactly `1/D_w≤1/(k+1)`; at `k=1` the payoff for immediate success is already 1 and birth drift is negative. Stop and apply the nonnegative supermartingale bound. For `c_w=0` a simpler upper bound is `1/max(1,k-1)` for the time until the marked leaf is promoted: the same drift test with `f(k)=1/(k-1)` for `k≥2` gives at most `-1/(k-1)+1/k≤0` at `k≥3` (and a trivial promotion bound at `k=2`). In fact at `k≥3` the combined promotion term for this `f` is `(1-f(k))/k + (k-1)(f(k-1)-f(k))/k`, which is at most `1/(k-1)`; the birth term is exactly `-1/(k-1)`.

Sharper bound `F_c(k)` for starting with at least `c≥1` core neighbours: freeze `c`, use actual leaf birth rate lower bound `k`, and assign rates `1/(k+c)` for marked promotion and `(k-1)/(k+c)` for unmarked promotion. On `1..K` solve

```
(k(k+c)+k) F_c(k) = k(k+c) F_c(k+1) + (k-1) F_c(k-1) + 1,
F_c(K+1)=1/(K+1).
```

`exact.py` solves this at `K=60` by rational elimination and verifies monotonicity `F_c(k)≥F_c(k+1)` exactly for every table. This makes it a supersolution for the actual generator: extra births reduce it, and the nonnegative promotion increment is multiplied by a smaller rate if the true core degree is larger. Evaluation can additionally take the minimum with `1/k`. If `c_w=0`, the whole tree is a star; the first unmarked promotion exits the star and goes to the bound with `c=1`; marked promotion succeeds; a birth moves to the next larger star. Backward substitution of these first-step probabilities with boundary `1/40` at `k=41` is `uleaf` in `exact.py`.

For a core with `l=0`, to become a hub it first has to gain a leaf. At any neighbour `w` of `v`, until `v` gains a leaf its first selection of `v` is bounded by `1/max(1,L_w-1)` using the same supermartingale as for a marked leaf but with the marked target now **core**: before success, there are `k` unmarked losses at rate `1/D_w` each and one success at rate `1/D_w`, with `D_w≥k+1`. At `k≥3` with `f(k)=1/(k-1)`, birth drift is at most `-1/(k-1)` and the combined positive term is `(1-f(k)+k(f(k-1)-f(k)))/D_w = (k^2-3k+4)/[(k-1)(k-2)D_w] ≤ 1/(k-1)` for **`k≥4`**. The formula at `k=3` is also valid: it is `2/D_w≤1/2` since `D_w≥4`. At `k≤2` the bound is trivial. Union over `w`, clipping at 1. Finally the edge theorem gives the additional lower bound `1-upper(w)` for every neighbour `w` of `v`.

**Moments and positive-threshold bounds.** Use `b=max(0,L_v-1)/m`, `a=D_v/m`, so `b≤E[Z_v|F_m]≤a`. The one-step increment bound for `D_v(D_v+1)` yields `E[Z_v²|F_m]≤s=D_v(D_v+1)/[m(m+1)]` (multiply by `(j+2)/j` conditionally at step `j`). Cauchy–Schwarz gives `P(Z_v>t|F_m)≥(b-t)_+²/s` for `t≥0`; Markov gives `≤a/t` for `t>0`. Also `E[(Z_v-a)^2|F_m]≤s-2ab+a²=:M`, so if `t<a` the lower bound can be increased to `1-M/(a-t)^2` and if `t>a` the upper bound decreased to `M/(t-a)^2`. The `t=0` upper/lower classifiers generally beat moments dramatically.

## PROVED — symmetry and first vertices

Given a rooted (unlabelled) finite tree, the root's entire future `Z` law depends only on this rooted isomorphism class. Consequently all the following linear relations hold for **laws**, or for their evaluations on any common measurable subset of `[0,1]`.

* At size 3 let `C,L` denote centre and endpoint laws of the path. Then `law(Z_1)=law(Z_2)=(C+L)/2` and `law(Z_3)=L`. A centre–endpoint comparison is not supplied by symmetry.
* At size 4 let `S_c,S_l` denote star centre/leaf; let `I,E` denote path interior/end. The star occurs with probability `2/3`; conditional on a path, either old endpoint was selected with probability `1/2`. Thus `C=2 S_c/3+I/3`, `L=2 S_l/3+(I+E)/6`, and `law(Z_4)=2 S_l/3+E/3`. In particular the comparison `3 versus 4` is exactly `I versus E`.
* At size 5 put `A_l` for the 5-star leaf; `B_i,B_o,B_h` for, respectively, the degree-two vertex, its leaf and either of the two short-arm leaves on the (3,2,1,1,1)-degree tree; `B_c` denotes its degree-three vertex. On the 5-path write `Q_m,Q_i,Q_e` (middle, next-to-end, end). The transitions of the 4-star are to the 5-star with weight `3/4` or to `B` with weight `1/4`. The 4-path goes to `B` with weight `3/4` or to the 5-path with weight `1/4`. Hence `law(Z_5)= A_l/2+B_o/6+B_h/4+Q_e/12`. The statically equivalent role formulas are `S_c=3A_c/4+B_c/4`, `S_l=3A_l/4+B_i/12+B_h/6`, `I=3(B_c+B_i)/8+(Q_m+Q_i)/8`, `E=3(B_h+B_o)/8+(Q_i+Q_e)/8`. Beyond this, root-isomorphism grouping is more convenient.

**Two short strict comparisons at the hub atom.** On the 3-path the centre has lower bound `B(2,0)=2/e>0.735`, while the endpoint has upper bound `11/18` by three first-step cases: any birth at the centre (probability `2/3`) leaves a star of three leaves with bound `1/2`; the marked promotion has probability `1/6`; the unmarked promotion leaves one marked leaf and one core with bound `F_1(1)≤2/3` (use cutoff `K=1`). The resulting upper bound is `1/3+1/6+1/9=11/18`. Thus `p_1=p_2>p_3` without enumeration.

On the 4-path use a **one jump** calculation. The two interior attachment sites each have probability `3/8` and the endpoint sites each have probability `1/8`. When attachment is at the marked interior root, it has `(l,c)=(2,1)`. At the opposite interior or at the far endpoint the marked root has `(1,1)`. At its own adjacent endpoint the marked root loses its only leaf; use 0. Hence the interior lower bound is `3 B(2,1)/8+B(1,1)/2 = 17/[16(e-1)]>0.618`. For the marked endpoint upper: at its own interior use `F_1(2)`; at the other interior use `F_1(1)`; at its own endpoint use 1; at the opposite endpoint use `F_1(1)`. This gives `3 F_1(2)/8+F_1(1)/2+1/8 = 239/392 <0.610` at cutoff `K=3`, where `F_1(2)=22/49`, `F_1(1)=31/49` (the bound `e<87/32` suffices for the strict comparison). Thus `p_3>p_4`. The size-5 conditional enumeration in `orbit.log` gives the slightly sharper endpoint upper `146190742515/240000000000`.

**Small-label result via exact enumeration: `p_1=p_2>p_3>p_4>p_5>p_6>p_7>p_8`.** For each size `n`, each possible attachment site gives one distinct parent history extension. The integer weight numerator is `sum_{x~w} lcm(1,...,n-1)/D_x`; all histories have common denominator `prod_{n=2}^{m-1} n*lcm(1,...,n-1)`. At `m=9,10,11` the checked counts were 40,320, 362,880 and 3,628,800 histories and 286, 719 and 1,842 rooted shapes. The Python rational bounds are rounded outward on grid `10^-10` before the exact integer calculation (`__int128`, sanitizer checked). For differences, group histories by marked rooted unlabelled tree, sort the balanced-parentheses child encodings to obtain a canonical key, and cancel the two exact rooted masses before multiplying by bounds: `p_i-p_j = sum_s (q_i(s)-q_j(s))p(s)`. Use a lower endpoint for a positive coefficient and an upper endpoint for a negative coefficient (reverse for upper differences). This step matters: separate marginal intervals overlap from label 3 onward. Size 9 already suffices for every adjacent comparison through 8. At size 11 the intervals below are rounded outward (full numerators, denominators, and all pairings in `exact.log`):

| Label | Hub lower | Hub upper | Mean lower | Mean upper |
|---:|---:|---:|---:|---:|
|1, 2|.586344|.744183|.230418|.393950|
|3|.264247|.489775|.043874|.178002|
|4|.178031|.418222|.019765|.138821|
|5|.127678|.375341|.010057|.120794|
|6|.092930|.345899|.005209|.110428|
|7|.066696|.324104|.002563|.103718|
|8|.045716|.307159|.001090|.099039|

| Adjacent difference | Hub lower | Hub upper |
|---:|---:|---:|
|2−3|.246544|.329961|
|3−4|.043227|.114543|
|4−5|.021662|.071572|
|5−6|.012190|.051999|
|6−7|.007086|.040943|
|7−8|.003997|.033929|

The same grouped computation also proves strict ordering of the means through label 5 (adjacent lower bounds `.1499925, .0106009, .0021528`); subsequent adjacent mean comparisons are unresolved by this certificate.

## PROVED — decay and aggregate scale

For `u≥4` the newborn is a leaf attached to the selected parent, which is core *after* birth. In the conditioning tree at `u-1`, a core attachment site `w` has rate `L_w+rho_w` where `rho_w=sum_{core x~w}1/D_x`. If selected, its post-birth leaf count is `L_w+1`. The marked promotion bound is at most `1/(L_w+1)` for `c_w≥1`, and at most `1/L_w` in the star case. Summing gives a contribution at most `C_{u-1}+sum_w rho_w` in units `1/(u-1)`. Attachments to leaves have total rate `sum_w L_w/D_w`; bound success by 1 in this case. Since `sum_w rho_w=sum_w c_w/D_w`, the credit and leaf-attachment terms sum to `C_{u-1}`. Thus `p_u≤2 E[C_{u-1}]/(u-1)`.

Conversely, condition on the newborn's parent `w` having degree `d`. Suppress `w`'s rings and run the rest of the process; this agrees with the true tree until the first ring of `w`, and the newborn stays a leaf up to that ring. The suppressed process satisfies `E[D_w(s)]≤d e^s` by its generator (all neighbour rates sum to at most current degree). The first parent ring is independent of this suppressed process, with density `e^{-s}`. Uniform choice at that ring selects the newborn with probability at least `∫_0^∞ e^{-s} (d e^s)^-1 ds = 1/(2d)` by Jensen. At promotion the newborn has one core neighbour and one leaf, so its chance to become a hub is at least `B(1,1)=1/(e-1)`. It remains to lower-bound `sum_w sum_{x~w} 1/[D_x(D_w+1)]` in the size-`u-1` tree, coming from averaging the attachment probability of `w` times the reciprocal of its post-birth degree. The following bound still holds if we count only core attachment sites `w`:

* For a core tree on `C≥2` vertices, write `A=#(c=1)` and `B=#(c≥3)`. Then `A≥B+2`. Every vertex of `c≤2` with `L≥1` contributes by `L/(D+1)≥1/4`; call this set `J`, so `A≤|J|` and `B≤|J|`.
* The remaining set `H` of `(c,L)=(2,0)` has both neighbours in the core, and `|H|≥C-2|J|`. Since the sum of core degrees outside `H` is `2(C-1)-2|H|`, at least `4|H|-2C+2` directed edges run from `H` to `H`. Each such directed edge contributes `1/(2*3)=1/6`. The combined contributions are at least `|J|/4 + max(0,2C-8|J|)/6 ≥ C/16`, by splitting at `|J|=C/4`. The case `C=1` is a star with `L/(D+1)≥2/3`.

It follows that, for every `u≥4`,

```
 E[C_{u-1}]/[32(e-1)(u-1)] ≤ p_u ≤ 2 E[C_{u-1}]/(u-1).
```

In particular `p_u≥1/[32(e-1)(u-1)]`. The existing mean-exponent proof and certified bracket in `proofs/existence.md` and `bracket-motifs.md` give `p_u=u^{mu-1+o(1)}` with `mu` the **mean** core exponent (bracket `[.5650,.5665]`); the a.s. exponent alone should not be substituted into an expectation without justification. Even just the a.s. sublinear core bound gives `p_u→0` by dominated convergence. These bounds do not give adjacent inequalities at large labels.

**Degree dichotomy and consistency of classifiers.** On `{Z_v=0}` the leaf count is bounded almost surely: stop on first reaching `K` and apply `B(K,0)→1`, then let `K→∞`. If a core vertex has a leaf at arbitrarily late times, it reaches every fixed `K` almost surely on this event. To see this, try at successive stopping times separated by at least one time unit, each with a leaf present. In the subsequent unit, the probability that the core has no rings and the chosen leaf rings at least `K` times is a fixed strictly positive number; then the leaf remains unpromoted and produces `K` leaves at the core. Iterating the conditional failure probability shows that infinitely many available trials without a single success have probability zero. Therefore every nonhub eventually has zero leaf neighbours (if core), so its degree is eventually constant: any future attachment at itself would create a leaf. A nonhub that never joined the core already has constant degree. Conversely a vertex of bounded degree is a nonhub. By the edge theorem all neighbours of a nonhub are hubs; once its degree has stabilized it has finitely many neighbours, each with leaf count tending to infinity by `L_w/n→Z_w` (Theorem 3.1). The simple upper classifiers therefore converge to 0 on nonhubs. On hubs the lower classifier converges to 1. Thus the terminal probability intervals have width tending to zero a.s., even around the zero atom.

**A general adjacent-label reduction.** Write `g(T,x)` for the expectation of a fixed bounded measurable test function of the limiting root degree fraction from root `x`. Conditional on `T_{n-1}`, expand `E[g(T_{n+1},n)-g(T_{n+1},n+1)|T_{n-1}]` by the successive attachment sites `x,y`. If `y=x`, the two births are exchangeable in the final tree. If `y` is an old vertex distinct from `x`, compare orders `x,y` and `y,x` on the same unlabeled extension with marks swapped: the surviving coefficient, multiplying the difference of the two leaf expectations in order `x,y`, is `r_x r'_y-r_y r'_x`, where `r_z=sum_{a~z}1/D_a` before the attachments and `r'_y` is the rate at `y` after the first attachment at `x` (likewise for `r'_x`). If `x` and `y` are not adjacent, `r'_y=r_y` and `r'_x=r_x` so these terms cancel exactly. For adjacent `x,y`, this coefficient is `r_y/[D_y(D_y+1)]-r_x/[D_x(D_x+1)]`. (CORRECTED by the lead, 24 Nov 12:10 UTC: as first written the two denominators were crossed, `r_y/[D_x(D_x+1)]-r_x/[D_y(D_y+1)]`; found by a child auditing Astra's late-label results and confirmed by the lead: `r'_y = r_y - 1/[D_x(D_x+1)]` after an attachment at the neighbour `x`, so `r_x r'_y - r_y r'_x = r_y/[D_y(D_y+1)] - r_x/[D_x(D_x+1)]`. Nothing in the certificate or in the bounds uses this formula.) The unmatched case is attachment at the first newcomer, with coefficient `r_x/(D_x+1)` and a comparison of the new degree-two vertex against its leaf. All terms have common denominator `n(n-1)`. This restricts a prospective proof to adjacent attachment sites and promotion of the first newcomer, but neither surviving term has a sign from symmetry alone.

## NUMERICAL — 62,125 standard-chain replicates, cutoff 1,000,000

`sim.cpp` samples a uniform vertex and then a uniform neighbour using a leaf-to-parent array and core adjacency lists; no extrapolation of the chains. The two endpoint bounds at the cutoff come from the classifiers above, using the simple upper bounds for speed and the truncated lower bound (large `l` or `c`: bound the omitted series geometrically). Average each endpoint to obtain a lower/upper sampling estimate for the true probability; the width estimates **conditional classification bias**, not sampling error. `stats.log` gives both endpoints and both Monte Carlo SEs for labels `1..60,100,300,1000`, and same-tree paired estimates for all adjacent differences. Seven pairs have negative point estimates even for the upper endpoint; the most negative is 59–60 (`[-.0034183,-.0025082]`, paired SE `.0016660`). No reversal exceeds two paired SEs on the upper endpoint and none is significant after adjustment for looking at many labels. The early order and overall decay are clear. Selected endpoints:

| Label | Average lower | Average upper | SE (approx.) |
|---:|---:|---:|---:|
|1|.698246|.698443|.001842|
|2|.697487|.697655|.001843|
|3|.424695|.424892|.001983|
|8|.247058|.247304|.001731|
|59|.096943|.097352|.001188|
|60|.099860|.100362|.001205|
|100|.077387|.077848|.001074|
|300|.050211|.050842|.000880|
|1000|.027770|.028614|.000667|

Under the usual independent sampling model, one-sided Hoeffding plus a union bound over the 126 endpoints gives a conservative simultaneous 95% half-width of `.00794`, in addition to the classification band. Positive-threshold paired comparisons use the conditional bounds above for 12 thresholds from `1e-5` to `.75`. `edf.tsv` reports each of the 708 comparisons, including the empirical tail difference of `D_v/m`, rigorous *conditional-on-terminal-tree* population-bound estimators, and paired SEs. No adjusted evidence against stochastic domination. At the least favourable tested points the samples are sparse (e.g. `47–48` at `.1`, upper endpoint `−.0000508±.0000646` SE). These are grid diagnostics, not a uniform-in-threshold proof.

## CONJECTURE / DEAD ENDS

* Global monotonicity of `p_u` and first-order stochastic domination remain open here. The simulation cannot distinguish a small genuine local reversal from sampling noise; the matched-pair standard error exceeds plausible successive differences for many later labels.
* Comparing birth and conversion rates of the two vertices gives no universal pointwise ordering of degrees or leaf counts. An old core can have no leaf and negligible rates from very large neighbours, while a newly promoted vertex has one leaf and birth rate at least one. Coupling on current degree alone loses the opposite effects of core adjacency and leaf reservoirs. In the fixed path cases I found the one-dimensional comparison and finite cancellation proofs above, not a two-vertex order-preserving coupling.
* Exact individual probability intervals are too wide for the middle adjacent labels; passing to the signed common rooted-state measure resolves this through 8. The same trick yields no automatic all-size sign: the coefficients depend on the label and some rooted types carry negative coefficients. First moment/second moment by themselves do not decide whether a tiny limiting mass is exactly zero.

## Reproduction and logs

All commands run in this directory; CPU pin and nice level were used for compilation and each run.

```
taskset -c 12 nice -n 5 python3 exact.py
taskset -c 12 nice -n 5 g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=undefined -fno-sanitize-recover=undefined exact.cpp -o exact
taskset -c 12 nice -n 5 bash -c 'for m in 3 4 5 9 10 11; do ./exact "$m" < "exact-$m.in"; done' > exact.log 2> orbit.log
taskset -c 12 nice -n 5 python3 tails.py > exact-summary.log
taskset -c 12 nice -n 5 g++ -std=c++17 -O3 -Wall -Wextra -pedantic sim.cpp -o sim
taskset -c 12 nice -n 5 timeout 560 ./sim 1000000 12000 81273 > sim-a.log
taskset -c 12 nice -n 5 timeout 560 ./sim 1000000 40000 2791840184 > sim-b.log
taskset -c 12 nice -n 5 timeout 540 ./sim 1000000 22000 61869110453 > sim-c.log
taskset -c 12 nice -n 5 round2/.venv/bin/python analyze.py sim-a.log sim-b.log sim-c.log > stats.log
```

Exact run exited successfully with mass/count assertions and no sanitizer diagnostics. Simulation runs A/B/C used 1m50s, 4m58s, and 3m19s user CPU, respectively (elapsed 2m11s, 9m20s—the timeout, and 3m28s); B contributed its 28,125 complete rows. The digital PRNG simulations and floating-point summaries are numerical. Source, rational certificate evidence, and numerical summaries are separate files named above.
