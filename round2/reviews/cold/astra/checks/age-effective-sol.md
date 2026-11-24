# Cold audit: a computable cutoff for adjacent hub ordering

**Verdict: the terminating search is justified conditional on the repository's physical renewal and local-limit inputs.** The author's proposal (`open/age-order/effective/README.md`) does not evaluate a physical starting label. The arithmetic script alone cannot certify one. Below I give an explicit analytic spectral decay rate `eta=2^-348` and an explicit positive bound on `h(T3)`, without a winding search. This still leaves the finite-cylinder search unevaluated; the spectral rate is not a numerical label cutoff.

## What the algorithm takes as input

The inputs are the nineteen-label residual generator and its weighted drift inequalities in finite-kernel (`fringe/finite-kernel.md`), the exact real Perron-root characterization and bracket, the physical immigration identity and depth reconstruction in spatial-limit (`fringe/spatial-limit.md`), and the two-birth signed score identity, local hub dichotomy and positive stationary gap from the adjacent (`open/age-order/adjacent/README.md`) and late (`open/age-order/late/README.md`) results. In particular `pi G=(1−mu)pi F≥beta`, since the discounted promotion flux implies the undiscounted hub score is at least `mu/2`. No effective almost-sure empirical convergence, convergence rate for a martingale, or lower bound on a time-to-hub uniform over environments is an input. Nor does the proof require evaluating stationary probabilities of infinite cylinders by an oracle.

## Spectral data: search and an explicit alternative

A complex point evaluation is a finite Gaussian-rational killed resolvent; it is invertible for `Re z>rho` by domination in `W`. Decomposition at first exit bounds the modulus of the error by the real exit account, uniformly on a half-plane as in finite-kernel (5)–(6). In particular for `s≥l` one can use `delta=1/24`; below `l` take any positive rational `delta≤1` with `rho+delta(1+q/2)≤s`. Rational power enclosures suffice. The strict real-Perron comparisons used to bracket `mu` are semidecidable: for `r>1` search for a positive rational Collatz witness on the lower matrix, and for `r<1` on the upper matrix. Dovetail rational brackets of shrinking width; never require a comparison at `mu` itself. The strict margins exist by positivity and the Perron characterization. Derivative bounds and elementary determinant product bounds give rational disks covering each short contour segment; once these avoid zero, compare with a rational polygon in the disks to obtain the winding by signed crossings of a ray chosen off the vertices. A sufficiently fine subdivision exists at any admissible rational `sigma` because the determinant has a positive minimum modulus there. One must dovetail over `sigma`, subdivisions and accuracy; waiting for any fixed wrong contour can hang. The author’s constants for inverse norm, residue and convolution are upper bounds: for the remainder of the renewal density invert `(I−K)^−1 K²`, split at height `H`, use the `2B1²/y²` high-frequency bound and integrate, then add the `k(t)` term. The `4H` residue multiplier is more than contour length divided by `2π`. The weighted planted error propagates through `A+2e^(qt)`; the shallow-layer sum is bounded by `2^k e^(t/2)` and is absorbed by `4^k e^(sigma t)`. These steps use the finite-label nonlattice/simple-pole proof for existence of the contour; plain convergence of the renewal density would not certify this search.

Here is an independent explicit margin, using only the same nineteen-label kernel inputs. Use the notation `l,u,rho,A` of the proposal and `b=(l−rho)/2`. Set `e` to the rational minimum in finite-kernel (81), `g_*=e(l−rho)/(19A)`, and choose Perron right vector `g` at `mu` with `max g=1`, so `min g≥g_*`. Conjugate by `diag(g)` and call the resulting kernel transform `P(z)`. The matrix `P(mu)` is row stochastic and every row contains at least `e g_*` in column 2. More strongly, the label-2 contribution from the original unchanged newborn gives a common measure in that column at `mu+iy` of density `g_* exp[−(mu+2)t]/2`. The original newborn's internal rate is `3/2` and its death is at most `1/2`; this is a submeasure of the time kernel, not an assumption of an independent time law. Its mass is at least `a=g_*/64`. Subtracting the submeasure leaves a nonnegative time kernel. Thus for `|y|≥r/2`, `0<r≤1/4`,

```
||P(mu+iy)|| <= 1 - a*r²/65536.
```

Indeed the removed exponential mass loses fraction `1−c/sqrt(c²+y²)`, `c=mu+2<3`. At `|y|≥r/2` this is at least `r²/144` (use `1−1/sqrt(1+x)≥x/4` at `x=(r/2)²/c²≤1/9`). This is stronger than the displayed bound. With `L=A/(g_* b²)`, `L2=A/(g_* b³)`, for `|z−mu|≤b` there are bounds `||P(z)−P(mu)||≤L|z−mu|` and `||P(z)−P(mu)−P'(mu)(z−mu)||≤L2|z−mu|²`. The first also bounds variation along the horizontal shift from `mu+iy` as long as the real part is at least `rho+b`.

For small `w=z−mu`, let `lambda` be the stationary row vector of `P(mu)`, project to constants with `1 lambda`, and let `E=I−1 lambda`, with norm at most two. The contraction of row oscillation from the common submeasure gives `||P(mu)^n−1 lambda||≤2(1−a)^n`: subtract the shared constant-row transition and iterate the resulting contraction in total variation on probability rows. On the subspace `lambda x=0`, `(I−P(mu))^−1` thus has norm at most `2/a`. For `8L|w|≤a` the lower-right inverse of `I−P(z)` in this splitting has norm at most `4/a` by Neumann. Its scalar Schur complement has linear term `−lambda P'(mu)1 · w`, with coefficient at least `a` by the same exponential submeasure's first moment (`g_*/[2(u+2)²]`). The first-order Taylor error plus the off-diagonal Schur correction is bounded by `(L2+8L²/a)|w|²`. Therefore take

```
r   = min(b,1/4,a/(8L),a/(2(L2+8L²/a)))
eta = min(r/2, a r²/(131072 L)),    sigma = mu−eta.
```

In the disk the Schur complement is nonzero except for its simple zero at zero, and has modulus at least `a|w|/2`. On the vertical line through `sigma`, for `|y|≤r/2` the disk estimate applies and `|w|≥eta`; on the rest the phase deficit minus the horizontal perturbation leaves at least `a r²/131072`. The same arguments exclude any other pole in `sigma≤Re z≤mu`, and the real comparison excludes poles to the right. In original coordinates, on the new vertical line use

```
M = max(2,(8/a+16/(a eta))/g_*, 131072/(g_* a r²)).
```

For the small disk bound, invert the lower-right block, then correct by the scalar Schur inverse; both off-diagonal factors have norm at most two when `8Lr≤a`, and the constant-free projection costs at most two. To bound the residue, take the `w→0` limit of `w(I−P(z))^−1`: it is `1 lambda/(-lambda P'(mu)1)`, so `D_* = 1/(a g_*)` suffices in original coordinates. This avoids needing an explicit pole location for the analytic bound. The author’s `C_v,C_J,C_P,C` now work with `M` on the vertical line, and with each `sigma` in a denominator replaced by the rational lower bound `rho+b`; choose `H>2B1+1` with that same lower bound. The fixture gives `g_*≥2^-23`, `r≥2^-134`, `eta≥2^-348`, `M≤2^404`. In particular one can use the rational rate `2^-348` in the final envelopes even though `sigma` in the inversion argument was written using `mu`.

One can also remove the search for a positive lower bound on `h(T3)`. With `max g=1`, normalize `ell g=1`; then `1≤ell 1≤1/g_*`. For the leading term write `d_2=ell∫e^(−mu t)JR_t1 dt / (−ell K'(mu)g)`. Unchanged newborns give a numerator at least `1/(u+2)` and the denominator is at most `A/[g_*(l−rho)²]`. The newborn value `H(J_i)=(K(mu)g)_i=g_i≥g_*`. In the physical immigration integral from `T3`, keep only the first direct birth term while the full tree has not changed. The probability of no jump to time `t` is `e^(−3t)` and the direct promotion rate then is one. Consequently

```
h(T3) >= g_*²(l−rho)² / [A(u+2)(u+3)] >= 2^-60.
```

An upper bound follows from the weighted leading operator and the immigration integral with `A+2e^(qt)`; no eigenvector approximation is required for these two bounds. Substituting the exact rational margins into the convolution constants gives `C≤2^453` and `h(T3)≤2^68`. Thus the explicit normalized spatial rate is `2^453 · 4^k exp(−2^-348 t)` times the test norm; these constants are deliberately conservative.

## Effective local convergence versus a limit without a modulus

For a placed tip or middle, use the lower hit test with cutoff `d≥3` so a hit always includes a fresh leaf; use the upper test with hit flag plus the capped first-send charges of neighbours when no hit has occurred. At fixed `d` one obtains limiting endpoints `P(hit d ever)` and `(1−1/d)P(hit d ever)`, because any hub hits and, on a nonhub path with no hit, the fixed eventual neighbours have diverging leaf stocks. Thus the width tends to at most `1/d`. The per-core uniform coefficient bound, applied after summing both endpoints per contrast, permits dominated convergence under `pi` despite arbitrary local degree. It establishes existence of a finite choice of `d,T` for a desired *stationary* tolerance; on its own it does not give a computable `T` from a pointwise rate.

For each candidate take rational `T` and use the old-cell coupling after insertion. The live queried cells are initially at distance at most five from the score root: the old source is within one ordinary edge of its charged core root, its old neighbour within two, and a new tagged vertex or its pre-existing anchor is within a further bounded number of edges. A neighbour newly created after insertion belongs to the relevant old cell. Promoting an old leaf in an insertion can create a new core but cannot enlarge old core distances except by adding its own edge. Keeping all cut slots inert preserves denominators. A path from an omitted cell has an ordered suffix of length at least `R−5`. Reversing it against random-walk probabilities bounds the total of its proposal weights by the target-cell initial core degree. The proposed sum `V` is therefore sufficient, including the length-zero/singleton issue by taking `R≥6`. The error is added unconditionally. Bound the size-cut probability by `Me^T/B`; the finite marked state space includes both the tested tags and the hit flag. Uniformization has nonnegative Poisson weights and a computable tail (for example the geometric ratio bound past an integer greater than twice its parameter). Use exact rational rates, bounded payoffs and outward exponential enclosures. Boundary total degrees may require one extra input level; the radius used in all subsequent estimates must include it.

For computability it is enough that for fixed `d,T` the coupling error goes to zero under `pi` as `R` grows; for chosen `R` only then raise the initial degree cap. The excluded part tends to zero by local finiteness. The optional moment bound also checks out: the old-root term vanishes after discounting a bounded truncation, promotion from the cell has intensity at most `e^(mu s)`, and the two age moment exponents are at most `1/4`. `c²≤40(5/4)^c` follows from the zeroth, first and second binomial terms; Cauchy–Schwarz gives the stated `J^(-1/8)` bound, and mass transport reverses each bad path. The moment bound by itself does not give an explicit large-`T` rate for the no-hit charge.

For any fixed rational cylinder, approximate `Q_t f(T3)` and `Q_t1(T3)` from a finite Yule-size kill. On cutoff the final size is at least `B` and `E N_t²≤12e^(2t)`, giving the stated error times the sup norm. The certified spatial error bounds the two leading coefficients before division; the denominator `h(T3)` has the positive lower bound above. As `t` is increased, the finite cutoff and uniformization precision must also be increased to approximate each leading coefficient to a requested error. These form finite evaluations, with a known stopping error, of `pi f`; there is no test of an equality of computable reals. Dovetail cylinder candidates and tests of the *strict* inequality `pi L>beta/2`. Some candidate satisfies it since `L≤G≤U` and widths can be arbitrarily small. That is a terminating search with a verifiable certificate, not just a convergence claim. The same reasoning semidecides a strict finite-label comparison of two true hub probabilities; shrinking enclosures is not a decision method at equality.

## Deterministic sizes and scope of the cutoff

For a supported cylinder, a single birth can only change its contribution at roots whose before or after ball contains a changed core mark or the promoted core. In the latter case, use the old parent on the side where the new core did not exist if a before-ball contribution is being lost. Counting backwards along paths in a *capped* contributing ball bounds affected supported roots by the stated `16(k+1)(D+1)^(k+1)` times the sup norm, including crossings of the cap. For the constant offset only promotions change the sum. Include boundary degree reads in `k`; if a library represents the cylinder using an extra uncapped level, use the jump count with the original capped degree-read radius and the semigroup bound with the enlarged cylinder radius. Taking the larger radius in both constants is safe.

From size `m=floor(n^(19/20))`, the full population at clock time `log(n/m)` is Yule and independent of the **discrete** jump chain because holding rates are exactly the current size. The signed score itself need not be monotone. The jump bound times `E|N−n|≤n/sqrt(m)` compares it on the same discrete history; the cell expectation is the full expected sum because the cells partition all cores. Centering at `pi f` cancels every initial `h` coefficient. For the core-count bounds use only monotonicity and `a_{j+1}≤(1+1/j)a_j`. From `T3` at time `log n`, `E a_N≤4 a_n` by `E N=3n`, and `E a_N≥a_n/4` by the geometric Yule-family tail. This gives both ratios in the proposal; in particular the floor loss is bounded by two and `l−21/40=1/25`. The two error envelopes are decreasing for all large `n`, so testing at one cutoff is sound for the entire subsequent tail. For `pi L>beta/2` and total error less than `beta/4`, one obtains the stated gap lower bound with factor `beta h_*/32`. Neither this result nor the early prefix decides the intervening finite set without further computations.

**Numerical status.** The analytic rate and `h_*` above are evaluated; no actual accepted cylinder (`k,D,L`) has been evaluated, hence no numerical physical `N0` is given. Once a cylinder passes, one explicit substitution (no huge integer expansion needed) is

```
J = 96(k+1)(D+1)^(k+1)+6
N0 = 2^s,  s = max(515·2^348, 20·2^348(982+2k), 25(71+ceil(log2 J))).
```

The three envelopes use coefficients at most `2^514`, `2^(975+2k)`, and `2^64 J`, respectively, and `beta/8>2^-7`. This fixes the spectral inputs in the label formula but does not fix the cylinder search result. The printed `2^525` in the smoke test is for illustrative inputs only, not the friend-tree cutoff. The `delta` from finite-kernel is an exit modulus for evaluating the kernel; it is not itself an exponent for spatial or deterministic-size convergence.

## Replay

`python3 -B round2/checks/age-effective-sol.py` returned:

```text
rational arithmetic checks passed; fixture N0 = 2^525
g_* >= 2^-23
r >= 2^-134
eta >= 2^-348
h_* >= 2^-60
inverse bound M <= 2^404; 505 clock floor controls
C <= 2^453; h(T3) <= 2^68
symbolic cutoff formula: 3 radius/cap controls
```

The independent fixture evaluates the rational margins and checks their inequalities; the imported author self-check is solely for the final cutoff arithmetic. It does not run a finite kernel matrix, a cylinder search, or an empirical convergence test. Author files were not edited.
