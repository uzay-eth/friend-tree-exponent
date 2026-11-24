RESULT: NO CANDIDATE ABOVE CHANCE

## Friend tree: preliminary input

`VALUE=0.56572366`, `SIGMA=0.000000007`, closed window `VALUE ± 2.1e-8`. **Zero hits in all 16 families** among 125,782,496 attempted descriptions (including overlaps). Counts in the target and controls are distinct numerical roots *within a family*. The three nearest values and formulas in each family are in the complete stdout (`run-friend.txt`); the nearest distances (signed) are also below. No unresolved root-isolation subintervals in this run.

| Family | Descriptions tried | Target | Control mean | Control max | Three nearest, in sigma |
|---|---:|---:|---:|---:|---|
| Rational | 30,397,485 | 0 | 1.316406 | 4 | −5.18, −6.77, −8.58 |
| Expressions in constants | 6,193,330 | 0 | 0.062500 | 1 | −102.08, −118.66, −120.80 |
| Algebraic, d=2 | 1,426,433 | 0 | 0.017578 | 1 | −976.69, −1185.58, −1420.83 |
| Algebraic, d=3 | 2,595,743 | 0 | 0.066406 | 1 | −26.05, −135.56, +147.08 |
| Algebraic, d=4 | 2,996,593 | 0 | 0.054688 | 1 | −94.66, +100.97, +252.05 |
| Algebraic, d=5 | 5,258,711 | 0 | 0.076172 | 2 | +7.52, +78.75, −96.94 |
| Algebraic, d=6 | 8,825,097 | 0 | 0.132812 | 2 | −43.65, −51.35, +61.91 |
| Multiquadratic | 12,122,400 | 0 | 0.072266 | 2 | −71.63, −94.54, +216.72 |
| HKRZ F[1] | 5,626,871 | 0 | 0.099609 | 2 | +27.14, +45.46, −114.41 |
| HKRZ other z | 3,506,652 | 0 | 0.041016 | 1 | −124.95, +166.34, +212.82 |
| Gamma kernels | 6,136,641 | 0 | 0.076172 | 2 | +37.49, −80.87, +92.00 |
| Kummer U / hypergeom. / Mellin | 3,344,315 | 0 | 0.058594 | 1 | +76.30, +97.23, +110.06 |
| Sparse elementary / special | 7,622,607 | 0 | 0.142578 | 3 | +46.28, −56.06, +59.58 |
| Affine two constants | 25,119,864 | 0 | 0.115234 | 2 | +9.07, +17.66, +22.36 |
| Self powers and Lambert | 2,006,244 | 0 | 0.019531 | 1 | +18.69, −304.83, +504.41 |
| Affine fixed points | 2,603,510 | 0 | 0.023438 | 1 | −76.16, +415.39, −424.89 |

Closest rational: `5539/9791 = 0.56572362373608415892`, at −5.18σ. Closest nonrational: `0.56572371265405470019`, root of `3 x^5 − 5 x^4 − 4 x^3 − 7 x^2 − 3 x + 5 = 0`, at +7.52σ. Closest HKRZ F[1]: `0.56572385000898219546`, root of `−7 + 7x + (4+x+7x²)F(x) = 0`, at +27.14σ.

The 512 control centers are independent uniform draws on `[VALUE−0.045, VALUE+0.045]` after removing centers within `6 SIGMA` of `VALUE` and keeping the windows inside `(0,1)` (seed `68133021`, reused across families). Control windows can intersect one another. The reported upper-tail rank `(1 + #controls with count ≥ target count)/513` equals 1 in every row. This is a local-density baseline, not a formal p-value under a proved probability distribution for the exponent. The nominal Bonferroni cutoff for 16 families is 0.003125; none approaches it. Heights have independently reasonable empirical densities: the control means are all between 0.01 and 10 at the supplied preliminary sigma. At other sigmas the heights remain fixed so that only centers/window widths change. The catalogue was not chosen to make a target hit stand out; the initial constant-expression height of 12 had zero control hits and was raised to 84 to bring it into the requested density range.

## IR-0

The same exhaustive *bounded catalogues* give **zero roots in all 16 families**, within `±1e-50` around `0.4822669183233687182688976229135491851852409825893293466110915`, the midpoint of the certified bracket in `../../checks/cold-estimator/REPORT2.md`. No unresolved subintervals; candidate confirmation used 95 decimal digits.

Separately, PSLQ (95 decimal working digits, `maxsteps=20000`, bounds strict) was run at the bracket midpoint. These are relation-search settings, not exhaustive nonexistence certificates:

| Basis | Coefficient bound | Internal residual tolerance | Outcome |
|---|---:|---:|---|
| powers 0–2 | 10¹⁸ | 10⁻⁴⁰ | none returned |
| powers 0–4 | 10¹¹ | 10⁻⁴⁷ | none returned |
| powers 0–6 | 10⁸ | 10⁻⁵⁰ | approximate relation, rejected |
| powers 0–8 | 10⁶ | 10⁻⁵² | none returned |
| 1, x, x², F, xF | 10¹¹ | 10⁻⁴⁷ | none returned |
| 1, x, x², F, xF, x²F | 10⁹ | 10⁻⁴⁹ | none returned |
| 1, x, e^x, xe^x, e^−x | 10¹¹ | 10⁻⁴⁷ | none returned |
| 1, x, log x, x log x | 10¹⁴ | 10⁻⁴⁴ | none returned |
| 1, x, Gamma(x), Gamma′(x) | 10¹⁴ | 10⁻⁴⁴ | none returned |

**Important distinction:** the loose, high-height degree-six PSLQ search returned coefficients `[34600289, −95405394, 60930212, −26222804, −36140252, 52814813, 60233216]` in ascending order. Its root **does** lie in the coarse `±1e-50` acceptance window, but is `8.4657035250946378535e-60` below the midpoint, outside the certified bracket of radius `5e-61`. At the midpoint the polynomial is `−4.0156951199114e-52`; it can change by at most `2.37174329812507e-53` across one bracket radius. `verify.py` rejects the polynomial throughout the bracket with exact rational arithmetic. Thus the coarse window does not imply chance matches are impossible for a search at height `10^8`. For comparison the elementary bounded degree-six catalogue uses height 5; for height `H=10^8` there are `H(2H+1)^6 ≈ 6.4×10^57` integer degree-six descriptions with positive leading coefficient before primitivity. PSLQ did not enumerate these descriptions or give their number actually considered internally.

## Catalogue and bounds

`search.py` is independent; the other ISC implementation was never opened. For equations, the catalogue count is exact **before testing where roots lie** and accounts for primitive integer coefficients and one normalized sign; Möbius inversion counts constant terms screened implicitly. In algebraic rows this denominator counts *all primitive polynomials tried*, including reducible ones screened out from the roots by exact irreducibility testing. The expression count includes formulas outside `(0,1)` and different descriptions of the same value. The tables therefore count candidate descriptions tried, not an estimate of the global number of distinct eligible real roots. Cross-family and algebraically equivalent formulas overlap.

| Catalogue | Heights and scope |
|---|---|
| Rational | All reduced `p/q` in `(0,1)`, `q≤10000` |
| Algebraic | Exact degrees 2,3,4,5,6 with coefficient sup-norm bounds 75,24,11,7,5 respectively |
| Multiquadratic | `d x + a+b√r+c√s+e√(rs)=0`, `1≤d≤16`, other absolute coefficients ≤9; `(r,s)=(2,13),(2,3),(2,5),(3,5),(3,7),(5,13)`; 2,020,400 descriptions per field |
| HKRZ F[1] | `P+Q F[1]`, degree ≤2 each, `Q≠0`, height 7 |
| HKRZ other z | Same with `F[z](x)=∫₀¹ t^x exp(z(t−1)) dt`, `z=−1,1/2,2,3`; height 5, 876,663 per kernel |
| Gamma kernels | Same with `Gamma(1+x)` and lower/upper incomplete gamma of shape `1+x` at cutoff `1/2,1,2`; height 5, 876,663 per kernel |
| Kummer / U / hypergeom. / Mellin | Same at height 4, 257,255 per kernel: `1F1(1+x;b;z)` at `(b,z)=(2,−1),(2,1),(3,2)`; `U(1,b+x,z)` at `(1,1),(2,1/2),(2,2)`; `2F1(1,1+x;2;z)` at `−1/2,1/2,3/4`; Mellin kernels `(1−zt)^−1` at `z=±1/2` and `1_{t>z}` at `z=1/2,3/4` |
| Elementary / special | All supports from `x, e^x, e^−x, log x, 2^x,3^x,π^x,Gamma(x),psi(x),cos(πx),sin(πx)` with optional constant; ≤4 nonconstant terms at height 4, 5 at height 2, 6–11 at height 1 |
| Affine two constants | `d x+a+b A+c B=0` over unordered pairs of the 22 constants below, `1≤d≤14`, `|a|,|b|,|c|≤10`, `bc≠0` |
| Expressions in constants | For ordered distinct `A,B`, `A^B,A/B,log_A B` multiplied by `p/q` in lowest terms; also `A^(±p/q)`; `1≤p,q≤84` |
| Self powers / Lambert | `a+d x+b f(x)=0`, `db≠0`, height 32; `f=x log x`, `x^(k x)`, `x exp(k x)`, `k=±1/2,±1,±2`; also `W_0(c x)` for `c=1/2,1,2,−1/4` and `W_−1(−x/4)` |
| Affine fixed points | Same affine form at height 28 over 35 kernels: all 25 HKRZ/Gamma/Kummer/hypergeom./Mellin kernels and 10 elementary/special functions other than `x` |

Constants: π, e, log 2, log 3, Euler's gamma, ζ(3), Catalan, 1/π, W(1), golden ratio, and square roots of `2,3,5,6,7,10,11,13,14,15,17,19`.

The integer-window sieve evaluates the nonconstant part at each window center and uses the sum of coefficient-weighted derivative bounds times the radius, with a floating-point pad, to test feasible constant integers. For survivors it reevaluates in 95-digit arithmetic. An endpoint sign change with monotonicity permits a unique root; a midpoint Taylor bound or, failing that, bisection of the interval checks for hidden root pairs. The derivative bounds come from differentiated integrals for the kernels (for example, multiplying a Mellin integrand by powers of `log t`), standard integral or series formulas for Gamma/psi, the Lambert W derivative formula, and direct differentiation elsewhere. Root isolation here is numerical rather than interval-certified; the IR-0 PSLQ polynomial rejection is exact. Roots that differ by at most `1e-54` within a window are treated as identical. For nearest candidates the code widens a separate window until it has three roots, without changing the target/control window.

## Reproduction and actual output

From `round2/children/isc-indep`:

```sh
taskset -c 14 nice -n 5 /usr/bin/time -p round2/.venv/bin/python -B search.py 0.56572366 0.000000007 > run-friend.txt
```

```text
real 63.32
user 58.54
sys 4.76
```

The full stdout, including exact formulas and distances for the three nearest in every family, 512-window histograms, IR results and timings, is in run-friend.txt. It ends `Tests=16; Bonferroni alpha per row=0.00312500; wall=63.091s`.

```sh
taskset -c 14 nice -n 5 round2/.venv/bin/python -B verify.py
```

```text
poly: exhaustive height 2 primitive count 223: OK
field: exhaustive height 2 primitive count 1169: OK
kernel: exhaustive height 2 primitive count 7399: OK
sparse: exhaustive height 2 primitive count 148: OK
F[z] versus independent quadrature at all 5 z: OK
Independent height-1 HKRZ enumeration; PAN, DAN and friend centers: [3, 1, 0] : OK
Independent rational sieve: 3043 reduced rational candidates, 5 hits: OK
Exact rational Taylor bound rules out the PSLQ degree-6 near relation throughout the certified IR-0 bracket: OK
```

**Not done:** Final friend-tree value was not yet provided. No general formula grammar, unconstrained special parameters/higher degrees/heights, model-derived spectral identity, or proof of transcendence or nonexistence of a closed form was attempted. The program accepts a replacement `VALUE SIGMA` with the target interval strictly in `(0,1)`; the runtime and calibration shown here pertain to the supplied small-uncertainty estimate.
