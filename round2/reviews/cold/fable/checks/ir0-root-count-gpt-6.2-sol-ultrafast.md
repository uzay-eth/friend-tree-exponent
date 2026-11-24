RESULT: CERTIFIED for Re s >= -0.551315

The unique zero there is the real mu_0 from the question (counted with multiplicity). In fact the computed left endpoint is the slightly more negative dyadic -606177253068/2^40. In the larger half-plane Re s >= -9/16 there are exactly three simple zeros: the real zero and one conjugate pair. The upper zero lies within a closed complex disk of radius 2^-17 about (-606178799255 + 551435146100 i)/2^40, with the lower enclosure obtained by conjugation. The method is an independent interval argument-principle checker, `check.cpp`, and a disk Krawczyk check for the upper root. It does not use the exploratory scan for the winding calculation, and it uses no external numerical library.

### Continuation and exclusion at infinity

Let G be the single-core generator, L=d-c, and b=L/d (zero at absorption). Define r=Gb. Direct evaluation gives

    r(d,c) = b(d,c) [c/(d+1) - 1/d].

The chain is nonexplosive, by comparison of its jump rates with a linear pure birth process. Kill when it is absorbed. For d >= D=512 use the weight W=(c+17)/(L-1/2) at live states, W=0 at absorption. We have |r|<=W and GW<=-(2/3)W. In fact the degree jump contributes -L/(L+1/2) to GW/W. At L=1 the birth jump kills and contributes -b. At L>=2 the birth jump contributes b times

    (1 + 1/(c+17)) 2/(2L-3) + 1/(c+17).

The checker uses b<=L/D for 2<=L<=15 and b<=1 otherwise, with the decreasing positive factor bounded at L=16 for the latter range. With beta=2/3, Dynkin's inequality (localize to finite sets, then use Fatou for the nonnegative supermartingale) gives E W(X_t)<=exp(-beta t) W(X_0). For d<D, the first-exit calculation extends integrability to every initial state: jump rates at live states are at least 1 and the strictly increasing d+c gives a finite triangular system.

For Re s>0 integration by parts yields R_s:=s H_s=b+integral_0^infty exp(-st) E r(X_t) dt. The right-hand side is holomorphic for Re s>-2/3, by locally uniform absolute convergence; below D this follows by the finite first-jump recursion. Thus F=(R_s-s)/s is meromorphic in this strip and G0(s):=R_s(2,1)-s is the holomorphic function used in the count. The checker proves that R_0(2,1)>0; hence F(s)=R_0(2,1)/s+O(1) and clearing this pole adds no zero. Equivalently the residue is the limit of E b(X_t): its derivative is E r(X_t), which is integrable. On surviving paths d increases without bound, since only finitely many transitions can occur in bounded degree and the nonabsorbing holding times are almost surely finite. Apply the weight estimate from first entrance to d=D and dominated convergence (the entrance has only finitely many possible states). This gives E[(1-b(X_t)) 1_live] -> 0 because c/d <= W after entrance, and the probability of being live before entrance tends to zero. Thus the residue is the survival probability.

For a=9/16, the checker bounds the absolute discounted r integral by W/(beta-a) at d=D, by zero on the diagonal, and below D by first-jump recursion with positive denominator d+1-ad/L and source |c/(d+1)-1/d|. At (2,1) it obtains the upper bound U=32172867510368/2^40. Therefore |H_s(2,1)| <= (1/2+U)/|s| < 1 for Re s>=-9/16 and |Im s|>=30. For Re s>=1, the original integral has modulus strictly less than 1 since b<=1 and b(2,1)=1/2 with positive holding probability for each finite time. This leaves rectangles with right side 1 and imaginary endpoints ±30.

### Explicit boundary error

For R_s use the approximate boundary Q_s(d,c)=sum_{i=0}^6 q_i(c)d^(-i) with q_0=1. Expand (s-G)Q_s-sb in powers of x=1/d using the degree substitution x -> x/(1+x), and the birth substitution c -> c+1. The coefficient of x^n is (s+n-Delta)q_n + B_n, where Delta p(c)=p(c+1)-p(c); B_n is the coefficient computed from previous q_i, including sc at n=1. The checker constructs q_n as the negative of the finite polynomial inverse sum sum_{j=0}^n Delta^j B_n/(s+n)^(j+1). Thus the residual vanishes through degree six algebraically, including on complex intervals. Multiplying the residual by (1+x)^6 gives a polynomial N in x of degrees 7 through 13, with coefficients polynomial in c of degree at most 7 (the multiplier clears the x/(1+x) substitution). The checker obtains its terms by binomial convolution and uses

    |(s-G)Q_s-sb| <= M ((c+20)/d)^7,
    M = sum_{j>=7,k<=7} |N[j,k]| / (D^(j-7) 20^(7-k)),

uniformly on each checked contour segment. It bounds complex coefficient modulus by the sum of coordinate moduli.

Stop in d>=D at the first exit from c/d<=2/3. For w_{p,A0}=((c+A0)/d)^p with (p,A0)=(7,20) or (12,40), we have Gw<=-(3/4)w before exit. Indeed the degree term divided by w is -L[1-(1+1/d)^(-p)], at most -p/3+p(p+1)/(2D) by Taylor's inequality and d/3<=L<=d. The birth term is at most ((A0+1)/A0)^p-1. Both bounds are checked by exact intervals. For either w, optional stopping yields E[exp(a tau) w(X_tau);tau<infinity]<=w(X_0), and the time integral of the discounted pre-exit expectation of w is bounded by w(X_0)/(3/4-a).

On exit, 2/3<c/d<=2/3+1/D (a birth jump). The checker bounds |R_s-Q_s| there by 1+W/(beta-a) plus the interval polynomial absolute majorant for Q_s with this ratio. Call this uniform bound B. Before exit, |R_s-1| <= c/d + W/(beta-a) and |Q_s-1| is O((c+20)/d) by its degree bounds, uniformly on each compact contour segment.

At degree cutoff states 1<=c<C=36 the resulting tail bound is

    |R_s-Q_s| <= [M/(3/4-a)] w_{7,20} + [B/(2/3)^12] w_{12,40}.

For clarity, the vanishing time-terminal bound used here is obtained using the baseline absorbing weight W, not the first-power drift pre-exit. Both |R_s-1| and |Q_s-1| are bounded by a fixed multiple of W before exit (because c/d<=2/3 and the polynomial degrees are <=i). Therefore the exponentially discounted expected time-terminal difference on no exit is O(exp(-(beta-a)t)) by the global killed W bound. This proves the displayed cutoff error by applying Dynkin to R_s-Q_s; its generator equation holds by analytic continuation of the original finite-step identity, so localization is legitimate with bounded functions and residual on the stopped region. At the c cutoff below D, the checker instead sets its approximation to zero and bounds |R_s|<=1+U(d,C) from the same positive majorant recurrence. The diagonal is exact zero.

### Segment enclosure and winding

The finite triangular approximation has source s and denominator Z=d+1+sd/L, obtained directly from (s-G)R_s=sb. Propagate boundary error with the exact recurrence in absolute value, using the rational lower bound max(Re Z, |Im Z|, 7/10 (Re Z+|Im Z|)) for |Z|, uniform over each segment. For midpoint rounding, every finite state is rounded to a dyadic complex representative; the difference from its local interval calculation plus propagated input errors gives an explicit numerical-error radius for the exact finite approximation. For variation within a segment, subtract the finite recurrence at its midpoint from the recurrence at s, divide by s-midpoint, and bound the result by the same step factors with additional forcing 1+(d/L)|midpoint finite value|. Cauchy's estimate on radius 1/8 bounds the derivative of the degree-boundary polynomial, evaluated on the enlarged segment rectangle (all possible polynomial coefficient poles are at the negative integers). Add the unit slope for subtracting s in G0. These bounds together give a coordinate rectangle containing the full image of each segment. The checker accepts a segment only if the hull of this rectangle and both dyadic polygon endpoint representatives excludes zero. Thus the polygon and true boundary have an explicit straight-line homotopy avoiding zero. The polygon winding is computed by signed positive-ray crossings with exact 128-bit integer cross products.

The planted test replaces G0 by (s-(1/4+3i/8))G0 and runs the same outer rectangle with no winding shortcuts. The only new zero is explicitly inside. Overflow causes subdivision or failure; it cannot accept a box. Fixed-point endpoints are in units of 2^-40 and all elementary operations round outward by integer division.

### Krawczyk disk

At the upper disk center z, take the complex preconditioner A to be the reciprocal of the centered difference of the two *dyadic* midpoint representatives at z±h, h=2^-11. Their interval errors divided by 2h bound the error due to using representatives. On the circle |u|=R=2^-7, the checker bounds |G0(z+u)-v| by M=err+2R*slope, with v the center representative and with residual and variation bounds taken on the containing square. Cauchy's coefficient bound M/R^n gives the centered-difference Taylor error (M/R)(h/R)^2/(1-(h/R)^2). It also gives the variation of G0' on |u|<=r=2^-17 as at most (M/R)(r/R)(2-r/R)/(1-r/R)^2. Thus the checker bounds |1-A G0'| by a strict contraction constant k and checks |A G0(z)|+rk<r. The map s -> s-A G0(s) sends the disk strictly into itself and is a contraction there, giving a unique zero. Its containing square was also counted independently by winding; the exact upper square endpoints are real [-606187187863,-606170410647]/2^40 and imaginary [551426757492,551443534708]/2^40. The Krawczyk image bound alone tightens both coordinates to within 4904260/2^40 of the center. Neither calculation asserts additional digits beyond its enclosure.

### Checker command and results

Run from `round2/children/ir0-roots/`:

```sh
g++ -O3 -std=c++17 -Wall -Wextra -Wno-missing-field-initializers -o check check.cpp && taskset -c 14 nice -n 5 timeout 1150 ./check
```

Exit status 0; actual stdout:

```text
global majorant upper endpoint (units 2^-40): 32172867510368
high frequency cutoff 30
right half-plane Re >= -1/2: winding=1, segments=508
right half-plane Re >= -9/16: winding=3, segments=754
scan point, midpoint residual and error units 2^-40: -1121725 -1354860 5916342
upper root: winding=1, segments=171
Krawczyk disk contraction, image radius upper (units 2^-40): 455040577127 4904260
planted extra zero: winding=2, segments=586
pole residue real endpoints, units 2^-40: 674214593914 674214679584
right half-plane Re >= -551315/1000000 (left endpoint -606177253068 * 2^-40): winding=1, segments=1580
```

The known positive real root identifies the counted root without needing to reproduce its 60-digit bracket. Alternatively existence follows because H is positive and strictly decreasing on the positive real axis, has the proved positive pole at zero, and is strictly below one at s=1. The winding counts exclude boundary roots. Together with conjugation and the half-plane count of three, they ensure there are no unaccounted roots between the two nonreal isolations and the real root anywhere in Re s>=-9/16. No zero count farther left than -9/16 was attempted; nor was the optimal real cutoff determined to the last digit. The smallest executed, successful uniqueness cutoff was -0.551315.
