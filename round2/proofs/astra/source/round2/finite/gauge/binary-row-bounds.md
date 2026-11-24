# The two true newborn rows and their excess

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## 1. Controlled resolvents on root degree and core degree

Use the certified four-bin weight, $`q=2/5`$, $`d=7/50`$, $`\rho=0.54`$, from [the binary construction](binary-and-return.md). Let $`g\ge0`$ be a test on the two physical labels. Denote $`g(a(D))`$ by $`g_D`$, the newborn weight credit by $`\chi(D)`$, and $`\max(1,\sup W(Y^+)/W(Y))`$ by $`\kappa`$. At a root with permanent bin $`i`$ let $`h=h_i(D,c)`$ and $`b=1-c/D`$. Its silent intrinsic weighted drift from its own local row is at most
```math
T=(D-c)\Delta_Dh+[(c-1)/2](\Delta_Dh)_++b[\Delta_ch+\chi(D)].
```
The other weighted vertices have nonnegative slack for $`s\ge\rho`$. Therefore
```math
 (s-R)W\ge h+\frac{(s-d)h-T}{a_i+d},\quad
 U_sg\le\alpha W,\quad
 \alpha=\max_{i,D,c}\frac{b g_D}{h+[(s-d)h-T]/(a_i+d)} .
```
An upper bound on the maximum suffices. The degree tail has constant $`h,\chi,g_D`$, and its maximum is bounded by setting $`b=1`$. The remaining maximum is finite rational arithmetic.

Write $`\mathcal H_{\lambda,e}^\pm r`$ for the lower and upper bounds on the discounted integral of a running reward $`r(D,c)`$. The controlled root has gain rate $`D-c+\xi`$ with $`0\le\xi\le(c-1)/2+e`$ and promotion rate $`b`$. For $`D<N`$ the two bounds take the minimum or maximum over the two endpoints of $`\xi`$ of
```math
 \frac{r+(D-c+\xi)v(D+1,c)+b\,v(D,c+1)}{\lambda+D-c+\xi+b}.
 \qquad\text{(1)}
```
The quotient is monotone or constant as a function of $`\xi`$. Stopped jump comparison gives the inequalities also for a predictable time-dependent input. Direct calculation yields $`-\eta b\le G b\le1-b`$ for such a root when $`e\le1/2`$, with $`\eta=\mathbf1_{c=1}/[D(D+1)]`$ nonincreasing along paths. On the degree tail $`D\ge N=32`$, $`\chi`$ and $`g_D`$ are constant. For a nonnegative constant multiplier $`m`$, the lower terminal for the reward $`mb`$ is $`mb/(\lambda+\eta)`$ and the upper is $`m[1/\lambda-(1-b)/(\lambda+1)]`$.

## 2. One old child and two intrinsic residence times

Put $`X(Y)=\sum_{Z\ {\rm child\ of}\ Y}a_ZW(Z)`$, $`v=s+\ell`$, and $`c_\ell=\kappa\alpha/(v-q)`$. For a root of mark $`\ell`$ define the first two tables
```math
t_\ell^-=\mathcal H^-_{v,0}(b g_D),\qquad
 t_\ell^+=\mathcal H^+_{v,0}\big(b[g_D+c_\ell\chi(D)]\big). 
```
The actual original-particle equation implies
$`t_\ell^-(D,c)\le U_\ell(Y)\le t_\ell^+(D,c)+c_\ell X(Y)`$. Indeed current copies have total score at most $`\alpha\kappa X(Y)`$. Each old child's expected weight grows at most $`e^{qt}W`$ while attached. New children begin with $`a_Z W(J)=\chi(D)`$ at rate $`b`$. Integrating these two contributions gives the stated coefficient.

For an original root of mark $`a`$ let $`\lambda=s+a`$. A child of mark $`\ell`$ being copied from it can receive a plus leaf, so set $`\widehat t^-=\min(t^-_\ell(D,c),t^-_\ell(D+1,c))`$ and $`\widehat t^+=\max(t^+_\ell(D,c),t^+_\ell(D+1,c))`$. Set
```math
 H_{a,\ell}^-=\mathcal H^-_{\lambda,\ell}\widehat t^-,
 \quad
 H_{a,\ell}^+=\mathcal H^+_{\lambda,\ell}
 \left[\widehat t^++c_\ell b\chi(D)/(\lambda-q)\right].
```
Before copying, the input from the current parent is at most $`\ell`$. An extra leaf at the eventual copy leaves $`X`$ unchanged. Old grandchildren contribute at most $`c_\ell X/(\lambda-q)`$ after integration, and this term vanishes at a newborn. These two tables thus bound the integrated copy score per unit rate at a newborn child. A valid lower terminal at $`D=N`$ is $`b g_N/[(v+\eta)(\lambda+\eta)]`$ since a plus leaf raises $`b`$ and decreases $`\eta`$. A valid upper terminal is $`\{[g_N+c_\ell\chi(N)]/v+c_\ell\chi(N)/(\lambda-q)\}/\lambda`$, using $`b\le1`$. In each case the boundary amounts apply to all higher degrees as well.

Finally apply (1) at the newborn root with silent parent and running reward $`b[g_D+a(D)H_{a,a(D)}^\pm(2,1)]`$. This gives true row bounds for $`(K_sg)(a)`$. To bound the row difference directly, use the original-root lifetime resolvent identity
```math
(K_sg)(\beta)-(K_sg)(A)=(A-\beta)(V_{s+\beta}^0U_A)(J).
```
Use precisely the preceding $`H_{\beta,A}^\pm`$ with $`e=0`$, and omit the min and max with the plus state. The same substitution with $`(a,\beta)`$ and prefactor $`\beta`$ bounds the derivative or divided difference $`P_{s+a}(J_\beta)`$ in the high-degree continuation. All the tables can round lower values down and upper values up to integer multiples of $`10^{-12}`$. [Implementation](binary_rows.py)

## 3. A uniform positive denominator on the known bracket

Write $`k_j=K_s(A,j)`$ and $`d_j=K_s(\beta,j)-k_j`$. The measures defining $`k_j,d_j`$ are nonnegative chronological Laplace transforms. At $`u=0.5665`$, the exact lower units divided by $`10^{12}`$ for $`k_A,k_\beta,d_A`$ are respectively $`471314072924`$, $`49408249826`$, $`193161058005`$. They are also lower bounds throughout $`[l,u]`$ for $`l=0.565`$. At $`s_0=0.54001`$ a coupled test $`g=(1/5,1)`$ gives $`r(K_{s_0})\le5.434538335003`$. The bound applies the controlled equations directly to the test vector.

Log convexity of $`r=r(K_s)`$, monotonicity and $`r(K_\mu)=1`$ imply that $`r\le5.435^{(u-s)/(u-s_0)}`$ on the prior interval. At its lower end the exponent is $`50/883`$, giving $`r<1.10061`$. Normalise the Perron vector to $`(x,1)`$, where $`x=k_\beta/(r-k_A)>0`$. Subtracting the two row equations gives
```math
d_\beta=r-(r+d_A)x.
```
Using the three positive lower bounds and replacing $`r`$ by its upper bound is valid because $`R-(R+d_A)k_\beta/(R-k_A)`$ increases in $`R>k_A`$. Rational power comparisons give

| Point | Upper bound on $`d_\beta`$ |
|---|---:|
| Every $`s`$ in the known bracket | 0.9990328 |
| $`s=0.565615`$ | 0.9529624 |
| $`s=0.565809`$ | 0.9386957 |
| $`s=\mu`$ | 0.8884949 |

Thus the scalar first-return formula $`k_A+k_\beta(1+d_A)/(1-d_\beta)`$ has a strictly positive denominator throughout the prior bracket. All its ingredients have positive path expansions except for the formal subtraction used to *name* the excess rows, whose positive lifetime representation is given in section 2. These inequalities use the existing exponent bracket and give no new two-sided signs at the two interior endpoints. [Exact comparisons](return_constants.py) · Input and replay

There is a sharper interpolation using $`d_\beta`$ itself. Its row-difference upper bound at $`0.54001`$ is $`4.999797019077<5`$. Since $`d_\beta`$ is a positive Laplace transform and $`d_\beta(\mu)<0.888495`$, putting $`p=(u-s)/(u-0.54001)`$ gives
```math
d_\beta(s)\le 5^p(0.888495)^{1-p}\qquad(0.54001\le s\le\mu).
```
For $`s\in[\mu,u]`$ monotonicity gives the same upper bound. The controls therefore give $`d_\beta<0.97982`$ throughout the old bracket, $`d_\beta<0.94130`$ at $`0.565615`$ and $`d_\beta<0.92947`$ at $`0.565809`$. Slightly below the bracket, $`d_\beta(0.56470)<0.99917`$ and $`d_\beta(0.56480)<0.99268`$. These estimates give [a uniform time-domain gap](dormant-return.md).
