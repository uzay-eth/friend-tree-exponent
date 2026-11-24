# Intrinsic transfer and an obstruction to a rational newborn response

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## 1. The deletion process on isolated fringes

For a finite planted tree $`S`$, let $`E f(S)=\sum_Y f(Y)`$, summed over its immediate core children. Put $`T f(S)=\sum_{x\ {\rm core}}f(S_x)`$. On the algebra of scores on finite trees, $`T^{-1}=I-E`$. The plain intrinsic generator $`L^0`$ runs the physical tree with its parent slot silent. Write $`M_a f(S)=a_Sf(S)`$, where $`a_S`$ is the permanent mark of the root. The full positive fringe generator is $`P_a=R_a+B`$ from [the renewal construction](../../fringe/renewal.md).

Under $`L^0`$, each old strict descendant receives one input from its parent missing from its standalone generator. Its signed increment is $`p[f(Y^+)-f(Y)]`$. Each actual promotion creates one new summand. These are respectively the signed part of the copy terms and the direct birth $`B`$. The remaining $`a_Y f(Y)`$ in each parent cancels the killing of that child, leaving just root killing. It follows that
```math
 TP_a=L^0 T-M_a,\qquad TP_aT^{-1}=Q_a=L^0+M_a(E-I). \qquad\text{(1)}
```
The latter has a direct independent-particle interpretation. Run each isolated fringe intrinsically. At rate $`a_S`$ remove its root and release all the old core-child subtrees as new isolated particles, retaining the parent slots and their entire old descendant trees. Leaf children are discarded. To compare gauges $`a,a'`$, record both rules on one common decorated tree so that the same $`L^0`$ implements all intrinsic changes. Write $`\mathcal G_a=(s-Q_a)^{-1}`$. In the norm $`\|f\|_n=\sup_S|f(S)|/n(S)`$, $`E`$ has norm at most one, and the total forest size under $`Q_a`$ has expected growth at most $`e^t`$. For real $`s>1`$, the bounded perturbation $`(M_a-M_{a'})(E-I)`$ therefore gives the explicit comparison
```math
 \mathcal G_a-\mathcal G_{a'}=\mathcal G_a(M_a-M_{a'})(E-I)\mathcal G_{a'},\qquad
 \|\mathcal G_a-\mathcal G_{a'}\|_n\le
 \frac{2\sup_S|a_S-a'_S|}{(s-1)^2}.
```
One obtains the identity by conditioning between perturbation times, or by the localised time-dependent jump formula followed by Laplace integration. These formulas compare distinct deletion generators.

For a geometric depth factor $`z>0`$, let $`T_z f=\sum_x z^{|x|}f(S_x)`$ and write $`C_af=\sum_Y[(a_Y-p)f(Y)+pf(Y^+)]`$. The weighted identity is
```math
T_z\{L^0-M_a+z(C_a+B)\}=L^0T_z-M_a,\qquad T_z^{-1}=I-zE.
```
Indeed the intrinsic signed input and the birth of a new summand are now charged at the child's depth. These identities are pointwise. $`T`$ need not be bounded in the original residual weight norm.

For a single constant mark $`A`$, $`B=b\,{\rm ev}_J`$ has rank one. Since $`({\rm ev}_J)T^{-1}={\rm ev}_J`$, (1) gives $`TR_AT^{-1}=Q_A-(Tb){\rm ev}_J`$. In particular, with $`V_\lambda=(\lambda-L^0)^{-1}`$, define the intrinsic deletion occupation
```math
G_A(s)={\rm ev}_J(s-Q_A)^{-1}Tb,\qquad
 (s-Q_A)^{-1}=\sum_{m\ge0}(V_{s+A}A E)^m V_{s+A}. \qquad\text{(2)} 
```
Each term conditions on successive root deletions and carries all the actual child shapes. In the size norm and for $`s>1`$, $`\|V_{s+A}A E\|_n\le\theta=A/(s+A-1)<1`$. The error from keeping terms $`0,\ldots,m-1`$ in (2) is at most $`\theta^m/(s-1)`$ in operator norm. This large-$`s`$ estimate does not reach the physical exponent. If $`K_A`$ is the residual newborn row applied to the unit reward,
```math
G_A=K_A+K_A G_A,\qquad K_A=G_A/(1+G_A). \qquad\text{(3)}
```
Here $`K_A(s)={\rm ev}_J(s-R_A)^{-1}b`$. Formula (3) follows by grouping full descendants by generations of direct births.

For details of integrability, take $`A=1/2`$ and $`n(S)`$ the number of vertices including leaves. We have $`L^0 n=n-p`$, and $`\sum_Y n(Y)=n-1-L`$. Hence $`R_A n=n-p-A(1+L)+p\#Y\le n+1`$ and $`P_A n\le n+3\le(5/2)n`$. The killed-and-replaced forest satisfies $`Q_A n\le n`$. Localisation of these drift inequalities controls total size and prevents explosion, since all event intensities are bounded by a multiple of population size. It also bounds full or residual particle counts on finite time intervals. To pass (1) to first moments for $`f=b`$, start with bounded nonnegative $`f`$ and the time-dependent full mean $`P_tf`$. On finite times it and its derivative are bounded in absolute value by constants times $`n`$ and $`n^2`$, respectively, by the first-event equation. Their subtree sums have polynomial growth. The total size of the $`Q_A`$ forest increases at rate at most total size and decreases at each deletion, so it is dominated by Yule size with finite moments of every order at bounded times. The stopped jump identity gives $`T P_tf=e^{tQ_A}T f`$. The nonnegative deletion construction gives (2) by monotone convergence. For example, all Laplace and renewal manipulations above hold for real $`s>3`$. The formula also holds wherever the individual positive integrals in it are finite.

The same algebra with $`A=0`$ gives an exact transfer in terms of an ordinary isolated particle. Its copy generator can have signed entries, so this special case is an algebraic resolvent identity. On the space of linearly size-bounded scores, $`Q_0=L^0`$, $`Tb\le n`$, and evaluation at $`J`$ is bounded. The perturbation $`-(Tb){\rm ev}_J`$ is bounded there. Its resolvent series therefore converges for all sufficiently large real $`s`$. On each fixed finite tree the transformed $`R_0`$ response agrees with this construction. With $`C`$ the physical core count,
```math
 G_0(s)=\int_0^\infty e^{-st}\mathbb E_J^0\sum_x b_x(t)\,dt
       =s\int_0^\infty e^{-st}\mathbb E_J^0 C(t)\,dt-1,
 \qquad K_0(s)=1-\left[s\int_0^\infty e^{-st}\mathbb E_J^0 C(t)\,dt\right]^{-1}.
```
The second equality in $`G_0`$ uses the actual intensity of intrinsic promotions and $`C(0)=1`$. It removes the unknown fringe response in favour of a single observable of an ordinary isolated tree. Evaluating that observable is still an infinite-state problem.

## 2. Algebraicity of the scalar response as a function of the Laplace parameter

Let $`J`$ be the planted newborn of degree 2 and let $`b=L/D`$. Fix $`A=1/2`$. The derivatives $`d_k=(R_A^k b)(J)`$ are rational. These are indeed the right time derivatives at zero. Each fixed generator iterate is polynomially bounded in size because the rates and number-weighted copy contributions grow at most linearly with $`n`$. Intrinsic changes increase size by one and copies have size at most the old size. For $`m\ge1`$, the exponential first-moment upper bound for each fixed residual size power follows from $`\sum_Y n(Y)^m\le n^m`$ and $`\sum_Y p(n(Y)+1)^m\le 2^m n^m`$ and the $`O(n^{m-1})`$ one-step increment of $`n^m`$ for intrinsic changes. Localisation and one more moment justify successive jump identities at zero.

For every odd prime $`p`$, put $`k=p-2`$. After at most $`j`$ successive generator factors starting at $`J`$, the current size is at most $`2+j`$. All diagonal and off-diagonal coefficients before the last action therefore have denominators prime to $`p`$. A terminal reward can have denominator divisible by $`p`$ only if its root has degree $`p`$. Size is then exactly $`p`$ and the root is a star. A core child would have an additional descendant and would require size at least $`p+1`$. The only contributing path consists of $`p-2`$ intrinsic leaf gains from the root's own leaves. Every gain increases size, and an intrinsic promotion would permanently create a core child. Consequently
```math
 d_{p-2}=(p-2)!\,(p-1)/p+\hbox{a rational number with denominator prime to }p.
```
Thus the reduced denominators of these derivatives have infinitely many distinct prime divisors.

A rational Laplace transform for an exponentially bounded response would give a constant-coefficient homogeneous linear differential equation in time and hence a constant-coefficient recurrence on its derivatives, possibly starting after finitely many terms. If a real recurrence exists, the equations for its monic coefficients have rational entries. Solving a maximal independent finite subset by rational elimination gives a rational solution of all the equations. A rational recurrence introduces only the finitely many denominator primes in its coefficients and starting block. This contradicts the preceding calculation. In particular the obstruction applies to the fixed scalar newborn output, even if only this one initial condition is of interest.

The same unique star term works for the response $`b`$ under $`L^0`$, $`P_A`$, or $`Q_A`$. Their other transitions cannot increase size by one and reach the terminal star except through the same intrinsic gains. For any fixed finite rational birth-mark scheme, exclude the finitely many primes dividing a mark denominator and take the newborn output summed over all labels. The proof still applies to that row sum under $`R`$. The same argument applies to the scalar return of the [active/dormant construction](dormant-return.md). A dormant emission or activation cannot give a size increase of one on a path to the final star. Its running reward is $`b`$ on active states and zero on dormant states. The result excludes rational Laplace transforms and finite time-homogeneous linear systems with constant coefficients. It leaves open general differential equations with variable coefficients, nonlinear finite systems, and integrals of classical functions.

In fact the scalar Laplace function itself is transcendental over $`\mathbb C(s)`$. To see this, let $`H(z)=\sum_{k\ge0}d_k z^k`$ formally. A complete large-$`s`$ asymptotic expansion of the scalar response is $`s^{-1}H(1/s)`$, by its successive right derivatives at zero. For each fixed order, Taylor expansion on a short initial interval and the exponential tail bound justify the Laplace expansion. A polynomial relation over $`\mathbb C(s)`$ would thus give a nonzero polynomial relation for $`H`$ as a formal series. Since its coefficients are rational, linear dependence of the finitely many participating monomials over $`\mathbb C`$ implies their dependence over $`\mathbb Q`$ by ordinary rational elimination. Clear denominators and take an irreducible polynomial $`P(z,H)=0`$ over $`\mathbb Q`$. In characteristic zero $`P_y(z,H)`$ has finite order $`r`$ with initial coefficient $`c\ne0`$. For $`n>r`$, set $`H_{<n}=\sum_{k<n}d_kz^k`$. The coefficient of degree $`n+r`$ in the polynomial identity gives
```math
c\,d_n=-[z^{n+r}]P(z,H_{<n}).
```
Higher Taylor powers of $`H-H_{<n}`$ start at degree $`2n>n+r`$. Consequently the denominators of every $`d_n`$ can involve only primes from $`c`$ and the finite initial block, contradicting the star path.

**This is a statement about the scalar response as a function of $`s`$. It gives no implication for the algebraicity or transcendence of the critical root $`\mu`$.** The proof does not decide the existence of a classical integral or an ODE with variable coefficients for the complete response.

Two related obstructions concern universal state compression. At a star of arbitrary large degree $`D`$, $`R_A^k b`$ is a rational function of $`D`$. Its pole at $`D=-k`$ is simple with residue $`(-1)^{k+1}(k+1)!`$. Only $`k`$ consecutive root gains introduce the final factor $`(D+k)^{-1}`$, since intermediate rates use root degrees at most $`D+k-1`$ and all other degrees remain bounded independently of $`D`$. Lower iterates lack this pole. The iterates are linearly independent. Also, on the sector consisting of a root of degree $`d>2`$ with one virgin-or-older star child of degree $`e`$, the immediate first-event reward is
```math
\frac{(d-2)/d}{s+A+d+e-1-1/d}.
```
Any square matrix of these values on distinct $`d`$ and distinct $`e`$ has nonzero Cauchy determinant for large positive $`s`$. Thus even this first-event term has unbounded separated rank in the two degree variables. This does not give a rank lower bound on its sum with later terms.

The checker evaluates (1) and its weighted version on 388 finite cases and checks the unique-prime contribution at $`p=3,5,7`$. [Arithmetic replay](transfer_check.py)
