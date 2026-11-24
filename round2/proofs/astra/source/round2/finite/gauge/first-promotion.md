# An exact first-promotion integral and a summable boundary error

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## 1. All star degrees in one scalar series

Start the original residual particle from $`J_a`$, write $`\lambda=s+a`$ and let $`g`$ be a nonnegative newborn test with $`\|g\|_\infty\le1`$. The direct reward and the intrinsic root promotion have the same predictable rate. Before its first intrinsic promotion this particle has no core children, so it cannot emit residual copies. If its degree is $`k`$, its next leaf gain has rate $`k-1`$ and its promotion has rate $`(k-1)/k`$. Let $`q_k`$ be the expected factor $`e^{-\lambda\tau}`$ on the event that its first promotion is at degree $`k`$ at time $`\tau`$, under the intrinsic process. Writing $`a(k)`$ for the birth mark and $`V_k`$ for the postpromotion tree with one virgin core child, conditioning gives
```math
U_sg(J_a)=\sum_{k\ge2}q_k\,[g(a(k))+U_sg(V_k)]. \qquad\text{(1)}
```
Here $`V_k`$ retains root mark $`a`$ with $`k-2`$ old leaves and one $`J_{a(k)}`$. The contribution from any trajectory still waiting for its first promotion vanishes in this decomposition at positive $`\lambda`$. The factors are
```math
q_2=\frac1{3+2\lambda},\qquad
 \frac{q_{k+1}}{q_k}=\frac{k^2}{k^2+(2+\lambda)k+\lambda}. \qquad\text{(2)} 
```
Put $`u+v=\lambda`$, $`uv=-1`$. The discounted mass $`T_K=\sum_{k\ge K}q_k`$ is the single classical function
```math
T_K=q_K\,{}_3F_2(1,K,K;K+1+u,K+1+v;1). \qquad\text{(3)}
```
For a real ordinary integral of elementary functions, take the positive root $`u`$ and the negative root $`v\in(-1,0)`$. Then
```math
\frac{T_K}{q_K}=
 \frac{1}{\mathrm B(K,1+u)\mathrm B(K,1+v)}
 \int_0^1\!\int_0^1
 \frac{(xy)^{K-1}(1-x)^u(1-y)^v}{1-xy}\,dx\,dy. 
```
Expanding the geometric factor and applying Tonelli gives the coefficient ratio (2) term by term. Equivalently, integrating one variable leaves a one-dimensional Euler integral of $`\,{}_2F_1`$. Thus the direct part in (1) is one tail in the constant gauge, and two tails in the binary gauge. The continuations in (1) retain the physical interaction of the root and its new child. Formula (3) does not evaluate those continuations.

The following rational bounds give certified values without a special-function library. Put $`\gamma=\lambda+2/(\lambda+2)`$. Substitution of (2) gives
```math
 (k+\gamma)q_k-(k+1+\gamma)q_{k+1}
 =(1+\lambda)q_k-\frac{\lambda^2q_{k+1}}{(\lambda+2)k^2}.
```
The ratio product in (2), or elementary summation of its logarithms, gives $`q_k=O(k^{-(2+\lambda)})`$ at each $`\lambda>0`$. Telescoping and $`q_{k+1}\le q_k`$ now give
```math
 \frac{(K+\gamma)q_K}{1+\lambda}\le T_K\le
 \frac{(K+\gamma)q_K}{1+\lambda-\lambda^2/[(\lambda+2)K^2]}.
 \qquad\text{(4)}
```
An exact prefix can be added before applying (4). At $`s=0.565615`$ and $`\lambda=s+1/2`$, $`T_2`$ is between 0.35490752 and 0.35490753 and $`T_{16}>0.00547946`$. For $`\lambda=s+1/16`$, $`T_2`$ is between 0.49208220 and 0.49208222 and $`T_{16}>0.01727635`$. The $`q_k`$ are discounted masses, so these are values of the first-promotion contribution only. [Rational replay](transfer_check.py)

## 2. A common continuation for all large first-promotion degrees

Use the binary gauge of [binary-and-return.md](binary-and-return.md), with the basic certificate $`q=49/125`$, $`d=4/25`$, $`\beta=1/16`$, $`\alpha=1/h_{\min}`$, and $`\kappa=\max(1,\sup W(Y^+)/W(Y))`$. For $`s>q+d`$ we have $`0\le U_sg\le\alpha W`$. Let $`V_\lambda^0=(\lambda-L^0)^{-1}`$ and put $`P_\lambda(Y)=a_Y V_\lambda^0 U_sg(Y)`$. Write $`H(x)`$ for the newborn response with its original root killing rate changed to $`x`$, keeping every descendant mark fixed. The resolvent identity gives
```math
P_\lambda(J_\beta)=\beta[H(\beta)-H(a)]/(a-\beta),\qquad
 P_{s+\beta}(J_\beta)=-\beta H'(\beta). \qquad\text{(5)} 
```
For the last derivative one can differentiate the original-particle lifetime integral or enclose the one-parameter integral as in the [binary row tables](binary-row-bounds.md).

At a tree whose root already has degree at least 16, put
```math
F(S)=[g(\beta)+P_\lambda(J_\beta)]/\lambda+\sum_{Y\ {\rm old}}P_\lambda(Y).
```
The parent sends at actual rate $`p=1/D`$. Subtract the intrinsic resolvent equation for $`F`$ from that for $`U=U_sg`$. At each old child $`Y`$, its own $`L^0_Y`$ contribution cancels $`a_Y U(Y)`$. The remaining terms are
```math
(\lambda-L^0)(U-F)
  =(b-1)[g(\beta)+P_\lambda(J_\beta)]
    +p\sum_Y[(U+P_\lambda)(Y^+)-(U+P_\lambda)(Y)]. \qquad\text{(6)} 
```
Here the root degree is nondecreasing, so all newly born children have mark $`\beta`$. The descendant degree bound for $`W`$ includes a child's actual incoming rate from its parent. Therefore $`\mathbb E W(Y_t)\le e^{qt}W(Y)`$ before a copy, and
```math
|\Delta(U+P_\lambda)(Y)|\le
 \alpha\kappa[1+a_Y/(\lambda-q)]W(Y). 
```
Indeed both nonnegative values at $`Y,Y^+`$ are bounded by the right side. Equation (6) follows under localisation. The old-child sum at a bounded time is polynomially bounded by the Yule size, and its discounted mean tends to zero at infinity for $`\lambda>q`$.

Put $`w=1/(\beta+d)`$, $`A_*=1+\alpha\beta w/(\lambda-q)`$. At $`V_k`$ for $`k\ge16`$ there is one old child of weight $`w`$, $`c_0=2`$, and the parent emits new children at rate at most one. We have $`|b-1|\le c(t)/k`$, $`\mathbb E c(t)\le2+t`$. The integral of the expected weight of all children is at most $`w(1+1/\lambda)/(\lambda-q)`$. Thus
```math
 \left|U_sg(V_k)-\frac{g(\beta)+P_\lambda(J_\beta)}{\lambda}
                    -P_\lambda(J_\beta)\right|
 \le B_\lambda/k,
```
```math
B_\lambda=A_*(2/\lambda+1/\lambda^2)
 +\alpha\kappa[1+\beta/(\lambda-q)]w(1+1/\lambda)/(\lambda-q).
```
Replace every continuation at $`k\ge K\ge16`$ in (1) by this common value and keep the entire mass $`T_K`$. The resulting absolute error is at most $`B_\lambda T_K/K=O(K^{-(2+\lambda)})`$. In particular the summable bound applies after removing all arbitrarily large degrees in this first phase. It does not truncate later descendant histories. Concrete one-sided errors are in the replay.
