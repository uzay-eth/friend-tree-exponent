# Moments of the two successive placements

In the original edge-seeded random friend tree, fix the history through size \(n\ge2\). For **every real \(p\ge3\)**, the older of the next two births has strictly greater conditional \(p\)-th moment of its limiting degree share. The first and second moments have separate sharper estimates using the same serial-middle calculation. The bounds here compare the two future labels on every old tree, with no label cutoff. [First- and second-moment edge estimates](../mean-ledger/README.md)

## A common clock for comparing the two labels

Write \(m=n+2\) for the size after both births. Restart the rate-one clocks there and define \(W_i=\lim_{t\to\infty}e^{-t}D_i(t)\). Conditional on this deterministic-size completion, the total-size amplitude \(S\) has law \(\operatorname{Gamma}(m,1)\). All holding times, which have rates depending only on total size, are independent of the entire embedded tree chain. Thus \(W_i=S Z_i\) and \(\mathbb E W_i^p=(m)_p\,\mathbb E Z_i^p\) at each fixed completion, where \((a)_p=\Gamma(a+p)/\Gamma(a)\). We use
\[
 h_p(i)=\frac{\mathbb E W_i^p}{\Gamma(1+p)} .
\]
The normalization is common to every two-birth placement. At a later random event we restart the clocks and discount its remaining \(W\)-moment by \(e^{-p\tau}\). No size-amplitude factorization at such an event is required. [Physical clock and conditional hub comparison](../README.md)

In the exact two-birth identity substitute \(h_p\) for hub probability. With \(w_u=r_u/(D_u+1)\), \(s_u=w_u/D_u\), and \(A=\sum_u w_u>0\), each serial middle gets its weight \(w_u\). A serial tip has a degree-two parent. On an old-edge term only the adverse sign needs an upper bound. The sum \(X\) of the absolute coefficients whose adverse leaf attaches to an old leaf is at most \(A/3\). To see this, at a core of degree \(d\) with \(L\) old leaves, their total serial weight is \(A_l=L/(2d)\) and the core's own weight is at least \(L/(d+1)\). For \(L>0\), their adverse coefficient is at most \(A_l(1-2L/(d+1))_+\le A_l(d-1)/(d+1)\). Sum \(3X\le A_l+w_u\) over cores; on the edge seed \(X=0\). The sum of *all* absolute edge coefficients is at most \(\sum_u D_u s_u=A\). A new leaf attached to an old core has parent degree at least three. It follows that a uniform middle lower \(M_p\) and two tip uppers \(U_{2,p}\ge U_{3,p}\) give the normalized signed numerator lower
\[
 A\left(M_p-\tfrac43 U_{2,p}-\tfrac23 U_{3,p}\right).
\]
The factor from this numerator to the conditional degree-share moment gap is \(\Gamma(1+p)/[n(n+1)(m)_p]\). [Signed two-birth identity](../../adjacent/README.md) · [Paired-edge coefficients](../../global-transport/annealed/extension/README.md)

## The serial middle

The middle of a serial placement starts with degree two and one tracked leaf. At true degree \(d\), track \(j\) leaves. They produce fresh leaves at rate \(j\). Their actual losses by promotion occur at rate \(j/d\). Additional physical gains at the focal vertex advance both \(d\) and \(j\). At a terminal degree 30 use the following normalized lower bound, valid when \(p\ge1\):
\[
 \phi_p(0)=\phi_p(1)=0,\qquad
 \phi_p(j)=\frac{(j-1)_p}{\Gamma(1+p)}=\prod_{i=1}^{j-2}(1+p/i),\quad j\ge2 .
\]
In all subsequent time the tracked stock gains at least \(j\) and loses at most one per unit time. For \(j\ge2\), its lower generator applied to \(\phi_p\), less \(p\phi_p\), is
\[
 p\phi_p(j)\Big(\frac1{j-1}-\frac1{j+p-2}\Big)\ge0.
\]
The boundary rows at zero and one also have nonnegative drift. At an arbitrary fixed age, the process \(e^{-pt}\phi_p(j(t))\) is a submartingale after localization. The rescaled total-size Yule process has bounded moments of every fixed positive order, including orders above \(p\). It supplies uniform integrability for the localization and the large-age limit. Since \(j(t)\le D(t)\) and the latter has physical limit \(W\) after rescaling, \(\mathbb EW^p/\Gamma(1+p)\ge\phi_p(j(0))\).

When the last tracked leaf is actually promoted, its resulting degree-two core can send a fresh leaf back. Before that return, a neighbour of degree \(k\) sends to the tag at rate \(1/k\), and its own gain rate is at most \(k-1+1/d\). Stop tracking at \(k=40\). The discounted return lower bound is
\[
 B_d(40)=0,\qquad
 B_d(k)=\frac{1/k+(k-1+1/d)B_d(k+1)}{p+1/k+k-1+1/d}.
\]
These values decrease as \(k\) increases. Indeed, even a permanent degree \(k+1\) gives a discounted return of at most \(1/[1+p(k+1)]\). Backward induction gives \(B_d(k+1)\le1/(1+pk)\), which implies \(B_d(k)\ge B_d(k+1)\) and \(B_d(k)\le1/(1+pk)\). Extra focal gains before the return also supply a fresh leaf. The true process therefore earns at least this comparison at the **actual** loss of its last tracked leaf.

Work backwards from \(V_p(30,j)=\phi_p(j)\), using
\[
 V_p(d,0)=B_d(2)V_p(d+1,1),\qquad
 V_p(d,j)=\frac{jV_p(d+1,j+1)+(j/d)V_p(d,j-1)}{p+j+j/d}.
\]
At zero the displayed value applies in the presence of the newly promoted neighbour. The value is nondecreasing in \(j\) and \(V_p(d,j)\le V_p(d+1,j+1)\). For the induction, write \(v_j=V_p(d,j)\), \(g_j=V_p(d+1,j+1)\); the next degree row is increasing by induction. Initially \((1+p)v_0\le V_p(d+1,1)\le g_1\), using \(B_d(2)\le1/(1+2p)\). At \(j\ge2\), the recurrence in the preceding column and \(v_{j-2}\le v_{j-1}\) imply \(g_{j-1}\ge(1+p/(j-1))v_{j-1}\). Thus \(g_j\ge(1+p/j)v_{j-1}\), which gives \(v_{j-1}\le v_j\le g_j\). Every additional physical gain has nonnegative drift in the comparison. Each monitored stage ends or advances before a finite number of bounded cutoffs. Discounted comparison starting at the serial middle gives \(M_p=V_p(2,1)\). The exact directed lower values at \(p=1,2,3\) are respectively \(0.563194756241683,\ 0.508336964997524,\ 0.480548462305765\). [Source of the physical last-leaf return](../../global-transport/hostile/README.md)

## Real powers from three through seven

At a parallel or serial tip the new leaf sends to its own parent at rate one up to the first promotion. Each such send increases the parent's degree. When its parent has degree \(k\), first promotion occurs at rate \(1/k\). At that event the tag has degree two and its permanent parent will have degree at least \(k\) thereafter. The tag's future gain rate is at most \(D-1+1/k\). A pure birth process at this rate has a Gamma limit of shape \(1+1/k\) when started from tag degree two, giving a normalized upper \(A_p(1+1/k)\) for \(A_p(r)=\Gamma(p+r)/[\Gamma(r)\Gamma(1+p)]\). This is an upper comparison even at zero leaf stock. [Physical first-promotion proof](../tip-moments/README.md)

The checker partitions \([3,7]\) into 2,048 closed intervals \([l,h]\). In the middle calculation use the increasing terminal value at \(l\), and use \(h\) in every discount denominator, rounding down. For the tips, \(A_p(r)\) with \(0<r<1\) is the \(p\)-th moment of \(\operatorname{Beta}(r,1-r)\), hence decreasing and convex in \(p\). On each interval the linear interpolation at \(l\) between the enclosing integer powers is an upper bound. Multiply it by \(1+h/r\) to bound \(A_p(1+r)\). At \(k\ge48\) the first promotion alone gives the uniform ceiling \(A_p(1+1/48)/(1+48p)\). For \(k<48\) take the larger of the bound at \(k+1\) and
\[
 \frac{A_p(1+1/k)/k+U(k+1)}{l+1+1/k},
\]
using the stated upper in its numerator and rounding up. Every additional gain by the parent can only reduce this decreasing upper bound. The resulting two classes and serial lower make \(M_p-(4U_{2,p}+2U_{3,p})/3>0\) on every interval. The accompanying 256 intervals on \([1,3]\) export serial lower values for use with the more detailed old-edge ledger.

## All powers above seven

There is also a direct, grid-free estimate. Couple the true serial middle's tracked stock from below with \(J(0)=1\), births at rate \(J\), deaths at rate one while positive, and absorption at zero. For \(q>0\) put \(P_q(j)=\Gamma(j+q)/[\Gamma(j)\Gamma(1+q)]\) for \(j\ge1\), \(P_q(0)=0\); take \(P_0(j)=\mathbf1_{\{j>0\}}\). Then
\[
 (\mathcal L-q)P_q=-P_{q-1}\quad(q\ge1).
\]
For \(u=1-e^{-t}\) and \(r_q(u)=e^{-qt}\mathbb E P_q(J_t)\), this gives \(r_q'(u)=-r_{q-1}(u)\). All \(r_q(0)=1\). For \(q\ge0\), \(0\le r_q\le1\) by comparison with the pure Yule process started at one, for which the normalized expectation equals one. Iterating the integral for real \(p\ge7\) and bounding the final remainder by \(1/7!\) gives \(M_p^{\rm true}\ge \sum_{j=0}^7(-1)^j/j!\). The physical limit and moment domination justify letting \(u\uparrow1\).

The simpler locked upper for a leaf at parent degree at least \(k\) is \(A_p(1/k)\) before promotion, by the same affine gain comparison from degree one. Since \(A_p(1/k)\) decreases in \(p\), the exact cost at 7 bounds every larger exponent. The normalized signed margin per unit \(A\) is greater than \(0.021\) throughout this range. This estimate uses neither a threshold cutoff nor a finite enumeration of old trees.

## Exact replay

Run `python3 -IB round2/open/age-order/hub-ranking/mellin/check.py`. It verifies the three integer-power inputs, reproduces the 256 lower-power rows byte for byte, proves the 2,048 signed interval inequalities, and checks the alternating constant and beta debit as exact fractions. The parameter cap of 40 is on one return attempt in the lower comparison. All physical degrees and all old tree sizes are covered.
