#!/usr/bin/env python3
"""Independent integer corner verification and unbounded current-gauge probes."""
from fractions import Fraction as Q
from functools import cache
from math import lcm
from pathlib import Path
import sys

from current_gauge_sol import J, degree, gain, moves


def main():
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / 'current-mark' / 'weight'))
    from uncapped import exact_inputs  # only the proposed edge factors

    M, d, rho, A = 400, Q(73,200), Q(53,100), Q(1,2)
    gamma, C = rho-d, (d+A)/d**2
    ticks = exact_inputs(M,d)
    P = lambda i: ticks[i] if i < M else C/i
    z = lambda i: 1/(1+d*i)
    k = lambda i: z(i)*z(i+1)
    assert all(1 >= P(i) >= max(P(i+1),z(i)) for i in range(2,M+1))

    # (P_j, deltaP_j, k_j) on a common integer denominator; the tail
    # normalizes by the true P_j, never assuming P_j >= P_M.
    kids = []
    for j in range(2,M+1):
        abc = (P(j),P(j)-P(j+1) if j<M else Q(0),k(j))
        den = lcm(*(x.denominator for x in abc))
        kids.append((*map(lambda x:int(x*den),abc),j))
    K = max(k(j)/P(j) for j in range(2,M+1))
    assert d*M >= 1 and C-1/d >= K and C*d*d == d+A
    n = 0
    for i in range(2,M):
        for a in (Q(0),A):
            u = a+d
            x, y = u*(P(i)-z(i)), P(i)-P(i+1)-u*k(i)
            assert x >= 0 and x+(i-Q(3,2))*y >= 0
            n += 2
            span, paid, charge = x+(i-2)*y, P(i)/i, u/i
            den = lcm(*(t.denominator for t in (x,span,y,paid,charge)))
            X,S,Y,B,D = (int(t*den) for t in (x,span,y,paid,charge))
            for pv,slip,drop,j in kids:
                lo = pv*(X*j+Y)+j*(B*slip-D*drop)
                hi = lo + pv*(S-X)*j
                assert min(lo,hi) >= 0, (i,j,a)
                n += 2
                if j == M:
                    assert min(lo,hi)-pv*Y >= 0, (i,j,a,'infinity')
                    n += 2
    assert n == 638392

    H = {M:[P(M)/gamma]*(M+1)}
    m = 0
    for i in reversed(range(2,M)):
        row = H[i] = [None]*(i+1)
        for c in reversed(range(1,i+1)):
            L = i-c
            b = Q(L,i)
            next_promo = b*(row[c+1]+P(i)) if L else Q(0)
            nxt = H[i+1][c]
            gains = (L,L+Q(c,2)+A+d)
            need = max(Q(1,100),*((r*nxt+next_promo)/(gamma+r+b) for r in gains))
            row[c] = Q(-(-need*10**9//1),10**9)
            if (i,c)==(2,1):
                assert row[c] == Q(999609761,10**9)
                row[c] = Q(1)
            for r in gains:
                assert r*nxt+next_promo <= (gamma+r+b)*row[c]
                m += 1
    hmin = min(x for row in H.values() for x in row[1:])
    assert hmin == Q(86271979,10**9) and max(x for row in H.values() for x in row[1:]) == 1
    assert m == 159598
    print('independent exact edges/local endpoints',n,m,'hmin',hmin)

    @cache
    def V(S):
        i=degree(S)
        return H[min(M,i)][min(M,1+len(S[1]))]+P(i)*sum(map(V,S[1]))

    def wp(S,l,u):
        return V(S)/(l+d) if l==u else V(S)/((l+d)*(u+d))

    # Use individual physical arrows, sparse arity, and fresh label-zero
    # states; separately test very large degrees without materializing stars.
    for i,j in ((3,4),(23,19),(M-1,M),(M,M-1),(M+1,M+2)):
        S = (i-2,((j-2,(J,)),))
        for a in (Q(0),Q(1,1000003),A):
            w=wp(S,a,a)
            continuation = sum(r*(wp(Y,a,a)-w) for r,Y in moves(S))-a*w
            # Use the true degrees here; no cap argument is consulted.
            child=S[1][0]; p, pj = Q(1,i),Q(1,j)
            continuation += p*wp(gain(child),p,p)
            continuation += (i-2+pj)*k(i)*V(child)+p*k(j)*V(J)
            assert w >= hmin/(A+d) and continuation <= rho*w
            assert wp(S,Q(1,j+1),pj) >= hmin/(A+d)**2
            rawdiff=1/(Q(1,j+1)+d)-1/(pj+d)
            assert rawdiff*V(S)==(pj-Q(1,j+1))*wp(S,Q(1,j+1),pj)

    def margin(i,j,r,a,wrong=False):
        u=a+d
        p=Q(1,i)
        return (u*(P(i)-z(i))+r*(P(i)-P(i+1)-u*k(i)))*P(j) + p*P(i)*(P(j)-P(j+1)) - (Q(1,j) if wrong else p)*u*k(j)

    for i in (M-1,M,M+1,10007,1000003):
        for j in (2,19,M-1,M,M+1,200003):
            for a in (Q(0),A):
                for r in (Q(1,j),i-2+Q(1,j)):
                    assert margin(i,j,r,a)>=0
    assert margin(1000003,2,Q(1,2),A,True)<0
    assert (C/4-1/d)<K  # a weakened tail cannot meet the analytic payment
    old=Q(565662790716854568223,1000008252270910000000)
    assert gamma<rho<old and Q(11299,20000)<old
    # Exact Taylor tail bound, including a=0.
    s=Q('.56572365')
    L=(A+d)/(hmin*(s-rho))
    B=1/hmin+L*(1/hmin+(2+A/d)/d)
    bound=B/(s-gamma)*(1/(20*(s-gamma)))**16
    assert bound<Q('5.044e-11')
    print('zero and near-zero marks; two-sided degree seams; true/weak tails; Taylor bound',float(bound))


if __name__=='__main__':
    main()
