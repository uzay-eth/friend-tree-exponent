"""Exact path and collar controls against Astra's actual permanent 19-mark rows."""

from fractions import Fraction as F
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'fringe'))
from check import frac_mark, gain, h, intrinsic_transitions, newborn, node, weight
from finite_rows import LABELS, cap_marks, degree, size, step


def intrinsic(x):
    return {cap_marks(y): rate for y, rate in intrinsic_transitions(x).items()}


def residence(s, x):
    return s + frac_mark(x[1]) + size(x) - F(1,degree(x,1))


def direct(x):
    return F(x[0],degree(x,1))


def main():
    s = F(1133,2000)  # denominators at the top of the physical bracket
    common = F(1,2)/(s+22)**19
    paths = 0
    for i in (0,)+LABELS:
        for j in LABELS:
            x = newborn(i)
            mass = F(1)
            while degree(x,1)<j:
                y = gain(x)
                assert step(x)[y] == intrinsic(x)[y] == x[0]
                mass *= intrinsic(x)[y]/residence(s,x)
                x = y
            assert min(degree(x,1),20)==j
            assert mass*direct(x)/residence(s,x) >= common
            paths += 1
        x = newborn(i)
        y = node(0,i,(newborn(2),))
        u = gain(newborn(2))
        assert intrinsic(x)[y]==F(1,2) and step(y)[u]==F(1,2) and u not in intrinsic(y)
        assert tuple(size(t)-F(1,degree(t,1)) for t in (x,y,u))==(F(3,2),F(5,2),F(8,3))
        l2 = F(1,6)/(residence(s,x)*residence(s,y)*residence(s,u))
        l1 = intrinsic(x)[gain(x)]/residence(s,x)*direct(gain(x))/residence(s,gain(x))
        assert min(l1,l2) >= F(1,6)/(s+4)**3
        assert weight(gain(x)) == weight(x)
    assert residence(s,newborn(0)) == s+F(3,2)

    # Above the cap, the same row must use both true p and capped a.
    at_cap = node(23,0)
    promoted = node(22,0,(newborn(20),))
    assert degree(at_cap,1)==24 and intrinsic(at_cap)[promoted]==F(23,24)
    grew = gain(promoted)
    for x in (promoted,grew):
        p = F(1,degree(x,1))
        assert step(x)[newborn(20)]==F(1,20)-p>0
        assert step(x)[gain(newborn(20))]==p
        assert sum(intrinsic(x).values()) == size(x)-p

    # Independent capped Lyapunov arithmetic; the high-degree tail is uniform.
    q, d, rho = F(9,25), F(3,20), F(51,100)
    for D in range(2,20):
        for c in range(D+1):
            assert (1-F(c,D))*(h(c+1)-h(c)+1/(1+d*D)) <= q*h(c)
    assert 1/(1+20*d)<=q*F(7,10)
    # V_delta=n^delta W works on a complex product neighbourhood.
    r, delta, sigma = F(51,50), F(1,40), F(14,25)
    beta = rho+F(13,20)*(r-1)+delta*(1+q/2)
    assert beta==F(221,400)<sigma
    # Neumann germ uniform on Re s >= sigma, |z| <= 1/200.
    norm = F(1,200)/(d*(sigma-rho))
    assert norm==F(2,3)<1
    print(f'PASS: {paths} actual no-copy depth-one row/column paths (361 ordinary + 19 inert)')
    print('PASS: depth-two plus-copy path in all 20 rows; true-degree split above cap; zero inert killing')
    print(f'PASS: capped drift controls, V_delta exponent {beta} < {sigma}; early Neumann norm <= {norm}')


if __name__=='__main__':
    main()
