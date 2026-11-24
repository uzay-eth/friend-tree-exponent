"""Independent rational shared-future and moment fixtures."""

from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from math import comb, factorial


T3 = ((1, 2), (0,), (0,))


@cache
def steps(g):
    n = len(g)
    rates = defaultdict(Q)
    for adj in g:
        for v in adj:
            rates[v] += Q(1, len(adj))
    return tuple((q/n, tuple(adj+(n,) if i == v else adj for i, adj in enumerate(g))+((v,),))
                 for v,q in sorted(rates.items()))


def paths(g, end):
    if len(g) == end:
        yield Q(1), (g,)
    else:
        for p, h in steps(g):
            for q, trace in paths(h, end):
                yield p*q, (g,)+trace


def fixture(first, last, end, degree):
    labels = tuple(range(first, last+1))

    def hit(g, v):
        return int(len(g[v-1]) >= degree)

    @cache
    def endpoint(g, v):
        return sum(p*hit(trace[-1], v) for p,trace in paths(g,end))

    @cache
    def forecasts(g):
        n = len(g)
        terminal = remaining = Q(0)
        for p,trace in paths(g,end):
            terminal += p*sum(hit(trace[-1],v) for v in labels)
            for x in trace[:-1]:
                v = len(x)+1
                if v in labels:
                    remaining += p*sum(q*endpoint(y,v) for q,y in steps(x))
        direct = sum(p*sum(hit(trace[-1],v) for v in labels if v>n)
                     for p,trace in paths(g,end))
        assert direct == remaining
        past = sum(endpoint(g,v) for v in labels if v<=n)
        assert terminal-remaining == past
        return terminal, remaining

    def increments(trace, j, wrong=False):
        g,h = trace[j:j+2]
        m,p = forecasts(g)
        mm,pp = forecasts(h)
        rate = sum(q*endpoint(y,len(h)) for q,y in steps(g)) if len(h) in labels else 0
        inserted = endpoint(h,len(h)) if wrong and len(h) in labels else 0
        return mm-m, rate+pp+inserted-p

    square = bracket = Q(0)
    residuals = defaultdict(Q)
    cross = wrong_mean = Q(0)
    support = set()
    for p,trace in paths(T3,end):
        old = forecasts(trace[0])
        assert old[0] == old[1]
        a = sum(hit(trace[-1],v) for v in labels)
        marks = [hit(trace[-1],v)-endpoint(trace[v-3],v) for v in labels]
        cross += p*marks[0]*marks[1]
        residuals[trace[labels[1]-3]] += p*marks[0]
        pairs = [increments(trace,j) for j in range(len(trace)-1)]
        ds = [m-q for m,q in pairs]
        square += p*sum(ds)**2
        bracket += p*sum(z*z for z in ds)
        wrong_mean += p*sum(increments(trace,j,True)[1] for j in range(len(trace)-1))
        support.add(tuple(hit(trace[-1],v) for v in labels))
        assert sum(ds) == a-sum(pair[1] for pair in pairs)-old[1]
    assert square == bracket and wrong_mean > 0
    assert any(residuals.values())
    return square,cross,wrong_mean,len(support)


def integrate(poly, lo=0, hi=1):
    return sum(a*Q(hi**(i+1)-lo**(i+1), i+1) for i,a in enumerate(poly))


def beta_rows():
    n = 0
    for m in range(3,18):
        density = [0]+[(m-1)*(m-2)*(-1)**k*comb(m-3,k) for k in range(m-2)]
        for a in range(6):
            assert integrate([0]*a+density) == Q(factorial(a+1)*factorial(m-1),factorial(m+a-1))
            n += 1
        for x in (Q(0),Q(1,7),Q(2,3),Q(1)):
            soon = (1-x)**(m-2)*(1+(m-2)*x)
            later = (1-x)**(m-1)*(1+(m-1)*x)
            assert soon == integrate(density,x)
            assert soon-later == (m-1)*x*x*(1-x)**(m-2)
            n += 1
    return n


def mesh_rows():
    n = 0
    # Keep mass in log coordinates so all three threshold tests are rational.
    l,d,r = Q(-1),Q(2),Q(1,3)
    for h in (Q(1,2), Q(1,7)):
        for shift in (Q(0),h/2,h):
            for b in (l,l+shift,l+h,d,d+shift,d+h):
                for w in (r,r-shift,r-h,r+h):
                    low = l+h < b <= d and w > r
                    actual = l+shift < b <= d+shift and w > r-shift
                    high = l < b <= d+h and w > r-h
                    assert low <= actual <= high
                    n += 1
    return n


if __name__ == '__main__':
    for args in ((4,5,7,2),(4,5,8,3),(5,6,8,2)):
        result = fixture(*args)
        print('window, terminal, degree:', args, 'variance, birth-residual cross, wrong drift, support:', *result)
    print(f'{beta_rows()} Beta moment/tail controls; {mesh_rows()} strict translated-window controls')
