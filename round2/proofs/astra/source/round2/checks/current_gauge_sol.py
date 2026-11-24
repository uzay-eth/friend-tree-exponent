#!/usr/bin/env python3
"""Exact referee fixtures; run from the repository root with python3 -B."""
from fractions import Fraction as Q
from functools import cache
from itertools import product
from pathlib import Path
import sys


J = (1, ())
ZERO = Q(0)
HALF = Q(1, 2)


def q(D, cap=20):
    return Q(1, min(D, cap))


def degree(S):
    return 1 + S[0] + len(S[1])


def gain(S):
    return S[0] + 1, S[1]


@cache
def size(S):
    return 1 + S[0] + sum(map(size, S[1]))


def flat(S, parent=-1, out=None):
    if out is None:
        out = []
    i = len(out)
    out.append(parent)
    for Y in S[1]:
        flat(Y, i, out)
    out.extend([i] * S[0])
    return out


def views(par, slot):
    kids = [[] for _ in par]
    for i, p in enumerate(par[1:], 1):
        kids[p].append(i)
    D = [len(k) + (slot if i == 0 else 1) for i, k in enumerate(kids)]
    shapes = {}
    for i in reversed(range(len(par))):
        if kids[i]:
            shapes[i] = (sum(not kids[j] for j in kids[i]),
                         tuple(sorted(shapes[j] for j in kids[i] if kids[j])))
    return D, shapes


def arrows(par, slot):
    D, _ = views(par, slot)
    for i, j in enumerate(par[1:], 1):
        yield i, j, Q(1, D[i])
        yield j, i, Q(1, D[j])


@cache
def moves(S):
    par = flat(S)
    return tuple((r, views(par + [j], 1)[1][0]) for _, j, r in arrows(par, 1))


def emissions(S, cap=20, full=True):
    D = degree(S)
    a, a1, p = q(D, cap), q(D + 1, cap), Q(1, D)
    r = S[0] + sum(Q(1, degree(Y)) for Y in S[1])
    if full and S[0]:
        yield S[0] * p, J, a, a
    for Y in S[1]:
        if a > p:
            yield a - p, Y, a, a
        yield p, gain(Y), a, a
        if a > a1:
            yield r, Y, a1, a
        e = degree(Y)
        if q(e + 1, cap) < q(e, cap):
            for Z in Y[1]:
                yield p, Z, q(e + 1, cap), q(e, cap)


def fixture():
    # root D=21 planted / 20 centre; Y D=4; Z D=19; X has two leaves.
    # Fixed indices preserve identity even when canonical ordering changes.
    return [-1, 0, 1, 2] + [0] * 19 + [1] * 2 + [2] * 17 + [3] * 2


def event_test():
    def f(S, a):
        return (size(S) + 2 * S[0]) ** 2 + a * (5 + size(S) ** 2) + 31 * a * a

    par = fixture()
    failures = dict.fromkeys(('post-gain child', 'stale grandchild', 'receiver rate',
                              'omit exterior', 'uncapped mark', 'gained shape in drop'), 0)
    special = []
    for slot in (0, 1):
        D, old = views(par, slot)

        def total(deg, shapes):
            return sum(f(S, q(deg[par[v]])) for v, S in shapes.items() if v)

        for src, dst, rate in arrows(par, slot):
            D1, new = views(par + [dst], slot)
            lhs = rate * (total(D1, new) - total(D, old))

            def split(plant=None):
                ans = ZERO
                for v, S in old.items():
                    if not v:
                        continue
                    p = par[v]
                    a = q(D[p])
                    # Intrinsic shapes and the externally charged shape increment
                    # are both read at the pre-event exterior label.
                    if not (plant == 'omit exterior' and src == p and dst == v):
                        ans += rate * (f(new[v], a) - f(S, a))
                    if p == dst and dst in old:
                        Y = gain(S) if plant == 'gained shape in drop' else S
                        high = Q(1, D[p]) if plant == 'uncapped mark' else a
                        low = Q(1, D1[p]) if plant == 'uncapped mark' else q(D1[p])
                        rr = rate
                        if src == par[dst] and plant == 'receiver rate':
                            rr = Q(1, D[dst])
                        if not (src == par[dst] and plant == 'stale grandchild'):
                            ans += rr * (f(Y, low) - f(Y, high))
                if dst not in old:
                    ans += rate * f(new[dst], q(D[src]))
                    if plant == 'post-gain child':
                        # A promotion preserves the promoter degree. Instead
                        # incorrectly treat it as another leaf gain there.
                        ans += rate * (f(J, q(D[src]+1)) - f(J, q(D[src])))
                return ans

            assert split() == lhs, (slot, src, dst)
            for m in failures:
                failures[m] += split(m) != lhs
            if (src, dst) in ((0, 1), (1, 2), (2, 1)):
                special.append((slot, src, dst, str(lhs)))
    assert all(failures.values()), failures
    print('single-arrow exact fixtures:', special)
    print('single-arrow plants rejected:', failures)


@cache
def terminal(S):
    return Q((S[0] + 2) % 3, 1 + len(S[1]))


@cache
def mean(k, S, l, u, plant=None):
    if not k:
        return terminal(S) if l == u else ZERO
    def score(Y, a=l, b=u):
        return mean(k - 1, Y, a, b, plant)
    x = score(S)
    ans = sum(r * (score(Y) - x) for r, Y in moves(S))
    if l < u:
        return ans - (l if plant == 'dormant diagonal' else u)*x + (u-l)*score(S, l, l)
    for r,Y,a,b in emissions(S):
        if plant == 'grandchild shape' and a != b:
            Y = gain(Y)
        ans += r * score(Y, a, b)
    return ans - l*x


@cache
def physical(k, S):
    if not k:
        return sum(terminal(Y) for Y in views(flat(S), 1)[1].values())
    v = physical(k - 1, S)
    return sum(r * (physical(k - 1, Y) - v) for r, Y in moves(S))


def shell_test():
    # A smaller asymmetric fixture for derivatives; still has depth three.
    S = (0, ((2, ((0, (J,)),)), (3, ())))
    a = q(degree(S))
    plants = dict.fromkeys(('missing shell', 'zero all roots', 'dormant diagonal', 'grandchild shape'), 0)
    for k in range(4):
        U = lambda Y, b=ZERO, m=None: mean(k, Y, b, b, m)
        M = lambda Y: physical(k, Y)
        exact = M(S) - sum(map(M, S[1]))
        shell = lambda m=None: U(S, ZERO, m) + sum(U(Y,a,m)-U(Y,ZERO,m) for Y in S[1])
        assert exact == shell(), k
        D, ys = views(flat(S), 1)
        p = flat(S)
        assert M(S) == sum(U(Y, q(D[p[v]]) if v else ZERO) for v,Y in ys.items())
        plants['missing shell'] += M(S) != shell()
        plants['zero all roots'] += exact != U(S)
        for m in ('dormant diagonal', 'grandchild shape'):
            plants[m] += exact != shell(m)
        l, u = Q(1, 5), Q(1, 4)
        assert mean(k, S, l, u) == U(S,l)-U(S,u)
    assert all(plants.values()), plants
    print('top-zero transform, shell and dormant derivatives: orders 0..3; rejected', plants)


def contraction_test():
    d, alpha, sigma = Q(3,20), Q(9,25), Q(13,25)
    h = lambda c: (Q(1), Q(17,20), Q(39,50))[c-1] if c <= 3 else Q(7,10)
    for D in range(2,21):
        for c in range(1,D+1):
            assert (1-Q(c,D))*(h(c+1)-h(c)+q(D)/(q(D)+d)) <= alpha*h(c)
    assert q(20)/(q(20)+d) <= alpha*h(4)  # degree tail
    norm = lambda exponent: max((q(j)+d)/(exponent+q(j)-alpha) for j in range(2,21))
    assert norm(sigma) == Q(65,66)
    assert norm(alpha+d) == 1 and norm(alpha+d-Q(1,100)) > 1
    # With 100 identical children, an unweighted child-sum estimate diverges
    # from the proper ratio using the cap-tail path weight at the root.
    a = q(20)
    gain = a/(sigma+a-alpha)
    assert 100*gain > 1 > 100*gain/(h(101)+100*a/(a+d))
    print('cap20 intrinsic profile and contraction: theta(.52)=65/66; threshold mutation rejected')


def laplace_test():
    @cache
    def current(k, S, l, u):
        if not k:
            return Q(S[0],degree(S)) if l == u else ZERO
        f = lambda Y,a=l,b=u: current(k-1,Y,a,b)
        x = f(S)
        ans = sum(r*(f(Y)-x) for r,Y in moves(S)) - u*x
        if l != u:
            return ans + (u-l)*f(S,l,l)
        return ans + sum(r*f(Y,a,b) for r,Y,a,b in emissions(S,full=False))

    @cache
    def permanent(k, par, tags):
        D, shapes = views(par,1)
        if not k:
            return Q(shapes[0][0],D[0])
        f = lambda p,a: permanent(k-1,p,a)
        x = f(par,tags)
        ans = -tags[0]*x
        for src,dst,r in arrows(par,1):
            mark = q(D[src]) if dst not in shapes else tags[dst]
            ntags = tags[:dst]+(mark,)+tags[dst+1:]+(ZERO,)
            ans += r*(f(par+(dst,),ntags)-x)
        p = Q(1,D[0])
        for child in shapes:
            if par[child] != 0:
                continue
            ids = [child]
            for j in range(child+1,len(par)):
                if par[j] in ids:
                    ids.append(j)
            subtree = (-1,) + tuple(ids.index(par[j]) for j in ids[1:])
            a = tuple(tags[j] for j in ids)
            ans += (a[0]-p)*f(subtree,a) + p*f(subtree+(0,),a+(ZERO,))
        return ans

    diff = [current(k,J,HALF,HALF)-permanent(k,(-1,0),(HALF,ZERO)) for k in range(5)]
    assert diff == [0,0,0,0,Q(1,96)], diff
    print('independent newborn residual derivative differences:', ','.join(map(str,diff)))


def killed(N, s, plant=None):
    shapes, level = {J}, {J}
    for _ in range(N - 2):
        level = {Y for S in level for _,Y in moves(S)}
        shapes.update(level)
    idx = {S:i for i,S in enumerate(sorted(shapes))}

    def add(out, vec, r):
        for i,x in vec.items():
            out[i] = out.get(i,ZERO) + r*x

    @cache
    def response(S, l, u):
        # Map the unknown common active forcing vector to this row's response.
        ans = {idx[S]: Q(1)} if l == u else {}
        for r,Y in moves(S):
            if Y in idx:
                add(ans,response(Y,l,u),r)
        if l < u:
            add(ans,response(S,l,l),u-l)
        event_rate = size(S)-Q(1, degree(S))
        if plant == 'erase boundary loss' and size(S) == N:
            event_rate = ZERO
        den = s + event_rate + (l if plant == 'dormant diagonal' else u)
        return {i:x/den for i,x in ans.items()}

    m = len(idx)
    rows = []
    for S,i in idx.items():
        D = degree(S)
        j = min(D,20)
        row = {i:Q(1)}
        for r,Y,l,u in emissions(S,full=False):
            add(row,response(Y,l,u),-r)
        rows.append([row.get(k,ZERO) for k in range(m)] + [Q(S[0],D)*Q(j,j+4)])
    for i in range(m):
        j = next(j for j in range(i,m) if rows[j][i])
        rows[i], rows[j] = rows[j], rows[i]
        for j in range(i+1,m):
            if rows[j][i]:
                a = rows[j][i]/rows[i][i]
                for k in range(i+1,m+1):
                    rows[j][k] -= a*rows[i][k]
    x = [ZERO]*m
    for i in reversed(range(m)):
        x[i] = (rows[i][-1]-sum(rows[i][j]*x[j] for j in range(i+1,m)))/rows[i][i]
    return tuple(sum(v*x[i] for i,v in response(J,q(j),q(j)).items()) for j in range(2,21))


def finite_test(root):
    sys.path.insert(0, str(root / 'eval'))
    import gauge
    s = Q('.565728')
    N = 6
    sh, qi, G, C = gauge.matrices(N, s)
    vals, _ = gauge.approximate(G, C, 80)
    bounds, _, _ = gauge.exact_certificate(G, C, sh, qi, s, vals)
    exact = killed(N, s)
    assert all(lo <= x <= hi for x,(lo,hi) in zip(exact,bounds))
    for m in ('erase boundary loss','dormant diagonal'):
        assert all(not lo <= x <= hi for x,(lo,hi) in zip(killed(N,s,m),bounds)), m
    print('exact killed N=6: 19 enclosures agree, both plants miss all 19')


def weight_test(root):
    sys.path.insert(0, str(root / 'weight'))
    import cert
    cap, d, rho = 36, Q(8,25), Q(141,250)
    gamma = rho-d
    # Borrow the candidate edge ticks only; no use of cert.prove.
    ticks = cert.propose(cap,d,rho,True)
    P = lambda D: ticks[min(D,cap)]
    z = lambda D: q(D,cap)/(q(D,cap)+d)
    k = lambda D: (z(D)-z(D+1))/d
    assert all(1 >= P(D) >= max(z(D), P(D+1)) for D in range(2,cap+1))
    pairs = [(D, Q(1,D)) for D in range(2,cap+1)] + [(cap,ZERO)]
    ages = q(cap,cap), HALF
    first = lambda D,r,a: (a+d)*(P(D)-z(D)-r*k(D)) + r*(P(D)-P(D+1))
    second = lambda D,p,j,r,a,send: (first(D,r,a)*P(j)
                    + p*P(D)*(P(j)-P(j+1)) - (a+d)*send*k(j))
    for D,a in product(range(2,cap+1),ages):
        assert min(first(D,r,a) for r in (0,D-Q(3,2))) > 0
    corners = [(D,p,j,r,a) for (D,p),(j,pj),a in product(pairs,pairs,ages)
               for r in (pj,D-2+pj)]
    assert min(second(D,p,j,r,a,p) for D,p,j,r,a in corners) > 0
    assert min(second(D,p,j,r,a,HALF) for D,p,j,r,a in corners) < 0

    @cache
    def h(D,c):
        if D >= cap:
            return P(cap)/gamma
        if (D,c) == (2,1):
            return Q(1)
        L = D-c
        b = Q(L,D)
        promo = b*(h(D,c+1)+P(D)) if L else ZERO
        need = max(Q(1,50), *((r*h(D+1,c)+promo)/(gamma+r+b)
                              for r in (Q(L), L+Q(c,2)+HALF+d)))
        return Q(-(-need*10**10//1),10**10)

    for D in range(2,cap+1):
        for c in range(1,D+1):
            L = D-c
            v = h(D,c)
            assert Q(1578186331,5000000000) <= v <= 1
            promo = Q(L,D)*(h(D,c+1)-v+P(D)) if L else ZERO
            for r in (Q(L), L+Q(c,2)+HALF+d):
                assert r*(h(D+1,c)-v)+promo <= gamma*v

    @cache
    def V(S):
        D = degree(S)
        return h(D,1+len(S[1])) + P(D)*sum(map(V,S[1]))

    def W(S,l,u):
        return V(S) * (1/(l+d) if l == u else 1/(l+d)-1/(u+d))

    for n in (32,33,34,38):
        S = (n, ((n+1, ((0,(J,)),)), (0,(J,))))
        for a in ages:
            v = W(S,a,a)
            drift = sum(r*(W(Y,a,a)-v) for r,Y in moves(S))-a*v
            drift += sum(r*W(Y,l,u) for r,Y,l,u in emissions(S,cap,False))
            assert drift <= rho*v
        l,u = Q(1,36),Q(1,35)
        v = W(S,l,u)
        bare = sum(r*(W(Y,l,u)-v) for r,Y in moves(S))
        switch = -u*v+(u-l)*W(S,l,l)
        assert switch == d*v and bare+switch <= rho*v
        assert -l*v+(u-l)*W(S,l,l) > d*v
    print('cap36 exact corners/profile, four seam shapes, both row types; wrong send/diagonal rejected')


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1] / 'current-mark'
    event_test()
    shell_test()
    contraction_test()
    laplace_test()
    finite_test(root)
    weight_test(root)
