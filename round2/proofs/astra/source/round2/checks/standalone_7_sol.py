"""Fresh exact finite controls for standalone Claim 7; no writer imports."""

from collections import defaultdict
from fractions import Fraction as F
from functools import cache


@cache
def state(path):
    adj = [[] for _ in range(len(path)+2)]
    for w, v in enumerate((0,)+path, 1):
        adj[v].append(w)
        adj[w].append(v)
    d = tuple(map(len, adj))
    p = [F(0) for _ in d]
    for u, row in enumerate(adj):
        for v in row:
            p[v] += F(1, len(d)*d[u])
    return d, tuple(p), sum(x >= 2 for x in d)


@cache
def descendants(path, n):
    if len(path)+2 == n:
        return ((path, F(1)),)
    return tuple((end, p*q) for w,p in enumerate(state(path)[1]) if p
                 for end,q in descendants(path+(w,), n))


def check():
    horizon = 7
    levels = {n: descendants((), n) for n in range(2,horizon+1)}
    promo_checks = drift_checks = stop_checks = 0
    old_late = below_factor = None
    for n, rows in levels.items():
        assert sum(p for _,p in rows) == 1
        if n < 3:
            continue
        for path,_ in rows:
            d,p,c = state(path)
            q = tuple(F(x,n) if x>=2 else F(0) for x in d)
            pk = sum(w for x,w in zip(d,p) if x == 1)
            assert sum(q) == F(n+c-2,n)
            assert 0 < pk <= F(c,n)
            l1 = sum(abs(a-b) for a,b in zip(p,q))
            assert l1 == F(c-2,n)+2*pk <= F(3*c-2,n)
            degree_law = defaultdict(F)
            for x,w in zip(d,p):
                degree_law[F(x,n)] += w
            assert degree_law[F(1,n)] == pk
            assert len([x for x in d if x == 1])*F(1,n) == F(n-c,n)
            v = path[-1]
            if state(path[:-1])[0][v] == 1:
                assert d[v] == 2
                for end_n in range(n,horizon+1):
                    avg = sum(p*F(state(end)[0][v],end_n)
                              for end,p in descendants(path,end_n))
                    assert avg <= F(2,n)
                    promo_checks += 1
                if n==5 and v==3:
                    below_factor = F(d[v],n)
        for m in range(3,n):
            tail = budget = F(0)
            for path,p in rows:
                d = state(path)[0]
                tail += p*sum((F(d[v],n) for v in range(m,n) if d[v]>=2),F(0))
                for r in range(m+1,n+1):
                    prev = path[:r-3]
                    v = path[r-3]
                    if state(prev)[0][v] == 1:
                        budget += p*F(2,r)
                        if v<m:
                            old_late = (m,r,v+1)
            assert tail <= budget

    # Include threshold stops, an early path-dependent stop, and a stop
    # triggered by the target selected at the very step under audit.
    rules = (
        lambda n,path,d,c: c*c>n,                  # K=1, beta=1/2
        lambda n,path,d,c: 4*c*c>9*n,              # K=3/2, beta=1/2
        lambda n,path,d,c: n>=3 and path[0]==0,   # before the cutoff on half the paths
        lambda n,path,d,c: n>=5 and path[-1]!=1,  # adaptively selects the old-term gain
        lambda n,path,d,c: False,
    )
    hits = set()
    for m in range(3,horizon):
        for rule in rules:
            # Carry (frozen size if already stopped, tracked labels, X, paid).
            def walk(path, mass, stop, tracked, x, paid):
                nonlocal drift_checks, stop_checks
                r = len(path)+2
                if r == horizon:
                    return mass*x, mass*paid
                d,p,c = state(path)
                if not stop and r>=3 and rule(r,path,d,c):
                    stop = r
                sx = sp = avg = F(0)
                for w,pw in enumerate(p):
                    nxt = path+(w,)
                    nd,_,nc = state(nxt)
                    size = r+1
                    live = not stop and size>m
                    add = live and d[w]==1
                    ntx = tracked+((w,) if add else ())
                    injection = F(2,size) if add else F(0)
                    y = sum((F(nd[v],size) for v in ntx),F(0)) if live else x
                    assert injection == (F(2*(nc-c),size) if live else 0)
                    if live and rule(size,nxt,nd,nc):
                        hits.add(size)
                    if live:
                        assert y >= sum(F(nd[v],size) for v in range(m,size) if nd[v]>=2)
                    if rule is rules[0]:
                        b = nc if not stop else state(path[:stop-2])[2]
                        assert (b-1)**2 <= size
                    avg += pw*(y-injection)
                    u,v = walk(nxt,mass*pw,stop,ntx,y,paid+injection)
                    sx += u
                    sp += v
                assert avg<=x, (path,m,stop,tracked)
                drift_checks += 1
                return sx,sp

            x,budget = walk((), F(1), 0, (), F(0), F(0))
            assert x<=budget
            stop_checks += 1

    assert hits and below_factor==F(2,5) and old_late
    # Prefix (1,2) makes path 3--1--2--4. At size 4, label 2
    # was promoted after m=3. The crossing branch can appear at r=5.
    path = (0,1)
    d,p,c = state(path)
    assert d[1]==2 and p[1]==F(3,8)
    selected_rule = rules[3]
    assert not selected_rule(4,path,d,c)
    survived = strict_next = survival_mass = F(0)
    for w,pw in enumerate(p):
        nxt = path+(w,)
        nd,_,nc = state(nxt)
        if selected_rule(5,nxt,nd,nc):
            strict_next += pw*F(d[1],4)
        else:
            strict_next += pw*F(nd[1],5)
            survived += pw*F(nd[1],5)
            survival_mass += pw
    assert strict_next-F(2,4)==F(3,80)>0
    # If the crossing r=5 promotes label 3, the tracked numerator is 4.
    d_hit,_,c_hit = state(path+(2,))
    assert c_hit==3 and c_hit**2>5 and (c_hit-1)**2<=5
    old_promotions = (1,2)
    assert all(v+1<=3 for v in old_promotions)
    assert not any(x>=2 for x in d_hit[3:])
    numerator = sum(d_hit[v] for v in old_promotions)
    correct_frozen, moving_denominator = F(numerator,5), F(numerator,7)
    assert correct_frozen>moving_denominator
    # Conditioning on a realized next promotion is not pre-step drift:
    # the recipient already known to be promoted has degree 2, not 1.
    assert F(2,5)>F(1,4) and below_factor>F(1,5)
    # Conditioning on survival through a stop that rejects other recipients
    # gives the tracked vertex the entire next gain.
    assert survived/survival_mass==F(3,5)>F(2,4)
    # Double star: the C/n budget cannot replace the label L1 budget.
    d,p,c = state((0,0,1,1))
    assert c==2 and len(d)==6
    leaf_recipient = sum(p[v] for v,x in enumerate(d) if x==1)
    q = tuple(F(x,6) if x>=2 else F(0) for x in d)
    assert leaf_recipient==F(2,9) != F(2,3)  # uniform vertex would have mass 2/3
    assert sum(abs(a-b) for a,b in zip(p,q))==F(4,9)>F(c,6)
    print(f'PASS: {promo_checks} conditional promotion bounds, {drift_checks} prefix drifts, '
          f'{stop_checks} horizon budgets; live crossing sizes {sorted(hits)}')
    print('PASS: exact label/degree/TV identities and mean tail at sizes 3..7')
    print('MUTANTS: anticipative strict stop drift +3/80; moving denominator 4/7 != 4/5; '
          'entry 1/5 < 2/5; postselected survival gain 3/5 > 2/4; late promotion may be old-label')
    print('DOUBLE STAR: P(leaf recipient)=2/9, uniform leaf mass=2/3, label L1=4/9 > C/n=1/3')


if __name__ == '__main__':
    check()
