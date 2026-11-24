#!/usr/bin/env python3
"""Sparse Gold1 top rows and finite offspring-Palm lifetime bounds."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'higher'))
from screen import Rows, Q, LO, HI, J, bud, change, gain, layout
from price import GRID, solve_row


KEYS = ('p5', 'b5', 's4')
ADJ = {k: tuple(map(frozenset, layout.neighbors(5, layout.SHAPES[k]))) for k in KEYS}
# Tip, attachment, envelope, retained internal roles, retained endpoint roles.
LEGS = {
    'p5': ((0,1,'45p',(1,2,3),(4,)), (4,3,'45p',(3,2,1),(0,))),
    'b5': ((1,0,'45bp',(0,3),(2,4)), (2,0,'45bp',(0,3),(1,4)),
           (4,3,'45bs',(0,3),(1,2))),
    's4': tuple((i,0,'45s4',(0,),tuple(j for j in range(1,5) if j!=i)) for i in range(1,5)),
}
# Parent indices in chronological order, then permutation into the row roles.
MASKS = {'p5': ((0,1,2,3),(0,1,2,3,4)),
         'b5': ((0,0,0,3),(0,1,2,3,4)),
         's4': ((0,0,0,0),(0,1,2,3,4))}


def pick(p, roles):
    return tuple(p[i] for i in roles)


class Top(Rows):
    @lru_cache(None)
    def value(self, key, p):
        i, e = layout.ROLES[key]
        return self.s.value(key, pick(p,i), pick(p,e))

    def support(self, key, p):
        return self.s.support(key, pick(p,layout.ROLES[key][0]))

    def row(self, key, p, projected=False):
        """Signed exact intercept, gain slopes, and separate incoming arrows."""
        assert self.support(key,p)
        s, sig = self.s, self.sig
        val = self.value(key,p)
        a = -sig*self.lam*val
        slopes = []
        for i,(D,c) in enumerate(p):
            L = D-c
            assert L >= 0 and c >= len(ADJ[key][i])
            dz = sig*(self.value(key,change(p,i,gain(p[i])))-val)
            a += (L+sum(Q(1,p[j][0]) for j in ADJ[key][i]))*dz
            if L:
                a += sig*Q(L,D)*(self.value(key,change(p,i,bud(p[i])))-val)
            slopes.append(dz)
        terms = []
        for w,x,f,ii,ee in LEGS[key]:
            inner, ends = pick(p,ii), pick(p,ee)
            delta = layout.source_delta(s,f,inner,ends)
            z = sig*(Q(1,p[w][0])-s.r1(p[x]))*delta
            terms.append(z)
            a += self.project(f,inner,(p[w],),ends,delta,p[x]) if projected else z
        return a,tuple(slopes),tuple(terms)


def literal(r):
    s = r.s
    counts = Counter()
    # Test the exact index predicates with isolated boundary transitions.
    # P5's inherited SP4=21 pair test is redundant under SP5=18.
    seams = {'p5': (((9,2),(4,2),(2,2)),0),
             'b5': (((8,3),(5,2)),0),
             's4': (((8,4),),0)}
    for key,(inner,i) in seams.items():
        assert s.support(key,inner) and not s.support(key,change(inner,i,gain(inner[i])))
    assert s.support('p5',((6,2),)*3) and not s.support('p5',((7,2),(6,2),(6,2)))
    assert s.support('b5',((6,3),(8,2))) and not s.support('b5',((6,3),(9,2)))
    for degree,arity in ((10,1),(10,5),(10,6),(11,1),(12,1),(22,1)):
        co = s.G[5]
        dist = co.dist((degree,arity))
        assert sum(dist.values())==1 and sum(co.inv(i)*w for i,w in dist.items())==Q(1,degree)
        if degree>11:
            assert dist[co.inf]==1-Q(11,degree)
    # The spoke fixture has physically distinct equal-labelled arms. Other
    # fixtures cross internal degree/sum/arity and endpoint tail seams.
    fixtures = [([(0,i) for i in range(1,5)]+[(i,i+4) for i in range(1,5)],
                 (8,4,4,4,4,2,2,2,2))]
    for key, Ds in [('p5', ((10,9,5,3,12),(11,6,6,6,22))),
                    ('b5', ((8,10,12,6,11),(6,11,11,8,22))),
                    ('s4', ((8,10,11,12,22),))]:
        edges = list(layout.SHAPES[key])
        for D in Ds:
            fixtures.append((edges+[(0,5)],D+(2,)))
        # Raise a retained internal vertex to its last supported arity.
        # Star keeps four equal spokes to exercise four separate transfers.
        i = 2 if key == 'p5' else 0
        extra = 6-len(ADJ[key][i])
        fixtures.append((edges+[(i,5+j) for j in range(extra)],
                         change((5,)*5,i,7)+(2,)*extra))

    for edges,D in fixtures:
        n = len(D)
        adj = list(map(frozenset,layout.neighbors(n,edges)))
        labels = tuple(zip(D,map(len,adj)))
        pre = layout.census(s,adj,D)
        plus = [layout.census(s,adj,change(D,x,D[x]+1)) for x in range(n)]
        birth = [layout.census(s,layout.neighbors(n+1,edges+[(x,n)]),D+(2,))
                 if D[x]>len(adj[x]) else {} for x in range(n)]
        raw = sum(Q(d-c,d)*s.f(*J) for d,c in labels)
        donor_ledger = defaultdict(Q)
        for S,(key,v,order) in pre.items():
            z = -r.lam*v
            for x in S:
                L = D[x]-len(adj[x])
                delta = plus[x].get(S,(None,0))[1]-v
                z += L*delta
                for w in adj[x]:
                    full = Q(1,D[w])*delta
                    moved = Q(0)
                    T = S|{w}
                    if w not in S and T in pre and labels[x] in s.G[2].ex:
                        moved = (Q(1,D[w])-s.r1(labels[x]))*delta
                        donor_ledger[S] -= moved
                        donor_ledger[T] += moved
                        counts['nonzero arrows'] += bool(moved)
                    assert (full-moved)+moved==full
                    z += full
                if L:
                    nb = birth[x].get(S|{n},(None,0))[1]
                    z += Q(L,D[x])*(birth[x].get(S,(None,0))[1]-v+nb)
                    counts['nonzero newborn'] += bool(nb)
            raw += z
            if key not in KEYS:
                continue
            p = pick(labels,order)
            a,b,t = r.row(key,p)
            rho = tuple(sum(Q(1,D[w]) for w in adj[x]-S) for x in order)
            actual = a+sum(vv*rr for vv,rr in zip(b,rho))
            assert len(t)==layout.INCOMING[key]
            # The incoming ledger may be populated later in iteration order.
            legs = []
            for w in S:
                if len(adj[w]&S)==1:
                    x, = adj[w]&S
                    T = S-{w}
                    legs.append(r.sig*(Q(1,D[w])-s.r1(labels[x]))*
                                (plus[x].get(T,(None,0))[1]-pre[T][1]))
            assert sorted(t)==sorted(legs)
            assert actual==r.sig*z+sum(legs)
            assert r.row(key,p,True)[0]>=a
            counts[key] += 1
            counts[key+' degree exit'] += sum(S not in plus[x] for x in S)
            counts[key+' arity exit'] += sum(bool(birth[x]) and S not in birth[x] for x in S)
            counts[key+' repeated leg'] += len(set(t))<len(t) and any(t)
            counts[key+' negative'] += actual<0
        transfer = sum(donor_ledger.values())
        assert transfer == 0
        # Independent full potential jump, each newborn occurrence counted once.
        total = lambda rows: sum(v for _,v,_ in rows.values())
        v = total(pre)
        full = -r.lam*v
        for x,(d,c) in enumerate(labels):
            full += (d-c+sum(Q(1,D[w]) for w in adj[x]))*(total(plus[x])-v)
            if birth[x]:
                full += Q(d-c,d)*(total(birth[x])-v)
        assert raw+transfer==full
    for key in KEYS:
        assert counts[key+' degree exit'] and counts[key+' arity exit']
    assert counts['b5 repeated leg'] and counts['s4 repeated leg']
    print(r.side,'literal',dict(sorted(counts.items())),flush=True)


def palm(r,key,caps,high,upper):
    """Strict nonfounder crown; first child born at root degree >= high."""
    parents, order = MASKS[key]
    discount = LO if upper else HI
    nxt = {}
    counts = []
    for n in range(5,0,-1):
        adj = layout.neighbors(n,tuple(zip(parents[:n-1],range(1,n))))
        c = tuple(len(a)+(i==0) for i,a in enumerate(adj))
        axes = [range(cap,max(2,k,high if i==0 and n>1 else 0)-1,-1)
                for i,(cap,k) in enumerate(zip(caps,c))]
        table = {}
        checks = positive = 0
        for ds in product(*axes):
            p = tuple(zip(ds,c))
            if n==5:
                row = pick(p,order)
                if not r.support(key,row):
                    continue
                a,b,_ = r.row(key,row)
                b = tuple(b[order.index(i)] for i in range(n))
            else:
                a,b = Q(0),(0,)*n
            # The root's live exterior parent contributes to reward and the
            # same gain hazard. All other participating core neighbors are
            # internal; an extra child exits the strict mask at any prefix.
            rho = ((Q(0),Q(1,2)),)+((Q(0),Q(0)),)*(n-1)
            den, num = discount, Q(0)
            added = parents[n-1] if n<5 else -1
            succ = []
            for i,(d,ci) in enumerate(p):
                L = d-ci
                born = Q(L,d)
                g = L+sum(Q(1,ds[j]) for j in adj[i])
                den += g+born
                y = table.get(change(ds,i,d+1),0)
                succ.append(y)
                num += g*y
                if i==added and (n>1 or d>=high):
                    num += born*nxt.get(ds+(2,),0)
            v,k = solve_row(a,b,rho,den,num,succ,upper)
            table[ds] = v
            checks += k
            positive += bool(v)
        counts.append((len(table),checks,positive))
        nxt = table
    value = (HI if upper else LO)*Q(nxt[(2,)],GRID)
    print(r.side,key,'caps',caps,'first birth D>=',high,'upper' if upper else 'lower',
          'ticks/core',float(value),'exact',value,'potential/core',float(value/r.s.h['S']),
          'prefix (states,corners,nonzero)',counts,flush=True)
    Top.value.cache_clear()
    return value


def registration():
    # Physical births including a forbidden early sibling and identical spoke
    # descriptors. Compare the irreversible prefix automaton with subsets.
    P = (None,0,1,1,1,1,2,3,6,8,2,2,2,0,13,14,15,16,0,18,18,18,21)
    counts = {}
    for key,(pa,_) in MASKS.items():
        for end in range(5,len(P)+1):
            prefixes = [set() for _ in range(5)]
            for x in range(1,end):
                active = [set() for _ in range(5)]
                for i,ps in enumerate(prefixes):
                    for p in ps:
                        if P[x] not in p:
                            active[i].add(p)
                        elif i<4 and P[x]==p[pa[i]]:
                            active[i+1].add(p+(x,))
                active[0].add((x,))
                prefixes = active
            expected = set()
            for p in combinations(range(1,end),5):
                if tuple(P[x] for x in p[1:])!=tuple(p[i] for i in pa):
                    continue
                if any(P[x] in p for x in range(1,end) if x not in p):
                    continue
                expected.add(p)
            assert prefixes[-1]==expected
            counts[key] = counts.get(key,0)+len(expected)
    assert all(counts.values()),counts
    print('chronological registration snapshots',counts,flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--side',choices=('lower','upper'),action='append')
    p.add_argument('--palm',action='store_true')
    p.add_argument('--direction',choices=('floor','ceiling'),action='append')
    a = p.parse_args()
    if a.palm:
        registration()
    for side in a.side or ('lower','upper'):
        r = Top(side)
        if a.palm:
            config = {'p5': ((20,5,5,5,20),10), 'b5': ((8,10,10,5,14),7),
                      's4': ((8,9,9,9,9),7)}
            for direction in a.direction or ('floor','ceiling'):
                upper = direction=='ceiling'
                result = sum(palm(r,k,caps,high,upper) for k,(caps,high) in config.items())
                if upper:
                    assert result < Q(r.s.h['S'],10**6)
                print(side,'disjoint total',direction,
                      'ticks/core',float(result),'exact',result,'potential/core',float(result/r.s.h['S']),
                      flush=True)
        else:
            literal(r)
        r.s.close()


if __name__=='__main__':
    main()
