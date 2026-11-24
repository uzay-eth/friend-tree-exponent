#!/usr/bin/env python3
"""Gold1 B5 lifetime with a genealogically tagged, evolving sixth core."""
import argparse
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from replay import Top, ADJ, GRID, HI, LO, Q, change, layout, solve_row


def structure(receiver):
    parents = (None,0,0,0,3,receiver)
    layers = [[(0,)]]
    next_birth = {}
    for _ in range(5):
        steps = set()
        for p in layers[-1]:
            allowed = tuple(j for j in range(1,6) if j not in p and parents[j] in p
                            and not (receiver!=1 and j==2 and 1 not in p))
            next_birth[p] = allowed
            steps.update(tuple(sorted((*p,j))) for j in allowed)
        layers.append(sorted(steps))
    return parents,layers,next_birth


def evaluate(r,receiver,p):
    a,b,_ = r.row('b5',p[:5])
    return a+b[receiver]/p[5][0],b[0]


def compare(r,receiver,caps,tip_min,upper,old=None):
    parents,layers,extensions = structure(receiver)
    tables = {}
    for layer in reversed(layers):
        for p in layer:
            n = len(p)
            pos = {v:i for i,v in enumerate(p)}
            adj = layout.neighbors(n,[(pos[parents[v]],pos[v]) for v in p[1:]])
            c = tuple(len(a)+(i==0) for i,a in enumerate(adj))
            axes = [range(caps[v],max(2,k)-1,-1) for v,k in zip(p,c)]
            table = tables[p] = {}
            checks = positive = 0
            children = [(pos[parents[j]],j,tables[tuple(sorted((*p,j)))])
                        for j in extensions.get(p,())]
            for D in product(*axes):
                a,b = Q(0),Q(0)
                if n==6:
                    labels = tuple(zip(D,c))
                    if not r.support('b5',labels[:5]):
                        continue
                    if (max(D[1],D[2],D[4]) if old is None else D[old])>=tip_min:
                        a,b = evaluate(r,receiver,labels)
                num,den = Q(0),LO if upper else HI
                y0 = 0
                births = []
                for i,(d,k) in enumerate(zip(D,c)):
                    L = d-k
                    B = Q(L,d)
                    births.append(B)
                    hazard = L+sum(Q(1,D[j]) for j in adj[i])
                    y = table.get(change(D,i,d+1),0)
                    num += hazard*y
                    den += hazard+B
                    if i==0:
                        y0 = y
                # For the same birth, each genuinely different admissible
                # role may seed a different prefix. The old one always exits.
                for i,j,dst in children:
                    k = sum(v<j for v in p)
                    num += births[i]*dst.get(D[:k]+(2,)+D[k:],0)
                v,k = solve_row(a,(b,),((Q(0),Q(1,2)),),den,num,(y0,),upper)
                table[D] = v
                checks += k
                positive += bool(v)
            if n==6:
                print('terminal',len(table),'corners',checks,'nonzero',positive,flush=True)
    value = (HI if upper else LO)*Q(tables[(0,)][(2,)],GRID)
    print(r.side,'receiver',receiver,'caps',caps,'tip_min',tip_min,'ceiling' if upper else 'floor',
          'ticks/core',float(value),'exact',value,'potential/core',float(value/r.s.h['S']),flush=True)
    Top.value.cache_clear()
    return value


def fixture(r):
    counts = Counter()
    for receiver in range(5):
        parents,layers,ext = structure(receiver)
        # A fixed six-vertex template has several chronological orders. The
        # recursive prefix construction and an independent literal map agree.
        edges = list(layout.SHAPES['b5'])+[(receiver,5),(0,6)]
        adj = layout.neighbors(7,edges)
        S = frozenset(range(5))
        for d in (2,3,5):
            D = (7,38,4,5,4,d,17)
            labels = tuple(zip(D,map(len,adj)))
            pre = layout.census(r.s,adj,D)
            assert pre[S][0]=='b5'
            value = pre[S][1]
            plus = [layout.census(r.s,adj,change(D,x,D[x]+1)) for x in S]
            birth = [layout.census(r.s,layout.neighbors(8,edges+[(x,7)]),D+(2,))
                     if D[x]>len(adj[x]) else {} for x in S]
            z = -r.lam*value
            for x in S:
                gain = D[x]-len(adj[x])+sum(Q(1,D[y]) for y in adj[x])
                z += gain*(plus[x].get(S,(None,0))[1]-value)
                z += Q(D[x]-len(adj[x]),D[x])*(birth[x].get(S,(None,0))[1]-value)
            for w in S:
                if len(adj[w]&S)==1:
                    x, = adj[w]&S
                    T = S-{w}
                    z += (Q(1,D[w])-r.s.r1(labels[x]))*(plus[x].get(T,(None,0))[1]-pre[T][1])
            a,b = evaluate(r,receiver,labels[:6])
            assert a+b/17 == r.sig*z
            # An actual gain at the donor updates the receiver reciprocal;
            # no B5 descriptor changes on that event.
            v,_ = evaluate(r,receiver,change(labels[:6],5,(d+1,labels[5][1])))
            assert v-a == r.row('b5',labels[:5])[1][receiver]*(Q(1,d+1)-Q(1,d))
            assert Q(D[receiver]-len(adj[receiver]))+sum(Q(1,D[y]) for y in adj[receiver]) == (
                Q(D[receiver]-len(adj[receiver]))+
                sum(Q(1,D[y]) for y in ADJ['b5'][receiver])+Q(1,d)+Q(receiver==0,17))
            counts['moving reward'] += v!=a
            counts['exact row'] += 1
        assert layers[-1]==[tuple(range(6))]
        allowed = set()
        def walk(order):
            if len(order)==6:
                allowed.add(order[1:])
                return
            for j in ext[tuple(sorted(order))]:
                walk(order+(j,))
        walk((0,))
        expected = {order for order in permutations(range(1,6))
                    if all(p==0 or order.index(p)<order.index(j)
                           for j,p in enumerate(parents) if j>0)
                    and (receiver==1 or order.index(1)<order.index(2))}
        assert allowed==expected
        counts['chronological orders'] += len(allowed)
    assert counts['moving reward']
    print(r.side,'fixtures',dict(counts),flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--side',choices=('lower','upper'),required=True)
    ap.add_argument('--receiver',type=int,action='append',choices=(0,1,3,4))
    ap.add_argument('--tip-cap',type=int,default=8)
    ap.add_argument('--tip-min',type=int,default=2)
    ap.add_argument('--old-role',type=int,choices=(1,2,4))
    ap.add_argument('--old-cap',type=int,default=48)
    ap.add_argument('--internal-cap',type=int,default=8)
    ap.add_argument('--ceiling',action='store_true')
    ap.add_argument('--fixture',action='store_true')
    opt = ap.parse_args()
    r = Top(opt.side)
    if opt.fixture:
        fixture(r)
    else:
        receivers = opt.receiver or [0,1,3,4]
        assert len(set(receivers))==len(receivers)
        total = Q(0)
        for x in receivers:
            ic,tc = opt.internal_cap,opt.tip_cap
            caps = change((ic,tc,tc,ic,tc,5),x,min(5,ic if x in (0,3) else tc))
            if opt.old_role is not None:
                assert opt.old_role!=x
                caps = change(caps,opt.old_role,opt.old_cap)
            total += compare(r,x,caps,opt.tip_min,opt.ceiling,opt.old_role)
        print('union ticks/core',float(total),'exact',total,'potential/core',float(total/r.s.h['S']),flush=True)
    r.s.close()


if __name__=='__main__':
    main()
