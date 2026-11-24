#!/usr/bin/env python3
"""Paired exact physical witnesses; no independent charge box per sign."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'higher'))
from extract import Cert, Extract, NB, SPEC, embedding, physical, shape, up, bud, label


def tree(leaves, *children):
    return (0, leaves, tuple(sorted(children)))


J = tree(1)


class Pair:
    def __init__(self, lower, upper):
        self.c = lower, upper
        assert tuple(c.mode for c in self.c) == (2, 1)
        assert tuple(c.lam for c in self.c) == (Q('.565673'), Q('.565775'))
        assert all(c.stored == (3495091968-96)//8 and c.S == 10**11 for c in self.c)
        self.e = tuple(Extract(c) for c in self.c)
        self.weight = lru_cache(maxsize=40000)(self.weight)

    def weight(self, kind, states):
        # Contract both witnesses on *one* Bernstein embedding. Nonlinear
        # generator arithmetic comes after evaluating the exact physical rates.
        if kind == 'f':
            return tuple(Q(e.weight(kind, states), e.c.S) for e in self.e)
        if kind == 'hV' and states[1][1]>self.c[0].K or kind == 'hC' and states[0][1]>self.c[0].K:
            return Q(0), Q(0)
        if kind in ('P4','Yt','Yb','S3') and not self.e[0].supported(kind,states):
            return Q(0), Q(0)
        w = [Q(0), Q(0)]
        for row in product(*(self.c[0].embed(s) for s in states)):
            xs, ws = zip(*row)
            p = Q(1)
            for v in ws:
                p *= v
            for i, e in enumerate(self.e):
                w[i] += p * Q(e.weight(kind, xs), e.c.S)
        return tuple(w)

    def local(self, p, x):
        """Exact all-eight *physical*, uncharged drift, pooled per receiver."""
        ne = [[] for _ in x]
        children = [[] for _ in x]
        for i, parent in enumerate(p[1:], 1):
            ne[parent].append(i)
            ne[i].append(parent)
            children[parent].append(i)
        assert all(sum(s)>=2 and s[1]==len(children[i]) and 0<=s[2]<=1
                   for i,s in enumerate(x))
        value = [Q(0), Q(0)]
        constant = [Q(0), Q(0)]
        slopes = [[Q(0), Q(0)] for _ in x]

        def add(kind, xs, rate, out):
            for e, w in enumerate(self.weight(kind, xs)):
                out[e] += rate*w

        for kind, ids in embedding(p):
            xs = tuple(x[i] for i in ids)
            add(kind, xs, 1, value)
            for j, i in enumerate(ids):
                for sign, new in ((1, label(xs, j, up(x[i]))), (-1, xs)):
                    add(kind, new, sign, slopes[i])
                if x[i][0]:
                    beta = Q(x[i][0], sum(x[i]))
                    add(kind, label(xs, j, bud(x[i])), beta, constant)
                    add(kind, xs, -beta, constant)
        for i, (l, k, h) in enumerate(x):
            if l:
                post = label(x, i, bud(x[i])) + (NB,)
                new = len(x)
                for kind, ids in embedding(p+(i,)):
                    if new in ids:
                        add(kind, tuple(post[j] for j in ids), Q(l, l+k+h), constant)
            gain = l+sum(Q(1, sum(x[j])) for j in ne[i])
            for e in range(2):
                constant[e] += gain*slopes[i][e]

        drift = tuple(constant[e] - c.lam*value[e] for e, c in enumerate(self.c))
        return (drift[0], -drift[1]), (slopes[0][0], -slopes[0][1]), tuple(value)

    def full(self, z, verify=False):
        p, x = shape(z, 0)
        d, _, v = self.local(p, x)
        assert min(d) >= 0 and min(v) > 0
        if verify:
            kids = [[i for i, par in enumerate(p) if par == j] for j in range(len(p))]
            for i,c in enumerate(self.c):
                phi = lambda xs, kk: c.potential(kk, xs)/c.S
                old = phi(x, kids)
                direct = -c.lam*old
                for j,(l,k,h) in enumerate(x):
                    donors = kids[j] + ([p[j]] if j else [])
                    gain = l+sum(Q(1,sum(x[t])) for t in donors)
                    direct += gain*(phi(label(x,j,up(x[j])),kids)-old)
                    if l:
                        after = [row[:] for row in kids]+[[]]
                        after[j].append(len(x))
                        direct += Q(l,sum(x[j]))*(phi(label(x,j,bud(x[j]))+(NB,),after)-old)
                assert (old, direct*(-1)**i)==(v[i],d[i])
        return d,v

    def charged(self, z, h=1, apex=True):
        """All eight charged inner rows with exactly one shared parent port."""
        p,x=shape(z,h)
        assert max(map(sum,x))<=16
        rows=[]
        for kind, ids in embedding(p):
            if apex and ids[0]!=0:
                continue
            states=tuple(x[i] for i in ids)
            both=[]
            for e in self.e:
                a,free=e.row(kind,states)
                s=Q(0)
                for i,upward,n,delta in free:
                    j=ids[i]
                    if upward and p[j]<0:
                        assert n==1
                        s+=delta
                        continue
                    donors=([p[j]] if upward else
                            [u for u in range(len(p)) if p[u]==j and u not in ids])
                    assert len(donors)==n
                    a+=delta*sum(Q(1,sum(x[u])) for u in donors)
                both.append((a/e.c.S,s/e.c.S))
            rows.append((kind,ids,both))
        return rows


def paired_stop(pair, cap, ratios):
    # Reuse the independently generated physical arrow/crown transitions.
    sys.path.insert(0,str(HERE.parent/'check'))
    from crowns import moves, n, gain, degree

    def from_crown(z):
        l,kids=z
        return tree(l,*(from_crown(y) for y in kids))

    @lru_cache(None)
    def rows(z):
        return pair.charged(from_crown(z))

    lower,upper=Q('.565674722414'),Q('.565773666188')
    grid=10**19
    def solve(ratio, parent_cap=0):
        @lru_cache(None)
        def run(z, parent):
            if n(z)>cap:
                return Q(0),Q(0)
            joint=[]
            for _,_,((dm,sm),(dp,sp)) in rows(z):
                joint.append((dm-ratio*dp,sm-ratio*sp))
            v=sum(x for x,_ in joint); slope=sum(s for _,s in joint)
            pts=moves(z)
            r0=sum(r for _,r in pts)
            edge=run(gain(z),parent)
            if 2<=parent<=parent_cap:
                influx=(Q(1,parent),)
                # The focal sends to the parent. The remaining parent
                # gain hazard can only increase that lower bound. No
                # independence from crown events is imposed.
                exterior=(Q(1,degree(z)),Q(parent))
                lifted=run(z,parent+1)
            else:
                influx=(Q(0),Q(1,max(2,parent)))
                exterior=(Q(0),)
                lifted=(Q(0),Q(0))
            cand=[[],[]]
            for i in (0,1):
                jump=sum(r*run(y,parent)[i] for y,r in pts)
                for discount,p,t in product((lower,upper),influx,exterior):
                    cand[i].append((v+jump+p*(slope+edge[i])+t*lifted[i])/(discount+r0+p+t))
            lo=min(cand[0]); hi=max(cand[1]); assert lo<=hi
            return Q(lo*grid//1,grid),Q(-(-hi*grid//1),grid)

        if not parent_cap:
            return run((1,()),0)
        vs=[run((1,()),d) for d in range(2,parent_cap+2)]
        return (min(v[0] for v in vs),max(v[1] for v in vs)), vs

    for ratio in ratios:
        result=solve(ratio)
        tight, bins=solve(ratio,6)
        assert result[0]<=tight[0]<=tight[1]<=result[1]
        if cap==8 and ratio==Q(6,5):
            assert result[0]<0<result[1] and tight[1]<Q('-6.72e-7')
        if cap==8 and ratio==Q(1,4):
            assert tight[0]>Q('2.52e-6')
        print('paired killed lifetime',cap,'ratio',ratio,'free and tracked bounds',
              *(show(q) for q in result+tight),'birth degree 2',*(show(q) for q in bins[0]),flush=True)
    # All 8 are present across apex rows of the all-inner fixture.
    fixture=tree(1,tree(1,tree(1,J,J),tree(1,J),J),tree(1,J),J)
    rr=pair.charged(fixture,0,False)
    assert {f for f,_,_ in rr}==set(SPEC)
    for i in (0,1):
        full=pair.full(fixture,True)[0][i]
        assert sum(both[i][0] for _,_,both in rr)==full
        assert not any(both[i][1] for _,_,both in rr)
    # Same occurrence counts and participant degrees, with the actual
    # parent send moving from 1/2 to 1/8. Both have whole-tree realizations.
    rho=Q(6,5)
    ps=[Q(1,8),Q(1,2)]
    # Use a deliberately nontrivial common-port fixture, with all charges.
    planted=tree(1,tree(1,J,J),tree(1,J),J)
    rr=pair.charged(planted)
    terms=[(dm-rho*dp,sm-rho*sp) for _,_,((dm,sm),(dp,sp)) in rr]
    each=[min(d+p*s for p in ps) for d,s in terms]
    pooled=min(sum(d+p*s for d,s in terms) for p in ps)
    assert pooled>=sum(each)
    for d in (2,8):
        pair.full(tree(d-1,planted),True)
    assert pooled-sum(each)==Q('0.00006098845275')
    print('shared actual-parent port advantage',pooled-sum(each),show(pooled-sum(each)),flush=True)


def train(m,k):
    # Every interior spine vertex has one spine child, optionally one hub.
    # The k children of a hub have one leaf each; the hub retains no leaves.
    z=J
    hubs=(tree(0,*([J]*k)),) if k else ()
    for _ in range(m):
        z=tree(int(bool(k)),z,*hubs)
    return tree(1,z)


def obstruction(pair):
    a,b=(c.lam for c in pair.c)
    l,r=Q('.565674722414'),Q('.565773666188')
    lo,hi=(l-a)/(b-l),(r-a)/(b-r)
    for k, compare in ((0,lambda t:t<lo),(16,lambda t:t>hi)):
        vals=[]
        for m in (5,6,7):
            dd,vv=pair.full(train(m,k),m==6)
            vals.append(dd+vv)
        d=tuple(y-x for x,y in zip(vals[0],vals[1]))
        assert d==tuple(y-x for x,y in zip(vals[1],vals[2]))
        dm,dp,vm,vp=d
        dd,pp,mm,vv=vals[1]
        assert all(y>0 for y in vals[1]+d)
        assert compare(dd*vv/(pp*mm)) and compare(dm*vp/(dp*vm))
        assert compare(dd/pp) and compare(dm/dp)
        print('reachable train',k,'six-cell exact (slack-,slack+,phi-,phi+)',*(str(q) for q in vals[1]),flush=True)
        print('reachable train',k,'exact per inserted cell',*(str(q) for q in d),flush=True)
        print('reachable train',k,'normalized finite odds and bulk odds',show(dd*vv/(pp*mm)),show(dm*vp/(dp*vm)),flush=True)


def odds(pair, z, verify=False):
    (dm,dp),(vm,vp) = pair.full(z, verify)
    return dm/dp, vp/vm, dm*vp/(dp*vm)


def show(q):
    return f'{float(q):.12g}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('lower', type=Path)
    ap.add_argument('upper', type=Path)
    ap.add_argument('--stop-cap',type=int,default=8)
    args = ap.parse_args()
    pair = Pair(Cert(args.lower), Cert(args.upper))
    candidates = []
    for n in (2,16,17,18,14000,14001):
        candidates.append((f'star L={n}', tree(n)))
    for d,k in ((4,3),(17,3),(18,3),(19,3)):
        for child in (J, tree(13999), tree(1000000)):
            candidates.append((f'crown D={d} k={k} leafstock={child[1]}', tree(d-k,*([child]*k))))
    branch = tree(1, tree(1,J,J),tree(1,J),J)
    mixed = tree(1,branch,tree(16,J),tree(0,tree(1,J)))
    all8 = {kind for kind, ids in embedding(shape(mixed,0)[0])}
    assert all8 == set(SPEC)
    candidates.append(('all8 mixed near tail',mixed))
    bounds = [[],[],[]]
    for desc,z in candidates:
        result = odds(pair,z,True)
        for j,r in enumerate(result):
            bounds[j].append((r,desc))
    print('direct all8 physical drift/paired matches', len(candidates))
    for name, rows in zip(('slack ratio','potential ratio','local normalized odds'),bounds):
        for side,ext in (('min',min),('max',max)):
            q, desc = ext(rows)
            print(name,side,show(q),desc,'exact',q)
    a,b = (c.lam for c in pair.c)
    l,r = Q('.565674722414'),Q('.565773666188')
    lo, hi = (l-a)/(b-l), (r-a)/(b-r)
    print('incoming odds',show(lo),show(hi))
    obstruction(pair)
    paired_stop(pair,args.stop_cap,(Q(1,4),Q(1),Q(6,5),Q(5,4),Q(3,2),Q(2)))


if __name__ == '__main__':
    main()
