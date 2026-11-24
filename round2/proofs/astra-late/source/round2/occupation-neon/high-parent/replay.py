#!/usr/bin/env python3
"""Exact physical occupation of virgin arms on finite Neon D18 hub bands."""
from argparse import ArgumentParser
from functools import cache
from itertools import product
from math import comb, prod
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'higher'))
from extract import Extract, Cert, NB, Q, LO, HI, bud, up, label

GRID = 10**24
ANCHORS = (18, 32, 64, 128)
DEN = (Q('1.356833021311'), Q('1.367169995210'))
PAID = (Q('4.37936743978386489e-6'), Q('3.99314053622308967e-6'))


class Physical(Extract):
    @cache
    def value(self, family, x):
        return sum(self.weight(family, tuple(y for y,_ in row))*prod(w for _,w in row)
                   for row in product(*(self.c.embed(s) for s in x)))

    def center(self, x, parent):
        c = self.c
        if sum(x) < c.D:
            return Q(int((c.rpar if parent else c.rch)[c.state[x]]), c.Sr)
        return Q(int((c.trpar if parent else c.trch)[x[1]]), c.Sr)

    def transfer(self, x):
        p,w = x
        c = self.c
        t = Q(c.D, sum(p))
        slope = self.weight('g', (('inf',),w))-self.weight('g', (('T1',p[1],1),w))
        return 3*t*t*(1-t)*Q(slope,c.D)+t**3*int(c.yp[p[1],c.child[w]])


def check_star(e, p, edge, pair, dg, tr, dh):
    """Independently sum whole star drift, then undo external row transfers."""
    k, d = p[1], sum(p)
    x = (p,)+(NB,)*k
    children = [list(range(1,k+1))]+[[] for _ in range(k)]
    for ext in (Q(0),Q(1,2)):
        def phi(z, kids):
            return e.c.potential(kids,z)

        v = phi(x,children)
        drift = -e.c.lam*v
        for i in range(k+1):
            dd = d if i == 0 else 2
            gain = p[0]+Q(k,2)+ext if i == 0 else 1+Q(1,d)
            drift += gain*(phi(label(x,i,up(x[i])),children)-v)
            if x[i][0]:
                kids = [a[:] for a in children]+[[]]
                kids[i].append(k+1)
                drift += Q(x[i][0],dd)*(phi(label(x,i,bud(x[i]))+(NB,),kids)-v)
        first = Q(0)
        for i in range(k+1):
            s = x[i]
            a,kc,_ = s
            dd = sum(s)
            old = e.weight('f',(s,))
            inc = e.weight('f',(up(s),))-old
            rho = Q(k,2)+ext if i == 0 else Q(1,d)
            psi = Q(int(e.c.f2[k]),d*d)+Q(int(e.c.f3[k]),d**3) if i == 0 else inc
            center = e.center(s,True)+kc*e.center(s,False)
            first += (a+rho)*inc+(center-rho)*psi-e.c.lam*old
            if a:
                bb = bud(s)
                first += Q(a,dd)*(e.weight('f',(bb,))+e.weight('f',(NB,))+e.value('g',(bb,NB))-old)
        # There is no actual exterior edge/motif in phi. Add back the
        # corrections sent to that edge and its hV instances.
        omitted = (ext-e.center(p,True))*(psi_at(e,p)+k*tr)
        ep = edge+(ext+Q(k-1,2))*(dg-tr)
        hp = pair+(ext+Q(k-2,2))*dh
        assert drift == first+ep*k+hp*comb(k,2)+omitted, (p, ext)


def psi_at(e,p):
    d,k = sum(p),p[1]
    return Q(int(e.c.f2[k]),d*d)+Q(int(e.c.f3[k]),d**3)


def charges(e, d, k):
    """Signed floors for one virgin edge and one unordered virgin pair."""
    c = e.c
    p = (d-k-1, k, 1)
    l = p[0]
    b = Q(l, d)
    rp, rc = e.center(p, True), e.center(p, False)
    rw = e.center(NB, True)
    # f -> g. The child is inner; the parent sends only the tail psi.
    psi = psi_at(e,p)
    df = e.weight('f', (up(NB),)) - e.weight('f', (NB,))
    x = (p, NB)
    v = e.value('g', x)
    dg = e.value('g', (up(p), NB))-v
    # Remove only the prescribed transfer, retaining the physical remainder
    # multiplying the unknown sum of all exterior reciprocals at the hub.
    tr = e.transfer(x)
    edge = ((l+Q(1, 2))*dg + (rp+(k-1)*rc)*tr
            + (1+Q(1, d))*(e.value('g', (p, up(NB)))-v) - c.lam*v
            + Q(1, 2)*(e.value('g', (p, bud(NB))) + e.value('hV', (p, bud(NB), NB))-v)
            + (Q(1, 2)-rc)*psi + (Q(1, d)-rw)*df)
    if l:
        edge += b*(e.value('g', (bud(p), NB))+e.value('hC', (bud(p), NB, NB))-v)
    # j is not in the row table; using the full exterior box is safe even
    # when some other child reciprocals are known to be 1/2.
    edge_floor = e.sign*edge + min(0, e.sign*(dg-tr)*Q(k, 2))

    pair = Q(0)
    if k >= 4:
        # No Level-4 transfer/newborn at a tail hC apex. Both tips are virgin.
        x = (p, NB, NB)
        v = e.value('hC', x)
        dh = e.value('hC', (up(p), NB, NB))-v
        a = ((l+1)*dh - c.lam*v + (1+Q(1, d))*2*(e.value('hC', (p, up(NB), NB))-v)
             + e.value('hC', (p, bud(NB), NB))-v + 2*(Q(1, 2)-rc)*tr)
        if l:
            a += b*(e.value('hC', (bud(p), NB, NB))-v)
        pair = e.sign*a + min(0, e.sign*dh*Q(k-1, 2))
        if d in ANCHORS and k in (4,7):
            check_star(e,p,edge,a,dg,tr,dh)
    assert edge_floor >= 0 and pair >= 0, (d, k, edge_floor, pair)
    return edge_floor/c.S, pair/c.S


def solve(first, last, cap, reward):
    """Backward integer-grid subsolution at HI, returned in stationary units."""
    nxt = {}
    for d in range(last, 1, -1):
        here = {}
        for k in range(min(d-1, cap), -1, -1):
            l = d-k-1
            b = Q(l, d)
            for j in range(k+1):
                z = k, j
                r = reward(d, k, j)*GRID if d >= first else 0
                loss = j*(Q(3, 2)+Q(1, d))
                base = r + b*here.get((k+1, j+1), 0) + loss*here.get((k,j-1), 0)
                den = HI+b+loss
                v = nxt.get(z, 0)
                low, high = l+Q(j, 2), l+Q(k+1, 2)
                # This fractional-linear row is monotone in gain rate;
                # determine the minimizing end before rounding down.
                g = low if v*den >= base else high
                out = (base+g*v)//(den+g)
                assert 0 <= out and all((den+t)*out <= base+t*v for t in (low, high))
                here[z] = out
        nxt = here
    return LO*Q(nxt[0, 0], GRID)


def fmt(x):
    # Print a lower bound even when the nearest float would round upward.
    ticks = x*10**16//1
    return f'{ticks//10**16}.{ticks%10**16:016d}'


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('lower', type=Path)
    ap.add_argument('upper', type=Path)
    ap.add_argument('--cap', type=int, default=10)
    ap.add_argument('--last', type=int, default=191)
    args = ap.parse_args()
    assert 4 <= args.cap < 17 and args.last >= max(ANCHORS)
    for d in ANCHORS:
        floors = [solve(d, d, args.cap, r) for r in (
            lambda d,k,j: 1, lambda d,k,j: j, lambda d,k,j: comb(j, 2)*(k >= 4))]
        assert 0 < floors[0] and 0 < floors[1] and 0 < floors[2]
        print('exact D', d, 'hub mass, virgin edges, virgin pairs (k>=4):', *map(fmt, floors), flush=True)
    for i,path in enumerate((args.lower, args.upper)):
        c = Cert(path)
        assert (c.mode, c.D, c.K, c.S, c.lam) == ((2,1)[i],18,17,10**11,Q(('.565673','.565775')[i]))
        e = Physical(c)
        rows = {(d,k): charges(e,d,k) for d in range(18,args.last+1) for k in range(1,args.cap+1)}
        edge = lambda d,k,j: rows[d,k][0]*j if j else 0
        pair = lambda d,k,j: rows[d,k][1]*comb(j,2) if j >= 2 else 0
        re = solve(18,args.last,args.cap,edge)
        rp = solve(18,args.last,args.cap,pair)
        rr = solve(18,args.last,args.cap,lambda d,k,j: edge(d,k,j)+pair(d,k,j))
        assert re > 0 and rp > 0
        name = ('lower','upper')[i]
        for d in ANCHORS:
            band = solve(d,d,args.cap,lambda d,k,j: edge(d,k,j)+pair(d,k,j))
            assert 0 < band <= rr
            print(name, 'exact D',d,'joint reward',fmt(band),flush=True)
        add = max(rr, re+rp)
        end = c.lam+e.sign*(PAID[i]+add)/DEN[i]
        ticks = e.sign*((e.sign*end*10**12)//1)
        print(name, 'bands 18..'+str(args.last), 'rows',len(rows),
              'edge',fmt(re),'pair',fmt(rp),'joint',fmt(rr),'endpoint move',fmt(add/DEN[i]),
              'outward endpoint',f'{ticks//10**12}.{ticks%10**12:012d}',flush=True)


if __name__ == '__main__':
    main()
