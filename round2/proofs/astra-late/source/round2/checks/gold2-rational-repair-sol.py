#!/usr/bin/env python3
"""Independent exact Gold2 repair: original host weights, physical Bellman, scan."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha1, sha256
from itertools import product
from math import ceil, comb, gcd
from pathlib import Path
import json
import runpy


R2 = Path(__file__).resolve().parents[1]
GOLD = R2/'occupation-quintet-gold2'
FLOORS = (207990, 395261, 52341, 251635)
TARGET = Q(99, 175)
WRITER = {'README.md': '43f1bc6dd9574096b60cfad82936d26f44a43d60',
          'check.py': '5f40e60547dd4e4d07e368126748e45d1e16eec1'}


def require(ok, why):
    if not ok:
        raise ValueError(why)


def blob(path):
    raw = path.read_bytes()
    return sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()


def host_upper(w, cap):
    """Read selected original weight blocks; do not call host.prepare/cells."""
    p, raw = w.p, w.raw
    C = {(d,c): int(w.f(d,c)) for d in range(2,cap+1) for c in range(d+1)}
    require(p.h['d2'] > cap, 'selected hosts cross degree seam')
    max_w = {}

    def credit(f, state, x):
        if state in C:
            max_w[f,state] = max(max_w.get((f,state),0), x)

    def peak(start, count):
        if len(raw[start:start+count]) != count:
            raise ValueError('truncated selected slab')
        return int(raw[start:start+count].max())

    E = p.ends[2]
    host_ids = {}
    for d,c in C:
        if c:
            host_ids.setdefault(E.ids[d,min(c,p.h['K2']+1)],[]).append((d,c))
    # Triangular G2 rows: each reading gives an entire prefix of partner roles.
    # The child endpoint may occur as either index and zero allows a founder.
    start = p.offset['g']
    for i in range(E.n):
        row = raw[start:start+i+1]
        if i in host_ids:
            for host in host_ids[i]:
                credit('g',host,int(row.max()))
        for j, hosts in host_ids.items():
            if j < i:
                for host in hosts:
                    credit('g',host,int(row[j]))
        start += i+1
    require(start == p.offset['h3'], 'G2 extent')

    for f,t,end,k in (('h3','3m',3,2), ('s3','S3',4,3), ('s4','S4',5,4)):
        width = comb(p.ends[end].n+k-1,k)
        for state,i in p.ins[t].items():
            if state in C:
                credit(f,state,peak(p.offset[f]+i*width,width))

    labels = {t:list(i) for t,i in p.ins.items()}
    for f in ('p4','p5','b5'):
        for idx,start in p.block[f].items():
            if f == 'p4':
                a,b = [labels['4i'][i] for i in idx]
                roles = [(a,b[1]-1), (b,a[1]-1)]
                n = p.ends[4].n
                size = comb(n+1,2) if a == b else n*n
            elif f == 'p5':
                a,b,c = [labels['5i'][i] for i in idx]
                roles = [(b,(a[1]-1)*(c[1]-1))]
                n = p.ends[5].n
                size = comb(n+1,2) if a == c else n*n
            else:
                i,j = idx
                a,b = labels['5o'][i], labels['5c'][j]
                roles = [(a,b[1]-1)]
                n = p.ends[5].n
                size = comb(n+1,2)*n
            if any(host in C for host,_ in roles):
                mx = peak(start,size)
                for host,k in roles:
                    credit(f,host,k*mx)

    for d,c in C:
        multiplicities = dict(g=1,h3=comb(c,2),p4=comb(c,2),s3=comb(c,3),
                              p5=comb(c,2),b5=c*comb(c-1,2) if c else 0,s4=comb(c,4))
        C[d,c] += sum(k*max_w.get((f,(d,c)),0) for f,k in multiplicities.items())
    require(all(v <= w.global_high for v in C.values()), 'pointwise global bound')
    return C


def bellman(U, cut, ceiling, hi):
    v = {}
    rows = 0
    for d in range(cut,1,-1):
        for c in range(d,0,-1):
            L = d-c
            birth = Q(L,d)
            gain, promote = v.get((d+1,c),0), v.get((d,c+1),0)
            B = ceiling-U[d,c]
            rates = (Q(L), L+Q(c,2))
            q = min((B+r*gain+birth*promote)/(hi+r+birth) for r in rates)
            v[d,c] = t = q.numerator//q.denominator
            for r in rates:
                require(0 <= t and hi*t <= B+r*(gain-t)+birth*(promote-t), 'Bellman row')
                rows += 1
    return v[2,1], rows


def fractions(lo, hi, n):
    return sorted(Q(p,q) for q in range(1,n+1)
                  for p in range(-((-lo.numerator*q)//lo.denominator), hi.numerator*q//hi.denominator+1)
                  if gcd(p,q)==1)


def physical():
    # Enumerate all increasing ordinary parent trees through size seven and all
    # literal arrows. No appeal to the boxed rates to compute the observed drift.
    count = arrows = 0
    for n in range(3,8):
        for pa in product(*(range(i) for i in range(1,n))):
            g = [[] for _ in range(n)]
            for i,p in enumerate(pa,1):
                g[p].append(i)
                g[i].append(p)
            degree = list(map(len,g))
            arity = [sum(degree[y]>1 for y in adj) for adj in g]
            for x,d in enumerate(degree):
                if d < 2:
                    continue
                c = arity[x]
                if c == 0:  # solitary core is not a nonfounder newborn
                    continue
                gain = bud = Q(0)
                for adj in g:
                    for y in adj:
                        arrows += 1
                        dd = d+(y==x)
                        cc = sum(degree[z]+(y==z)>1 for z in g[x])
                        rate = Q(1,len(adj))
                        require((dd-d,cc-c) in ((0,0),(1,0),(0,1)), 'joint gain and promotion')
                        if dd>d:
                            gain += rate
                        if cc>c:
                            require(degree[y]==1 and y in g[x], 'physical newborn target')
                            bud += rate
                require(bud==Q(d-c,d) and d-c<=gain<=d-c+Q(c,2), 'literal projected rates')
                count += 1
    return count,arrows


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument('--writer', nargs=2, action='append', metavar=('FILE','GIT_BLOB'), default=[])
    args = p.parse_args()
    pins = dict(args.writer) if args.writer else WRITER
    for name,pin in pins.items():
        path = GOLD/'rationals/repaired'/name
        require(path.parent == GOLD/'rationals/repaired' and blob(path)==pin, f'writer {name}')
    print(f'PASS {len(pins)} writer release pins', flush=True)
    path = GOLD/'host/reader.py'
    require(blob(path)=='796d43b52d9f472513c77823e862b26fc580bcc4', 'original host reader')
    reader = runpy.run_path(str(path))
    conditioned = (GOLD/'host/conditioned24.json').read_bytes()
    require(sha256(conditioned).hexdigest()=='3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790',
            'prior conditioned receipt')
    check = json.loads(conditioned)['results']['lower']
    exact = (Q(check['f_floor']), Q(39526146756920663734311,10**17),
             Q(1046834122639903380009,2*10**16),
             Q(16363762149897743918379249111,65029541606195200000000))
    require(tuple(int(v) for v in exact)==FLOORS, 'four rounded floors')
    count, arrows = physical()
    print(f'PASS {count} literal core drift fixtures, {arrows} arrows', flush=True)
    with reader['open_host']('lower') as w:
        U = host_upper(w,24)
        ceiling = w.global_high
        lo, hi = w.prior
        E, lam = w.golden['err_per_core'], w.lam
        require((ceiling,E,lam)==(1808603193693,19217,Q(7071421,12500000)), 'Gold2 lower data')
        cells = {tuple(map(int,k.split(','))): int(v[1]) for k,v in check['cells'].items()}
        mismatch = [(d,c,U[d,c],cells[d,c-1]) for d,c in U if c and U[d,c]!=cells[d,c-1]]
        require(not mismatch, f'independent host cells vs prior frozen cells: {mismatch[:20]}')
    print(f'PASS 299 matched original U cells and 23 solitary-case singleton cells through degree 24', flush=True)
    net = sum(FLOORS)-E
    require(net > 0, 'sign for improved denominator')
    ends = []
    for cut, needed in ((6,4700),(18,25700),(24,27680)):
        v, rows = bellman(U,cut,ceiling,hi)
        bound = ceil(ceiling-lo*v)
        margin = net-(TARGET-lam)*bound
        require(0 < bound < 1459000000000 and margin > needed, 'repaired margin')
        lower = lam+Q(net,bound)
        require(TARGET < lower <= hi < Q(142,251) and 142*175-99*251==1, 'Farey interval')
        require(fractions(lower,hi,425)==[] and fractions(lower,hi,426)==[Q(241,426)], 'denominator scan')
        ends.append(lower)
        print(f'PASS D<={cut}: {rows} Bellman corners, V(2,1)={v}, M<={bound}, margin={margin}',flush=True)
        print(f'  repaired [{lower}, {hi}]; no rational through 425; only 241/426 through 426',flush=True)
    require(ends == sorted(ends), 'nested cutoffs')
    require(ends[-1] == Q(10043824822832044427,17754235186087500000), 'pinned repair endpoint')


if __name__ == '__main__':
    main()
