#!/usr/bin/env python3
"""Gold2 upper: live 3999->4000 seam and disjoint edge reserves."""
from fractions import Fraction as Q
from hashlib import file_digest, sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
import json
from math import ceil, floor
import os
from pathlib import Path
from subprocess import check_output
from types import SimpleNamespace


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
M, E, C, GRID, X = 4000, 32, 44, 10**12, 10**9
LO, HI = 565713486, 565736987


def load(name, path):
    spec = spec_from_file_location(name, ROOT / path)
    m = module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def bind(path, digest):
    raw = path.read_bytes()
    assert sha256(raw).hexdigest() == digest, path
    return json.loads(raw)


def extremum(n, gain, other, scale, upper):
    disc = LO if upper else HI
    fixed = sum(r for r, _ in other)
    fixed_n = n + sum(r*v for r, v in other)
    points = [(X*(fixed_n + sum(r*v for r, v in zip(rr, vs))), X*(fixed + sum(rr)) + disc*scale)
              for rr in product(*(rates for rates, _ in gain))
              for vs in [tuple(v for _, v in gain)]]
    if upper:
        out = max(-(-a//b) for a, b in points)
        assert all(out*b >= a for a, b in points)
    else:
        out = min(a//b for a, b in points)
        assert all(out*b <= a for a, b in points)
    return out


def host(cut, upper, source=None, arity_at_hit=False):
    """Either a physical first-hit transform or paid actual child registrations."""
    assert not arity_at_hit or (upper and source is None)
    zero = GRID if upper and source is None and not arity_at_hit else 0
    nxt = [GRID*(c if arity_at_hit else int(source is None)) for c in range(C+1)] + [zero]
    for d in range(cut-1, 1, -1):
        cur = [zero]*(C+2)
        for c in range(min(C,d), 0,-1):
            b = 2*(d-c)
            pay = 0 if source is None else source[d]
            cur[c] = extremum(0, [((2*d*(d-c), d*(2*d-c)), nxt[c])],
                              [(b, cur[c+1]+pay)], 2*d, upper)
        nxt = cur
    return Q(nxt[1], GRID)


def geom(k):
    return (1+Q(LO,X))**(-k)


def face(p):
    e = p.ends[2]
    tails = [e.ids[e.d,c] for c in range(1,e.k+2)]
    word = lambda i, j: p.word(p.offset['g']+max(i,j)*(max(i,j)+1)//2+min(i,j))
    assert all(word(i,j)==0 for i,j in product([e.inf,*tails],repeat=2))
    vals = {}
    signs = [Q(0), Q(0)]
    for j in range(len(e.labels)):
        w = [Q((M-e.d)*word(e.inf,j)+e.d*word(i,j), M) for i in tails]
        # The seam arity may be arbitrary. Do not assume the whole face is negative.
        x,y = min(w), max(w)
        signs = [min(signs[0],x), max(signs[1],y)]
        vals[e.labels[j]] = floor(max(0,-y)), ceil(max(0,-x))
    return vals, signs


def child(vals, upper):
    """Physical marked child, with exact reciprocal in each live direction."""
    idx = int(upper)
    nxt = {(d,c): vals[d,min(c,23)][idx]*GRID for d in range(2,E+1) for c in range(1,d+1)}
    carry = max(v[idx] for (d,_),v in vals.items() if d>E)*GRID if upper else 0
    born = [0]*(M+1)
    for d in range(M-1, 1, -1):
        tail = -(-carry*d*X//(d*X+LO)) if upper else 0
        cur = {}
        for q in range(E,1,-1):
            for c in range(q,0,-1):
                scale = 2*d*q
                high = ((0 if upper else d-min(d,C+1))*scale+2*d,
                        (d-2)*scale+2*d+d*q)
                low = ((q-c)*scale+2*q, (q-c)*scale+2*q+(c-1)*d*q)
                g = [(high,nxt[q,c]), (low,cur.get((q+1,c),tail))]
                b = [((q-c)*2*d,cur.get((q,c+1),0))]
                if not upper:
                    b.append((scale,0))  # worst-case kill on any later host promotion
                cur[q,c] = extremum(0,g,b,scale,upper)
        born[d] = cur[2,1]
        nxt, carry = cur, tail
    return born


def fixtures(p, wire, cap):
    ordinary = wire.old

    def selected(graph):
        adj, marks = ordinary.physical_core(graph)
        return {(u,v): p.g(marks[u],marks[v]) for u, ns in adj.items() for v in ns if u<v
                and max(marks[u][0],marks[v][0])>=M}

    slots = exterior = 0
    for d,q in ((M-2,2),(M-1,3),(M-1,251),(M-1,252),(M,2),(M,4000)):
        tree = ordinary.physical_tree([(0,1),(0,2),(0,3),(1,4)], [d,q,5,252,7])
        core, marks = ordinary.physical_core(tree)
        old = selected(tree)
        for u, nbrs in enumerate(tree):
            for v in nbrs:
                slots += 1
                new = list(tree)+[{v}]
                new[v] = tree[v]|{len(tree)}
                post = selected(new)
                old_jump = sum(post[k]-a for k,a in old.items())
                birth = v not in core
                if birth and marks[u][0]>=M:
                    b, i = post[tuple(sorted((u,v)))], 0
                elif not birth and marks[v][0]==M-1:
                    b, i = 0, sum(post[tuple(sorted((v,w)))] for w in core[v]
                                  if marks[w][0]<M)
                else:
                    b=i=0
                assert sum(post.values())-sum(old.values()) == old_jump+b+i
                if not birth and any(v in edge and u not in edge for edge in old):
                    exterior += bool(old_jump)
        assert all(abs(a)<=cap for a in old.values())
    assert exterior
    print('literal arrow ownership:', slots, 'nonzero exterior changes:', exterior, flush=True)


def edge_reserves(writer, p, doc):
    # Same generic exact edge recurrence, isolated in this process with upper-only inputs.
    with writer.Snapshot('upper') as snap:
        e = writer.Rows(snap,p,doc)
        price = type(e).__mro__[1].price
        price.__globals__.update(LO=Q(LO,X), HI=Q(HI,X))
        d18 = price(e,18)
        small = writer.gold
        large = load('gold2_seam_payer','occupation-quintet-gold2/payer/price.py')
        fy = []
        for module, opts in (
                (small, dict(host=251, stock=17, receiver=5, rstock=5, donor=12, tstock=5, mask='both',
                             exclude=18, absorb=True, bins=[6,8,11,15,18,24,32,48,64,90,120,150])),
                (large, dict(host=150, stock=10, receiver=34, rstock=7, donor=5, tstock=5,
                             mask='small', bins=list(range(6,35))))):
            if module is small:
                module.LO, module.HI = Q(LO,X), Q(HI,X)
                m = module.Price(e,SimpleNamespace(**opts))
            else:
                module.gold.LO, module.gold.HI = Q(LO,X), Q(HI,X)
                m = module.Price(e,SimpleNamespace(**opts),False)
            module.ext.LO, module.ext.HI = Q(LO,X), Q(HI,X)
            v = snap.scale*sum(m.run())
            print('upper FY',opts['mask'],v,float(v),'states',m.path_count,flush=True)
            assert v>0
            fy.append(v)
        snap.check_unchanged()
        return d18,*fy


def global_test(doc, reserves):
    profile = bind(HERE.parent/'host/h3-joint/receipt.json',
                   'ae6f5cde2ef710b4196c3c445a3786580e5ccc9dfd9206b2962f60628866c6ed')['results']['upper']
    old = Q(profile['final_endpoint'])
    uni = Q(profile['uniform_endpoint'])
    denom = Q(profile['d24_denominator'][1])
    assert Q(profile['edge_floor'])==0 and profile['payload_sha256']==doc['verify']['sha256']
    target = Q(241,426)
    improved = uni-sum(reserves)/denom
    assert target < improved < old
    display = Q(ceil(improved*10**15),10**15)
    extra = Q(2087,3689)
    assert target < display < extra < old
    # Old ceiling was below 142/251. These Farey neighbours show that
    # 2087/3689 is the only fraction newly excluded through order 3689.
    below, above = Q(1945,3438), Q(142,251)
    assert below < display and old < above
    assert above.numerator*below.denominator-below.numerator*above.denominator==1
    assert extra==Q(below.numerator+above.numerator,below.denominator+above.denominator)
    print('standalone edge rebates:',*map(float,reserves),'ticks/core')
    print('new upper (upward 15 decimals):',f'{float(display):.15f}',
          'gain vs best prior:',float(old-improved),'gap to 241/426:',float(improved-target))
    print('more ticks needed from new upper on same denominator:',float((improved-target)*denom))
    print('new rational exclusion:',extra,'through order',extra.denominator,'with 241/426 still inside')
    return display


def main():
    wire = load('gold2_seam_wire','attacks/neon2-quintet/gold2/physical/check.py')
    writer = load('gold2_seam_writer','occupation-quintet-gold2/post-fy/replay.py')
    source = wire.ROOT/'source'
    assert check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()==wire.PIN
    docs = bind(source/wire.CERT/'certificate_gold2.json',
                '5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284')
    assert Q(LO,X) <= Q(docs['lower']['certified_exact']) <= Q(docs['upper']['certified_exact']) <= Q(HI,X)
    doc = docs['upper']
    path = wire.ROOT/'quintet_upper.bin'
    p = wire.Gold(path)
    try:
        with os.fdopen(os.dup(p.fd),'rb') as f:
            assert file_digest(f,'sha256').hexdigest()==doc['verify']['sha256']==wire.DIGEST['upper']
        assert (p.h['mode'],p.h['DF'],p.ends[2].d,p.ends[2].k)==(1,M,252,22)
        vals, signs = face(p)
        cap = max(-doc['verify']['comp_fin']['g_far_min'],doc['verify']['comp_fin']['g_far_max'])
        assert max(map(abs,signs))<=cap
        fixtures(p,wire,cap)
        h0,h1 = (host(M,up) for up in (False,True))
        # P(c(A)>k) <= geom(k), for A independent Exp(mu), by rate-one birth domination.
        overflow = (C+1+Q(X,LO))*geom(C)
        # All children born after the registration sweep exits are charged here.
        clo = host(M,False,child(vals,False))
        chi = host(M,True,child(vals,True))+cap*overflow
        # First moment above seam; used for a global-cap comparator and opposite sign.
        c_mass = host(M,True,arity_at_hit=True)+overflow
        lo = Q(LO,X)*clo - Q(HI,X)*signs[1]*c_mass
        up = Q(HI,X)*(chi - signs[0]*h1)
        base = Q(HI,X)*cap*c_mass
        assert 0<lo<up<base
        print('seam face ticks min/max:',*map(float,signs))
        print('high mass:',float(h0),float(h1),'crossing arity upper:',float(c_mass))
        print('negative seam flux ticks/core:',float(lo),float(up),'vs cap upper',float(base),flush=True)
        assert lo>=155105227 and up<=3288652241 and base>=7124998973
        # Price a new edge joining two G2-tail vertices. Only pay histories
        # with the mother already high at the birth, then track the child.
        j = p.ends[2].d
        a,b = host(j,False), host(j,True)
        c_high = 10*b + geom(10)/(1-geom(1))
        high_bud = max(Q(0),a-c_high/j)
        double_tail = doc['cap_g']*high_bud*a/Q(HI,X)
        assert double_tail>0
        global_test(doc,(*edge_reserves(writer,p,doc),double_tail))
        stat = os.fstat(p.fd)
        assert (stat.st_dev, stat.st_ino, stat.st_size, stat.st_ctime_ns,stat.st_mtime_ns)==(
            p.stat.st_dev,p.stat.st_ino,p.stat.st_size,p.stat.st_ctime_ns,p.stat.st_mtime_ns)
        print('PASS exact upper-only seam corners and edge composition; sparse words:',len(p.reads))
    finally:
        os.close(p.fd)


if __name__=='__main__':
    main()
