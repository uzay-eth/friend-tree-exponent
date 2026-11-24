#!/usr/bin/env python3
"""Independent Fraction reference on every corner of a small, nonzero, asymmetric fixture."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement as multisets, permutations, product
from pathlib import Path
import random
import struct
import subprocess
import tempfile

from oracle import ARITY, CUTS, FAMILIES, HEADER, INF, NB, Certificate, bud, gain, inv, newborn_extensions, replace

PARENTS = {'V':(-1,0,1),'C':(-1,0,0),'P':(-1,0,1,2),'F':(-1,0,0,2),'Y':(-1,0,1,1),'T':(-1,0,0,0)}


def make_table(template, file, mode):
    c = Certificate(template)
    data = bytearray(Path(template).read_bytes())
    h = list(HEADER.unpack_from(data))
    h[1] = mode
    h[10] = c.Sr//4 if mode == 2 else 0
    HEADER.pack_into(data,0,*h)
    rng = random.Random(819)
    for name,(off,n) in c.sections.items():
        for i in range(n):
            if name in FAMILIES: v = rng.randrange(-100,101)
            elif name in ('rch','rpar','xch','xpar','xcc'): v = rng.randrange(c.Sr//2+1)
            elif name.startswith('e'):
                # +/-200 contains all fine differences. Deliberately asymmetric
                # means and radii test the symmetrization taken at Y/T sources.
                v = rng.randrange(311,601) if int(name[-1])%2 else rng.randrange(-90,91)
            else: continue
            if name == 'E':
                g = c.levels[2]
                if i//g.n not in g.parent or i//g.n == i%g.n == g.inf: v = 0
            struct.pack_into('<q',data,HEADER.size+(off+i)*8,v)
    file.write_bytes(data)


class Reference:
    def __init__(self, c, f):
        self.c, self.f = c,f
        self.g = c.levels[ARITY[f]]
        g = self.g
        self.xyz = list(g.internal)+[INF]+[(g.d0-k-1,k) for k in range(g.kt+2)]
        self.parents = PARENTS[f]
        self.ch = [[i for i,p in enumerate(self.parents) if p == j] for j in range(len(self.parents))]

    @lru_cache(maxsize=None)
    def value(self, labels):
        st = tuple(self.xyz[i] for i in labels)
        a = self.c.active(self.f,st)
        return self.c.sign*self.c.tables[self.f][self.c.offset(self.f,a,labels)] if a >= 0 else 0

    def source(self, row, st, active):
        c, f = self.c,self.f
        def entry(name,i):
            start,n = c.sections[name]
            assert 0 <= i < n
            return c.data[start+i]
        def raw(labels):
            val = Q(0)
            for k,j,_,_,cen in CUTS[f]:
                if f == 'V':
                    group, pair, index = 'e23', (0 if k == 0 else 2), active*self.g.n+labels[2 if k==0 else 0]
                elif f == 'C':
                    group, pair, index = 'e23', 0, c.act['V'][st[0]]*self.g.n+labels[3-k]
                elif f in ('P','F'):
                    group = 'ep' if f == 'P' else 'ef'
                    pair = 2 if k == 3 else 0
                    index = active*self.g.n+labels[(0 if f == 'P' else 1) if k==3 else 3]
                elif f == 'Y' and k:
                    group,pair = 'eyv',0
                    index = (active*len(self.g.parent)+self.g.parent[labels[0]])*self.g.n+labels[5-k]
                else:
                    group,pair = 'eyc',0
                    center = st[1] if f == 'Y' else st[0]
                    tips = (2,3) if f == 'Y' else tuple(i for i in (1,2,3) if i != k)
                    i,l = (labels[u] for u in tips)
                    index = (c.act['Y'][center]*self.g.n+i)*self.g.n+l
                a = inv(self.xyz[labels[k]])-c.center(cen,st[j])
                val += a*c.sign*entry(group+str(pair),index)+abs(a)*entry(group+str(pair+1),index)
            return val
        if f not in ('Y','T'): return raw(row)
        split = 2 if f == 'Y' else 1
        orders = tuple(permutations(row[split:]))
        return sum(raw(row[:split]+o) for o in orders)/len(orders)

    def evaluate(self, row):
        c, g, f = self.c, self.g, self.f
        st = tuple(self.xyz[i] for i in row)
        active = c.active(f,st)
        w = self.value(row)
        total = -c.lam*w+self.source(row,st,active)
        for i,lab in enumerate(row):
            if lab == g.inf: continue
            nb = self.ch[i]+([self.parents[i]] if i else [])
            if lab > g.inf:
                d = self.value(replace(row,i,g.inf))-w
                kap = lab-g.inf-1
                assert len(nb) == 1 and row[nb[0]] < g.inf
                gl = max(0,(g.d0-kap-1+inv(st[nb[0]]))/(g.d0+1)) if kap <= g.kt else 0
                total += max(d,gl*d)
                if kap <= g.kt:
                    b = self.value(replace(row,i,lab+1))-w
                    total += max(b,Q(max(0,g.d0-kap-1),g.d0)*b)
                continue
            s = st[i]
            next_id = g.mu(gain(s))[0][0]
            def jump(x): return self.value(replace(x,i,next_id))-self.value(x)
            d = jump(row)
            total += s[0]*d
            for y in nb:
                if row[y] < g.inf: total += inv(st[y])*d
                elif row[y] > g.inf: total += Q(max(d,jump(replace(row,y,g.inf))),g.d0)
            if i == 0:
                centered = (f == 'V' and c.active('P',(NB,)+st)>=0) or (f == 'C' and c.active('Y',(NB,)+st)>=0)
                total += c.center(0,s)*d if centered else Q(max(0,d),2)
            ext = s[1]-len(self.ch[i])
            if ext:
                centered = False
                if f in ('V','C'):
                    ff,ns = newborn_extensions(f,st,i)
                    centered = c.active(ff,ns)>=0
                total += ext*(c.center(1,s)*d if centered else Q(max(0,d),2))
            if s[0]:
                post = replace(st,i,bud(s))
                nbval = 0
                if f in ('V','C'):
                    ff,ns = newborn_extensions(f,post,i)
                    nbval = c.sign*c.weight(ff,ns)
                total += Q(s[0],s[0]+s[1]+1)*(self.value(replace(row,i,g.id[post[i]]))-w+nbval)
        return total

    def rows(self):
        g, f = self.g,self.f
        for center in self.c.act[f]:
            if f in ('C','T'):
                for ends in multisets(range(g.n),ARITY[f]-1):
                    yield (g.id[center],)+ends
            elif f in ('V','Y'):
                for a in g.parent:
                    for ends in multisets(range(g.n),ARITY[f]-2):
                        yield (a,g.id[center])+ends
            else:
                for a,b in product(g.parent if f == 'P' else range(g.n),range(g.n)):
                    x,y = map(g.id.get,center)
                    yield (a,x,y,b) if f == 'P' else (x,a,y,b)


def run(path,checker):
    c = Certificate(path)
    for f in FAMILIES[1:]:
        ref = Reference(c,f)
        cases = tuple(ref.rows())
        assert len(cases) == len(c.tables[f])
        inp = ''.join(f+' '+' '.join(map(str,x))+'\n' for x in cases)
        out = subprocess.run([str(checker),str(path),'--query'],input=inp,
                             text=True,capture_output=True,check=True)
        vals = out.stdout.splitlines()
        assert len(vals) == len(cases)
        peak = None
        for x,v in zip(cases,vals):
            want = ref.evaluate(x)
            if want != Q(v): raise AssertionError((f,x,want,v))
            if peak is None or want > peak: peak = want
        top = -(-peak.numerator//peak.denominator)
        print(f"mode={c.sign} family={f} independent_exact_rows={len(vals)} max_ceil={top}",flush=True)


if __name__ == '__main__':
    out = Path('round2/outputs')
    with tempfile.TemporaryDirectory(dir=Path(__file__).parent,prefix='corner-fixture-') as tmp:
        path = Path(tmp)/'witness.bin'
        for mode in (1,2):
            make_table(out/'mc-fixture.bin',path,mode)
            run(path,out/'mh-check')
    print('PASS',flush=True)
