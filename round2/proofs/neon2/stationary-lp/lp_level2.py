"""Chains and ordered sibling pairs in the physical occupation LP."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement, product

from lp_level1 import Model, acc, add, cli, v


class MotifModel(Model):
    def build(self):
        super().build()
        if not 2 <= self.D3 <= self.N:
            raise ValueError('need 2 <= motif cutoff <= degree cutoff')
        S = [s for s in self.states if 2 <= s[0] <= self.D3]
        C, V, omitted, credit = {}, {}, {}, defaultdict(dict)
        count = defaultdict(dict)
        for p,m,c in product(S, repeat=3):
            if p[1] < 2 or m[1] < 2:
                continue
            j = C[p,m,c] = self.var(('chain',p,m,c))
            acc(count,('ancestor',(m,c)),j,1)
            acc(count,('child',(p,m)),j,1)
            acc(credit,('ancestor',(m,c)),j,F(1,p[0]))
            acc(credit,('child',(p,m)),j,F(1,c[0]))
        for p in S:
            if p[1] < 3:
                continue
            for c,z in combinations_with_replacement(S,2):
                j = V[p,c,z] = self.var(('wedge',p,c,z),self.D3-2)
                for pivot,donor in ((c,z),(z,c)):
                    acc(count,('sibling',(p,pivot)),j,1)
                    acc(credit,('sibling',(p,pivot)),j,F(1,donor[0]))
        for (p,c), pop in self.em.items():
            if p not in S or c not in S:
                continue
            for mode,mult in (('ancestor',1),('child',c[1]-1),('sibling',p[1]-2)):
                if not mult:
                    continue
                key = mode,(p,c)
                j = omitted[key] = self.var(('omitted',key),mult)
                e = self.ext[key]
                self.eq(('motif-marginal',key),add((1,count[key]),(1,v(j)),(-mult,pop)))
                self.row(('motif-credit-min',key),add((1,credit[key]),(-1,v(e))))
                self.row(('motif-credit-max',key),add((1,v(e)),(-1,credit[key]),(-F(1,self.D3+1),v(j))))
        self.motif_count = len(C)+len(V)
        self.motif_station(C,((1,),(0,2),(1,)),'C')
        self.motif_station(V,((1,2),(0,),(0,)),'V')

    def motif_station(self, table, neighbors, tag):
        # Each unordered pair is counted once, even when its child states agree.
        low, high = defaultdict(dict), defaultdict(dict)
        def canon(t):
            return (t[0], *sorted(t[1:])) if tag=='V' else t
        def put(target, expr, scale=1, which=(low,high)):
            if target not in table:
                return
            for dst in which:
                for j,x in expr.items():
                    acc(dst,target,j,x*scale)
        for st,j in table.items():
            put(st,v(j,-1*sum(self.bval(s) for s in st)))
            for i,s in enumerate(st):
                base = s[0]-s[1]+sum((F(1,st[z][0]) for z in neighbors[i]),F(0))
                ext_cap = F(s[1]-len(neighbors[i]),2)
                ext = self.var(('motif-external',tag,st,i),ext_cap*max(1,self.D3-2)) if ext_cap else None
                if ext is not None:
                    self.row(('motif-external-cap',tag,st,i),{ext:1,j:-ext_cap})
                rate = add((1,v(j,base)),(1,v(ext) if ext is not None else {}))
                put(st,rate,-1)
                for dest,hazard in ((self.grow(s),rate),(self.promo(s),v(j,self.bval(s)))):
                    put(canon(st[:i]+(dest,)+st[i+1:]),hazard)
        birth = (2,1)
        for (p,c), e in self.em.items():
            if not p[0] or not c[0]:
                continue
            if tag=='C':
                put((p,self.promo(c),birth),e,self.bval(c))
            else:
                put(canon((self.promo(p),c,birth)),e,self.bval(p))
        for st,j in table.items():
            self.station((tag,st),v(j),low[st])


if __name__=='__main__':
    cli(MotifModel, default_motif=6)
