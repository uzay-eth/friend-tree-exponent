"""A slab-local strengthening of the unchanged reciprocal ChainModel."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb, prod
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from chains import ChainModel
from model import Cell


def basis(n,x):
    return tuple(comb(n,j)*x**j*(1-x)**(n-j) for j in range(n+1))


class BandModel(ChainModel):
    """Two slabs: D>=18 (including equality), and H<=D<=17."""
    def __init__(self, threshold, arity, order, *, local=True, **options):
        if not 3<=threshold<=16:
            raise ValueError('threshold must leave a distinct D<=16 remainder')
        if not all(options.get(k, True) for k in ('edge_balances','balances')):
            raise ValueError('bands require the original root, edge and chain balances')
        super().__init__(threshold,arity,order,breaks=(Q(threshold,18),),**options)
        self.local_rows={}
        self.seams={}
        self.empty=set()
        self.local_added=False
        self.original_counts=(len(self.variables),len(self.equalities),len(self.upper_rows))
        if local:
            self.add_local()

    def impossible(self,c):
        return not c.degree and c.slab==1 and c.arity>=17

    def degree_bounds(self,c):
        if c.degree: return Q(1,c.degree),Q(1,c.degree)
        if c.slab==0: return Q(0),Q(1,max(18,c.arity+1))
        return Q(1,17),Q(1,max(self.threshold,c.arity+1))

    def rate_bounds(self,cells,i):
        lo=hi=Q(16-cells[i].arity)
        for j in (i-1,i+1):
            if 0<=j<len(cells):
                a,b=self.degree_bounds(cells[j]); lo+=a; hi+=b
        hi+=sum(Q(self.slots(cells,kind),2) for kind in self.ports(cells) if kind[0]==i)
        return lo,hi

    def add_seam(self,cells,i):
        H,N=self.threshold,self.n
        lo,hi=self.rate_bounds(cells,i)
        left,right=self.boundaries[1:]
        t=(Q(H,17)-left)/(right-left)
        cut=(Q(H,16)-left)/(right-left)
        weights=tuple(v/(1+Q(1,17))**self.m for v in basis(N,t))
        other=tuple(j for j in range(len(cells)) if j!=i)
        ids={}
        # Leaf-stock never enters these global caps; D at the seam is 17.
        budget=Q(1)/self.low
        group=(len(cells),i,'D17 seam')
        for idx in product(*(range(N+1) if not cells[j].degree else (0,) for j in other)):
            a=self.variable(('alpha',cells,i,idx),group+('alpha',),budget)
            g=self.variable(('Gamma',cells,i,idx),group+('Gamma',),20*budget)
            ids[idx]=(a,g)
            self.add_ineq({a:lo,g:-1}); self.add_ineq({a:-hi,g:1})
            mass=[]
            for j,w in enumerate(weights):
                index=idx[:i]+(j,)+idx[i:]
                cols=self.blocks[cells,'mass'][index]
                mass.append(cols)
                self.add_ineq({a:w,**{u:-1 for u in cols}})
                if cells[i].arity==16:
                    # A nonfounder with k=16 already has D>=17. On
                    # the upper reciprocal slab its entire mass is D=17.
                    self.add_eq({a:-w,**{u:1 for u in cols}})
            # In the normalized slab coordinate t, (t-cut)B_{N-1,j}
            # = A B_{N,j}+B B_{N,j+1}. Positive affine scale is omitted.
            for j in range(N):
                A=-cut*Q(N-j,N)
                B=(1-cut)*Q(j+1,N)
                self.add_ineq({a:A*weights[j]+B*weights[j+1],
                               **{u:-A for u in mass[j]},**{u:-B for u in mass[j+1]}})
        self.seams[cells,i]=ids

    def crossing(self,row,cells,e,i,sign):
        other=tuple(j for j in range(len(cells)) if j!=i)
        axes=[self.axis(cells[j],e[j]) for j in other]
        at=Q(self.threshold,18)**e[i]
        for idx,(_,g) in self.seams[cells,i].items():
            row[g]+=sign*at*prod(a[k] for a,k in zip(axes,idx))

    def local(self,cells,e):
        return self.local_rows.setdefault((cells,e),defaultdict(Q))

    def gain(self,row,cells,e,i,sign,*,difference=False):
        if len(cells)==3:
            self.chain_gain(row,cells,e,i,sign,difference=difference)
        else:
            self.add_gain(row,cells,e,i,sign,difference=difference)

    def add_local(self):
        if self.local_added:
            raise ValueError('slab balances already installed')
        self.local_added=True
        support=[{(c,) for c in self.cells},set(self.edges),set(self.chains)]
        for family in support:
            for cells in sorted(family):
                if any(self.impossible(c) for c in cells):
                    self.empty.add(cells)
                    for kind in ('mass',*self.ports(cells)):
                        for cols in self.blocks[cells,kind].values():
                            self.add_ineq(dict.fromkeys(cols,Q(1)))
                else:
                    for i,c in enumerate(cells):
                        if not c.degree and c.slab==1:
                            self.add_seam(cells,i)
        valid=[family-self.empty for family in support]
        for family in valid:
            for cells in sorted(family):
                for e in self.tests(tuple(c.label for c in cells)):
                    row=self.local(cells,e)
                    self.add_observation(row,cells,e,discount=True)
                    for i,c in enumerate(cells):
                        self.gain(row,cells,e,i,1,difference=not c.degree)
                        if (cells,i) in self.seams:
                            self.crossing(row,cells,e,i,1)
                        self.add_observation(row,cells,e,factor=(i,'birth'))
                for i,c in enumerate(cells):
                    def dest(new): return cells[:i]+(new,)+cells[i+1:]
                    if c.degree:
                        d=c.degree+1
                        to=dest(Cell(d,c.arity) if d<self.threshold else Cell(0,c.arity,1))
                        if to in family:
                            for e in self.tests(tuple(z.label for z in to)):
                                self.gain(self.local(to,e),cells,e,i,-1)
                    elif (cells,i) in self.seams:
                        to=dest(Cell(0,c.arity,0))
                        if to in family:
                            for e in self.tests(tuple(z.label for z in to)):
                                self.crossing(self.local(to,e),cells,e,i,-1)
                    to=dest(c._replace(arity=c.arity+1))
                    if to in family:
                        for e in self.tests(tuple(z.label for z in to)):
                            self.add_observation(self.local(to,e),cells,e,-1,factor=(i,'birth'))
                if len(cells)<3:
                    i=len(cells)-1
                    to=cells[:i]+(cells[i]._replace(arity=cells[i].arity+1),self.virgin)
                    if to in valid[len(cells)]:
                        for e in self.tests(tuple(z.label for z in to)):
                            self.add_observation(self.local(to,e),cells,e[:-1],-1,factor=(i,'birth'))
        # An order-zero slab sum must reproduce every global row. More
        # generally the same statement holds exponent by exponent.
        for (cells,e),row in self.local_rows.items():
            self.add_eq(row,int(cells==(self.virgin,) and e==(0,)))
