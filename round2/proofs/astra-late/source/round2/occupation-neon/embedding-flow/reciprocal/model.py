"""All-degree stationary root and edge tests polynomial in inverse degree."""

from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import comb, prod
from pathlib import Path
import sys
from typing import NamedTuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from solve import LO, HI, optimize


class Cell(NamedTuple):
    degree: int
    arity: int
    slab: int = -1

    @property
    def label(self):
        return self.degree, self.arity


def add_scaled(row, values, scale):
    if scale:
        for index, coefficient in values.items():
            row[index] += scale * coefficient


@cache
def bernstein_on_slab(power, degree, slab, subdivisions):
    """Convert a univariate polynomial exactly on a subinterval of [0,1]."""
    if len(power) > degree+1:
        raise ArithmeticError("polynomial exceeds stored degree")
    left, width = Q(slab,subdivisions), Q(1,subdivisions)
    shifted = [sum((power[j]*comb(j,i)*left**(j-i)*width**i
                    for j in range(i,len(power))), Q(0)) for i in range(len(power))]
    return tuple(sum((v*Q(comb(j,i),comb(degree,i))
                      for i,v in enumerate(shifted) if j>=i), Q(0))
                 for j in range(degree+1))


class Model:
    def __init__(self, threshold, arity, test_order, subdivisions=1, *,
                 edge_balances=True, share=True, low=LO, high=HI):
        if threshold < 3 or arity < 0 or test_order < 0 or subdivisions < 1 or not 0 < low <= high:
            raise ValueError("invalid reciprocal census or discount interval")
        self.threshold, self.arity, self.m, self.slabs = threshold, arity, test_order, subdivisions
        self.low, self.high = low, high
        self.n = 2*test_order+2
        # All small degrees are explicit. The only missing cores have both D>=H and k>K.
        self.cells = tuple(Cell(d,k) for d in range(2,threshold) for k in range(d)) + tuple(
            Cell(0,k,h) for k in range(arity+1) for h in range(subdivisions))
        self.labels = {c.label for c in self.cells}
        self.virgin = Cell(2,0)
        self.equalities, self.upper_rows, self.variables, self.group_bounds = [], [], [], {}
        self.blocks = {}
        self.balance = {}
        self.total = self.variable("total","total",1/low)
        self.add_ineq({self.total:1},1/low)
        self.add_ineq({self.total:-1},-1/high)
        self.roots = [(c,) for c in self.cells]
        self.edges = [(p,x) for p in self.cells if p.arity for x in self.cells] if share or edge_balances else []
        for cells in self.roots+self.edges:
            self.make_block(cells,"mass")
            if len(cells)==1 or edge_balances:
                for i,c in enumerate(cells):
                    if i==0:
                        self.make_block(cells,(i,"mother"))
                    if self.spare(cells,i):
                        self.make_block(cells,(i,"children"))
        self.global_constraints()
        for cells in self.roots:
            self.root_projections(cells,share)
            self.existing_balance(cells)
            if edge_balances:
                self.edge_entrance(cells)
        if edge_balances:
            for cells in self.edges:
                self.existing_balance(cells)
                for i in range(2):
                    if i==0:
                        self.box_port(cells,(i,"mother"),1)
                    if self.spare(cells,i):
                        self.box_port(cells,(i,"children"),self.spare(cells,i))
        for key,row in self.balance.items():
            self.add_eq(row, int(key == ((self.virgin.label,),(0,))))

    def variable(self, name, group, upper):
        self.group_bounds[group] = upper
        self.variables.append((name,group))
        return len(self.variables)-1

    def add_eq(self,row,rhs=0):
        self.equalities.append(({i:v for i,v in row.items() if v},Q(rhs)))

    def add_ineq(self,row,rhs=0):
        self.upper_rows.append(({i:v for i,v in row.items() if v},Q(rhs)))

    @staticmethod
    def spare(cells,i):
        return cells[i].arity-int(i==0 and len(cells)==2)

    def make_block(self,cells,kind):
        indices = tuple(product(*(range(self.n+1) if not c.degree else (0,) for c in cells)))
        group = len(cells),kind
        # Each core has one mother in the limit. Total edge and root masses are both 1/mu.
        factor = 1
        if kind!="mass":
            factor = Q(1,2) if kind[1]=="mother" else Q(max(self.arity,self.threshold-2),2)
        by_index = {}
        for idx in indices:
            lanes = (self.low,self.high) if kind=="mass" and self.low<self.high else (self.low,)
            by_index[idx] = tuple(self.variable((cells,kind,idx,s),group,factor/self.low) for s in lanes)
        self.blocks[cells,kind] = by_index

    @cache
    def axis(self,c,e,form="old",factor="one"):
        """Coefficients of W times the test or the within-tail old-minus-new difference."""
        d,k,slab = c
        if d:
            if form!="old":
                raise ValueError("small-degree gain must enter its next census")
            x=Q(1,d)
            return ({"one":Q(1),"p":x,"birth":1-(k+1)*x,"leaves":Q(d-k-1)}[factor],)
        H, m = self.threshold,self.m
        if form=="difference":
            if e==0:
                return (Q(0),)*(self.n+1)
            if e>m:
                raise ValueError("gain-difference test exceeds generator degree")
        length=e+m+1
        polynomial=[Q(0)]*length
        for j in range(m+1):
            polynomial[e+j]=Q(comb(m,j),H**j)
        if form=="difference":
            # W(t)(t^e - (t/(1+t/H))^e) has a double zero at t=0.
            for j in range(m-e+1):
                polynomial[e+j]-=Q(comb(m-e,j),H**j)
        if factor=="leaves":
            if polynomial[0]:
                raise ArithmeticError("uncompensated high-degree leaf rate")
            polynomial=[H*a-(k+1)*b for a,b in zip(polynomial[1:]+[0],polynomial)]
        elif factor=="p":
            polynomial=[Q(0)]+[a/H for a in polynomial]
        elif factor=="birth":
            polynomial=[b-Q(k+1,H)*a for a,b in zip([0]+polynomial,polynomial+[0])]
        elif factor!="one":
            raise ValueError(factor)
        return bernstein_on_slab(tuple(polynomial),self.n,slab,self.slabs)

    def add_product(self,row,cells,kind,axes,sign=1,discount=False):
        if not sign or any(not any(a) for a in axes):
            return
        for idx,cols in self.blocks[cells,kind].items():
            coefficient=sign*prod(a[j] for a,j in zip(axes,idx))
            if coefficient:
                for i,s in zip(cols,(self.low,self.high)):
                    row[i] += coefficient*s if discount else coefficient

    def tests(self,labels):
        return product(*(range(self.m+1) if not d else (0,) for d,k in labels))

    def row(self,labels,exps):
        return self.balance.setdefault((tuple(labels),tuple(exps)),defaultdict(Q))

    def add_observation(self,row,cells,e,sign=1,kind="mass",factor=None,discount=False):
        axes=[self.axis(c,j,factor=(factor[1] if factor and i==factor[0] else "one"))
              for i,(c,j) in enumerate(zip(cells,e))]
        self.add_product(row,cells,kind,axes,sign,discount)

    def add_gain(self,row,cells,e,receiver,sign,difference=False):
        form = "difference" if difference else "old"
        x = cells[receiver]
        axes=[self.axis(c,j) for c,j in zip(cells,e)]
        axes[receiver]=self.axis(x,e[receiver],form,"leaves")
        self.add_product(row,cells,"mass",axes,sign)
        axes[receiver]=self.axis(x,e[receiver],form)
        if len(cells)==2:
            other=1-receiver
            internal=list(axes)
            internal[other]=self.axis(cells[other],e[other],factor="p")
            self.add_product(row,cells,"mass",internal,sign)
        if receiver==0:
            self.add_product(row,cells,(receiver,"mother"),axes,sign)
        if self.spare(cells,receiver):
            self.add_product(row,cells,(receiver,"children"),axes,sign)

    def existing_balance(self,cells):
        labels=tuple(c.label for c in cells)
        for e in self.tests(labels):
            row=self.row(labels,e)
            self.add_observation(row,cells,e,discount=True)
            for i,c in enumerate(cells):
                self.add_gain(row,cells,e,i,1,difference=not c.degree)
                self.add_observation(row,cells,e,factor=(i,"birth"))
        for i,c in enumerate(cells):
            d,k,_=c
            if d:
                new=(d+1 if d+1<self.threshold else 0,k)
                if new in self.labels:
                    dest=labels[:i]+(new,)+labels[i+1:]
                    # At the first tail degree H the incoming test argument is exactly one.
                    for e in self.tests(dest):
                        self.add_gain(self.row(dest,e),cells,e,i,-1)
            after=(d,k+1)
            if after in self.labels:
                dest=labels[:i]+(after,)+labels[i+1:]
                for e in self.tests(dest):
                    self.add_observation(self.row(dest,e),cells,e,-1,factor=(i,"birth"))

    def edge_entrance(self,cells):
        c,=cells
        after=c.degree,c.arity+1
        if after in self.labels:
            dest=after,self.virgin.label
            for e in self.tests(dest):
                self.add_observation(self.row(dest,e),cells,(e[0],),-1,factor=(0,"birth"))

    def box_port(self,cells,kind,count):
        for idx,port in self.blocks[cells,kind].items():
            row={port[0]:1}
            for col in self.blocks[cells,"mass"][idx]:
                row[col]=-Q(count,2)
            self.add_ineq(row)

    def root_projections(self,cells,share):
        c,=cells
        for mother in (True,False):
            count=1 if mother else c.arity
            if not count:
                continue
            kind=0,"mother" if mother else "children"
            if not share:
                self.box_port(cells,kind,count)
                continue
            for idx,port in self.blocks[cells,kind].items():
                lower=defaultdict(Q,{port[0]:-1})
                upper=defaultdict(Q,{port[0]:1})
                for i in self.blocks[cells,"mass"][idx]:
                    upper[i]-=Q(count,self.threshold)
                for donor in self.cells:
                    if mother and not donor.arity:
                        continue
                    edge=(donor,c) if mother else (c,donor)
                    weight=self.axis(donor,0)
                    recip=self.axis(donor,0,factor="p")
                    for j,(w,p) in enumerate(zip(weight,recip)):
                        pair=(j,idx[0]) if mother else (idx[0],j)
                        for i in self.blocks[edge,"mass"][pair]:
                            lower[i]+=p
                            upper[i]+=w/self.threshold-p
                self.add_ineq(lower)
                self.add_ineq(upper)

    def global_constraints(self):
        # All omitted roots have k>K and D>=H. Count both ends of every omitted edge.
        omitted=defaultdict(Q,{self.total:1})
        missing_k=defaultdict(Q,{self.total:1})
        missing_b=defaultdict(Q)
        missing_g=defaultdict(Q,{self.total:1})
        for cells in self.roots:
            k=cells[0].arity
            self.add_observation(omitted,cells,(0,),-1)
            self.add_observation(missing_k,cells,(0,),-k)
            self.add_observation(missing_b,cells,(0,),-1,factor=(0,"birth"))
            self.add_observation(missing_g,cells,(0,),-1,kind=(0,"mother"))
            if k:
                self.add_observation(missing_g,cells,(0,),-1,kind=(0,"children"))

        # Missing birth is 1+missing_b, and missing incoming pressure is missing_g-1.
        for v,constant in ((omitted,0),(missing_k,0),(missing_b,1),(missing_g,-1)):
            self.add_ineq({i:-x for i,x in v.items()},constant)
        def linear(parts,rhs=0):
            v=defaultdict(Q)
            for vector,scale in parts:
                add_scaled(v,vector,scale)
            self.add_ineq(v,rhs)

        ratio=1/(1+self.low)
        tail=ratio**(self.arity+1)
        self.add_ineq({**omitted,self.total:omitted[self.total]-tail})
        self.add_ineq({**missing_k,self.total:missing_k[self.total]-tail*(self.arity+1+1/self.low)})
        linear([(omitted,self.arity+1),(missing_k,-1)])
        linear([(missing_b,1),(omitted,-1)],-1)
        # M-B = integral (1+k)/D on omitted roots, bounded by (M+K)/H.
        linear([(omitted,1-Q(1,self.threshold)),(missing_b,-1),
                (missing_k,-Q(1,self.threshold))],1)
        linear([(missing_g,1),(omitted,-Q(1,2)),(missing_k,-Q(1,2))],1)
        if self.edges:
            unmapped=defaultdict(Q,{self.total:1})
            for cells in self.edges:
                self.add_observation(unmapped,cells,(0,0),-1)
            linear([(unmapped,1),(omitted,-1),(missing_k,-1)])

    def objective(self,name):
        cost=defaultdict(Q)
        if name=="birth-virgin":
            self.add_observation(cost,(self.virgin,),(0,),factor=(0,"birth"))
        elif name=="old-mother-birth":
            for cells in self.roots:
                c,=cells
                if not c.degree and c.arity>=2:
                    self.add_observation(cost,cells,(0,),factor=(0,"birth"))
        elif name=="mixed-edge":
            for cells in self.edges:
                p,x=cells
                if not p.degree and x==self.virgin:
                    self.add_observation(cost,cells,(0,0),factor=(0,"p"))
        else:
            raise ValueError(name)
        return cost
