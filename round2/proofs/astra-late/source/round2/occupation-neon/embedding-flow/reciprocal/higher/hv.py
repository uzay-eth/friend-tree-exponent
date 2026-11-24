"""Original D18 hV allocated row, in exact reciprocal polynomials (no LP spend)."""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import comb, prod
from pathlib import Path
import sys

from sympy import Poly, QQ, symbols

OCC = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(OCC/'higher'))
from extract import Cert, Extract, PAYLOAD_SHA

H = 18
Y = symbols('y0:3')
B = tuple(comb(3,j)*Y[0]**j*(1-Y[0])**(3-j) for j in range(4))
NB = (1,0,1)


def P(x=0):
    return Poly(x, *Y, domain=QQ)


ZERO = P()
ONE = P(1)
AXIS = tuple(P(y) for y in Y)
SHIFT = tuple((ONE+y.mul_ground(Q(1,H)))**3 for y in AXIS)


@dataclass(frozen=True)
class Role:
    """degree=None is a physical tail coordinate, at fixed actual arity."""
    degree: int | None
    k: int


@dataclass(frozen=True)
class Slot:
    role: Role
    axis: int

    @property
    def tail(self): return self.role.degree is None

    @property
    def idx(self):
        d,k=self.role.degree,self.role.k
        return (d-k-1,k,1) if d is not None else None

    @property
    def p(self): return AXIS[self.axis].mul_ground(Q(1,H)) if self.tail else P(Q(1,self.role.degree))

    @property
    def b(self): return ONE-(self.role.k+1)*self.p

    def bud(self): return Slot(Role(self.role.degree,self.role.k+1),self.axis)

    def gain(self): return Slot(Role(self.role.degree+1,self.role.k),self.axis)


@dataclass
class Row:
    denominator: Poly
    base: Poly
    ports: dict
    terms: dict
    roles: tuple

    def at(self, degrees, ports):
        args=[]
        for r,d in zip(self.roles,degrees,strict=True):
            if type(d) is not int or (r.degree is None and d<max(H,r.k+1)) or (r.degree is not None and d!=r.degree):
                raise ValueError('degree is outside this physical cell')
            args.append(Q(H,d) if r.degree is None else Q(0))
        if set(ports)-self.ports.keys():
            raise ValueError('allocated or absent port supplied')
        def ev(p): return Q(p.eval(tuple(args)))
        n=ev(self.base)
        for key,(count,slope) in self.ports.items():
            rate=Q(ports[key])
            if not 0<=rate<=Q(count,2):
                raise ValueError('impossible live aggregate receive rate')
            n+=rate*ev(slope)
        return n/ev(self.denominator)


def opener(path, witness, *, authenticate=False):
    """Sparse mmap by default; caller is responsible for inherited blob authentication."""
    from hashlib import file_digest
    if witness not in PAYLOAD_SHA:
        raise ValueError('unknown original witness')
    path=Path(path)
    if authenticate:
        with path.open('rb') as f:
            if file_digest(f,'sha256').hexdigest()!=PAYLOAD_SHA[witness]:
                raise ValueError('original D18 digest mismatch')
    c=Cert(path)
    if (c.D,c.K,c.S,c.Sr,c.mode,c.lam)!=(18,17,10**11,10**9,
            2 if witness=='lower' else 1,Q('.565673' if witness=='lower' else '.565775')):
        raise ValueError('unexpected D18 witness header')
    return c


class Reader:
    def __init__(self,c):
        self.c=c
        self.raw=Extract(c)

    @cache
    def embed(self,s):
        if not s.tail:
            if s.role.degree<H:
                return ((s.idx,ONE),)
            # This is only the actual degree-17 gain successor.
            if s.role.degree!=H:
                raise ValueError('unexpected fixed degree in a tail successor')
            return ((('T3',min(s.role.k,self.c.K+1),1),ONE),)
        k=min(s.role.k,self.c.K+1)
        states=(('inf',),*((f'T{j}',k,1) for j in (1,2,3)))
        return tuple(zip(states,(P(b.subs(Y[0],Y[s.axis])) for b in B)))

    def weight(self,f,x):
        return sum((prod(ws)*self.raw.weight(f,ids) for entry in product(*(self.embed(s) for s in x))
                    for ids,ws in [zip(*entry)]),ZERO)

    def edge_transfer(self,m,neighbor,from_above):
        """Original g→hV transfer coefficient before multiplying by p_donor−center."""
        c=self.c
        if not m.tail:
            if from_above: # donor is apex; old edge is (m,tip)
                return self.weight('g',(m.gain(),neighbor))-self.weight('g',(m,neighbor))
            return self.weight('g',(neighbor,m.gain()))-self.weight('g',(neighbor,m))
        k=m.role.k
        t1=('T1',k,1)
        y=AXIS[m.axis]
        # T1 and INF have no incoming transfer, T2 has the divided
        # first difference, T3 uses its own original reference array.
        idx=c.state[t1]
        if from_above:
            def controls(a):
                return (Q(int(c.g[c.state[('inf',)],c.state[a]])-int(c.g[idx,c.state[a]]),H),
                        int(c.yp[k,c.child[a]]))
        else:
            def controls(a):
                return (Q(int(c.g[c.state[a],c.state[('inf',)]])-int(c.g[c.state[a],idx]),H),
                        int(c.yw[c.par[a],k]))
        return sum((w*(3*y*y*(ONE-y)*t2+y**3*t3)
                    for a,w in self.embed(neighbor) for t2,t3 in [controls(a)]),ZERO)

    def center(self,s,above):
        c=self.c
        arr=(c.trpar if above else c.trch) if s.tail else (c.rpar if above else c.rch)
        return Q(int(arr[s.role.k if s.tail else c.state[s.idx]]),c.Sr)

    def row(self,roles):
        """Return unsigned raw witness drift divided by S; no clipping at zero."""
        roles=tuple(roles)
        if len(roles)!=3 or any(type(r.k) is not int or r.k<int(i<2) or
                  (r.degree is not None and (type(r.degree) is not int or not 2<=r.degree<H or r.k>=r.degree))
                  for i,r in enumerate(roles)):
            raise ValueError('not a valid oriented D18 cell')
        x=tuple(Slot(r,i) for i,r in enumerate(roles))
        den=prod((SHIFT[i] for i,s in enumerate(x) if s.tail),start=ONE)
        if roles[1].k>self.c.K:
            return Row(den,ZERO,{}, {},roles)
        own=self.weight('hV',x)
        terms={'discount':-self.c.lam*own*den}
        ports={}
        cut=lambda i,arr: not x[i].tail and x[i].idx in arr
        p4=lambda i: cut(i,self.c.p4)
        routes={(0,'mother'):p4(0) and p4(1),
                (0,'children'):cut(0,self.c.ybu) and cut(1,self.c.ybw),
                (1,'children'):cut(1,self.c.yt), (2,'children'):p4(1) and p4(2)}
        virgin=Slot(Role(2,0),0)
        for i,s in enumerate(x):
            replace=lambda v:x[:i]+(v,)+x[i+1:]
            if s.tail:
                # Clear only the receiver's own shift denominator. The
                # resulting difference is divisible by y_i^2, not just y_i.
                y=AXIS[i]
                head=sum((Poly.from_dict({tuple(0 if j==i else powr[j] for j in range(3)):c},Y,domain=QQ)
                          *y**powr[i]*(ONE+y.mul_ground(Q(1,H)))**(3-powr[i]) for powr,c in own.terms()),ZERO)
                delta=head-own*SHIFT[i]
                quotient=delta.exquo(y)
                quotient.exquo(y) # certify the removable origin
                rest=den.exquo(SHIFT[i])
                leaf=rest*(H*quotient-(s.role.k+1)*delta)
                gain=rest*delta
            else:
                gain=den*(self.weight('hV',replace(s.gain()))-own)
                leaf=(s.role.degree-s.role.k-1)*gain
            terms[i,'leaves']=leaf
            terms[i,'internal']=gain*sum((x[j].p for j in range(3) if abs(i-j)==1),ZERO)
            for side,count in [('mother',int(i==0)),('children',s.role.k-int(i<2))]:
                if not count: continue
                key=i,side
                if routes[key]:
                    terms[i,side,'centered']=count*self.center(s,side=='mother')*gain
                else:
                    ports[key]=(count,gain.mul_ground(Q(1,self.c.S)))
            if s.b!=ZERO:
                post=replace(s.bud())
                w=self.weight('hV',post) if post[1].role.k<=self.c.K else ZERO
                # Four-body newborn ownership tests are POST-promotion. The
                # other old roles remain physical and unshifted.
                a,m,z=post
                if i==0 and not a.tail and not m.tail and a.idx in self.c.ybu and m.idx in self.c.ybw:
                    newborn=self.weight('Yb',(a,virgin,m,z))
                elif i==1 and not m.tail and m.idx in self.c.yt:
                    newborn=self.weight('Yt',(a,m,z,virgin))
                elif i==2 and not m.tail and not z.tail and m.idx in self.c.p4 and z.idx in self.c.p4:
                    newborn=self.weight('P4',(a,m,z,virgin))
                else: newborn=ZERO
                terms[i,'promotion']=den*s.b*(w-own)
                terms[i,'newborn']=den*s.b*newborn
        terms['g→hV','from apex']=den*(x[0].p-self.center(x[1],True))*self.edge_transfer(x[1],x[2],True)
        terms['g→hV','from tip']=den*(x[2].p-self.center(x[1],False))*self.edge_transfer(x[1],x[0],False)
        terms={k:v.mul_ground(Q(1,self.c.S)) for k,v in terms.items()}
        return Row(den,sum(terms.values(),ZERO),ports,terms,roles)


def contrast(lower,upper,roles,theta=Q(1)):
    """Sign-lower minus theta sign-upper = raw lower + theta raw upper."""
    a,b=(reader.row(roles) for reader in (lower,upper))
    if a.denominator!=b.denominator or {k:v[0] for k,v in a.ports.items()}!={k:v[0] for k,v in b.ports.items()}:
        raise ValueError('witness allocation mismatch')
    terms={k:a.terms.get(k,ZERO)+theta*b.terms.get(k,ZERO) for k in a.terms.keys()|b.terms.keys()}
    return Row(a.denominator,a.base+theta*b.base,
               {k:(n,v+theta*b.ports[k][1]) for k,(n,v) in a.ports.items()},terms,a.roles)


def bernstein(poly, bounds, degrees):
    """Exact tensor Bernstein coefficients over physical reciprocal intervals."""
    if any(l>=h for l,h in bounds):
        raise ValueError('empty Bernstein interval')
    power={}
    for exp,c in poly.terms():
        if any(e>d for e,d in zip(exp,degrees,strict=True)):
            raise ValueError('control degree too small')
        for i in product(*(range(e+1) for e in exp)):
            val=Q(c)*prod(comb(e,j)*l**(e-j)*(h-l)**j for e,j,(l,h) in zip(exp,i,bounds,strict=True))
            power[i]=power.get(i,0)+val
    return {j:sum((v*prod(Q(comb(a,b),comb(d,b)) for a,b,d in zip(j,i,degrees,strict=True))
                   for i,v in power.items() if all(a>=b for a,b in zip(j,i))),Q(0))
            for j in product(*(range(d+1) for d in degrees))}


def moment_controls(row, order, bounds=None):
    """Coefficients for moments of B_j(y)/W(y), W=prod(1+y_i/H)^order.

    The base uses the triple mass block; a slope uses the *same* triple's
    rate-weighted port block. Never replace that block by a marginal rate.
    """
    if type(order) is not int or order<3:
        raise ValueError('reciprocal test order must include the cubic gain denominator')
    domain=tuple((Q(0),min(Q(1),Q(H,r.k+1))) for r in row.roles)
    bounds=domain if bounds is None else tuple(bounds)
    if len(bounds)!=3 or any(not 0<=l<h<=top for (l,h),(_,top) in zip(bounds,domain)):
        raise ValueError('invalid reciprocal slab')
    degrees=tuple(2*order+2 if r.degree is None else 0 for r in row.roles)
    W=prod(((ONE+AXIS[i].mul_ground(Q(1,H)))**order for i,r in enumerate(row.roles)
            if r.degree is None), start=ONE)
    factor=W.exquo(row.denominator)
    return (bernstein(row.base*factor,bounds,degrees),
            {key:(count,bernstein(slope*factor,bounds,degrees)) for key,(count,slope) in row.ports.items()})
