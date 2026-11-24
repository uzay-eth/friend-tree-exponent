"""Independent exact small-crown occupation and sparse Neon f-table prices."""

from fractions import Fraction as Q
from functools import cache
from pathlib import Path
from struct import Struct
from collections import defaultdict
import argparse
import hashlib
import mmap

LOW = Q('.565673570322')
HIGH = Q('.565775')
J = (1, ())
TICKS = 10**16
PINS = {
    1: '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41',
    2: '9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
}


class Price:
    def __init__(self, path):
        self.file = Path(path).open('rb')
        self.data = mmap.mmap(self.file.fileno(), 0, access=mmap.ACCESS_READ)
        if len(self.data) != 3495091968:
            raise ValueError(f'{path}: wrong production payload size')
        magic, mode, version, p, q, self.S, self.Sr, *caps = Struct('<Qii4q12i').unpack_from(self.data)
        if mode not in PINS or hashlib.sha256(self.data).hexdigest() != PINS[mode]:
            raise ValueError(f'{path}: SHA-256 does not match the pinned production payload')
        print(f'{path}: SHA-256 verified ({"upper" if mode == 1 else "lower"})', flush=True)
        assert magic == 0x4e454f4e42334648 and version == 1
        assert caps == [14000, 180, 18, 17, 10, 7, 9, 6, 10, 7, 8, 6]
        self.lam, self.sign = Q(p, q), (1 if mode == 1 else -1)
        assert (mode, self.lam) in ((1, HIGH), (2, Q('.565673')))
        self.M = Q(195820675182 if mode == 1 else 184124847315, self.S)
        self.states = {s: i for i,s in enumerate((l-k-h,k,h) for h in (1,0)
                       for l in range(2,18) for k in range(l-h+1))}
        i = len(self.states)
        self.states['inf'] = i
        for j in range(1,4):
            for h in (1,0):
                for k in range(19):
                    i += 1
                    self.states[f'T{j}',k,h] = i
        self.ng = i+1
        inner = 320
        assert self.states['inf'] == inner
        self.rch = inner
        self.tail = 2*inner
        children = sum(s == 'inf' or s[-1]==1 for s in self.states)
        parents = sum(s == 'inf' or s[1]>=1 for s in self.states)
        self.f0 = self.tail+4*38+38*children+19*parents
        self.g0 = self.f0+2*181*14001+2

    def get(self, index):
        return Struct('<q').unpack_from(self.data,96+8*index)[0]

    def f(self, l,k):
        return self.get(self.f0+min(l,14000)+k*14001)

    def g(self, s):
        return self.get(self.g0+self.states[s]*self.ng+self.states[1,0,1])

    def bud_edge(self,l,k):
        d=l+k+1
        if d<18:
            return self.g((l,k,1))
        t=Q(18,d); j=min(k,18)
        return sum(w*self.g(s) for w,s in zip(((1-t)**3,3*t*(1-t)**2,3*t*t*(1-t),t**3),
                                               ('inf',('T1',j,1),('T2',j,1),('T3',j,1))))

    @cache
    def row(self,l,k):
        d=l+k+1; v=self.f(l,k); inc=self.f(l+1,k)-v
        b=Q(l,d)*(self.f(l-1,k+1)+self.f(1,0)-v+self.bud_edge(l-1,k+1)) if l else 0
        base=l*inc+b-self.lam*v
        coef=0
        if d<18:
            i=self.states[l,k,1]
            base+=Q(self.get(i)+k*self.get(self.rch+i),self.Sr)*inc
        elif k<=17:
            cp=self.tail+k
            center=Q(self.get(cp)+k*self.get(cp+38),self.Sr)
            psi=Q(self.get(cp+76),d*d)+Q(self.get(cp+114),d**3)
            base+=center*psi
            coef=inc-psi
        else:
            coef=inc
        return -self.sign*base/self.S,-self.sign*Q(coef,self.S)


@cache
def n(x):
    return 1+x[0]+sum(map(n,x[1]))


def degree(x):
    return 1+x[0]+len(x[1])


def gain(x):
    return x[0]+1,x[1]


@cache
def moves(x,ext=Q(0)):
    # Target hazards from a labelled adjacency graph, then canonicalize.
    tree = []
    def flatten(y, parent):
        i=len(tree); tree.append(parent)
        tree.extend([i]*y[0])
        for w in y[1]:
            flatten(w,i)
    flatten(x,-1)
    adj=[[] for _ in tree]
    for i,p in enumerate(tree[1:],1):
        adj[i].append(p);adj[p].append(i)
    dd=[len(a)+(i==0) for i,a in enumerate(adj)]
    out=defaultdict(Q)
    for w,nb in enumerate(adj):
        rate=sum(Q(1,dd[i]) for i in nb)+(ext if w==0 else 0)
        def encode(i):
            kids=[]; leaf=0
            for j in adj[i]:
                if tree[j]!=i:
                    continue
                if dd[j]+(j==w)==1:
                    leaf+=1
                else:
                    kids.append(encode(j))
            return leaf+(i==w),tuple(sorted(kids))
        out[encode(0)]+=rate
    assert sum(out.values())==n(x)-Q(1,degree(x))+ext
    assert all(n(y)==n(x)+1 for y in out)
    return tuple(out.items())


def value(c, cap, exterior, initial=J):
    @cache
    def rec(x):
        if n(x)>cap:
            return Q(0)
        l,ch=x
        a,b=c.row(l,len(ch))
        inner=sum(Q(1,degree(y)) for y in ch)
        d=degree(x)
        assert all(a+b*r>=0 for r in (0,Q(len(ch)+1,2)))
        continuation=sum(r*rec(y) for y,r in moves(x))
        ans=min((a+b*(inner+p)+continuation+p*rec(gain(x))) /
                (HIGH+n(x)-Q(1,d)+p) for p in exterior)
        return Q(ans*TICKS//1,TICKS)
    return rec(initial)


def rich_birth_mass(d):
    # At a stable-age newborn gain d-2 leaves, then promote once at D=d.
    # Parent incoming may vary arbitrarily in [0,1/2].
    assert d>=2
    path=Q(1)
    for l in range(1,d):
        want=Q(l) if l<d-1 else Q(l,d)
        path*=want/(HIGH+l+Q(l,l+1)+Q(1,2))
    # Compensating that one selected birth per initial newborn gives a
    # fraction of the stationary newborn Palm intensity, not of all cores.
    return path


def rho_star(c,l,k):
    # Not an extra credit: an improved floor on a tail row already counted
    # by its affine value in value(). Core-neighbour contribution only.
    a,b=c.row(l,k)
    return min(a+b*r for r in (Q(0),Q(k+1,2))), min(a+b*r for r in (Q(1,6),Q(k+1,2)))


def short_promotion(c, d, maximum):
    # Newborn of a parent in [d,maximum]. During the next short window, require one
    # focal gain (ext lower rate counts), then one focal promotion, with
    # no more fringe events. Freeze parent degree by prohibiting all gains.
    # Use slot/2 for each ordered event to keep an exact simple bound.
    rate=maximum+5
    h=Q(1,2*(rate+1))
    segment=h/2
    # The exterior degree event and crown event have disjoint targets;
    # only gain into the exterior parent has been forbidden. Crown total
    # rate <= 4+.5, and forbidden exterior hazard <= maximum.
    # After the two jumps, the crown is (one leaf, one virgin core child).
    anchor=(1,(J,))
    a,b=c.row(1,1)
    floor=min(a+b*r for r in (Q(1,2)+Q(1,maximum), Q(1,2)+Q(1,d)))
    pending=(1+Q(1,maximum))*Q(2,3)*segment**2*(1-rate*h)
    # Persist another h for an occupation lower bound.
    conditional=pending*h*(1-(rate+HIGH)*h)*(1-HIGH*h)
    return anchor,floor,conditional


def root_value(c, cap, parent_lower_degree):
    # Entire crown collapsed to the root's two evolving counts: each
    # child's inverse degree remains between 0 and 1/2, even adaptively.
    @cache
    def rec(l,k):
        d=l+k+1
        if d>cap:
            return Q(0)
        a,b=c.row(l,k)
        boot=Q(l,d)
        next_k=rec(l-1,k+1) if l else 0
        next_d=rec(l+1,k)
        hi=Q(1,parent_lower_degree)+Q(k,2)
        answer=min((a+b*r+(l+r)*next_d+boot*next_k)/(HIGH+l+r+boot)
                   for r in (0,hi))
        assert answer>=0
        return Q(answer*TICKS//1,TICKS)
    return rec(1,0)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('certs', nargs=2, metavar='PAYLOAD'); parser.add_argument('--cap',type=int,default=12)
    args=parser.parse_args()
    prices = [Price(path) for path in args.certs]
    if {c.sign for c in prices} != {-1, 1}:
        parser.error('supply one original upper and one original lower payload')
    for c in prices:
        print('endpoint',c.lam,'signed stars',[str(c.row(l,0)[0]) for l in (1,2,3)])
        for cap in (2,4,7,args.cap):
            v=value(c,cap,(Q(0),Q(1,2)))
            print(f'cap={cap} reward={float(v):.12g} improvement={float(LOW*v/c.M):.12g}')
        for d,m in ((2,3),(3,4),(5,6),(7,10)):
            anchor,floor,conditional=short_promotion(c,d,m)
            mass=LOW*rich_birth_mass(d)*conditional
            print('parent bin',d,m,'rich crown',anchor,'slack',f'{float(floor):.12g}','pi floor',f'{float(mass):.12g}')
        baseline=root_value(c,64,2)
        refined=baseline+sum(rich_birth_mass(d)*(root_value(c,64,d)-baseline) for d in (3,4,5,6,7,8))
        assert refined>=baseline
        improved=LOW*refined/c.M
        endpoint=c.lam-c.sign*improved
        outward=(endpoint*10**12)//1 if c.sign==-1 else -((-endpoint*10**12)//1)
        assert outward <= 565773985336 if c.sign==1 else outward >= 565674338794
        print('root cap=64 baseline',f'{float(baseline):.12g}','disjoint rich-bin bonus',f'{float(refined-baseline):.12g}',
              'gap',f'{float(improved):.12g}','outward endpoint',Q(outward,10**12))
        print('rho samples',[(l,*(float(q) for q in rho_star(c,l,0))) for l in (17,18,24,39)])
