"""Certified active rows and a complete account of boundary width by first exit."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache
from math import ceil, floor
from pathlib import Path
import sys
from time import monotonic

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'hub'))
import hybrid_rows as hy
from transfer_check import A, BETA, info, mark, newborn, physical, plus

UNIT = 10**23
RAY = 32
TAGS = (2, 16)


@cache
def geometry(tree):
    degree, _ = info(tree)
    parts = [geometry(y) for y in tree[2]]
    return (max([degree]+[x[0] for x in parts]), 1+sum(x[1] for x in parts),
            int(degree>6)+sum(x[2] for x in parts))


def retag(tree, tag):
    return (tag,tree[1],tree[2])


def category(tree):
    if not tree[2]:
        assert info(tree)[0]==RAY
        return 'star_tail'
    maxd, cores, high = geometry(tree)
    if not high:
        name='low'
    elif high==1:
        name='one_root' if info(tree)[0]>6 else 'one_below'
    else:
        name='several'
    return f'{name}_{cores}'


def enumerate_shapes(low, high_cores, high_degree, high_position='all'):
    def keep(t):
        if not t[2]:
            return info(t)[0]<RAY
        degree, cores, high = geometry(t)
        eligible=(high_position=='all' or (info(t)[0]>6)==(high_position=='root'))
        return (info(t)[1]<=9 or cores==2 and degree<16 or cores<=low and not high
                or cores<=high_cores and high==1 and degree<=high_degree and eligible)
    seen={newborn()}
    pending=list(seen)
    while pending:
        tree=pending.pop()
        for nxt,_ in physical(tree):
            if nxt not in seen and keep(nxt):
                assert all(y[0]==2 for y in nxt[2])
                seen.add(nxt)
                pending.append(nxt)
    return sorted(seen,key=lambda t:(info(t)[1],t))


class Sector:
    def __init__(self, s, low, high_cores, high_degree, high_position='all'):
        hy.ex.S=hy.rows.S=hy.tc.S=s
        hy.rows.N=RAY
        hy.tc.ORDER=2
        hy.rows.direct_pair=hy.tc.improved_direct_pair
        shapes=enumerate_shapes(low,high_cores,high_degree,high_position)
        self.shapes=shapes
        size=len(shapes)
        self.size=size
        index={t:i for i,t in enumerate(shapes)}
        self.start=index[newborn()]
        self.trees=[retag(t,tag) for tag in TAGS for t in shapes]
        self.weights=[hy.ex.weight(t) for t in self.trees]
        self.moves=[None]*(2*size)
        self.exits=[None]*(2*size)
        self.slack=[None]*(2*size)
        self.direct=[None]*(2*size)
        for i,t in enumerate(shapes):
            d,n=info(t)
            p=F(1,d)
            intr=physical(t)
            copies=[]
            for y in t[2]:
                assert mark(y)==A
                if A>p:
                    copies.append((y,A-p))
                copies.append((plus(y),p))
            for z,tag in enumerate(TAGS):
                k=z*size+i
                den=s+F(1,tag)+n-p
                self.direct[k]=F(t[1],d)/den
                moves=[]
                exits=[]
                slack=self.weights[k]
                for is_copy,edges in ((False,intr),(True,copies)):
                    for nxt,rate in edges:
                        prob=rate/den
                        j=index.get(nxt)
                        if j is None:
                            out=retag(nxt,2 if is_copy else tag)
                            exits.append((out,prob))
                        else:
                            j=j if is_copy else z*size+j
                            assert is_copy or j>k
                            moves.append((is_copy,j,prob.numerator,prob.denominator))
                            slack-=prob*self.weights[j]
                assert slack>0
                self.moves[k]=moves
                self.exits[k]=exits
                self.slack[k]=slack
        simplex=Counter(geometry(t)[1] for t in shapes if geometry(t)[0]<=6)
        print('label-erased low shapes by core count',dict(sorted(simplex.items())))
        print('distinct shapes',size,'labelled multicore rows',2*(size-(RAY-2)),
              'labelled stars',2*(RAY-2),'live edges',sum(map(len,self.moves)),
              'exiting edges',sum(map(len,self.exits)),flush=True)

    def occupation(self, tag, generations):
        """Positive directed entrance masses. Intrinsic order is increasing size."""
        length=2*self.size
        seed=(tag==16)*self.size+self.start
        lower=[0]*length
        upper=lower.copy()
        lower[seed]=upper[seed]=UNIT
        visits_lo=[0]*length
        visits_hi=[0]*length
        for _ in range(generations):
            next_lo=[0]*length
            next_hi=[0]*length
            for k,edges in enumerate(self.moves):
                hi=upper[k]
                if not hi:
                    continue
                lo=lower[k]
                visits_lo[k]+=lo
                visits_hi[k]+=hi
                for is_copy,j,n,d in edges:
                    a,b=(next_lo,next_hi) if is_copy else (lower,upper)
                    a[j]+=lo*n//d
                    b[j]+=(hi*n+d-1)//d
            lower,upper=next_lo,next_hi
        tail=sum(m*w for m,w in zip(upper,self.weights) if m)/UNIT
        return visits_lo,visits_hi,tail

    def payments(self,col):
        @cache
        def outside(t):
            if t[2]:
                lo,hi=hy.ex.boundary(t,col)
            else:
                l,u=hy.tc.star_terminal(t[0],RAY,col)
                lo,hi=F(l,hy.tc.UNIT),F(u,hy.tc.UNIT)
            assert 0<=lo<=hi
            return lo,hi-lo,category(t)
        points=[]
        bound=F(0)
        for t,direct,exits,slack in zip(self.trees,self.direct,self.exits,self.slack):
            d,_=info(t)
            low=direct*hy.rows.g(col,d)
            diff=defaultdict(F)
            for out,p in exits:
                l,w,c=outside(out)
                low+=p*l
                diff[c]+=p*w
            bound=max(bound,(low+sum(diff.values()))/slack)
            points.append([('base',low)]+list(diff.items()))
        return points,F(ceil(1000*bound),1000)


def bracket(a,b,digits=12):
    scale=10**digits
    def fmt(x):
        return f'{x//scale}.{x%scale:0{digits}d}'
    return f'[{fmt(floor(a*scale))},{fmt(ceil(b*scale))}]'


def evaluate(points, bound, visits):
    vlo,vhi,tail=visits
    lo=defaultdict(int)
    hi=defaultdict(int)
    for m,M,source in zip(vlo,vhi,points):
        if not M:
            continue
        for label,prob in source:
            n,d=prob.numerator,prob.denominator
            lo[label]+=m*n//d
            hi[label]+=(M*n+d-1)//d
    extra=tail*bound
    lower={k:F(v,UNIT) for k,v in lo.items()}
    upper={k:F(v,UNIT) for k,v in hi.items()}
    low=lower.pop('base')
    high=upper.pop('base')
    print('lower endpoint',bracket(low,high+extra),
          'upper endpoint',bracket(low+sum(lower.values()),high+sum(upper.values())+extra))
    print('total width',bracket(sum(lower.values()),sum(upper.values())+extra),
          'unvisited upper',ceil(extra*10**16),'/ 1e16')
    for label in sorted(lower):
        print('width',label,bracket(lower[label],upper[label]+extra))
    return low,high+sum(upper.values())+extra


if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('--s',type=F,default=F('.565809'))
    p.add_argument('--low',type=int,default=4)
    p.add_argument('--high-cores',type=int,default=0)
    p.add_argument('--high-degree',type=int,default=7)
    p.add_argument('--high-position',choices=('all','root','below'),default='all')
    p.add_argument('--columns',default='B')
    p.add_argument('--seeds',default='B')
    p.add_argument('--generations',type=int,default=27)
    args=p.parse_args()
    assert 2<=args.low<=7 and 0<=args.high_cores<=7 and 7<=args.high_degree<=15
    assert args.s>hy.rows.Q+hy.rows.DAMP
    stamp=monotonic()
    print('s',args.s,'full size 9, two-core degree <16, low cores',args.low,
          'one-high extension',args.high_cores,args.high_degree,args.high_position,'generations',args.generations,
          flush=True)
    model=Sector(args.s,args.low,args.high_cores,args.high_degree,args.high_position)
    print('compiled seconds',round(monotonic()-stamp,3),flush=True)
    occ={se:model.occupation({'A':2,'B':16}[se],args.generations)
         for se in args.seeds.split(',')}
    print('compiled plus forward seconds',round(monotonic()-stamp,3),flush=True)
    for name in args.columns.split(','):
        col={'A':A,'B':BETA,'half':(F(1,2),F(1))}.get(name)
        if col is None:
            col=F(name),F(1)
        points,bound=model.payments(col)
        for seed,visits in occ.items():
            print('column',name,'seed',seed,'majorant',bound)
            evaluate(points,bound,visits)
        print('cumulative seconds',round(monotonic()-stamp,3),flush=True)
