"""Check every unsigned profile against the independent Python snapshot layout."""

from bisect import bisect_left
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement
import gzip
import sys

from run import resources
from v5 import CHILD, DEGREE, OUT, ROOT, group

sys.path.insert(0, str(ROOT/'round2/motif-high'))
from oracle import Certificate, HEADER, degree
from score_tail import build_lines


def sym_maximum(values, level, arity):
    """Maximum with one member fixed, and with two members fixed for stars."""
    positions = list(combinations_with_replacement(range(level.n), arity))
    sites = np.asarray(positions, dtype=np.int32).T
    one = np.full((len(values), level.n), np.iinfo(np.int64).min, dtype=np.int64)
    two = None
    if arity == 3:
        two = np.full((len(values), level.pairs), np.iinfo(np.int64).min, dtype=np.int64)
        pairs = [[level.pair(t[i], t[j]) for t in positions] for i,j in ((0,1),(0,2),(1,2))]
    for row, v in enumerate(values):
        for idx in sites:
            np.maximum.at(one[row], idx, v)
        if two is not None:
            for idx in pairs:
                np.maximum.at(two[row], idx, v)
    return one, two


class Tables:
    def __init__(self, side):
        self.side = side
        c = self.cert = Certificate(OUT/f'mc-{side}.bin')
        L2, L3, L4 = c.levels.values()
        def raw(family, *shape):
            return np.frombuffer(c.tables[family], dtype='<i8').reshape(shape)
        self.edge = raw('E',L2.n,L2.n)
        self.vertical = raw('V',len(L3.parent),-1,L3.n)
        self.chain = raw('P',len(c.act['P']),len(L4.parent),L4.n)
        self.fork = raw('F',len(c.act['F']),L4.n,L4.n)
        self.pmax = self.chain.max(axis=1)
        self.fmax = self.fork.max(axis=1)
        for raw_data in (self.edge,self.vertical,self.chain,self.fork):
            assert max(abs(int(raw_data.min())),abs(int(raw_data.max()))) < 2**46
        for (a,_),i in c.act['F'].items():
            assert self.fmax[i].max() <= (1<<60)//max(1,a[1]-1)
            self.fmax[i] *= a[1]-1
        self.bottom = {
            'E': self.edge[list(L2.parent)].max(axis=0),
            'V': np.maximum(0,self.vertical.max(axis=(0,1))),
            'P': np.maximum(0,self.pmax.max(axis=0)),
            'F': np.maximum(0,self.fmax.max(axis=0)),
        }
        self.pairs = {
            'C': raw('C',len(c.act['C']),L3.pairs),
            'Y': raw('Y',len(c.act['Y']),len(L4.parent),L4.pairs).max(axis=1),
            'T': raw('T',len(c.act['T']),L4.tris),
        }
        for raw_data in self.pairs.values():
            assert max(abs(int(raw_data.min())),abs(int(raw_data.max()))) < 2**46
        self.apex = {f: a.max(axis=1) for f,a in self.pairs.items()}
        self.one = {}
        for f,rows in self.pairs.items():
            self.one[f],two = sym_maximum(rows,c.levels[c.level(f)],3 if f=='T' else 2)
            if f=='T': self.pair_star=two
        self.vmax = self.vertical.max(axis=0)
        # This bound also covers products by degree in the finite bottom interpolation.
        largest=max(max(abs(int(a.min())), abs(int(a.max()))) for a in self.bottom.values())
        assert largest*(8*c.par['DF']) < 2**63
        self.finite = self.check_comparison()

    def check_comparison(self):
        c=self.cert
        path=OUT/f'mc-mbar-{self.side}.i64'
        target=np.fromfile(path,dtype='<i8')
        off,length=c.sections['f']
        singleton=np.frombuffer(c.data[off:off+length],dtype='<i8')
        assert len(target)==len(singleton)
        apex_terms=defaultdict(list)
        for family,rows in self.apex.items():
            for state,i in c.act[family].items():
                k=state[1]
                copies=k*(k-1)//2 if family!='T' else k*(k-1)*(k-2)//6
                apex_terms[degree(state)].append((k,int(rows[i])*copies))
        for d in range(2,c.par['DF']):
            start=d*(d-1)//2-1
            k=np.arange(d)
            v=singleton[start:start+d].copy()
            for f,vals in self.bottom.items():
                lvl=c.levels[c.level(f)]
                if d<lvl.d0:
                    v+=vals[start:start+d]
                else:
                    corner=lvl.inf+1+np.minimum(k,lvl.kt+1)
                    num=lvl.d0*vals[corner]+(d-lvl.d0)*vals[lvl.inf]
                    v+=-np.floor_divide(-num,d)
            for child,value in apex_terms[d]:
                v[child]+=value
            assert np.array_equal(v,target[start:start+d]),(self.side,'mbar',d)
        tail=OUT/f'mc-mbar-{self.side}-tail.tsv'
        rows=[tuple(map(int,z.split())) for z in tail.read_text().splitlines()[1:]]
        assert len(rows)==c.levels[2].kt+2
        for k,intercept,slope in rows:
            fixed=0
            delta=0
            for f,a in self.bottom.items():
                l=c.levels[c.level(f)]
                fixed+=int(a[l.inf])
                delta+=l.d0*(int(a[l.inf+1+min(k,l.kt+1)])-int(a[l.inf]))
            assert (fixed,delta)==(intercept,slope)
        return target

    def classes(self,l,K,parent,coarse=False):
        corners=defaultdict(set)
        # The largest child band starts at 64. All levels clamp below 73.
        assert l.d0<73 and l.kt+1<73
        for d in range(2,74):
            for k in range(int(parent),d):
                g=group(d,k,K) if not coarse or d<=K else -1
                corners[g].update(i for i,_ in l.mu((d-1-k,k)))
        return {g:np.asarray(sorted(a),dtype=np.int32) for g,a in corners.items()}

    def inputs(self,K,M):
        c=self.cert
        path=OUT/f'moments-v5-{self.side}-k{K}.tsv.gz'
        extrema=np.full((len(DEGREE),len(CHILD)+1,16),np.iinfo(np.int64).min,dtype=np.int64)
        for d in range(M+1,c.par['DF']):
            i=d*(d-1)//2-1
            k=np.arange(d)
            np.maximum.at(extrema[bisect_left(DEGREE,d)],(np.searchsorted(CHILD,k),16*(d-1-k)//d),self.finite[i:i+d])
        parent=self.classes(c.levels[2],K,True)
        child=self.classes(c.levels[2],K,False)
        counts=defaultdict(int)
        with gzip.open(path,'rt') as st:
            next(st)
            for line in st:
                tag,*z=line.split()
                a=list(map(int,z))
                if tag=='v':
                    d,k,_,price,_=a
                    assert price==self.finite[d*(d-1)//2-1+k]
                elif tag=='b':
                    db,kb,q=a[:3]
                    assert a[8]==extrema[db,kb,q]
                elif tag=='e':
                    p,w,_,price=a
                    actual=self.edge[np.ix_(parent[p],child[w])]
                    assert price<=(self.bottom['E'][child[w]]-actual).min()
                else: raise ValueError(tag)
                counts[tag]+=1
        assert counts['v']==M*(M+1)//2-1 and counts['b']==np.count_nonzero(extrema!=np.iinfo(np.int64).min)
        assert counts['e']==len(parent)*len(child)
        return counts

    def tail_lines(self,K,M):
        path=OUT/f'moments-v5-{self.side}-k{K}.tsv.gz'
        cells={}
        with gzip.open(path,'rt') as st:
            header=list(map(int,next(st).split()[1:]))
            for line in st:
                if line.startswith('b '):
                    z=list(map(int,line.split()[1:]))
                    cells[tuple(z[:3])]=z[3:]
        side,_,_,unit,df,p0,p1,p2,kt=header
        frozen=HEADER.unpack_from(self.cert.mapping)
        assert (side,unit,df,p0,p1,p2,kt)==(frozen[1],frozen[5],self.cert.par['DF'],*frozen[7:10],self.cert.par['K2'])
        far={}
        for j in range(16):
            for kb in range(len(CHILD)+1):
                kl=CHILD[kb-1]+1 if kb else 0
                kh=CHILD[kb] if kb<len(CHILD) else None
                q0,q1=F(j,16),F(j+1,16)
                if kh is not None: q0=max(q0,1-F(kh+1,df))
                if q0<=q1:
                    key=len(DEGREE),kb,j
                    far[key]=(q0,q1)
                    cells[key]=df,0,kl,kh
        tail=[tuple(map(int,z.split())) for z in (OUT/f'mc-mbar-{self.side}-tail.tsv').read_text().splitlines()[1:]]
        lines=build_lines(OUT/f'mc-mbar-{self.side}.i64',header,M,DEGREE,CHILD,cells,tail)
        by_slope={}
        for (db,kb,j),opts in lines.items():
            if db==len(DEGREE): continue
            for slope,intercept in opts:
                if slope not in by_slope:
                    by_slope[slope]=np.full((len(DEGREE),len(CHILD)+1,16),2**50,dtype=np.int64)
                by_slope[slope][db,kb,j]=intercept
        for d in range(M+1,df):
            start=d*(d-1)//2-1
            k=np.arange(d); leaves=d-k-1
            keys=bisect_left(DEGREE,d),np.searchsorted(CHILD,k),16*leaves//d
            score=self.finite[start:start+d]
            for slope,intercepts in by_slope.items():
                assert np.all((score-intercepts[keys])*d<=slope*leaves),(self.side,K,d,slope)
        for key,(q0,q1) in far.items():
            _,_,kl,kh=cells[key]
            for slope,a in lines[key]:
                for k,base,delta in tail:
                    if k<min(kl,kt+1) or k>min(kt+1,kt+1 if kh is None else kh): continue
                    for x in (0,F(1,df)):
                        A,B,C=(0,slope-p1,a-p0-base-delta*x-4) if side==1 else (
                            p1,slope-p0+p2*x,a-base-delta*x-4)
                        assert min(A*q*q+B*q+C for q in (q0,q1))>=0
                        if A and 2*A*q0 < -B < 2*A*q1:
                            assert 4*A*C>=B*B
        return sum(map(len,lines.values()))

    def overlap(self,K):
        c=self.cert
        l3,l4=c.levels[3],c.levels[4]
        dom=c.act
        known={group(d,k,K):(d-k-1,k) for d in range(2,K+1) for k in range(d)}
        outer={(n,p):self.classes(c.levels[n],K,p,True) for n in (3,4) for p in (False,True)}
        sub3={g:np.asarray([l3.parent[i] for i in xs]) for g,xs in outer[3,True].items()}
        sub4={g:np.asarray([l4.parent[i] for i in xs]) for g,xs in outer[4,True].items()}
        # Each lookup has the symmetric packed storage order of the Python oracle.
        pair={n:np.asarray([[l.pair(i,j) for j in range(l.n)] for i in range(l.n)])
              for n,l in ((3,l3),(4,l4))}
        path=OUT/f'moments-v5-{self.side}-score-k{K}.tsv.gz'
        counts=defaultdict(int)
        with gzip.open(path,'rt') as st:
            for line in st:
                tag,*fields=line.split()
                x,y,w,price=map(int,fields)
                b3,b4=outer[3,False][w],outer[4,False][w]
                if tag=='t':
                    mid=known[y]
                    true=self.vertical[:,dom['V'][mid],:][np.ix_(sub3[x],b3)]
                    room=int((self.bottom['V'][b3]-true).min())
                    if x in known:
                        centers=(known[x],mid)
                        p,f=(dom[a].get(centers) for a in ('P','F'))
                        d=self.bottom['P'][b4]+self.bottom['F'][b4]
                        if p is not None: d=d-self.pmax[p,b4]
                        if f is not None: d=d-self.fmax[f,b4]
                        room+=int(d.min())
                elif tag=='q':
                    i=dom['P'][known[y],(0,1)]
                    room=int((self.pmax[i,b4]-self.chain[i][np.ix_(sub4[x],b4)]).min())
                elif tag=='c':
                    mid=known[x]
                    i=dom['C'][mid]
                    a=pair[3][np.ix_(outer[3,False][y],b3)]
                    room=int(self.apex['C'][i]-self.pairs['C'][i,a].max())
                    a=pair[4][np.ix_(outer[4,False][y],b4)]
                    d=np.zeros(a.shape,dtype=np.int64)
                    i=dom['Y'].get(mid)
                    if i is not None: d+=self.apex['Y'][i]-self.pairs['Y'][i,a]
                    i=dom['T'].get(mid)
                    if i is not None: d+=(mid[1]-2)*(self.apex['T'][i]-self.pair_star[i,a])//3
                    room+=int(d.min())
                else: raise ValueError(tag)
                assert 0<=price<=room,(self.side,K,tag,x,y,w,price,room)
                counts[tag]+=1
        return counts

    def conditioned(self,K):
        c=self.cert
        w3=self.classes(c.levels[3],K,False)
        w4=self.classes(c.levels[4],K,False)
        coarse=self.classes(c.levels[4],K,False,True)
        path=OUT/f'moments-v5-{self.side}-edge-score-k{K}.tsv.gz'
        with gzip.open(path,'rt') as st:
            prices={(t,int(a),int(b)):int(p) for t,a,b,p in map(str.split,st)}
        visited=set()
        bg=self.bottom
        for d in range(2,74):
            for k in range(1,d):
                m=(d-1-k,k)
                g=group(d,k,K)
                above=d<=K
                vv=self.vmax[c.act['V'][m]] if m in c.act['V'] else 0
                pf={}
                for f,raw in (('P',self.pmax),('F',self.fmax)):
                    choices=[i for (a,h),i in c.act[f].items() if h==m and (not above or degree(a)>K)]
                    pf[f]=np.maximum(0,raw[choices].max(axis=0)) if choices else 0
                b4=bg['P']+bg['F']-pf['P']-pf['F']
                if above:
                    for b,cc in coarse.items():
                        key='t',g,b; visited.add(key)
                        assert prices.get(key,0)<=b4[cc].min(),(self.side,key,m)
                    continue
                b3=bg['V']-np.maximum(0,vv)
                for f in ('C','Y','T'):
                    i=c.act[f].get(m)
                    if i is None: continue
                    gap=self.apex[f][i]-self.one[f][i]
                    if f=='C': b3=b3+gap*(k-1)//2
                    elif f=='Y': b4=b4+gap*(k-1)//2
                    else: b4=b4+gap*((k-1)*(k-2))//6
                for b,cc in w3.items():
                    key='e',g,b; visited.add(key)
                    room=int(b3[cc].min())+int(b4[w4[b]].min())
                    assert 0<=prices.get(key,0)<=room,(self.side,key,m,room)
        assert set(prices)<=visited
        return len(prices)


if __name__=='__main__':
    resources()
    import numpy as np
    for side in ('upper','lower'):
        tables=Tables(side)
        print(f'{side}: all {len(tables.finite)} comparison maxima and all analytic tail coefficients agree.',flush=True)
        for K,M in ((4,32),(5,80),(7,144)):
            groups=tables.inputs(K,M)
            overlaps=tables.overlap(K)
            edges=tables.conditioned(K)
            affine=tables.tail_lines(K,M)
            print(f'{side} K={K}: {groups["v"]} singleton states, {groups["b"]} finite tail cells, '
                  f'{groups["e"]} edge rows, {sum(overlaps.values())} overlap prices and {edges} conditional '
                  f'prices passed the independent integer comparison. {affine} birth-conditioned score '
                  'lines passed finite cross-multiplication and analytic tail checks.',flush=True)
