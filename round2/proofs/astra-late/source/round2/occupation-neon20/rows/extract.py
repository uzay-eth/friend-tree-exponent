#!/usr/bin/env python3
"""Authenticated sparse D20 affine charges and grouped potential comparison."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
import json
from math import comb, prod
import mmap
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROUND2 = HERE.parents[1]
AUDIT_PATH = ROUND2/'attacks/neon1/audit.py'
assert hashlib.sha256(AUDIT_PATH.read_bytes()).hexdigest() == 'd548a9aafa899ecab2b00a7f6ad0a4e6f828600088780f0a15f08927f17f0bfe'
_spec = spec_from_file_location('neon_layout_audit', AUDIT_PATH)
audit = module_from_spec(_spec)
_spec.loader.exec_module(audit)

COMMIT = '51f4cfadac675fdf1cef84d44935eebb5456aabb'
SOURCE = {'round2/proof-b3f-motif4.md':'22c61b0f858fc22d9f1eddae24275676df0ea77f',
          'round2/cert/cert_b3f_format.h':'eebd8b61c0a75bd85dc6c55854346883e6c96720',
          'round2/cert/cert_b3f_verify_gmp.c':'9df09928b75ea74d27e92ef885c7010331b84aeb',
          'round2/outputs/cert_b3f_20_lo.log':'8067ff1c9948902ab49b85822c4d9c5ea3986480',
          'round2/outputs/cert_b3f_20_hi.log':'8e0c496aeae4a7b8ee470caf555bbb5dfa91ad7c'}
HASH = {'lower':'eb93985b2c90a516b65e85934f4ac847779f42b3fdbfcd2af4dc07fd30c37469',
        'upper':'55dd6015eca85ca2fba6267d9444ec8807554c456999803d8f9abcc695c090e2'}
ROOT = {'lower':(113472636536,184219674916,146041684821),
        'upper':(115454594604,200854054323,146044513514)}  # 4Sc0, SC0, SΦ(T3)
LO, HI = Q('0.565680'), Q('0.565768')
NB = (1,0,1)
SPEC = {'f':(-1,), 'g':(-1,0), 'hV':(-1,0,1), 'hC':(-1,0,0),
        'P4':(-1,0,1,2),'Yt':(-1,0,1,1),'Yb':(-1,0,0,2),'S3':(-1,0,0,0)}
SIZE = 11190881056


def pinned_sources(path):
    import subprocess
    for name, expected in SOURCE.items():
        out = subprocess.check_output(['git','-C',str(path),'rev-parse',f'{COMMIT}:{name}'], text=True).strip()
        if out != expected:
            raise ValueError(f'wrong source blob for {name}')


def stamp(p):
    s = p.stat()
    return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)


class Cert(audit.Cert):
    def __init__(self,path,side,*,verify=True):
        path = Path(path)
        self.path,self.side = path,side
        mode,p = {'lower':(2,565680), 'upper':(1,565768)}[side]
        self.stamp = stamp(path)
        if self.stamp[2] != SIZE:
            raise ValueError('wrong D20 size')
        with path.open('rb') as stream:
            header = audit.HEADER.unpack(stream.read(audit.HEADER.size))
            expected = (0x4E454F4E42334648,mode,1,p,10**6,10**11,10**9,16000,200,20,19,12,9,11,8,12,9,10,7)
            if header != expected:
                raise ValueError(f'unexpected D20 {side} header: {header}')
            if verify:
                stream.seek(0)
                if hashlib.file_digest(stream,'sha256').hexdigest() != HASH[side]:
                    raise ValueError(f'{side} SHA-256 mismatch')
        super().__init__(path)
        self.check_unchanged()
        self.authenticated = verify
        self.sign = 1 if mode == 2 else -1
        self.c0,self.C0,self.phi_T3 = (Q(v,self.S*(4 if i==0 else 1)) for i,v in enumerate(ROOT[side]))
        assert self.potentials == 1398838433 and self.stored == 1398860120
        assert Q(int(self.f[self.vi(2,0,0)]), self.S) == self.phi_T3

    def check_unchanged(self):
        if self.stamp != stamp(self.path):
            raise ValueError('payload replaced or mutated during replay')


def children(p):
    ch=[[] for _ in p]
    for j,i in enumerate(p[1:],1): ch[i].append(j)
    return ch


@lru_cache(maxsize=10000)
def embeddings(p):
    ch=children(p); out=[]
    for a,kids in enumerate(ch):
        out.append(('f',(a,)))
        for b in kids:
            out.append(('g',(a,b)))
            for d in ch[b]:
                out.append(('hV',(a,b,d)))
                out.extend(('P4',(a,b,d,e)) for e in ch[d])
                out.extend(('Yb',(a,w,b,d)) for w in kids if w!=b)
            out.extend(('Yt',(a,b,d,e)) for d,e in combinations(ch[b],2))
        out.extend(('hC',(a,*ids)) for ids in combinations(kids,2))
        out.extend(('S3',(a,*ids)) for ids in combinations(kids,3))
    return tuple(out)


def pick(x,ids): return tuple(x[i] for i in ids)
def replace(x,i,s): return x[:i]+(s,)+x[i+1:]
def gain(s): return s[0]+1,s[1],s[2]
def bud(s): return s[0]-1,s[1]+1,s[2]


class Rows:
    def __init__(self,c): self.c=c

    def live(self,f,x):
        c=self.c
        return (f in ('f','g') or f=='hV' and 1<=x[1][1]<=c.K
                or f=='hC' and 2<=x[0][1]<=c.K
                or f=='P4' and x[1] in c.p4 and x[2] in c.p4
                or f=='Yt' and x[1] in c.yt
                or f=='Yb' and x[0] in c.ybu and x[2] in c.ybw
                or f=='S3' and x[0] in c.s3)

    @lru_cache(maxsize=40000)
    def value(self,f,x):
        """Physical potential of one motif in ticks, with no charge transfers."""
        c=self.c
        if f=='f': return Q(int(c.f[c.vi(*x[0])]))
        if not self.live(f,x): return Q(0)
        # Restricted slots must use their actual physical census. Only the
        # other slots are Bernstein-embedded.
        hold={'P4':(1,2),'Yt':(1,),'Yb':(0,2),'S3':(0,)}.get(f,())
        total=Q(0)
        for terms in product(*(c.embed(s) if i not in hold else [(s,Q(1))] for i,s in enumerate(x))):
            y,weights=zip(*terms); a,b,*rest=y
            match f:
                case 'g': u=c.g[c.state[a],c.state[b]]
                case 'hV': u=c.hvtab[c.par[a],c.hv[b],c.child[rest[0]]] if b in c.hv else 0
                case 'hC': u=c.hctab[c.hc[a],c.tp(b,rest[0])] if a in c.hc else 0
                case 'P4': u=c.p4tab[c.par[a],c.p4[b],c.p4[rest[0]],c.child[rest[1]]]
                case 'Yt': u=c.yttab[c.par[a],c.yt[b],c.tp(*rest)]
                case 'Yb': u=c.ybtab[c.ybu[a],c.ybw[rest[0]],c.child[b],c.child[rest[1]]]
                case 'S3': u=c.s3tab[c.s3[a],c.tt(b,*rest)]
                case _: raise ValueError(f)
            total+=int(u)*prod(weights)
        return total

    def routed(self,f,x,i,parent):
        c=self.c; k=x[i][1]
        if f=='f': return sum(x[i])<c.D or k<=c.K
        if f=='g': return k<=c.K
        if f=='hV':
            if i==0: return x[0] in c.p4 and x[1] in c.p4 if parent else x[0] in c.ybu and x[1] in c.ybw
            return x[1] in c.yt if i==1 else x[1] in c.p4 and x[2] in c.p4
        if f=='hC': return x[0] in (c.yt if parent else c.s3) if i==0 else x[0] in c.ybu and x[i] in c.ybw
        return False

    def center(self,s,parent):
        c=self.c; _,k,h=s
        if s in c.inner_set: return Q(int((c.rpar if parent else c.rch)[c.state[s]]),c.Sr)
        return Q(int((c.trpar if parent else c.trch)[(1-h)*(c.K+2)+min(k,c.K+1)]),c.Sr)

    def transfer(self,f,x,i):
        c=self.c; _,k,h=x[i]; d=sum(x[i]); j=(1-h)*(c.K+2)+k
        if d<c.D or f not in ('f','g'):
            return self.value(f,replace(x,i,gain(x[i])))-self.value(f,x)
        if f=='f': return Q(int(c.f2[j])*d+int(c.f3[j]),d**3)
        t=Q(c.D,d); out=Q(0); tail=('T1',k,h); inf=('inf',)
        for v,w in c.embed(x[1-i]):
            if i:
                a=Q(int(c.g[c.state[v],c.state[inf]])-int(c.g[c.state[v],c.state[tail]])); y=c.yw[c.par[v],k]
            else:
                a=Q(int(c.g[c.state[inf],c.state[v]])-int(c.g[c.state[tail],c.state[v]])); y=c.yp[j,c.child[v]]
            out+=w*(3*t*t*(1-t)*a/c.D+t**3*int(y))
        return out

    @lru_cache(maxsize=35000)
    def affine(self,f,x):
        """Return (signed intercept, (slot,is_parent,count,signed slope), ...)."""
        c=self.c; p=SPEC[f]; n=len(p)
        if not self.live(f,x): return Q(0),()
        old=self.value(f,x); a=-c.lam*old; ports=[]
        for i,(l,k,h) in enumerate(x):
            dp=self.value(f,replace(x,i,gain(x[i])))-old
            a+=(l+sum((Q(1,sum(x[j])) for j in range(n) if p[j]==i or p[i]==j),Q(0)))*dp
            if l:
                new=replace(x,i,bud(x[i])); y=new+(NB,)
                born=sum(self.value(g,pick(y,ids)) for g,ids in embeddings(p+(i,)) if len(ids)==n+1)
                if f=='f': born+=self.value('f',(NB,))
                a+=Q(l,l+k+h)*(self.value(f,new)-old+born)
            for parent,count in ((True,h-int(i>0)),(False,k-p.count(i))):
                if count:
                    shifted=self.transfer(f,x,i) if self.routed(f,x,i,parent) else Q(0)
                    if shifted: a+=count*self.center(x[i],parent)*shifted
                    ports.append((i,parent,count,c.sign*(dp-shifted)/c.S))
        for g,ids in embeddings(p):
            if len(ids)!=n-1: continue
            (donor,)=set(range(n))-set(ids)
            (j,)=(j for j in ids if p[j]==donor or p[donor]==j)
            y=pick(x,ids); i=ids.index(j); parent=p[j]==donor
            if self.routed(g,y,i,parent):
                a+=(Q(1,sum(x[donor]))-self.center(x[j],parent))*self.transfer(g,y,i)
        return c.sign*a/c.S,tuple(ports)

    def first_body(self,state):
        a,ports=self.affine('f',(state,)); slopes={s for _,_,_,s in ports}
        assert len(slopes)<=1
        return a,next(iter(slopes),Q(0))


def first_body(c,state): return Rows(c).first_body(state)


def stock(c,k,h=1):
    """Signed corners at L=LF and a uniform envelope for L>=LF+1."""
    e=Rows(c); L=c.LF; d=L+k+h; s=(L,k,h); old=e.value('f',(s,)); nb=e.value('f',(NB,))
    base=Q(0); knots=(Q(0),)
    if k<=c.K:
        i=(1-h)*(c.K+2)+k
        a,b=Q(int(c.f2[i]),d*d),Q(int(c.f3[i]),d**3)
        knots=(Q(0),a/3,a+b)
        base=h*e.center(s,True)+k*e.center(s,False)
    adjustment=[r*p for r,p in product((base,base-Q(k+h,2)),knots)]
    y=bud(s)
    at=Q(L,d)*(e.value('f',(y,))+nb-old+e.value('g',(y,NB)))-c.lam*old
    controls=[Q(int(c.g[c.state[z],c.state[NB]])) for z in c.tail(*y)]
    t=Q(c.D,d+1); restricted=[controls[0]]
    while len(controls)>1:
        controls=[(1-t)*a+t*b for a,b in zip(controls,controls[1:])]
        restricted.append(controls[0])
    f_new=e.value('f',((L,k+1,h),))+nb-old
    return (tuple(c.sign*(at+a)/c.S for a in adjustment),
            tuple(c.sign*(-c.lam*old+b*(f_new+g)+a)/c.S
                  for b,g,a in product((Q(L+1,d+1),Q(1)),restricted,adjustment)))


def bands(c,kmax=4):
    """Nonfounder floors for all legal leaf stocks at each k<=kmax."""
    e=Rows(c); out={}
    for k in range(kmax+1):
        best=None
        for l in range(max(0,1-k),c.LF):
            a,b=e.first_body((l,k,1)); v=a+min(0,b*Q(k+1,2))
            assert v>=0,(l,k,v)
            if best is None or (v,l)<best: best=(v,l)
        at,above=map(min,stock(c,k)); assert min(at,above)>=0
        out[k]={'finite_min':best[0],'finite_argmin_L':best[1],'at_LF':at,'above_LF':above,
                'floor':min(best[0],at,above)}
        print(c.side,'band',k,'floor',out[k]['floor'],flush=True)
    return out


def physical_deficit(c,q,t,degree):
    """Discounted lower bound on nonnegative recharge deficit, per core."""
    if q.c is not c: raise ValueError('charge is from another witness')
    cap=q.cap(t); v={}
    for d in range(degree,1,-1):
        for k in range(d-1,-1,-1):
            l=d-k-1; b=Q(l,d)
            r=Q(cap-q.exact(l,k,1)-t*(k-1),c.S); assert r>=0
            rhs=r+b*v.get((d,k+1),0); adv=v.get((d+1,k),0)
            rates=(Q(l),l+Q(k+1,2)); scale=10**21
            w=Q(min((rhs+g*adv)/(HI+b+g) for g in rates)*scale//1,scale)
            assert all(w*(HI+b+g)<=rhs+g*adv for g in rates)
            v[d,k]=w
    return LO*v[2,0]


def recharge(c,q,degree):
    """Optimize a small exact recharge grid with fresh D20 physical deficits."""
    best=None
    for j in range(51):
        t=-Q(j*c.S,1000); cap=q.cap(t); deficit=physical_deficit(c,q,t,degree)
        candidate=(Q(cap,c.S)-deficit,t,cap,deficit)
        if best is None or candidate<best: best=candidate
    value,t,cap,deficit=best
    return dict(degree=degree,ceiling=value,recharge_ticks=t,cap_ticks=cap,stationary_deficit=deficit)


class Charge:
    """All-eight upper envelope (ticks), with optional zero-sum recharge."""
    def __init__(self,c):
        self.c=c
        def reduce(a,axes):
            if 0 not in axes: raise ValueError('streamed reduction requires axis 0')
            step=max(1,(16<<20)//(a[0].size*a.itemsize))
            answer=None
            for i in range(0,len(a),step):
                part=a[i:i+step].max(axis=axes)
                answer=part if answer is None else np.maximum(answer,part)
                backing=getattr(a,'_mmap',None)
                if backing is not None and hasattr(backing,'madvise'):
                    backing.madvise(mmap.MADV_DONTNEED)
            return answer

        def bottom(a): return reduce(a.reshape(-1,a.shape[-1]),(0,))

        def first_axis_max(a):
            out=[]
            for block in a:
                out.append(int(block.max()))
                a._mmap.madvise(mmap.MADV_DONTNEED)
            return np.array(out,dtype=np.int64)

        g=c.g[np.ix_([c.state[x] for x in c.par],[c.state[x] for x in c.child])]
        hm=bottom(c.hvtab); pm=bottom(c.p4tab)
        self.hc=first_axis_max(c.hctab); self.yt=reduce(c.yttab,(0,2)); self.s3=first_axis_max(c.s3tab)
        bm=[bottom(c.ybtab[i]) for i in range(len(c.ybu))]
        largest=max(max(abs(int(a.min())),abs(int(a.max()))) for a in
                    (c.f,bottom(g),hm,pm,self.hc,self.yt,self.s3,*bm))
        assert largest*(4+max(s[1] for s in c.ybu)+2*comb(c.K,2)+comb(c.K,3))<2**63
        bend=np.maximum.reduce([(s[1]-1)*bm[i] for s,i in c.ybu.items()])
        self.edge=bottom(g)+np.maximum(0,hm)+np.maximum(0,pm)+np.maximum(0,bend)
        self.by_k={}
        for h,k in product((0,1),range(c.KF+2)):
            start=max(0,2-k-h); end=max(0,c.D-k-h)
            arr=np.array(c.f[c.vi(0,k,h):c.vi(c.LF,k,h)+1])
            if h:
                arr[end:]+=max(self.edge[c.child[t]] for t in c.tail(0,k,h))
                for l in range(start,end): arr[l]+=self.edge[c.child[l,k,h]]
            if 2<=k<=c.K:
                arr[end:]+=comb(k,2)*max(self.hc[c.hc[t]] if t in c.hc else 0 for t in c.tail(0,k,h))
                for l in range(start,end): arr[l]+=comb(k,2)*self.hc[c.hc[l,k,h]]
            for index,v,m in ((c.yt,self.yt,2),(c.s3,self.s3,3)):
                for (l,j,bh),i in index.items():
                    if (j,bh)==(k,h): arr[l]+=comb(k,m)*v[i]
            self.by_k[h,k]=int(arr[start:].max())
        assert self.cap(0)<=ROOT[c.side][1]

    @lru_cache(maxsize=8000)
    def exact(self,l,k,h):
        c=self.c; s=(l,k,h); out=Q(int(c.f[c.vi(*s)]))
        for t,w in c.embed(s):
            if h: out+=w*int(self.edge[c.child[t]])
            if t in c.hc: out+=comb(k,2)*w*int(self.hc[c.hc[t]])
        if s in c.yt: out+=comb(k,2)*int(self.yt[c.yt[s]])
        if s in c.s3: out+=comb(k,3)*int(self.s3[c.s3[s]])
        assert out<=self.by_k[h,min(k,c.KF+1)]
        return out

    def cap(self,t):
        if t>0: raise ValueError('positive recharge has no uniform overflow bound')
        return max(a+t*(k-h) for (h,k),a in self.by_k.items())


def finite_tables(c,degree):
    e=Rows(c); out=[]
    for d in range(2,degree+1):
        for k in range(d):
            s=(d-k-1,k,1); a,b=e.first_body(s); z=a+b*Q(k+1,2)
            assert min(a,z)>=0,(s,a,z)
            out.append([*s,a,b,a,z])
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('path',type=Path); p.add_argument('side',choices=HASH)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--degree',type=int,default=64)
    p.add_argument('--kmax',type=int,default=4);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--charge',action='store_true');p.add_argument('--recharge',action='store_true');a=p.parse_args()
    if not __debug__: p.error('assertions required')
    if not (2<=a.degree<=201 and 0<=a.kmax<=19): p.error('degree or band cut outside verified export scope')
    pinned_sources(a.source); c=Cert(a.path,a.side); print('authenticated',a.side,HASH[a.side],flush=True)
    data={'source_commit':COMMIT,'source_blobs':SOURCE,'payload_sha256':HASH[a.side],'side':a.side,
          'lambda':c.lam,'root_ticks':ROOT[a.side],'parameters':{'D0':c.D,'KT':c.K,'LF':c.LF,'KF':c.KF,'S':c.S,'Sr':c.Sr},
          'row_columns':['L','k','h','intercept','total_receive_slope','value_at_0','value_at_(k+h)/2'],
          'bands':bands(c,a.kmax),'f_rows':finite_tables(c,a.degree)}
    if a.charge or a.recharge:
        q=Charge(c); data['charge_by_h_k_ticks']=[[h,k,v] for (h,k),v in sorted(q.by_k.items())]
        print(c.side,'unrecharged stock charge ticks',q.cap(0),flush=True)
        if a.recharge:
            data['recharge']=recharge(c,q,a.degree)
            print(c.side,'physical recharge',data['recharge'],flush=True)
    c.check_unchanged()
    with a.output.open('x') as out: json.dump(data,out,default=str,separators=(',',':'),sort_keys=True);out.write('\n')
    print(c.side,len(data['f_rows']),'finite rows;',len(data['bands']),'all-stock floors',flush=True)


if __name__=='__main__':main()
