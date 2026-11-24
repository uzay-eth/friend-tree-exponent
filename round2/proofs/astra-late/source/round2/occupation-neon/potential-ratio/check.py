#!/usr/bin/env python3
"""Exact direct comparison of the pinned D18 potential tables."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'attacks'/'neon1'))
from audit import Cert

PINS = ('9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
        '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41')
NAMES = ('f','g','hvtab','hctab','p4tab','yttab','ybtab','s3tab')


def pin(path):
    with path.open('rb') as f:
        return sha256_file(f)


def sha256_file(f):
    h=sha256()
    for b in iter(lambda:f.read(8<<20),b''):
        h.update(b)
    return h.hexdigest()


def differences(num,den):
    """A blockwise signed comparison, from two frozen payloads."""
    assert num and den
    def delta(a,b):
        # Frozen payload entries are range-checked in both multiplicands
        # before any machine integer multiplication or subtraction.
        bound=2**62//(abs(num)+abs(den))
        assert max(int(a.max()),int(b.max()))<=bound
        assert min(int(a.min()),int(b.min()))>=-bound
        return den*b-num*a
    return delta


def reductions(lo,hi,num,den):
    c=lo;delta=differences(num,den)
    def fold(name,axis):
        a,b=getattr(lo,name),getattr(hi,name)
        # The retained first axis is kept, if requested; other axes can
        # be folded incrementally over slabs bounded in RAM.
        out=[]
        step=max(1, (1<<19)//a[0].size)
        for i in range(0,len(a),step):
            d=delta(a[i:i+step],b[i:i+step])
            out.append((d.min(axis=axis),d.max(axis=axis)))
        x,y=(np.concatenate(v,axis=0) for v in zip(*out))
        return x,y
    def ends(name):
        a=getattr(c,name)
        x,y=fold(name,tuple(range(1,a.ndim-1)))
        return x.min(axis=0),y.max(axis=0)
    aa=np.ix_([c.state[x] for x in c.par],[c.state[x] for x in c.child])
    g=delta(lo.g[aa],hi.g[aa]);g=(g.min(axis=0),g.max(axis=0))
    hv=ends('hvtab');p4=ends('p4tab')
    hc=fold('hctab',1);s3=fold('s3tab',1)
    yt=fold('yttab',2);yt=(yt[0].min(axis=0),yt[1].max(axis=0))
    yb=fold('ybtab',(1,2)); factor=[x[1]-1 for x in c.ybu]
    assert max(max(abs(int(x.min())),abs(int(x.max())))*max(factor) for x in yb)<2**63
    yb=( (yb[0]*np.array(factor)[:,None]).min(axis=0),
         (yb[1]*np.array(factor)[:,None]).max(axis=0))
    f=delta(lo.f,hi.f)
    return f,g,hv,p4,yb,hc,yt,s3


def local(c,red,side,groups=False):
    """Extrema over leaves for every (core-child count, parent flag)."""
    f,g,hv,p4,yb,hc,yt,s3=red
    ext=min if side==0 else max
    np_ext=np.min if side==0 else np.max
    opt=lambda a: (np.minimum if side==0 else np.maximum)(a,0)
    bottom=[sum(map(int,x)) for x in zip(g[side],opt(hv[side]),opt(p4[side]),opt(yb[side]))]
    def controls(vertices,table,index,zero=False):
        return ext([int(table[index[v]]) for v in vertices if v in index]+([0] if zero else []))
    out={}
    arrays=(f,hc[side],yt[side],s3[side]); largest=max(
        max(abs(int(a.min())),abs(int(a.max()))) for a in arrays)
    assert largest*(1+2*comb(c.K,2)+max(comb(x[1],3) for x in c.s3))+max(map(abs,bottom))<2**63
    for h,k in product((0,1),range(c.KF+2)):
        start=max(0,2-k-h);end=c.LF if k<=c.KF else start
        tail=max(0,c.D-k-h)
        # All subsequent combined operations have a checked common safe
        # integer magnitude; recharges are done separately in Python Q.
        a=f[c.vi(0,k,h):c.vi(end,k,h)+1].copy()
        if h:
            a[tail:]+=ext(bottom[c.child[t]] for t in c.tail(0,k,h))
            for l in range(start,tail): a[l]+=bottom[c.child[l,k,h]]
        if 2<=k<=c.K:
            a[tail:]+=comb(k,2)*controls(c.tail(0,k,h),hc[side],c.hc,True)
            for l in range(start,tail): a[l]+=comb(k,2)*int(hc[side][c.hc[l,k,h]])
        for index,arr,m in ((c.yt,yt[side],2),(c.s3,s3[side],3)):
            for (l,j,bh),i in index.items():
                if (j,bh)==(k,h):a[l]+=comb(k,m)*int(arr[i])
        if not groups:
            out[k,h]=int(np_ext(a[start:]))
        else:
            for l in range(start,tail):out[(l,k,h)]=int(a[l])
            out[('tail',k,h)]=int(np_ext(a[max(start,tail):]))
    return out


def edge_messages(lo,hi,ratio,red):
    """Derive and verify a type-dependent zero-sum edge recharge."""
    f,g,hv,p4,yb,hc,yt,s3=red
    z=lambda a:tuple(np.zeros_like(v) for v in a)
    base=local(lo,(f,z(g),z(hv),p4,yb,hc,yt,s3),1,True)
    states=list(base)
    children=[j for j,s in enumerate(states) if s[2]]
    delta=differences(ratio.numerator,ratio.denominator)
    def labels(s):return lo.tail(0,s[1],s[2]) if s[0]=='tail' else (s,)
    mapped=[[lo.state[v] for v in labels(s)] for s in states]
    full=delta(lo.g,hi.g)
    G=np.array([[max(int(full[i,j]) for i in a for j in mapped[b]) for b in children] for a in mapped],dtype=np.int64)
    # Move P3 from its bottom up to its middle. For a given incoming
    # edge the middle now owns one charge per child. Keep it in that
    # edge's row, using the exact jointly indexed parent and middle.
    tab=lo.hvtab; hit=hi.hvtab
    vals=np.empty(tab.shape[:2],dtype=np.int64)
    for i in range(len(tab)):vals[i]=delta(tab[i],hit[i]).max(axis=1)
    H=np.empty_like(G)
    sl=list(lo.state)
    for i,a in enumerate(mapped):
        p=[lo.par[sl[u]] for u in a if sl[u] in lo.par]
        for j,b in enumerate(children):
            x=[lo.hv[sl[u]] for u in mapped[b] if sl[u] in lo.hv]
            y=[int(vals[u,v]) for u in p for v in x]
            if len(x)!=len(mapped[b]) or not y:y.append(0)
            term=states[b][1]*max(y)
            assert abs(term)<2**61
            H[i,j]=term
    assert max(abs(int(G.min())),abs(int(G.max())))+max(abs(int(H.min())),abs(int(H.max())))<2**61
    G+=H
    A=np.array(list(base.values()),dtype=np.int64)
    k=np.array([s[1] for s in states],dtype=np.int64)
    cap=[j for j,s in enumerate(states) if s[1]==lo.KF+1]
    child=np.array(children)
    # Each child's recharge subtracts the message of its actual parent,
    # while each parent pays once per child. Local recursion below is
    # max-plus; the final verification, not convergence, is the proof.
    def evaluate(m):
        mx=max(map(abs,map(int,m)))
        assert mx*int(k.max())+max(map(abs,map(int,A)))+max(abs(int(G.min())),abs(int(G.max())))<2**62
        rows=A+k*m
        return rows,(G+rows[child]).max(axis=1)
    m=np.full(len(states),-4*10**10*ratio.denominator,dtype=np.int64)
    for step in range(128):
        rows,want=evaluate(m)
        root=max(int(v) for v,s in zip(rows,states) if not s[2])
        if root<=0 and max(m[cap])<=0 and np.all(want<=m):
            print('  PASS type-dependent upper recharge',ratio,'iterations',step,'root reserve',-root,
                  'edge reserve',int((m-want).min()),'overflow slopes',list(map(int,m[cap])),flush=True)
            return True
        m=np.maximum(m,want)
        if max(map(abs,map(int,m)))>=10**15:break
    print('  no type-dependent recharge found',ratio,flush=True)
    return False


def brute(lo,hi,endpoints):
    """Enumerate actual rooted trees, using the exact full motif products."""
    fixtures=[]
    # All chronological trees from T3; the founder stays vertex zero.
    for n in range(3,9):
        for parents in product(*[range(j) for j in range(3,n)]):
            p=(-1,0,0)+parents; deg=[0]*n
            for j in range(1,n):deg[j]+=1;deg[p[j]]+=1
            ids={j:i for i,j in enumerate(j for j in range(n) if deg[j]>=2)}
            ch=[[] for _ in ids]
            for j,i in ids.items():
                if j:ch[ids[p[j]]].append(i)
            fixtures.append((ch,[(deg[j]-len(ch[i])-bool(j),len(ch[i]),int(bool(j))) for j,i in ids.items()]))
    def add(par,leaves):
        ch=[[] for _ in par]
        for j,p in enumerate(par[1:],1):ch[p].append(j)
        assert all(l+len(ch[j])+bool(j)>=2 for j,l in enumerate(leaves))
        fixtures.append((ch,[(l,len(ch[j]),int(bool(j))) for j,l in enumerate(leaves)]))
    for d in (2,17,18,19,200,14001,10**30):
        add([-1],[d]);add([-1,0,1,2,3],[d,0,d,1,d])
    for k in (2,3,4,17,18,180,181,1000):
        add([-1]+[0]*k,[0]+[1]*k)
        if k<=18:
            add([-1]+[0]*k+list(range(1,k+1)),[0]+[1]*k+[30]*k)
    add([-1,0,1,2,2,1,0,0],[0,0,0,1,1,1,1,1])
    low=high=None;where=None;used=set()
    for ch,st in fixtures:
        a=lo.potential(ch,st,parts=True); b=hi.potential(ch,st,parts=True)
        used.update(i for i,(x,y) in enumerate(zip(a,b)) if x or y)
        av,bv=sum(a),sum(b)
        assert av>0 and bv>0
        q=bv/av
        assert endpoints[0]<=q<=endpoints[1]
        low=q if low is None else min(low,q)
        if high is None or q>high:
            high=q;where=(len(st),st[0],max(map(sum,st)))
    assert used==set(range(8)),used
    print('PASS whole-tree brute force',len(fixtures),'trees, all 8 motif families; observed ratio',float(low),float(high),flush=True)
    print('  max fixture (cores, founder state, max degree)',where,flush=True)


def feasible(rows,side):
    """Exact interval for t with sign*(row + t*(k-h)) >= 0."""
    sign=1 if side==0 else -1
    low=high=None
    for (k,h),v in rows.items():
        p=sign*(k-h);q=sign*v
        if not p:
            if q<0:return None
            continue
        t=Q(-q,p)
        if p>0:low=t if low is None else max(low,t)
        else:high=t if high is None else min(high,t)
    return (low,high) if low<=high else None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lo',type=Path);p.add_argument('hi',type=Path)
    p.add_argument('ratios',nargs='*',type=Q,default=[])
    p.add_argument('--verify-sha',action='store_true')
    p.add_argument('--diagnostic',action='store_true')
    p.add_argument('--edge-messages',action='store_true')
    a=p.parse_args()
    paths=(a.lo,a.hi)
    stamps=lambda: [(s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns) for s in (x.stat() for x in paths)]
    stamp=stamps()
    if a.verify_sha or not a.ratios:
        assert tuple(map(pin,paths))==PINS
        print('PASS both original SHA256 pins',flush=True)
    lo,hi=map(Cert,paths)
    assert (lo.mode,hi.mode)==(2,1) and (lo.lam,hi.lam)==(Q(565673,10**6),Q(22631,40000))
    for key in ('S','D','K','LF','KF','state','par','child','hv','hc','p4','yt','ybu','ybw','s3'):
        assert getattr(lo,key)==getattr(hi,key)
    endpoints=a.ratios or [Q('0.999'),Q('1.12')]
    for ratio in endpoints:
        side=0 if ratio<=1 else 1
        red=reductions(lo,hi,ratio.numerator,ratio.denominator)
        rows=local(lo,red,side)
        band=feasible(rows,side)
        print('ratio',ratio,'side',side,'raw extrema',min(rows.items(),key=lambda x:x[1]),max(rows.items(),key=lambda x:x[1]),'recharge',band,flush=True)
        if a.diagnostic:
            print('  base',[(s,rows[s]) for s in ((0,0),(0,1),(1,1),(181,0),(181,1))],flush=True)
            for cond,ex in ((lambda k,h:k>h,max),(lambda k,h:k<h,min)):
                print('  slope',ex((Q(v,k-h),k,h) for (k,h),v in rows.items() if cond(k,h)),flush=True)
            for s in ((0,1,1),(1,1,1),(0,17,1),(1,0,1)):
                for t,_ in lo.embed(s):
                    i=lo.child[t];f,g,hv,p4,yb,hc,yt,s3=red
                    print('  terms',s,t,[int(z) for z in (f[lo.vi(*s)],g[side][i],hv[side][i],p4[side][i],yb[side][i], hc[side][lo.hc[t]]*comb(s[1],2) if t in lo.hc else 0)],flush=True)
        if not a.ratios and side==0:
            assert band is not None
            t=Q(50*10**9)
            assert band[0]<=t<=band[1]
            slack=min((1 if side==0 else -1)*(v+t*(k-h)) for (k,h),v in rows.items())
            print('  certified recharge',t,'minimum signed row reserve',slack,flush=True)
            assert slack>0 and min(rows.values())<0
        if side==1 and (not a.ratios or a.edge_messages):
            ok=edge_messages(lo,hi,ratio,red)
            if not a.ratios:assert ok and band is None
    if not a.ratios:
        brute(lo,hi,endpoints)
    assert stamps()==stamp


if __name__=='__main__':main()
