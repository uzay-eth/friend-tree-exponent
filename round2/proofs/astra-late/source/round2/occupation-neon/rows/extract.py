"""Sparse exact replay of charged f and inner–inner g rows for pinned D18."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import sys
import json
import hashlib

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'attacks'/'neon1'))
from audit import Cert
from itertools import product


def f(c,s): return int(c.f[c.vi(*s)])
def g(c,p,w): return int(c.g[c.state[p],c.state[w]])
NB=(1,0,1)
HASH={'lower':'9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
      'upper':'1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41'}


def recoil(c,s,role):
    return F(int((c.rpar if role=='p' else c.rch)[c.state[s]]),c.Sr)


def frow(c,s,rho=None):
    L,k,h=s; d=sum(s); assert 2<=d and L<c.LF and k<=c.KF
    v=f(c,s); up=f(c,(L+1,k,h))-v
    q=L*up-c.lam*v
    if L:
        b=(L-1,k+1,h)
        q+=F(L,d)*(f(c,b)+f(c,NB)-v+sum(w*g(c,x,NB) for x,w in c.embed(b)))
    if d<c.D:
        q+=(h*recoil(c,s,'p')+k*recoil(c,s,'c'))*up
        return q/c.S
    p=0; base=0
    if k<=c.K:
        i=(1-h)*(c.K+2)+k
        p=F(int(c.f2[i]),d*d)+F(int(c.f3[i]),d**3)
        base=F(h*int(c.trpar[i])+k*int(c.trch[i]),c.Sr)
    if rho is None:
        return tuple((q+base*p+r*(up-p))/c.S for r in (0,F(k+h,2)))
    return (q+base*p+rho*(up-p))/c.S


def capped(c,k,h):
    # Verifier's independent corners at L=LF and L>=LF+1, with exact
    # de Casteljau restriction of the newborn g edge to the latter interval.
    L=c.LF; d=L+k+h; v=f(c,(L,k,h)); base=F(0); a=b=F(0)
    if k<=c.K:
        i=(1-h)*(c.K+2)+k
        base=F(h*int(c.trpar[i])+k*int(c.trch[i]),c.Sr)
        a,b=F(int(c.f2[i]),d*d),F(int(c.f3[i]),d**3)
    extra=[r*x for r,x in product((base,base-F(k+h,2)),(0,a/3,a+b))]
    nb=f(c,NB)
    bud=(L-1,k+1,h)
    gex=sum(w*g(c,s,NB) for s,w in c.embed(bud))
    exact=[-c.lam*v+F(L,d)*(f(c,bud)+nb-v+gex)+x for x in extra]
    cpoints=[F(g(c,s,NB)) for s in c.tail(*bud)]
    t=F(c.D,d+1); restricted=[cpoints[0]]
    while len(cpoints)>1:
        cpoints=[(1-t)*x+t*y for x,y in zip(cpoints,cpoints[1:])]
        restricted.append(cpoints[0])
    cap=[-c.lam*v+prob*(f(c,(L,k+1,h))+nb-v+y)+x
         for prob,y,x in product((F(L+1,d+1),F(1)),restricted,extra)]
    return tuple(x/c.S for x in exact),tuple(x/c.S for x in cap)


def grow(c,s):
    L,k,h=s
    return (L+1,k,h) if sum(s)<c.D-1 else (f'T3',k,h)


def grow_f(c,s):
    L,k,h=s
    return f(c,(L+1,k,h))-f(c,s)


def grow_core(c,s):
    L,k,h=s
    return (L-1,k+1,h)


def grow_edge_row(c,p,w):
    # Both endpoints exact inner states: no physical external credit left over.
    assert p in c.inner_set and w in c.inner_set and p[1]>=1 and w[2]==1
    v=g(c,p,w); L,k,h=p; l,m,_=w
    q=-c.lam*v
    q+=(F(1,sum(w))-recoil(c,p,'c'))*grow_f(c,p)
    q+=(F(1,sum(p))-recoil(c,w,'p'))*grow_f(c,w)
    q+=(L+F(1,sum(w))+h*recoil(c,p,'p')+(k-1)*recoil(c,p,'c'))*(g(c,grow(c,p),w)-v)
    q+=(l+F(1,sum(p))+m*recoil(c,w,'c'))*(g(c,p,grow(c,w))-v)
    if L:
        b=grow_core(c,p); x=g(c,b,w)
        if b in c.hc: x+=int(c.hctab[c.hc[b],c.tp(w,NB)])
        q+=F(L,sum(p))*(x-v)
    if l:
        b=grow_core(c,w); x=g(c,p,b)
        if b in c.hv: x+=int(c.hvtab[c.par[p],c.hv[b],c.child[NB]])
        q+=F(l,sum(w))*(x-v)
    return q/c.S


def scan(c):
    sign=1 if c.mode==2 else -1
    for size in (4,7,17):
        inner=sorted([(sign*frow(c,s),s) for s in c.inner if sum(s)<=size and s[2]==1],reverse=True)
        print('f inner',size,[(s,str(v),float(v)) for v,s in inner[:8]])
    for size in (4,7,17):
        edges=sorted([(sign*grow_edge_row(c,p,w),(p,w)) for p in c.inner if p[1]>=1 and p[2] and sum(p)<=size
                      for w in c.inner if w[2] and sum(w)<=size],reverse=True)
        print('g inner inner',size,[(s,str(v),float(v)) for v,s in edges[:8]])
    # Tail near boundary with exact tightened interval [0,(k+h)/2].
    for size in (24,45,120):
        tail=sorted([(min(sign*x for x in frow(c,(L,k,1))),(L,k,1)) for L in range(size)
                     for k in range(min(4,size-L)) if c.D<=L+k+1<=size],reverse=True)
        print('f tail',size,[(s,str(v),float(v)) for v,s in tail[:8]])


def bands(c):
    sign=1 if c.mode==2 else -1
    for k in range(5):
        finite=[]
        for L in range(max(0,1-k),c.LF):
            a=frow(c,(L,k,1)); a=a if isinstance(a,tuple) else (a,)
            finite.append((min(sign*x for x in a),L))
        caps=[min(sign*x for x in a) for a in capped(c,k,1)]
        assert min(*[a for a,_ in finite],*caps)>=0
        print('band',c.mode,k,'finite minimum',min(finite),'caps',list(zip(caps,map(float,caps))),
              'low',[(L,str(x)) for x,L in sorted(finite)[:3]],flush=True)


def buffer(c):
    # Neighbor-free lower envelope for the signed f charge. Exact tail endpoint
    # values are included; possible correlations between endpoints not required.
    sign=1 if c.mode==2 else -1
    def inputs(s):
        L,k,h=s; d=sum(s)
        b=(L-1,k+1,h)
        v={'state':s,'f':f(c,s),'f_gain':f(c,(L+1,k,h)),
           'f_bud':f(c,b) if L else 0,'f_newborn':f(c,NB),
           'g_bud_newborn':[[int(g(c,x,NB)),str(w)] for x,w in c.embed(b)] if L else []}
        if d<c.D:
            v.update(mode='inner',center=[str(recoil(c,s,'p')),str(recoil(c,s,'c'))])
        elif k<=c.K:
            idx=(1-h)*(c.K+2)+k
            v.update(mode='tail',center=[str(F(int(c.trpar[idx]),c.Sr)),str(F(int(c.trch[idx]),c.Sr))],
                     c2=int(c.f2[idx]),c3=int(c.f3[idx]))
        else:
            v.update(mode='tail',center=['0','0'],c2=0,c3=0)
        a=frow(c,s); corners=a if isinstance(a,tuple) else (a,)
        v['signed_corners']=[str(sign*x) for x in corners]
        assert min(map(F,v['signed_corners']))>=0
        return v
    states=[(d-k-1,k,1) for d in range(2,65) for k in range(d)]
    return {'mode':'lower' if sign==1 else 'upper','sign':sign,'lambda':str(c.lam),'S':c.S,'Sr':c.Sr,
            'D0':c.D,'KT':c.K,'LF':c.LF,'KF':c.KF,'rows':[inputs(s) for s in states]}


def portable(c):
    # Raw integer slab needed for k=0..4 at all L, including one extra k for budding.
    ks=range(5)
    return dict(lam=str(c.lam),sign=1 if c.mode==2 else -1,S=c.S,Sr=c.Sr,D0=c.D,LF=c.LF,
                f=[[f(c,(L,k,1)) for L in range(c.LF+1)] for k in range(6)],
                r=[[[int(c.rpar[c.state[L,k,1]]),int(c.rch[c.state[L,k,1]])]
                    for L in range(max(0,1-k),c.D-1-k)] for k in ks],
                tr=[[int(a[k]) for a in (c.trpar,c.trch,c.f2,c.f3)] for k in ks],
                e=[[g(c,(L-1,k+1,1),NB) for L in range(1,c.D-1-k)] for k in ks],
                t=[[g(c,s,NB) for s in c.tail(0,k+1,1)] for k in ks])


def edges_buffer(c):
    rows=[]
    for p in c.inner:
        if p[1]<1: continue
        for w in c.inner:
            if not w[2]: continue
            value=g(c,p,w); pb=wb=0; hC=hV=0
            if p[0]:
                b=grow_core(c,p); pb=g(c,b,w)
                if b in c.hc: hC=int(c.hctab[c.hc[b],c.tp(w,NB)])
            if w[0]:
                b=grow_core(c,w); wb=g(c,p,b)
                if b in c.hv: hV=int(c.hvtab[c.par[p],c.hv[b],c.child[NB]])
            rr=lambda s: [int(a[c.state[s]]) for a in (c.rpar,c.rch)]
            rows.append([*p,*w,value,g(c,grow(c,p),w),g(c,p,grow(c,w)),pb,wb,hC,hV,
                         grow_f(c,p),grow_f(c,w),*rr(p),*rr(w)])
    return dict(lam=str(c.lam),sign=1 if c.mode==2 else -1,S=c.S,Sr=c.Sr,rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('cert',nargs='+',type=Path)
    p.add_argument('--buffer',action='store_true')
    p.add_argument('--bands',action='store_true')
    p.add_argument('--portable',action='store_true')
    p.add_argument('--edges',action='store_true')
    p.add_argument('--verify-sha',action='store_true')
    args=p.parse_args(); out={}
    for fn in args.cert:
        c=Cert(fn); assert c.D==18 and c.K==17 and c.S==10**11
        mode='lower' if c.mode==2 else 'upper'
        assert c.lam=={'lower':F(565673,10**6),'upper':F(22631,40000)}[mode]
        if args.verify_sha:
            before=fn.stat()
            with fn.open('rb') as stream: digest=hashlib.file_digest(stream,'sha256').hexdigest()
            assert digest==HASH[mode] and before.st_size==3495091968
            print(f'{fn} SHA256 {digest} size {before.st_size}',file=sys.stderr)
        if args.portable and args.edges: out[mode]=dict(f=portable(c),g=edges_buffer(c))
        elif args.portable: out[mode]=portable(c)
        elif args.edges: out[mode]=edges_buffer(c)
        elif args.buffer:
            b=buffer(c); out[b['mode']]=b
        elif args.bands: bands(c)
        else:
            print(fn,c.mode,c.lam,flush=True); scan(c)
        if args.verify_sha:
            after=fn.stat(); assert (after.st_ino,after.st_size,after.st_mtime_ns)==(before.st_ino,before.st_size,before.st_mtime_ns)
    if args.buffer or args.portable or args.edges:
        out['payload_sha256']=HASH
        assert set(out)=={'lower','upper','payload_sha256'}
        # Receipt-linked sparse extraction: source payload identities are from B,
        # not a new streaming hash of each 3.5 GB file.
        print(json.dumps(out,sort_keys=True,separators=(',',':')))
