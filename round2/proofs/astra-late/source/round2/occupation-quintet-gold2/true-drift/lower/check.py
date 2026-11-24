#!/usr/bin/env python3
"""Original Gold2 lower: disjoint large-payer FY and post-FY edge reserves."""
import argparse
import ast
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import file_digest, sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
import json
from os import close, dup, fdopen
from pathlib import Path
from runpy import run_path
from subprocess import check_output
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[3]
LO = Q(35083855339718511,62017003637500000)
HI = Q(127993535256830043,226242120050000000)
S = 2**40
RECEIPTS = {
    ('wide-payer',251,12,125,3,5): Q(58162652610249882843115293,6201700363750000000000),
    ('wide-host',251,18,34,8,5): Q(7127543113684417780424361,1240340072750000000000),
    ('old-donor',251,11,40,5,9): Q(25026790798544028667583523,6201700363750000000000),
}
POST = Q(6316495315166525735295397,78502536250000000000000)


def union():
    old = Q(2339594681813540639638101,10**22)+Q(39526146756920663734311,10**17)+Q(1046834122639903380009,2*10**16)
    add = POST+sum(RECEIPTS.values())
    assert add==Q(119326162107114155945908957363,6201700363750000000000000)
    assert 19000<add<30000
    assert old+add==Q(23173422990275479137666991749575647391,49613602910000000000000000000000)
    print('new union ticks:',add,float(add),'normalized:',add/S,float(add/S))
    print('original booked union plus extension ticks:',old+add,float(old+add))


def load(name, path):
    spec = spec_from_file_location(name, ROOT/path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fy_price(gold, e, args):
    class Residual(gold.Price):
        def rewards(self,h):
            d,c=h
            if c<2 or h in e.mids or h not in e.snap.inner: return {}
            ak = [e.fy(h,x) for x in e.knots]
            bands = {}
            for t,i in self.band_id.items():
                if args.mode=='old-donor' and t<6 or args.mode!='old-donor' and t>5: continue
                bands[i] = (Q(1,t),)
            if self.tail:
                cut=Q(1,self.tail[0]); bands[len(self.donor_bands)] = tuple(k for k in e.knots if k<cut)+(cut,)
            ss={i:tuple((x,e.fy(h,x)) for x in xs) for i,xs in bands.items()}
            out={}
            for i,r in enumerate(self.rs):
                if r[0]<=5 or max(d,r[0])<=18: continue
                if args.mode=='wide-payer' and r[0]<=34: continue
                if args.mode=='wide-host' and (r[0]>34 or (d<=150 and c<=10 and r[1]<=7)): continue
                z=e.g((d+1,c),r)-e.g(h,r)
                m=max(e.sig*(x*z-a) for x,a in zip(e.knots,ak))
                for j,xs in ss.items():
                    v=m-max(e.sig*(x*z-a) for x,a in xs)
                    assert v>=0
                    out[i,j]=self.round(v)
                    self.bands+=1
            return out

    gold.LO,gold.HI = LO,HI
    gold.ext.LO,gold.ext.HI=LO,HI
    a=SimpleNamespace(host=args.host,stock=args.stock,receiver=args.payer,rstock=args.pstock,
                      donor=args.donor,tstock=args.donor,bins=list(range(6,args.donor+1)),
                      mask='both',absorb=args.mode=='old-donor')
    r=Residual(e,a)
    vals=r.run()
    print(args.mode,'orientation ticks:',*[str(v*S) for v in vals],flush=True)
    total=sum(vals)
    key=(args.mode,args.host,args.stock,args.payer,args.pstock,args.donor)
    if key in RECEIPTS: assert total*S==RECEIPTS[key]
    print(args.mode,'normalized',total,'ticks',total*S,float(total*S),'states',r.path_count,'slots',r.bands,flush=True)
    return total*S


class Edges:
    def __init__(self,e):
        self.e,self.s=e,e.snap
        self.k,self.d=self.s.K2,self.s.d2
        self.sig,self.cap=e.sig,e.cap
        assert (self.k,self.d,self.sig,self.cap)==(22,252,-1,792)
        assert e.knots[:2]==(0,Q(1,self.d)) and e.h['d3e']<=self.d
        assert all(self.s.g(i,j)==0 for i,j in product(range(self.s.tail_base,self.s.inf+1),repeat=2))

    @cache
    def knots(self,y): return tuple((self.sig*x,self.sig*self.e.fy(y,x)) for x in self.e.knots)

    def maximum(self,y,jump): return max(x*jump-a for x,a in self.knots(y))

    def relaxed(self,x,y):
        out=-self.sig*self.s.lam*self.e.g(x,y)
        b=[]
        for h,p in ((x,y),(y,x)):
            a,z,fy=self.e.end(h,p)
            out+=self.sig*a+((h[1]-1)*self.maximum(h,z) if fy else 0)
            b.append(0 if fy else self.sig*z)
        return out,tuple(b)

    @cache
    def inf(self,y):
        d,c=y
        e,s=self.e,self.s
        def g(p): return s.g(e.weights(p,s.inner,self.d,self.k)[0][0],s.inf)
        old=g(y); j=g((d+1,c))-old; r=s.center(*y)
        a=(d-c)*j-s.lam*old
        if d>c:
            post=(d,c+1); m=e.mids.get(post)
            h3=0 if m is None else s.item(e.off_h3+m*e.hstride+(e.n3-1)*e.n3//2)
            a+=Q(d-c,d)*(g(post)-old+h3)
        if r is not None: a-=r*(s.f(d+1,c)-s.f(d,c))
        if y in e.mids: return self.sig*(a+(c-1)*r*j),0,j
        if r is not None and c>=2: return self.sig*(a+(c-1)*e.fy(y,0))+(c-1)*self.maximum(y,j),0,j
        return self.sig*a,self.sig*j,j

    @cache
    def knot(self,c,y):
        e,s=self.e,self.s
        d,k=y; x=self.d,c
        a,b,fy=e.end(y,x)
        old=e.g(x,y); j=e.g((d+1,k),x)-old
        yy=e.weights(y,s.inner,self.d,self.k)[0][0]; g0=s.g(s.inf,yy)
        base=self.sig*(a-s.lam*old)+Q(max(self.sig*j,self.sig*self.inf(y)[2])-self.sig*j,self.d)
        if fy: base+=(k-1)*self.maximum(y,j)
        return base,self.sig*(g0-old),self.sig*(e.g((self.d,min(c+1,self.k+1)),y)-old),0 if fy else self.sig*b

    def fixtures(self):
        counts=Counter()
        samples=((3,1),(18,4),(31,30),(79,10),(79,11),(80,2),(87,4),(88,3),
                 (150,6),(251,22),(251,30),(252,1),(252,22),(252,23),(391,1),(391,62),(4000,100))
        for x,y in product(samples,repeat=2):
            for u,v in product((Q(0),Q(1,252),Q(1,17),Q(1,2)),repeat=2):
                a,b=self.relaxed(x,y)
                row=self.e.row(x,y,[u]*(x[1]-1),[v]*(y[1]-1))
                gaps=0
                for h,p,z in ((x,y,u),(y,x,v)):
                    _,j,fy=self.e.end(h,p)
                    if fy: gaps+=(h[1]-1)*(self.maximum(h,j)-self.sig*(z*j-self.e.fy(h,z)))
                actual=a+b[0]*(x[1]-1)*u+b[1]*(y[1]-1)*v
                assert gaps>=0 and actual==row+gaps<=self.cap,(x,y,u,v,actual)
                if x[0]>=self.d>y[0]:
                    d,c=x; t=Q(self.d,d)
                    f,g,k,b=self.knot(min(c,self.k+1),y); a0,b0,_=self.inf(y)
                    ext=(y[1]-1)*v
                    gamma=Q(d-c+Q(1,y[0])+(c-1)*u,d+1)
                    beta=Q(d-c,d)
                    assert actual<=t*(f+gamma*g+beta*k+b*ext)+(1-t)*(a0+b0*ext)<=self.cap
                    counts['mixed lifted']+=1
                if min(x[0],y[0])>=self.d: assert actual==0
                counts['row and all-FY checks']+=1
        print(dict(counts),flush=True)


def post_fy(e):
    post=Edges(e)
    post.hi=HI
    post.tails=range(1,post.k+2)
    # Reuse the generic finite/bridge recurrences only, with this module's
    # independently checked mixed-tail row supplying all signed rates.
    path=ROOT/'occupation-quintet/edges/post-fy/replay.py'
    tr=ast.parse(path.read_bytes())
    cls=next(n for n in tr.body if isinstance(n,ast.ClassDef))
    ns=dict(Q=Q,product=product,GRID=10**8,J=(2,1))
    fs=[n for n in tr.body if isinstance(n,ast.FunctionDef) and
        n.name in {'solve','promoter','tail_promoter'}]
    fs += [n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name in {'finite','bridge'}]
    exec(compile(ast.Module(fs,type_ignores=[]),str(path),'exec'),ns)
    post.tail_inf=lambda y: (post.inf(y)[0],(y[1]-1)*post.inf(y)[1])
    post.relaxed_original=post.relaxed
    post.relaxed=lambda x,y: (lambda a,b:(a,tuple(z*(p[1]-1) for p,z in zip((x,y),b))))(*post.relaxed_original(x,y))
    enters,n=ns['finite'](post,25,3)
    finite=ns['promoter'](HI,LO,sorted(enters,reverse=True),lambda d,c:enters[d,c])
    print('post-FY finite',finite,float(finite),'edges',n,flush=True)
    u={}
    for y in [(d,c) for d in range(10,1,-1) for c in range(d,0,-1)]:
        f,k=y
        b=Q(f-k,f)
        a0,z0,_=post.inf(y)
        for c in reversed(post.tails):
            a,j,h,z=post.knot(c,y)
            gamma=(Q(post.d-c+Q(1,f),post.d+1),Q(1)) if c<=post.k else (0,1)
            beta=(Q(post.d-c,post.d),Q(1)) if c<=post.k else (0,)
            corners=[]
            for g,p,r,t in product(gamma,beta,(0,Q(k-1,2)),(0,1)):
                assert a0+z0*r<=post.cap
                reward=post.cap-a-j*g-h*p-z*r
                receive=f-k+r+Q(t,post.d)
                assert reward>=0
                n=reward+receive*u.get((c,(f+1,k)),0)+b*u.get((c,(f,k+1)),0)+p*u.get((c+1,y),0)
                corners.append((n,HI+receive+b+p+g))
            u[c,y]=ns['solve'](corners)
    enters,n=ns['bridge'](post,u,10,8)
    mixed=ns['tail_promoter'](post,LO,u,enters)
    helper=run_path(str(Path(__file__).resolve().parents[1]/'check.py'))
    a,b=[helper['host'](252,up) for up in (False,True)]
    g=helper['geom']; x=helper['X']; l=helper['LO']; h=helper['HI']
    assert Q(l,x)<=LO<=HI<=Q(h,x)
    double=post.cap*max(0,a-(10*b+g(10)/(1-g(1)))/252)*a/Q(h,x)
    # One outward rounding to keep the receipt compact.
    double=Q(double*10**12//1,10**12)
    assert finite+mixed+double==POST
    for group,v in (('mixed',mixed),('double',double),('total',finite+mixed+double)):
        print('post-FY',group,'ticks',v,float(v),'normalized',v/S,flush=True)


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--mode',choices=('old-donor','wide-payer','wide-host','cap','post-fy','union'),default='old-donor')
    for name,val in (('host',150),('stock',10),('payer',28),('pstock',5),('donor',8)):
        cli.add_argument('--'+name,type=int,default=val)
    a=cli.parse_args()
    if a.mode=='union':
        union()
        return
    assert __debug__ and 2<=a.host<252 and 6<=a.payer<=251
    assert a.mode in {'cap','post-fy'} or (a.donor>=7 if a.mode=='old-donor' else a.donor==5)
    wire=load('lower_wire','attacks/neon2-quintet/gold2/physical/check.py')
    w=load('lower_rows','occupation-quintet-gold2/post-fy/replay.py')
    src=wire.ROOT/'source'
    assert check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()==wire.PIN
    data=(src/wire.CERT/'certificate_gold2.json').read_bytes()
    assert sha256(data).hexdigest()=='5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'
    docs=json.loads(data); receipt=docs['lower']
    assert (Q(receipt['certified_exact']),Q(docs['upper']['certified_exact']))==(LO,HI)
    p=wire.Gold(wire.ROOT/'quintet_lower.bin')
    try:
        with fdopen(dup(p.fd),'rb') as f:
            assert file_digest(f,'sha256').hexdigest()==receipt['verify']['sha256']==wire.DIGEST['lower']
        with w.Snapshot('lower') as s:
            assert s.stat==p.stat and s.scale==S
            e=w.Rows(s,p,receipt)
            if a.mode=='cap': Edges(e).fixtures()
            elif a.mode=='post-fy': post_fy(e)
            else: fy_price(w.gold,e,a)
            s.check_unchanged()
    finally:
        close(p.fd)
    print('PASS original Gold2 lower pins and exact comparisons',flush=True)


if __name__=='__main__': main()
