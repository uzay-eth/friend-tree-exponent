#!/usr/bin/env python3
"""Exact sparse D18 higher-body and shared-direction charge extractor."""
import argparse
import json
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


def load_source(name,path):
    spec=spec_from_file_location(name,path.resolve(strict=True))
    module=module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if __name__ != '__main__':
    # The adjacent occupation experiments import these two helpers from here.
    # The standalone production CLI instead requires explicit source paths.
    R2=Path(__file__).resolve().parents[2]
    Cert=load_source('neon_audit',R2/'attacks'/'neon1'/'audit.py').Cert
    physical=load_source('physical_moves',R2/'finite'/'residual_entries.py')

NB=(1,0,1)
LO,HI=Q('.565673570322'),Q('.565775')
C0={'upper':Q('1.95820675182'),'lower':Q('1.84124847315')}
# Outward endpoints and denominator ceilings from the f-only renewal replay.
F_EDGE={'upper':Q('.565773666188'),'lower':Q('.565674722414')}
F_DEN={'upper':Q('1.518769126087491'),'lower':Q('1.458909469275164')}
PAYLOAD_SHA={
    'upper':'1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41',
    'lower':'9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
}
SPEC={'f':(-1,), 'g':(-1,0),'hV':(-1,0,1),'hC':(-1,0,0),
      'P4':(-1,0,1,2),'Yt':(-1,0,1,1),'Yb':(-1,0,0,2),'S3':(-1,0,0,0)}


def label(x,i,v):return x[:i]+(v,)+x[i+1:]
def up(x):return x[0]+1,x[1],x[2]
def bud(x):return x[0]-1,x[1]+1,x[2]


@lru_cache(maxsize=12000)
def embedding(p):
    ch=[[] for _ in p]
    for i,j in enumerate(p[1:],1):ch[j].append(i)
    out=[]
    for a in range(len(p)):
        out.append(('f',(a,)))
        for b in ch[a]:
            out.append(('g',(a,b)))
            for c in ch[b]:
                out.append(('hV',(a,b,c)))
                for d in ch[c]:out.append(('P4',(a,b,c,d)))
                for w in ch[a]:
                    if w!=b:out.append(('Yb',(a,w,b,c)))
            for c,d in combinations(ch[b],2):out.append(('Yt',(a,b,c,d)))
        for u in combinations(ch[a],2):out.append(('hC',(a,*u)))
        for u in combinations(ch[a],3):out.append(('S3',(a,*u)))
    return tuple(out)


@lru_cache(maxsize=40000)
def shape(z,h):
    p=[];x=[]
    def visit(z,par):
        _,l,ys=z;i=len(p);p.append(par);x.append((l,len(ys),h if i==0 else 1))
        for y in ys:visit(y,i)
    visit(z,-1)
    return tuple(p),tuple(x)


class Extract:
    def __init__(self,c):
        self.c=c;self.sign=1 if c.mode==2 else -1
        self.weight=lru_cache(maxsize=50000)(self.weight)
        self.row=lru_cache(maxsize=40000)(self.row)

    def supported(self,f,x):
        c=self.c
        return (f in ('f','g') or f=='hV' and x[1] in c.hv or f=='hC' and x[0] in c.hc
                or f=='P4' and x[1] in c.p4 and x[2] in c.p4 or f=='Yt' and x[1] in c.yt
                or f=='Yb' and x[0] in c.ybu and x[2] in c.ybw or f=='S3' and x[0] in c.s3)

    def weight(self,f,x):
        c=self.c
        if not self.supported(f,x):return 0
        a,*b=x
        if f=='f':return int(c.f[c.vi(*a)])
        if f=='g':return int(c.g[c.state[a],c.state[b[0]]])
        if f=='hV':return int(c.hvtab[c.par[a],c.hv[b[0]],c.child[b[1]]])
        if f=='hC':return int(c.hctab[c.hc[a],c.tp(*b)])
        if f=='P4':return int(c.p4tab[c.par[a],c.p4[b[0]],c.p4[b[1]],c.child[b[2]]])
        if f=='Yt':return int(c.yttab[c.par[a],c.yt[b[0]],c.tp(*b[1:])])
        if f=='Yb':return int(c.ybtab[c.ybu[a],c.ybw[b[1]],c.child[b[0]],c.child[b[2]]])
        return int(c.s3tab[c.s3[a],c.tt(*b)])

    def route(self,f,x,i,upward):
        # Return whether the missing donor is charged to the indicated
        # supported larger embedding. A false entry keeps its live rate.
        c=self.c
        if f=='f':return True
        if f=='g':return (x[0] in (c.hv if upward else c.hc)) if i==0 else x[1] in c.hv
        if f=='hV':
            return ((x[0] in c.p4 and x[1] in c.p4) if upward else (x[0] in c.ybu and x[1] in c.ybw)) if i==0 else (x[1] in c.yt if i==1 else x[1] in c.p4 and x[2] in c.p4)
        if f=='hC':
            return (x[0] in (c.yt if upward else c.s3)) if i==0 else x[0] in c.ybu and x[i] in c.ybw
        return False

    def row(self,f,x):
        """Exact signed baseline and free directional coefficients (ticks)."""
        c=self.c; p=SPEC[f];n=len(p)
        assert all(sum(y)+1<c.D for y in x)
        if not self.supported(f,x):return Q(0),()
        own=self.weight(f,x);base=-c.lam*own; free=[]
        for i,y in enumerate(x):
            L,k,h=y; D=sum(y)
            ne=[j for j in range(n) if p[j]==i]+([p[i]] if i else [])
            dp=self.weight(f,label(x,i,up(y)))-own
            base+=(L+sum(Q(1,sum(x[j])) for j in ne))*dp
            # New motifs are charged to their old embedding at the unique
            # promotion. Enumerating old-set plus newborn chooses that route.
            if L:
                xb=label(x,i,bud(y))+(NB,)
                newborn=sum(self.weight(f1,tuple(xb[j] for j in ids)) for f1,ids in embedding(p+(i,)) if len(ids)==n+1)
                if f=='f':newborn+=self.weight('f',(NB,))
                base+=Q(L,D)*(self.weight(f,xb[:-1])-own+newborn)
            for role,count in ((True,h-int(i!=0)),(False,k-p.count(i))):
                if count<=0:continue
                if self.route(f,x,i,role):
                    r=c.rpar if role else c.rch
                    base+=count*Q(int(r[c.state[y]]),c.Sr)*dp
                else:free.append((i,role,count,self.sign*dp))
        # Incoming centered transfers from each connected deletion.  The
        # deleted tip donates to its sole neighbor in the remaining motif.
        if n>1:
            for donor in range(n):
                for f0,ids in embedding(p):
                    if len(ids)!=n-1 or donor in ids:continue
                    contact=[j for j in ids if p[j]==donor or p[donor]==j]
                    assert len(contact)==1
                    j=contact[0];upward=p[j]==donor
                    v=tuple(x[i] for i in ids);idx=ids.index(j)
                    if not self.supported(f0,v):continue
                    assert self.route(f0,v,idx,upward),(f,f0,x,ids)
                    r=c.rpar if upward else c.rch
                    dp=self.weight(f0,label(v,idx,up(x[j])))-self.weight(f0,v)
                    base+=(Q(1,sum(x[donor]))-Q(int(r[c.state[x[j]]]),c.Sr))*dp
        return self.sign*base,tuple(free)

    def charges(self,z,h=1,owned=True,detail=None):
        p,x=shape(z,h); n=len(x)
        if any(sum(y)+1>=self.c.D for y in x):return {},Q(0)
        totals=defaultdict(Q); unknown=Q(0)
        for f,ids in embedding(p):
            if owned and (len(ids)<3 or ids[0]!=0):continue
            states=tuple(x[j] for j in ids)
            if not self.supported(f,states):continue
            a,free=self.row(f,states)
            floor=a+sum(count*min(0,Q(dp,2)) for _,_,count,dp in free)
            assert floor>=0,(f,states,floor)
            slope=Q(0)
            for i,role,count,dp in free:
                j=ids[i]; worst=min(0,Q(dp,2)); a+=count*worst
                # Add only the excess of the same row term over its GMP
                # directional floor; use actual crown donors or joint parent.
                if role:
                    if p[j]<0:
                        unknown+=dp;slope+=dp;a-=worst
                    else:a+=Q(dp,sum(x[p[j]]))-worst
                else:
                    donors=[v for v in range(n) if p[v]==j and v not in ids]
                    assert len(donors)==count
                    a+=sum(Q(dp,sum(x[v]))-worst for v in donors)
            totals[f]+=a/self.c.S
            if detail is not None:detail.append((f,ids,floor/self.c.S,a/self.c.S,slope/self.c.S))
        return dict(totals),unknown/self.c.S

    def credit(self,z):
        a,b=self.charges(z);return sum(a.values())+min(0,b/2)


def validate(E,cap):
    # Independently evaluate the complete drift of real whole trees (h=0),
    # where no exterior receiver or donor remains. Compare all charge classes.
    layer={(0,2,())};checked=0
    while layer:
        nxt=set()
        for z in layer:
            p,x=shape(z,0); kids=[[j for j in range(len(p)) if p[j]==i] for i in range(len(p))]
            score=lambda x,k:E.c.potential(k,x)
            val=score(x,kids); drift=-E.c.lam*val
            for i,(L,k,h) in enumerate(x):
                rate=L+sum(Q(1,sum(x[j])) for j in kids[i]+([p[i]] if i else []))
                drift+=rate*(score(label(x,i,up(x[i])),kids)-val)
                if L:
                    new=[a[:] for a in kids]+[[]];new[i].append(len(x))
                    drift+=Q(L,sum(x[i]))*(score(label(x,i,bud(x[i]))+(NB,),new)-val)
            a,b=E.charges(z,0,False)
            assert b==0 and sum(a.values())==E.sign*drift/E.c.S,(z,a,drift)
            checked+=1
            if physical.statistics(z)[0]<cap:
                def trees(i,xx,children):return (0,xx[i][0],tuple(sorted(trees(j,xx,children) for j in children[i])))
                for i in range(len(x)):
                    nxt.add(trees(0,label(x,i,up(x[i])),kids))
                    if x[i][0]:
                        new=[a[:] for a in kids]+[[]];new[i].append(len(x))
                        nxt.add(trees(0,label(x,i,bud(x[i]))+(NB,),new))
        layer=nxt
    print('independent whole-tree charge/drift checks',checked,flush=True)


def replay(E,cap,dump=False):
    # Suboccupation keeping only intrinsic arrivals. Before the next event,
    # the missing parent rate is at most 1/2 for every birth environment.
    mass={physical.newborn(2):Q(1)}; occ=defaultdict(Q)
    grid=10**23
    down=lambda q:Q(q*grid//1,grid)
    while mass:
        nxt=defaultdict(Q)
        for z,w in mass.items():
            events=physical.intrinsic_moves(z)
            den=HI+sum(a for _,a in events)+Q(1,2)
            r=down(w/den);occ[z]+=r
            for y,a in events:
                if a and physical.statistics(y)[0]<=cap:nxt[y]+=down(r*a)
        mass=nxt
    total=baseline=direction=joint=Q(0); bins=defaultdict(Q); families=defaultdict(Q); witness=[]; records=[]
    for z,w in occ.items():
        detail=[]; a,b=E.charges(z,detail=detail); g=sum(a.values())+min(0,b/2)
        assert g>=0,(z,a,b)
        if g:
            m=LO*w; total+=m*g
            b0=sum(v for _,_,v,_,_ in detail)
            split=sum(v+min(0,r/2) for _,_,_,v,r in detail)
            assert 0<=b0<=split<=g
            baseline+=m*b0;direction+=m*(split-b0);joint+=m*(g-split)
            for f,_,_,v,r in detail:families[f]+=m*(v+min(0,r/2))
            fs=tuple(sorted(a));bins[fs]+=m*g
            witness.append((m*g,z,g,m,a,b))
            if dump:records.append([z,str(m),[[f,ids,str(v),str(w),str(r)] for f,ids,v,w,r in detail]])
    assert total==baseline+direction+joint
    assert sum(families.values())+joint==total
    if dump:dump.write_text(json.dumps(records,separators=(',',':'))+'\n')
    print('intrinsic-only fringe states',len(occ),'paid',len(witness),'joint pi slack lower',str(total),float(total),flush=True)
    for name,x in [('baseline row floors',baseline),('observed donor direction',direction),('shared exterior nonattainment',joint)]:
        print(name,str(x),float(x),flush=True)
    for fs,x in sorted(bins.items()):print('disjoint family-set bin',fs,float(x))
    for f,x in sorted(families.items()):print('unique family rows, before joint bonus',f,str(x),float(x))
    for w,z,g,m,a,b in sorted(witness,reverse=True)[:4]:print('witness',z,'price',float(w),'slack',str(g),'mass',float(m),'charges',a,'external slope',b)
    return total


def birth(E,cap):
    """Direct physical lifetime from newborn, allowing an adapted parent send."""
    grid=10**22; down=lambda q:Q(q*grid//1,grid)
    @lru_cache(maxsize=None)
    def V(z):
        size,degree,_=physical.statistics(z)
        if size>cap:return Q(0)
        a,b=E.charges(z);g=sum(a.values())
        edges=physical.intrinsic_moves(z)
        x=sum(r*V(y) for y,r in edges); rate=sum(r for _,r in edges)
        assert rate==size-Q(1,degree)
        advance=V(physical.plus(z))
        assert g+min(0,b/2)>=0
        ends=(Q(0),Q(1,2))
        v=down(min((g+x+r*(b+advance))/(HI+rate+r) for r in ends))
        assert all((HI+rate+r)*v<=g+x+r*(b+advance) for r in ends)
        return v
    v=LO*V(physical.newborn(2))
    print('direct actual stable-age pi lower',str(v),float(v),'contexts',V.cache_info().currsize,flush=True)
    return v


def main():
    global physical
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('blob',type=Path)
    p.add_argument('--audit-source',type=Path,required=True,help='Astra audit.py defining Cert')
    p.add_argument('--physical-source',type=Path,required=True,help='Astra residual_entries.py defining physical moves')
    p.add_argument('--cap',type=int,choices=range(2,13),default=12,metavar='2..12')
    p.add_argument('--dump',type=Path)
    p.add_argument('--ledger',action='store_true',help='also replay the intrinsic-only suboccupation by family')
    a=p.parse_args()
    if not __debug__:p.error('assertions must be enabled (do not use -O or PYTHONOPTIMIZE)')
    if a.blob.stat().st_size!=3495091968:p.error('wrong Neon D18 payload size')
    with a.blob.open('rb') as f:
        digest=hashlib.file_digest(f,'sha256').hexdigest()
    side=next((s for s,d in PAYLOAD_SHA.items() if d==digest),None)
    if side is None:p.error(f'payload SHA-256 does not match either pinned Neon D18 input: {digest}')
    print(side,'pinned payload SHA-256',digest,flush=True)
    c=load_source('neon_audit',a.audit_source).Cert(a.blob)
    physical=load_source('physical_moves',a.physical_source)
    assert c.mode==(2 if side=='lower' else 1)
    assert (c.D,c.K,c.S,c.lam)==(18,17,10**11,Q('.565673' if side=='lower' else '.565775'))
    # Rates and charges do not use birth marks. Collapse all of them to a
    # single state label; 2 here puts no lower bound on the exterior input.
    physical.MARK_CAP=2
    E=Extract(c);validate(E,8);direct=birth(E,a.cap)
    if a.ledger or a.dump:
        assert direct>=replay(E,a.cap,a.dump)
    print(side,'additional lambda improvement lower',str(direct/C0[side]),float(direct/C0[side]),flush=True)
    # Reuse the existing numerator/denominator proof, not the endpoint as a
    # new pointwise certificate. Both rewards are for the original lambda.
    edge=F_EDGE[side]+E.sign*direct/F_DEN[side]
    ticks=edge*10**12
    outward=E.sign*((E.sign*ticks)//1)
    print(side,'outward endpoint with first-body renewal',f'{outward//10**12}.{outward%10**12:012d}')


if __name__=='__main__':main()
