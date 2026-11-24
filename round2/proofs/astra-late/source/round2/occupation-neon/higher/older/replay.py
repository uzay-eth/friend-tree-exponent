#!/usr/bin/env python3
"""D18 exact local drift, sampled omitted-charge census, and controlled crown prices."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import cache, lru_cache
import hashlib
from itertools import combinations, product
from math import prod
from pathlib import Path
from random import Random
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from extract import Cert, Extract, HI, LO, NB, PAYLOAD_SHA, SPEC, bud, embedding, label, up

FAMILIES=tuple(SPEC)
DEN={'lower':Q('1.356833021311'),'upper':Q('1.367169995210')}


@cache
def deletions(f):
    t=SPEC[f]
    out=[]
    for g,ids in embedding(t):
        if len(ids)!=len(t)-1:continue
        donor,=set(range(len(t)))-set(ids)
        j,=(v for v in ids if t[v]==donor or t[donor]==v)
        out.append((g,ids,donor,ids.index(j),t[j]==donor))
    return out


@cache
def extensions(f,i):
    t=SPEC[f]+(i,)
    return tuple((g,ids) for g,ids in embedding(t) if len(ids)==len(t))


def pick(x,ids):return tuple(x[i] for i in ids)


class Local(Extract):
    def live(self,f,x):
        if f=='hV':return 1<=x[1][1]<=self.c.K
        if f=='hC':return 2<=x[0][1]<=self.c.K
        return self.supported(f,x)

    @lru_cache(maxsize=45000)
    def value(self,f,x):
        if f=='f':return Q(self.weight(f,x))
        if not self.live(f,x):return Q(0)
        return sum((prod(w for _,w in terms)*self.weight(f,tuple(v for v,_ in terms))
                    for terms in product(*(self.c.embed(s) for s in x))), Q(0))

    def moved(self,f,x,i,parent):
        l,k,h=x[i]; D=l+k+h
        if f=='f':return D<self.c.D or k<=self.c.K
        if f=='g':return k<=self.c.K
        return self.route(f,x,i,parent)

    def center(self,x,parent):
        c=self.c
        if sum(x)<c.D:return Q(int((c.rpar if parent else c.rch)[c.state[x]]),c.Sr)
        _,k,h=x; j=(1-h)*(c.K+2)+min(k,c.K+1)
        return Q(int((c.trpar if parent else c.trch)[j]),c.Sr)

    def transfer(self,f,x,i):
        c=self.c;y=x[i];D=sum(y)
        if D<c.D or f not in ('f','g'):
            return self.value(f,label(x,i,up(y)))-self.value(f,x)
        _,k,h=y; j=(1-h)*(c.K+2)+k
        if f=='f':return Q(int(c.f2[j])*D+int(c.f3[j]),D**3)
        t=Q(c.D,D); inf=('inf',); t1=('T1',k,h)
        ans=Q(0)
        for v,w in c.embed(x[1-i]):
            if i:
                d=self.weight('g',(v,inf))-self.weight('g',(v,t1)); ref=c.yw[c.par[v],k]
            else:
                d=self.weight('g',(inf,v))-self.weight('g',(t1,v)); ref=c.yp[j,c.child[v]]
            ans+=w*(3*t*t*(1-t)*Q(d,c.D)+t**3*int(ref))
        return ans

    @lru_cache(maxsize=35000)
    def terms(self,f,x):
        """Per occurrence: killing term, existing gain/bud deltas, newborn rewards."""
        old=self.value(f,x)
        gain=tuple(self.value(f,label(x,i,up(s)))-old for i,s in enumerate(x))
        promote=[]
        for i,s in enumerate(x):
            if not s[0]:promote.append(Q(0));continue
            y=label(x,i,bud(s)); z=y+(NB,)
            born=sum(self.value(g,pick(z,ids)) for g,ids in extensions(f,i))
            if f=='f':born+=self.value('f',(NB,))
            promote.append(self.value(f,y)-old+born)
        return old,gain,tuple(promote)

    def charge(self,f,x,ext):
        """Raw row in normalized units. ext[i]=(exterior parent sum, child sum)."""
        old,gain,born=self.terms(f,x); t=SPEC[f];c=self.c
        ans=-c.lam*old
        for i,(l,k,h) in enumerate(x):
            inner=sum(Q(1,sum(x[j])) for j in range(len(x)) if t[j]==i or t[i]==j)
            ans+=(l+inner+sum(ext[i]))*gain[i]+Q(l,l+k+h)*born[i]
            for parent,num,rho in ((True,h-int(i!=0),ext[i][0]),(False,k-t.count(i),ext[i][1])):
                if num and self.moved(f,x,i,parent):
                    ans+=(num*self.center(x[i],parent)-rho)*self.transfer(f,x,i)
        for g,ids,donor,i,parent in deletions(f):
            y=pick(x,ids)
            if self.moved(g,y,i,parent):
                ans+=(Q(1,sum(x[donor]))-self.center(y[i],parent))*self.transfer(g,y,i)
        return self.sign*ans/c.S

    def extrema(self,f,x,host_half_load=False):
        z=[(Q(0),Q(0))]*len(x); a=self.charge(f,x,z)
        low=high=a; t=SPEC[f]; gain=self.terms(f,x)[1]; widths=[]
        for i,(_,k,h) in enumerate(x):
            for parent,n in ((True,h-int(i!=0)),(False,k-t.count(i))):
                if not n:continue
                d=gain[i]-(self.transfer(f,x,i) if self.moved(f,x,i,parent) else 0)
                shift=Q(n,2)*d*self.sign/self.c.S
                left=shift/2 if host_half_load and i==0 and not parent else Q(0)
                low+=min(left,shift); high+=max(left,shift)
                widths.append((abs(shift-left),i,'parent' if parent else 'children'))
        assert 0<=low<=high,(f,x,low,high)
        return low,high,sorted(widths,reverse=True)

    def occurrences(self,ch,x):
        c=self.c
        for a,kids in enumerate(ch):
            yield 'f',(a,)
            for b in kids:
                yield 'g',(a,b)
                for d in ch[b]:
                    if x[b][1]<=c.K:yield 'hV',(a,b,d)
                    if x[b] in c.p4 and x[d] in c.p4:
                        for e in ch[d]:yield 'P4',(a,b,d,e)
                    if x[a] in c.ybu and x[b] in c.ybw:
                        for w in kids:
                            if w!=b:yield 'Yb',(a,w,b,d)
                if x[b] in c.yt:
                    for d,e in combinations(ch[b],2):yield 'Yt',(a,b,d,e)
            if len(kids)<=c.K:
                for b,d in combinations(kids,2):yield 'hC',(a,b,d)
            if x[a] in c.s3:
                for b,d,e in combinations(kids,3):yield 'S3',(a,b,d,e)

    def tree(self,p,x,verify=False):
        ch=[[] for _ in x]
        for j,i in enumerate(p[1:],1):ch[i].append(j)
        d=[sum(s) for s in x]; r=[Q(1,s) for s in d]
        child=[sum((r[j] for j in kids),Q(0)) for kids in ch]
        gain=[s[0]+child[i]+(r[p[i]] if i else 0) for i,s in enumerate(x)]
        flow=defaultdict(Q); phi=defaultdict(Q); raw=[]
        for f,ids in self.occurrences(ch,x):
            t=SPEC[f]; states=pick(x,ids)
            if not self.live(f,states):continue
            ext=[(r[p[j]] if k==0 and j else Q(0),
                  child[j]-sum(r[ids[b]] for b in range(1,len(ids)) if t[b]==k))
                 for k,j in enumerate(ids)]
            v=self.charge(f,states,ext)
            assert v>=0,(f,states,ext,v)
            raw.append((f,ids,v))
            old,dl,db=self.terms(f,states)
            phi[f]+=old
            flow[f]+=sum(gain[j]*dl[k]+x[j][0]*r[j]*db[k] for k,j in enumerate(ids))
        physical={f:self.sign*(flow[f]-self.c.lam*phi[f])/self.c.S for f in FAMILIES}
        assert sum(physical.values())==sum(v for _,_,v in raw)
        if verify:
            c=self.c; total=c.potential(ch,x)
            assert total==sum(phi.values())
            direct=-c.lam*total
            for i,s in enumerate(x):
                direct+=gain[i]*(c.potential(ch,label(x,i,up(s)))-total)
                if s[0]:
                    y=[z[:] for z in ch]+[[]]; y[i].append(len(x))
                    direct+=s[0]*r[i]*(c.potential(y,label(x,i,bud(s))+(NB,))-total)
            assert self.sign*direct/c.S==sum(v for _,_,v in raw)
        return raw,ch


def sample(n,rng):
    p=[-1,0,0]; adj=[[1,2],[0],[0]]
    for i in range(3,n):
        near=adj[rng.randrange(i)]; j=near[rng.randrange(len(near))]
        adj[j].append(i); adj.append([j]); p.append(j)
    return compress(p,adj)


def compress(p,adj):
    cores=[i for i,a in enumerate(adj) if len(a)>1]; index={j:i for i,j in enumerate(cores)}
    ps=tuple(index[p[u]] if u else -1 for u in cores); m=Counter(ps)
    return ps,tuple((len(adj[u])-m[i]-int(i>0),m[i],int(i>0)) for i,u in enumerate(cores))


def fixture(p):
    adj=[[] for _ in p]
    for j,i in enumerate(p[1:],1):adj[j].append(i); adj[i].append(j)
    return compress(p,adj)


def validate(E):
    rng=Random(8098)
    cases=[sample(n,rng) for n in range(3,23)]
    # P4, both forks, 17/18 boundaries, churn at k=17/18, and multiple tails.
    spine=[-1,0,1,2,3,0,1,2,2,3,4,5,6,7,8]
    cases.extend(fixture(spine+[0]*a+[1]*b+[2]*c) for a,b,c in ((0,0,0),(13,0,0),(14,20,25),(30,40,50)))
    fan=[-1]+[0]*20+[i for i in range(1,21)]
    cases.extend(fixture(fan[:1+k]+list(range(1,k+1))) for k in (16,17,18,20))
    seen=set()
    for p,x in cases:
        rows,_=E.tree(p,x,verify=True)
        seen.update(f for f,_,_ in rows)
    assert seen==set(FAMILIES),seen
    print('exact all-family incident/whole-potential checks',len(cases),'families',len(seen),flush=True)


def bands(d):return 0 if d<8 else 8 if d<18 else 18 if d<32 else 32 if d<64 else 64


def census(E,side,n,reps,seed):
    rng=Random(seed); mass=defaultdict(Q); counts=Counter(); best={}; sectors=defaultdict(Q); cores=0
    for _ in range(reps):
        p,x=sample(n,rng); rows,ch=E.tree(p,x); cores+=len(x)
        sz=[s[0]+1 for s in x]; nc=[1]*len(x)
        for j in range(len(x)-1,0,-1):sz[p[j]]+=sz[j];nc[p[j]]+=nc[j]
        for f,ids,val in rows:
            a=ids[0]
            if f in ('f','g') or a==0 or sz[a]<=12:continue
            degrees=tuple(sum(x[j]) for j in ids)
            # Each occurrence belongs to exactly one bucket; earliest max
            # degree slot breaks ties. 5+ is disjoint even from unlimited arms.
            slot=max(range(len(ids)),key=degrees.__getitem__); host=ids[slot]
            key=(min(nc[a],5),f,slot,bands(degrees[slot]),min(x[host][1],8),
                 sum(x[j]==NB for j in ids))
            mass[key]+=val;counts[key]+=1
            if nc[a]>=5:
                sectors['apex tail' if degrees[0]>=18 else 'non-apex tail' if degrees[slot]>=18 else 'all inner']+=val
            if nc[a]>=5 and degrees[slot]>=18 and (f not in best or val>best[f][0]):
                best[f]=(val,pick(x,ids))
    print('sample only: trees',reps,'ordinary size',n,'seed',seed,'divisor (all cores)',cores,flush=True)
    for nc in (3,4,5):
        chosen=[k for k in mass if k[0]==nc]
        print('sample bucket cores',nc,'higher reward/core',float(sum(mass[k] for k in chosen)/cores),flush=True)
        for k in sorted(chosen,key=mass.get,reverse=True)[:10]:
            print(' rank',k,'reward/core',float(mass[k]/cores),'occ/core',counts[k]/cores,
                  'mean row',float(mass[k]/counts[k]),flush=True)
    print('sample >=5 core sector reward/core', {k:float(v/cores) for k,v in sorted(sectors.items())},flush=True)
    target=DEN[side]*Q(1,10**5)
    for f,(r,x) in sorted(best.items(),key=lambda item:item[1][0],reverse=True):
        low,high,widths=E.extrema(f,x); assert low<=r<=high
        print('high witness',f,x,'sample row',float(r),'uniform row interval',(float(low),float(high)),
              'incidences/core necessary for 1e-5',float(target/high) if high else None,
              'sufficient using floor',float(target/low) if low else None,
              'largest rate width',[(i,s,float(w)) for w,i,s in widths[:2]],flush=True)
        if f in ('Yt','P4','hV') and sum(x[0])>=18:
            floor,_,_=E.extrema(f,x,host_half_load=True)
            print(' conditional same descriptor, omitted host child load >= count/4: floor',float(floor),
                  'sufficient incidences/core',float(target/floor) if floor else None,flush=True)
    return mass


def tagged(E,side,host_cap,tip_cap,child_cap,kcap,tag_min):
    """Follow one actual child born on an older host; other subtrees remain live."""
    assert 13<=tag_min<=host_cap and 4<=kcap<=host_cap and tip_cap>=2 and child_cap>=2
    grid=10**21
    def down(x):return Q(x*grid//1,grid)

    def bellman(a,b,known,base,advance,rmax):
        assert a>=0 and a+b*rmax>=0
        v=down(min((a+known)/base,(a+known+rmax*(b+advance))/(base+rmax)))
        assert 0<=v and base*v<=a+known and (base+rmax)*v<=a+known+rmax*(b+advance)
        return v

    @lru_cache(maxsize=30000)
    def reward(host,child,tips):
        D,k=host; s=(D-k-1,k,1); c=(child,len(tips),1)
        # D>=13 and k>=4: entire fringe has size>12 and at least five cores.
        # Host is outside all Level-4 receiver domains, so its unobserved
        # parent & other child sums share one live (unrouted) coefficient.
        a=b=Q(0)
        cells=[('hV',(s,c,(l,0,1)),sum(Q(1,j+1) for q,j in enumerate(tips) if q!=i))
               for i,l in enumerate(tips)]
        if len(tips)==2 and c in E.c.yt:
            cells.append(('Yt',(s,c,*((l,0,1) for l in tips)),Q(0)))
        for f,states,siblings in cells:
            zero=[(Q(0),Q(0))]*len(states)
            zero[1]=(Q(0),siblings)
            a+=E.charge(f,states,zero)
            assert not E.moved(f,states,0,True) and not E.moved(f,states,0,False)
            b+=E.sign*E.terms(f,states)[1][0]/E.c.S
        return a,b

    @cache
    def V(D,k,l,tips):
        if D>host_cap or k>kcap or l+len(tips)+1>child_cap or any(j+1>tip_cap for j in tips):return Q(0)
        assert D>=tag_min>=13 and k>=4 and D>k
        dc=l+1+len(tips); hl=D-k-1; bh=Q(hl,D); bc=Q(l,dc)
        gt=tuple(Q(1,j+1) for j in tips)
        gh=hl+Q(1,dc); gc=l+Q(1,D)+sum(gt)
        a,b=reward((D,k),l,tips)
        adv=V(D+1,k,l,tips)
        total=gh*adv+bh*V(D,k+1,l,tips)+gc*V(D,k,l+1,tips)
        if l and len(tips)<2:total+=bc*V(D,k,l-1,tuple(sorted((*tips,1))))
        # Promotions by a tip, and a third middle-child promotion, kill this
        # projection; their rates remain in its diagonal.
        rates=gh+bh+gc+bc
        for i,j in enumerate(tips):
            gain=j+Q(1,dc); rates+=gain+Q(j,j+1)
            total+=gain*V(D,k,l,tuple(sorted(label(tips,i,j+1))))
        return bellman(a,b,total,HI+rates,adv,Q(k,2))

    # Root-only physical occupancy of all birth starts. One charge per
    # actual tagged promotion, no assumption on sibling count or lifetimes.
    w=defaultdict(Q); m=defaultdict(Q); best=(Q(0),None)
    for D in range(host_cap,1,-1):
        for k in range(min(kcap-1,D-1),-1,-1):
            l=D-k-1; bp=Q(l,D)
            selected=l and D>=tag_min and k>=3
            price=V(D,k+1,1,()) if selected else Q(0)
            if price>best[0]:best=(price,(D,k+1))
            rmax=Q(k+1,2)
            w[D,k]=bellman(bp*price,0,l*w[D+1,k]+bp*w[D,k+1],HI+l+bp,w[D+1,k],rmax)
            m[D,k]=bellman(bp*int(bool(selected)),0,l*m[D+1,k]+bp*m[D,k+1],
                           HI+l+bp,m[D+1,k],rmax)
    val=LO*w[2,0]
    print('tagged D,k,child D,tip D caps',host_cap,kcap,child_cap,tip_cap,'host min',tag_min,
          'contexts',V.cache_info().currsize,'birth intens floor',float(LO*m[2,0]),flush=True)
    print('tagged best conditional lifetime',float(best[0]),'post-birth D,k',best[1],flush=True)
    print('tagged certified higher reward',str(val),float(val),'lambda improvement using supplied denominator',
          float(val/DEN[side]),flush=True)
    return val


def older_path(E,side,host_cap,kcap,child_cap,tip_cap):
    """All host and middle birth starts; track a selected tip and its virgin census."""
    assert 13<=host_cap<=64 and 4<=kcap<=17 and 2<=child_cap<=10 and 2<=tip_cap<=12
    grid=10**22

    def solve(a,base,rates):
        vs=tuple(v for _,v in rates)
        corners=tuple((a+sum(r*x for r,x in zip(rs,vs)),base+sum(rs))
                      for rs in product(*(set(r) for r,_ in rates)))
        n,d=min(corners,key=lambda p:p[0]/p[1])
        v=Q((n*grid/d)//1,grid)
        assert v>=0 and all(y*v<=x for x,y in corners)
        return v

    @lru_cache(maxsize=40000)
    def price(D,k,d,j,e,m,v):
        # Reserve the D18--64,k>=5 parent band for old-flow; 65+ is
        # outside the host cut regardless of the other row coordinates.
        if D<13 or k<3 or e==0 or (D>=18 and k>=5):return Q(0)
        states=((D-k-1,k,1),(d-j-1,j,1),(e-m-1,m,1))
        # No virgin in the paid occurrence (host/middle are nonvirgin
        # structurally). Four cores below the host are guaranteed by k>=3
        # and the two-edge path: total >=5. No constraint on global size.
        if states[-1]==NB:return Q(0)
        # Condition on a directly counted sender population at the tip;
        # unresolved sends in all three slots cover their full rate ranges.
        a=E.charge('hV',states,[(Q(0),Q(0)),(Q(0),Q(0)),(Q(0),Q(v,2))])
        gains=E.terms('hV',states)[1]; n=(k-1,j-1,m-v)
        for i,hi in enumerate(n):
            for par,b in ((True,int(i==0)),(False,hi)):
                if not b:continue
                delta=gains[i]-(E.transfer('hV',states,i) if E.moved('hV',states,i,par) else 0)
                a+=min(0,E.sign*delta*Q(b,2)/E.c.S)
        assert a>=0
        return a

    @cache
    def V(D,k,d,j,e,m,v):
        if D>host_cap or k>kcap or d>child_cap or e>tip_cap:return Q(0)
        L=D-k-1; l=d-j-1; h=e-m-1 if e else 0
        assert min(L,l,h)>=0
        b=Q(L,D); c=Q(l,d); s=Q(h,e) if e else 0
        loss=v*(Q(3,2)+Q(1,e)) if e and v else 0
        a=price(D,k,d,j,e,m,v)
        fixed=b*V(D,k+1,d,j,e,m,v) if b else Q(0)
        if c:fixed+=c*(V(D,k,d,j+1,e,m,v) if e else V(D,k,d,1,2,0,0))
        if s:fixed+=s*V(D,k,d,j,e,m+1,v+1)
        if loss:fixed+=loss*V(D,k,d,j,e,m,v-1)
        g0=L+Q(1,d); g1=l+Q(1,D)+Q(1,e)
        rates=[((g0,g0+Q(k,2)),V(D+1,k,d,j,e,m,v)),
               ((g1,g1+Q(j-1,2)),V(D,k,d+1,j,e,m,v))]
        if e:
            g2=h+Q(1,d)+Q(v,2)
            rates.append(((g2,g2+Q(m-v,2)),V(D,k,d,j,e+1,m,v)))
        out=solve(a+fixed,HI+b+c+s+loss,rates)
        return out

    @cache
    def A(D,k,d,j):
        if D>host_cap or k>kcap or d>child_cap:return Q(0)
        L=D-k-1; l=d-j-1; b=Q(L,D); c=Q(l,d)
        f=b*A(D,k+1,d,j) if b else Q(0)
        if c:f+=c*(A(D,k,d,j+1)+V(D,k,d,j+1,2,0,0))
        g0=L+Q(1,d); g1=l+Q(1,D)
        return solve(f,HI+b+c,[((g0,g0+Q(k,2)),A(D+1,k,d,j)),
                              ((g1,g1+Q(j,2)),A(D,k,d+1,j))])

    w=defaultdict(Q); strongest=(Q(0),None)
    for D in range(host_cap,1,-1):
        for k in range(min(D-1,kcap-1),-1,-1):
            l=D-k-1; b=Q(l,D)
            c=A(D,k+1,2,0) if l else Q(0)
            if c>strongest[0]:strongest=c,(D,k+1)
            w[D,k]=solve(b*(c+w[D,k+1]),HI+b,[((Q(l),l+Q(k+1,2)),w[D+1,k])])
    ans=LO*w[2,0]
    print('older hV no virgin',str(ans),float(ans),'endpoint move',float(ans/DEN[side]),
          'contexts',V.cache_info().currsize,A.cache_info().currsize,
          'strongest conditional',float(strongest[0]),strongest[1],flush=True)
    return ans


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('blob',type=Path)
    p.add_argument('--n',type=int,default=1800);p.add_argument('--reps',type=int,default=4)
    p.add_argument('--seed',type=int,default=122419); p.add_argument('--verify-hash',action='store_true')
    p.add_argument('--tagged',action='store_true'); p.add_argument('--host-cap',type=int,default=36)
    p.add_argument('--host-min',type=int,default=13); p.add_argument('--tip-cap',type=int,default=6)
    p.add_argument('--child-cap',type=int,default=8); p.add_argument('--kcap',type=int,default=7)
    p.add_argument('--older-path',action='store_true')
    a=p.parse_args(); assert __debug__ and a.n>=3 and a.reps>=0
    c=Cert(a.blob); side='upper' if c.mode==1 else 'lower'
    assert (c.D,c.K,c.S,c.lam,c.stored)==(18,17,10**11,Q('.565775' if side=='upper' else '.565673'),436886484)
    if a.verify_hash:
        with a.blob.open('rb') as fh:digest=hashlib.file_digest(fh,'sha256').hexdigest()
        assert digest==PAYLOAD_SHA[side]
        print(side,'pinned SHA-256',digest,flush=True)
    E=Local(c);print(side,flush=True);validate(E)
    if a.reps:census(E,side,a.n,a.reps,a.seed)
    if a.tagged:tagged(E,side,a.host_cap,a.tip_cap,a.child_cap,a.kcap,a.host_min)
    if a.older_path:older_path(E,side,a.host_cap,a.kcap,a.child_cap,a.tip_cap)


if __name__=='__main__':main()
