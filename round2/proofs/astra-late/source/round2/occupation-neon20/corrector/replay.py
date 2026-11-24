#!/usr/bin/env python3
"""D20 signed census and existing-edge correctors with actual promotion sources."""
import argparse
from hashlib import sha256
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json
import sys
from importlib.util import spec_from_file_location,module_from_spec

HERE=Path(__file__).resolve().parent
PINS={
 'rows/extract.py':'f9352f511c7eb38f94aee4cf6bbce4a3286ae0a4dcb13f577578773ab4dbf1fe',
 'rows/lower.json':'b59a2c20b31e7b09f88df436e94e75921d22ac424a759a58c3caf118f7692d4d',
 'rows/upper.json':'eeb67cf94792f0a6372f72b02214fa07656a166a65935997d06703a8d9ed85ad',
 'renewal/replay.py':'c72fc793e1aaff19d86d22fc72694e301dfa163b0d311d484217e97cb8bf4257',
}
for name,digest in PINS.items():
    if sha256((HERE.parent/name).read_bytes()).hexdigest()!=digest:
        raise ValueError(f'D20 source/export changed: {name}')
spec=spec_from_file_location('d20_renewal',HERE.parent/'renewal'/'replay.py')
renewal=module_from_spec(spec);spec.loader.exec_module(renewal)
Fringe,root=renewal.Fringe,renewal.root
from extract import COMMIT, SOURCE, HASH, LO, HI, ROOT, Cert, Rows, Charge, children, gain, bud, replace, NB, embeddings, pick, audit

Z=Q(0)
J=(2,0)
FAMS={'hV','hC','P4','Yt','Yb','S3'}
GRID=10**21


def down(x):
    return Q(x*GRID//1,GRID)


def show(x,sgn=1):
    n=sgn*(sgn*x*10**15//1)
    k,r=divmod(abs(n),10**15)
    return f'{"-" if n<0 else ""}{k}.{r:015d}'


def states(degree,arity=None,kmin=0):
    return [(d,k) for d in range(degree,1,-1)
            for k in range(d-1 if arity is None else min(arity,d-1),kmin-1,-1)]


def label(x):
    d,k=x
    return d-k-1,k,1


def moves(x):
    d,k=x
    return (d+1,k),(d,k+1)


class Edge:
    def __init__(self,rows,degree,arity,preserve):
        self.data={}
        for p,w in product(states(degree,arity,1),states(degree,arity)):
            L,k,h=label(p); M,j,_=label(w)
            dp,dw=p[0],w[0]
            # Minimum planted size consistent with both adjacent marks:
            # two cores, their leaves, and two vertices per off-edge child.
            size=2+L+M+2*(k-1+j)
            pay=not preserve or size>16 or k+j>=4
            a=b=e=Z
            if pay:
                a,terms=rows.affine('g',(label(p),label(w)))
                slopes=[set(),set()]
                for i,_,_,s in terms:
                    slopes[i].add(s)
                assert max(map(len,slopes))<=1
                b,e=(next(iter(s),Z) for s in slopes)
            pb,wb=Q(L,dp),Q(M,dw)
            pg,wg=L+Q(1,dw),M+Q(1,dp)
            box=(Q(k,2),Q(j,2))
            pts=list(product(*(tuple(sorted({0,t})) for t in box)))
            # Nonnegative production g permits the hinge. Check all box
            # corners and every intersection with the hinge hyperplane.
            for axis,coeff in enumerate((b,e)):
                if coeff:
                    for t in (0,box[1-axis]):
                        u=-(a+(e,b)[axis]*t)/coeff
                        if 0<u<box[axis]:
                            pts.append((u,t) if axis==0 else (t,u))
            self.data[p,w]=(pb,wb,tuple((pg+u,wg+v,max(Z,a+b*u+e*v)) for u,v in pts))

    def solve(self,s,flat=False):
        V={}
        for (p,w),(a,b,pts) in self.data.items():
            pg,pp=moves(p); wg,wp=moves(w)
            z=a*V.get((pp,w),Z)+b*V.get((p,wp),Z)
            pval=V.get((pg,w),Z);wval=V.get((p,wg),Z)
            floor=min(r for g,h,r in pts)
            pairs=[(z+g*pval+h*wval+(floor if flat else r),s+a+b+g+h) for g,h,r in pts]
            v=down(min(n/d for n,d in pairs))
            assert v>=0 and all(v*d<=n for n,d in pairs)
            V[p,w]=v
        return V


class Census:
    def __init__(self,rows,data,charge=None):
        self.pay={}
        c=rows.c
        assert (data['source_commit'],data['source_blobs'],data['payload_sha256'],data['root_ticks'])==(
            COMMIT,SOURCE,HASH[c.side],list(ROOT[c.side]))
        assert Q(data['lambda'])==c.lam
        for l,k,h,a,b,u,v in data['f_rows']:
            key=(l+k+h,k)
            assert h==1 and key not in self.pay
            A,B=rows.first_body((l,k,h))
            assert (A,B,A,A+B*Q(k+1,2))==tuple(map(Q,(a,b,u,v))) and min(Q(u),Q(v))>=0
            self.pay[key]=(A,B)
        self.degree=max(d for d,k in self.pay)
        assert self.pay.keys()==set(states(self.degree))
        self.bands=[Q(data['bands'][str(k)]['floor']) for k in range(len(data['bands']))]
        assert min(self.bands)>=0
        rc=data['recharge']; self.den=Q(rc['ceiling']);self.M=Q(rc['cap_ticks'])/c.S
        assert self.M-Q(rc['stationary_deficit'])==self.den
        caps={(h,k):Q(v) for h,k,v in data['charge_by_h_k_ticks']}; t=Q(rc['recharge_ticks'])
        assert self.M==max(a+t*(k-h) for (h,k),a in caps.items())/c.S
        self.K=c.KF
        if charge is not None:
            assert charge.by_k==caps and charge.cap(t)==self.M*c.S
        self.tail=[self.M-(caps[1,k]+t*(k-1))/c.S for k in range(self.K+2)]
        self.deficit={(d,k):(self.M-(charge.exact(*label((d,k)))+t*(k-1))/c.S if charge else self.tail[k])
                      for d,k in self.pay}
        assert min(*self.tail,*self.deficit.values())>=0

    def solve(self,s,E,delta=Z,flat=False):
        W={};tail={};H=self.degree
        # Capped stock classes are monotone; overflow uses its conservative
        # k=KF+1 cap, and the nonpositive recharge only helps beyond it.
        for k in range(self.K+1,-1,-1):
            r=(self.bands[k] if k<len(self.bands) else 0)+delta*self.tail[k]
            rates=(max(Z,1-Q(k+1,H+1)),Q(1)) if k<=self.K else (Z,)
            vals=[(r+b*tail.get(k+1,Z),s+b) for b in rates]
            v=down(min(a/b for a,b in vals))
            assert all(v*b<=a for a,b in vals)
            tail[k]=v
        for x in states(H):
            d,k=x; l=d-k-1; p,g=(d,k+1),(d+1,k); a,e=self.pay[x]
            a+=delta*self.deficit[x]
            if flat:
                a+=min(0,e*Q(k+1,2));e=0
            promotion=Q(l,d)
            # Born edge starts at the post-promotion census; singleton
            # birth gets zero after recentering W by W(J).
            a+=promotion*(W.get(p,Z)+E.get((p,J),Z))
            w=W[g] if d<H else tail[k]
            candidates=[(a+r*e+(l+r)*w,s+l+r+promotion) for r in (0,Q(k+1,2))]
            v=down(min(a/b for a,b in candidates))
            assert v>=0 and all(v*b<=a for a,b in candidates)
            W[x]=v
        return W,tail


def controls(rows,cs,s,E,W,tail,delta,margin,replace_edges):
    c=rows.c; cases=[]
    front=[(-1,0,0)]
    for n in range(3,7):
        cases+=front
        front=[p+(j,) for p in front for j in range(n)] if n<6 else []
    for ls in ((1,3,1,1),(1,17,1,2),(1,18,19,1),(1,65,1,3),(1,16001,1,1)):
        cases.append(audit.build((-1,0,1,1),ls))
    births=0; detected=0; full=0
    def qval(x):
        d,k=sum(x),x[1]
        return W.get((d,k),tail[min(k,c.KF+1)])
    def eval_correction(p,x):
        census=[(sum(z),z[1]) for z in x]
        return sum(qval(z)-W[J] for z in x[1:])+sum(E.get((census[p[i]],census[i]),Z) for i in range(1,len(x)) if p[i]>0)
    for par in cases:
        ch,x,core=audit.physical_view(par); x=tuple(x); ids={v:i for i,v in enumerate(core)}
        p=(-1,*(ids[par[v]] for v in core[1:]))
        base=eval_correction(p,x); physical=-s*base; local=-s*base; born=Z
        for i,y in enumerate(x):
            l,k,h=y;d=sum(y);b=Q(l,d)
            g=l+sum(Q(1,sum(x[j])) for j in ch[i]+([p[i]] if i else []))
            physical+=g*(eval_correction(p,replace(x,i,gain(y)))-base)
            if l:
                physical+=b*(eval_correction(p+(i,),replace(x,i,bud(y))+(NB,))-base)
            if not i:
                continue
            local+=g*(qval(gain(y))-qval(y))
            if l:
                cost=b*E.get(((d,k+1),J),Z)
                local+=b*(qval(bud(y))-qval(y))+cost
                born+=cost;births+=bool(cost)
            for j in ch[i]+([p[i]] if p[i]>0 else []):
                def v(z):
                    a,b=(z,x[j]) if p[j]==i else (x[j],z)
                    return E.get(((sum(a),a[1]),(sum(b),b[1])),Z)
                local+=g*(v(gain(y))-v(y))
                if l:
                    local+=b*(v(bud(y))-v(y))
        assert physical==local
        detected+=born!=0 and physical!=local-born
        # Every component in the independent residual is inventoried;
        # eligible edge charges are evaluated on the labelled physical tree.
        occupied=Z
        for f,v in embeddings(p):
            if f=='f' and v[0]:
                j=v[0];l,k,h=x[j];d=l+k+h
                a,e=rows.first_body(x[j]);r=sum(Q(1,sum(x[z])) for z in ch[j]+[p[j]])
                occupied+=a+e*r
            if f=='g' and v[0]:
                i,j=v; l,k,_=x[i]; m,n,_=x[j]
                if replace_edges or 2+l+m+2*(k-1+n)>16 or k+n>=4:
                    a,ports=rows.affine(f,pick(x,v))
                    for at,direction,count,slope in ports:
                        t=v[at]
                        donors=[p[t]] if direction else [q for q in ch[t] if q not in v]
                        assert len(donors)==count
                        a+=slope*sum(Q(1,sum(x[q])) for q in donors)
                    assert a>=0
                    occupied+=a
        slack=occupied+physical
        if delta:
            # The conditional census deficit is the exact envelope when
            # supplied, and otherwise the valid all-stock cap relaxation.
            for y in x[1:]:
                d,k=sum(y),y[1]
                slack+=delta*cs.deficit.get((d,k),cs.tail[min(k,c.KF+1)])
        assert slack>=margin*(len(x)-1)
        old=c.potential(ch,x)
        drift=Z
        for i,y in enumerate(x):
            l,k,h=y;d=sum(y)
            g=l+sum(Q(1,sum(x[j])) for j in ch[i]+([p[i]] if i else []))
            drift+=g*(c.potential(ch,replace(x,i,gain(y)))-old)
            if l:
                kids=[v[:] for v in ch]+[[]];kids[i].append(len(x))
                drift+=Q(l,d)*(c.potential(kids,replace(x,i,bud(y))+(NB,))-old)
        production=c.sign*(drift-c.lam*old)/c.S
        assert production>=occupied
        full+=1
    print('labelled fixtures',len(cases),'full production checks',full,'nonzero births',births,
          'omitted-birth identity failures',detected,flush=True)
    assert births and detected


def run(side,args):
    c=Cert(getattr(args,side),side); rows=Rows(c)
    print('authenticated',side,HASH[side],flush=True)
    data=json.loads((HERE.parent/'rows'/f'{side}.json').read_text())
    cs=Census(rows,data,Charge(c) if args.exact_den else None)
    f0=root(lambda d,k:cs.pay[d,k],cs.degree,cs.bands)
    fringe=Fringe(rows,16,4)
    reserve=fringe.solve(FAMS|{'g'});higher=fringe.solve(FAMS)
    orig=(f0+reserve)/cs.den
    assert show(c.lam+c.sign*orig,c.sign)=={'lower':'0.565681716887067','upper':'0.565766354357632'}[side]
    print(side,'renewal offset',show(orig),'finite joint',show(reserve),'higher-only',show(higher),'den',show(cs.den,-1),flush=True)
    if args.stationary:
        from stationary import Stationary
        lp=Stationary(args.stationary)
        print('stationary incidence model',len(lp.names),'rows',len(lp.cols),'columns',flush=True)
        r,repair=lp.bound(cs,{},Z)
        print(side,'stationary LP first-body only',show(r),'raw offset',show((r+reserve)/cs.den),
              'repair',show(repair,-1),flush=True)
    for mode in args.mode or ('preserve','replace'):
        keep=reserve if mode=='preserve' else higher
        edge=Edge(rows,args.edge_degree,args.arity,mode=='preserve')
        E=edge.solve(HI); W,T=cs.solve(HI,E)
        stat=(LO*W[J]+keep)/cs.den
        Wf,_=cs.solve(HI,edge.solve(HI,flat=True),flat=True)
        print(mode,'edge states',len(E),'positive',sum(v>0 for v in E.values()),'max',show(max(E.values()),-1),
              'flat stationary offset',show((LO*Wf[J]+keep)/cs.den),
              'direct stationary offset',show(stat),'inward',show(stat-orig),flush=True)
        def margin(delta,kind):
            s=c.lam+c.sign*delta if kind=='point' else HI
            e=edge.solve(s) if kind=='point' else E
            w,t=cs.solve(s,e,delta)
            return (s if kind=='point' else LO)*w[J]-delta*cs.M+keep,e,w,t,s
        for kind in ('point','birth-bellman'):
            l=Q(stat*10**13//1,10**13)
            # Refresh the affordable trial from the signed residual;
            # every finally reported rational is independently checked.
            for _ in range(args.iterations):
                r,*_=margin(l,kind)
                l=Q((l+r/cs.M)*10**13//1,10**13)
            for _ in range(10):
                if margin(l,kind)[0]>=0:break
                l-=Q(1,10**13)
            else:
                raise ArithmeticError('trial failed the signed residual check')
            sl,e,w,t,s=margin(l,kind)
            assert sl>=0
            if kind=='point':
                assert max(abs(v-w[J]) for v in (*w.values(),*t.values()))+max(e.values())<c.c0
            print(side,mode,kind,'offset',show(l),'inward vs renewal',show(l-orig),
                  'outward endpoint',show(c.lam+c.sign*l,c.sign),'exact offset',l,
                  'margin',sl,flush=True)
            if kind=='point' and args.fixtures:
                controls(rows,cs,s,e,w,t,l,s*w[J],mode=='replace')
        if args.stationary:
            r,repair=lp.bound(cs,E,Z)
            plain=(r+keep)/cs.den
            target=orig+Q('1e-6')
            r1,er=lp.bound(cs,E,target)
            test=r1+keep-target*cs.M
            trial=max(plain,target) if test>=0 else plain
            trial=Q(trial*10**13//1,10**13)
            r2,er2=lp.bound(cs,E,trial)
            signed=r2+keep-trial*cs.M
            if signed>0:
                trial+=signed/cs.den
            best=max(plain,trial)
            assert signed>=0 or best==plain
            print(side,mode,'stationary LP direct offset',show(plain),'coupled offset',show(best),
                  'inward',show(best-orig),'outward endpoint',show(c.lam+c.sign*best,c.sign),
                  'target signed margin',show(test),'repair≤',show(max(repair,er,er2),-1),flush=True)
        del edge
    c.check_unchanged();fringe.charge.cache_clear();fringe.rates.cache_clear()
    rows.value.cache_clear();rows.affine.cache_clear()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lower',type=Path,default=Path('b3f20_lo.bin'))
    p.add_argument('--upper',type=Path,default=Path('b3f20_hi.bin'))
    p.add_argument('--side',choices=('lower','upper'),action='append')
    p.add_argument('--mode',choices=('preserve','replace'),action='append')
    p.add_argument('--edge-degree',type=int,default=18)
    p.add_argument('--arity',type=int,default=4)
    p.add_argument('--iterations',type=int,default=3)
    p.add_argument('--stationary',type=int,choices=(19,24,32,40,48,64))
    p.add_argument('--exact-den',action='store_true')
    p.add_argument('--fixtures',action='store_true')
    a=p.parse_args();assert __debug__ and 2<=a.arity<a.edge_degree<=64
    for side in a.side or ('lower','upper'):run(side,a)


if __name__=='__main__':main()
