#!/usr/bin/env python3
"""Gold2 live maternal stock, real birth marks and up to two latest donors."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from json import loads
from math import ceil, floor, lcm
from mmap import mmap, ACCESS_READ
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('gold2_life',HERE.parent/'check.py')
life=module_from_spec(spec);spec.loader.exec_module(life)
G=life.GRID


def div(n,d,up): return -(-n//d) if up else n//d


def better(up): return max if up else min


def dec(x,up,places=9):
    n=(ceil if up else floor)(x*10**places)
    a,b=divmod(abs(n),10**places)
    return f'{"-" if n<0 else ""}{a}.{b:0{places}d}'


def stock_split(t,rec,lo,hi,W,joint,up,overflow):
    """Price the one score after first entry into actual c>=4."""
    lim=joint['limit']; sat=joint['saturation']
    m=life.load('stock_lift',life.R2/'occupation-neon/pressure-memory/birth-law/replay.py')
    m.LO,m.HI=lo,hi; m.SAT=sat; lp=m.build(lim,sided=True)
    j=4
    def density(x):
        _,k,b,*_=x
        return b if k==j-2 else 0
    lower=lo*max(0,lp.lower(lambda k,x:density(x) if k=='occ' else 0))
    upper=-hi*lp.lower(lambda k,x:-density(x) if k=='occ' else 0)
    print('stationary first-entry flux to c=4',dec(lower,False),dec(upper,True),flush=True)
    select=better(up)
    far=[Q(z[j],G) for d,z in W.items() if d>lim]
    if up: far.append(Q(life.suffix(t,4000,70,lo,rec['cap_f'],True)[j],G))
    else: far.append(Q(0))
    end=select(far); s=-1 if up else 1
    def flux_reward(x):
        d,k,b,*_=x
        if k!=j-2: return 0
        return b*(Q(W[d][j],G) if d<=lim and d>=j else end if d>lim else 0)
    price=(hi if up else lo)*max(0,s*10**6*lp.lower(lambda k,x:s*flux_reward(x)/10**6 if k=='occ' else 0))
    print('stationary c>=4 singleton ticks via actual promotion entrance',dec(price+overflow,up),
          'overflow degree entrance value ticks',dec(end,up),flush=True)
    # The pooled overflow flux above can put its entire mass at the worst
    # degree. Instead, enforce every actual degree transition on the way in.
    disc=lo if up else hi; p,q=disc.numerator,disc.denominator
    nxt=[0]*j
    if up:
        # All-reward analytic suffix is also an upper continuation for this
        # masked reward; the discarded pre-entry reward stays nonnegative.
        nxt=life.suffix(t,4000,70,lo,rec['cap_f'],True)[:j]
    for d in range(max(W),1,-1):
        cur=[0]*j
        for c in range(min(j-1,d),0,-1):
            L=d-c; w=2*d*q
            B=(W[d][j] if c==j-1 else cur[c+1]) if L else 0
            num=L*(w*nxt[c]+2*q*B);den=2*d*p+L*(2*q+w);dx=c*d*q
            cur[c]=select(div(num,den,up),div(num+dx*nxt[c],den+dx,up))
            v=cur[c]
            assert (v*den>=num and v*(den+dx)>=num+dx*nxt[c]) if up else (
                v*den<=num and v*(den+dx)<=num+dx*nxt[c])
        nxt=cur
    z=(hi if up else lo)*Q(nxt[1],G)+overflow
    print('stationary c>=4 singleton ticks via exact degree entrance',dec(z,up),flush=True)


def live(t,cap,W,H,K,M,T,J,disc,up):
    p,q=disc.numerator,disc.denominator
    beyond=(0,0)
    mother=(beyond,)+tuple((m,h) for m in range(M,1,-1) for h in range(m,1,-1))
    tips=((T+1,0),)+tuple((d,c) for d in range(T,1,-1) for c in range(d,0,-1))
    tipid={v:i for i,v in enumerate(tips)}; born=tipid[2,1]
    # Ordered by birth, latest first; equal descriptors need no extra birth.
    tiprows={k:tuple(product(range(len(tips)),repeat=min(k,J))) for k in range(K)}
    parents={x:i for i,x in enumerate(mother)}
    nxt={}; nrows=checks=0

    for D in range(H,1,-1):
        cur={}
        for c in range(min(K,D),0,-1):
            k=c-1; L=D-c; a,z=life.row(t,D,c,cap)
            assert min(a,a+Q(c,2)*z)>=0
            a,z=[(ceil if up else floor)(x*G) for x in (a,z)]
            for im,(m,h) in enumerate(mother):
                for tt in tiprows[k]:
                    dd=lcm(2,D,m or M+1,*(tips[j][0] for j in tt)); scale=dd*q; f=scale//D
                    known=scale//m if m else 0
                    rem=(k-len(tt))*(scale//2)+(0 if m else scale//(M+1))
                    for j in tt:
                        d0,h0=tips[j]
                        if h0: known+=scale//d0
                        else: rem+=scale//d0
                    Wg=nxt[c,im,tt] if D<H else W[D+1][c]
                    ttnew=(born,*tt)[:J]
                    Wb=(cur[c+1,im,ttnew] if c<K else W[D][c+1]) if L else 0
                    rate=L*scale+known
                    num=a*scale+z*known+rate*Wg+L*f*Wb
                    den=dd*p+rate+L*f
                    ports=[(z+Wg,rem)]
                    if m:
                        g=cur[c,parents[(m+1,h) if m<M else beyond],tt]
                        r=(m-h)*scale+f
                        num+=r*g;den+=r;ports.append((g,(h-1)*(scale//2)))
                        if m>h:
                            r=(m-h)*(scale//m)
                            b=cur[c,parents[m,h+1],tt]; num+=r*b;den+=r
                    for ix,j in enumerate(tt):
                        d0,h0=tips[j]
                        if not h0: continue
                        def change(new): return tt[:ix]+(tipid[new],)+tt[ix+1:]
                        g=cur[c,im,change((d0+1,h0) if d0<T else (T+1,0))]
                        r=(d0-h0)*scale+f
                        num+=r*g;den+=r;ports.append((g,(h0-1)*(scale//2)))
                        if d0>h0:
                            r=(d0-h0)*(scale//d0)
                            b=cur[c,im,change((d0,h0+1))];num+=r*b;den+=r
                    n,d=num,den
                    for x,r in sorted(ports,reverse=up):
                        if (x*d>n if up else x*d<n): n+=x*r; d+=r
                        else: break
                    v=div(n,d,up)
                    # Check the worst signed residual over the entire box.
                    rs=num-den*v+sum((max if up else min)(0,r*(x-v)) for x,r in ports)
                    assert (rs<=0 if up else rs>=0) and v>=0
                    if D in (2,3,t.d,H) and c<=3:
                        for bits in product((0,1),repeat=len(ports)):
                            rs=num-den*v+sum(b*r*(x-v) for b,(x,r) in zip(bits,ports))
                            assert rs<=0 if up else rs>=0
                            checks+=1
                    cur[c,im,tt]=(min if up else max)(v,W[D][c]);nrows+=1
        nxt=cur
    # Mother at actual focal birth has post-arity c_pre+1; no post-bare row
    # is created by a nonpromotion, and D is unchanged by her promotion.
    return {x:Q(nxt[1,i,()],G) for i,x in enumerate(mother)},nrows,checks


def mark(M,lo,hi,roots,up):
    """Physical birth-law contraction using one source and its actual marks."""
    # beta K_mu = beta; only beta K_mu 1=1 is required. The expected
    # discounted number of children need not be one in each incoming class.
    labels=tuple(sorted(roots)); mats={}
    for label in labels:
        m,_=label; cap=Q(1,m or M+1)
        for target in labels:
            if target==(0,0): continue
            vals={}
            for d in range(M,1,-1):
                for c in range(d,0,-1):
                    L=d-c; b=Q(L,d)
                    num=b*(int(target==(d,c+1))+vals.get((d,c+1),0))
                    gg=vals.get((d+1,c),0)
                    z=better(up)((num+(L+r)*gg)/((lo if up else hi)+b+L+r)
                                 for r in ((cap if m else 0),cap+Q(c-1,2)))
                    vals[d,c]=Q((ceil if up else floor)(z*G),G)
            mats[label,target]=vals[2,1]
    return mats


def pooled(M,lo,hi,roots,up):
    lower=mark(M,lo,hi,roots,False)
    upper=mark(M,lo,hi,roots,True)
    # Stationary birth-law intervals by simple uniformity over incoming
    # marks. Use these as marginal bounds with common total mass one.
    labels=tuple(sorted(roots)); exact=labels[1:]
    low={x:min(lower[y,x] for y in labels) for x in exact}
    high={x:max(upper[y,x] for y in labels) for x in exact}
    tail_lo,tail_hi=max(0,1-sum(high.values())),1-sum(low.values())
    low[0,0],high[0,0]=tail_lo,tail_hi
    assert sum(low.values())<=1<=sum(high.values()) and all(low[x]<=high[x] for x in labels)
    avail=1-sum(low.values()); val=sum(low[x]*roots[x] for x in labels)
    for x in sorted(labels,key=roots.__getitem__,reverse=up):
        amt=min(avail,high[x]-low[x]); val+=amt*roots[x];avail-=amt
    assert avail==0
    return val,{x:(low[x],high[x]) for x in labels}


def joint_price(joint,t,rec,roots,W,lo,hi,overflow):
    """Birth-mark/host joint trial and exact sided finite singleton residual."""
    r=life.R2
    m=life.load('signed_pressure_lp',r/'occupation-neon/pressure-memory/birth-law/replay.py')
    m.LO,m.HI=lo,hi; m.SAT=joint['saturation']; limit=joint['limit']
    lp=m.build(limit,sided=True)
    cells={tuple(map(int,k.split(','))):tuple(map(int,v)) for k,v in joint['cells'].items()}
    shift=Q(joint['signed_shift']); assert shift>0 and shift==Q(99,175)-t.lam
    M=rec['m_max']; H=rec['higher_per_core']; cap=rec['cap_f']
    a={ (d,c):life.row(t,d,c,cap) for d in range(2,limit+1) for c in range(1,d+1)}
    unit=10**6
    def local(x):
        d,k,b,p,u,v=x; c=k+1
        if d>limit: return Q(0)
        val,slope=a[d,c]; rho=u+v
        # Stop on gain out of the finite box. This is the actual gain
        # using the same rho as the original signed singleton reward.
        return val+slope*rho+((d-c+rho)*Q(W[d+1][c],G) if d==limit else 0)
    def residual(x):
        d,k,b,p,u,v=x; c=k+1
        if d>limit: return Q(0)
        z=Q(W[d][c],G); gain=d-c+u+v
        zg=Q(W[d+1][c],G) if d<limit else 0
        zb=Q(W[d][c+1],G) if b else 0
        # Dynkin with the original reward and the live receive, not a
        # separate adversarial corner for the potential's gain increment.
        return local(x)+gain*(zg-z)+b*(zb-z)-(lo if up else hi)*z+b*Q(W[2][1],G)
    def price(fn,up=False):
        s=-1 if up else 1
        v=s*unit*lp.lower(lambda k,x:s*fn(x)/unit if k=='occ' else 0)
        return (hi if up else lo)*max(0,v)
    def birth(x):
        d,k,b,*_=x
        return b*roots.get((d,k+2),roots[0,0])
    def rebate(x):
        d,k,*_=x; _,U,h=cells[d,k]
        assert M>=U and H>=h
        return H-h+shift*(M-U)
    up=overflow>0
    v=price(local,up)+overflow
    print('sided finite singleton stationary ticks',life.fmt(v,up,9),flush=True)
    b=price(birth,up)+overflow
    print('sided birth-mark singleton ticks',life.fmt(b,up,9),flush=True)
    if not up:
        for name,fn in [('marked',birth),('finite',local),('Dynkin residual',residual)]:
            gap=rec['err_per_core']+shift*M-price(lambda x:fn(x)+rebate(x))
            print('joint',name,'trial gap before disjoint edges',dec(gap,True),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--live-cut',type=int,default=500);p.add_argument('--live-arity',type=int,default=12)
    p.add_argument('--mother',type=int,default=3);p.add_argument('--tip',type=int,default=2)
    p.add_argument('--donors',type=int,default=2);p.add_argument('--side',choices=('lower','upper'),default='lower')
    p.add_argument('--upper-control',action='store_true')
    p.add_argument('--trial-discount',action='store_true',help='conditional at mu=99/175, instead of uniform over Gold2 prior')
    p.add_argument('--joint',action='store_true',help='use sided D24 occupation and original trial cap terms')
    p.add_argument('--stock-flux',action='store_true',help='bound physical promotion flux into c=4 using the same occupation constraints')
    a=p.parse_args();H,K,M,T,J=a.live_cut,a.live_arity,a.mother,a.tip,a.donors
    assert 2<=H<3999 and 1<=K<=30 and 2<=M<=8 and 2<=T<=4 and 1<=J<=2
    assert (life.tail.BASE/'hashes.log').read_text().splitlines()==['quintet_lower.bin: OK','quintet_upper.bin: OK']
    gold=loads(life.tail.RECEIPT.read_text());lo,hi=[Q(gold[x]['certified_exact']) for x in ('lower','upper')]
    side=a.side;up=a.upper_control;rec=gold[side];target=Q(99,175)
    prior=lo,hi
    if a.trial_discount:
        assert lo<=target<=hi
        lo=hi=target
    disc=lo if up else hi
    with (life.tail.BASE/f'quintet_{side}.bin').open('rb') as fp,mmap(fp.fileno(),0,access=ACCESS_READ) as raw:
        assert len(raw)==life.tail.SIZE
        t=life.tail.Snapshot(raw,rec);W={};C=70 if up else 30
        stop=3999 if up else 7999
        v,n=life.price(t,stop,C,disc,rec['cap_f'],up,(stop if a.stock_flux else max(H+1,25),W))
        roots,rows,checks=live(t,rec['cap_f'],W,H,K,M,T,J,disc,up)
        bound,marg=pooled(M,lo,hi,roots,up)
        assert bound<=v if up else bound>=v
        overflow=130000*2**40/(1+lo)**C if up else Q(0)
        reward=(hi if up else lo)*bound+overflow
        print(side,'ceiling' if up else 'floor','rows/checks',rows,checks,'roots',[(x,life.fmt(z,up,6)) for x,z in roots.items()],flush=True)
        print('birth marginals',[(x,life.fmt(l,False,9),life.fmt(h,True,9)) for x,(l,h) in marg.items()])
        print('mixed lifetime exact',bound,'stationary ticks',life.fmt(reward,up,9),'broad rows',n,flush=True)
        if side=='lower':
            required=rec['err_per_core']+rec['m_max']*(target-t.lam)
            potential=rec['err_per_core']+rec['m_min']*(target-t.lam)
            print('required original cap ticks',life.fmt(required,True,9),
                  'required minus singleton',dec(required-reward,not up),
                  'required at minimum potential ticks',life.fmt(potential,False,9),flush=True)
            raw_joint=(life.R2/'occupation-quintet-gold2/host/joint24.json').read_bytes()
            assert sha256(raw_joint).hexdigest()=='104ad923d984e68a185ac5ca8678a2305bb9aa31fb510470a5fea27df514c99a'
            joint=loads(raw_joint)['results']['lower']
            assert (joint['target'],joint['nominal'],joint['cap_error'],tuple(map(Q,joint['prior'])))==(
                str(target),str(t.lam),rec['err_per_core'],prior)
            assert Q(joint['edge_floor'])==0 and joint['reduced_f']
            # Joint24's rebate includes a different, broad f residual. Use
            # its potential envelope alone with the new full f replacement.
            m_upper=Q(joint['denominator'][1]); threshold=rec['err_per_core']+m_upper*(target-t.lam)
            print('D24 M upper only, threshold ticks',life.fmt(threshold,True,9),
                  'threshold minus singleton',dec(threshold-reward,not up),flush=True)
            if a.joint:
                joint_price(joint,t,rec,roots,W,lo,hi,overflow)
            if a.stock_flux:
                stock_split(t,rec,lo,hi,W,joint,up,overflow)


if __name__=='__main__': main()
