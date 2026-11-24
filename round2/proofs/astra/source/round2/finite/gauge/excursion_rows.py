"""Schur elimination at every star visit with two analytic tail moments."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import ceil, floor

import binary_rows as rows
from transfer_check import A, BETA, b, info, mark, newborn, physical, plus

UNIT = 10**14
N = 5
END = 96
SWEEPS = 28
S = rows.S


def star(parent, degree):
    return parent,degree-1,()


def promoted(parent, degree):
    return parent,degree-2,(newborn(2 if degree<16 else 16),)


def successors(tree):
    rates=defaultdict(F,physical(tree))
    p=1/F(info(tree)[0])
    for child in tree[2]:
        rates[child]+=mark(child)-p
        rates[plus(child)]+=p
    return {t:r for t,r in rates.items() if r}


def retained():
    seen={newborn(),newborn(16)}
    todo=list(seen)
    while todo:
        t=todo.pop()
        for nxt in successors(t):
            if info(nxt)[1]<=N and nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    return sorted(t for t in seen if t[2])


@cache
def weight(tree):
    D,n=info(tree)
    c=1+len(tree[2])
    if rows.VT is None:
        h=rows.H[D,c] if D<64 else rows.H[64,1]
    else:
        cut=rows.virtual_weight.CUTOFF
        allowed=[t for i,t in rows.VT.items() if rows.virtual_weight.actual(i)==mark(tree)]
        h=max(t[D,c] if D<cut else t[cut,1] for t in allowed)
    return (h+sum(mark(y)*weight(y) for y in tree[2]))/(mark(tree)+rows.DAMP)


@cache
def boundary(tree,col):
    """Two controlled root residence times, retaining each old grandchild weight."""
    a=mark(tree)
    D,n=info(tree)
    c=1+len(tree[2])
    assert D<=rows.N
    lo,hi=rows.second(a,col)[D,c]
    lo,hi=F(lo,rows.UNIT),F(hi,rows.UNIT)
    for y in tree[2]:
        ell=mark(y)
        dy,_=info(y)
        cy=1+len(y[2])
        lv,uv=rows.child(a,ell,col,True)[dy,cy]
        X=sum(mark(z)*weight(z) for z in y[2])
        coeff=rows.KAPPA*rows.alpha(col)/((S+ell-rows.Q)*(S+a-rows.Q))
        lo+=ell*F(lv,rows.UNIT)
        hi+=ell*(F(uv,rows.UNIT)+coeff*X)
    return lo,min(hi,rows.alpha(col)*weight(tree))


def ratio(lam,k):
    return F(k*k)/(k*k+(2+lam)*k+lam)


@cache
def moments(lam,K):
    """T and M relative to unit mass at the first tail degree K."""
    prefix_t=F(0)
    prefix_m=F(0)
    mass=F(1)
    for k in range(K,END):
        prefix_t+=mass
        prefix_m+=mass/k
        mass*=ratio(lam,k)
    gamma=lam+2/(lam+2)
    numer=(END+gamma)*mass
    tlo=numer/(1+lam)
    thi=numer/(1+lam-lam**2/((lam+2)*END**2))
    # Sum the second moment by telescoping a constant multiple of q_k.
    # (q_k-q_{k+1})/(q_k/k) = k*((2+lam)*k+lam)/(k*k+(2+lam)*k+lam).
    # The ratio increases in k to 2+lam.
    mlo=mass/(2+lam)
    mhi=mass/(END*(1-ratio(lam,END)))
    return tuple(F(fn(UNIT*x),UNIT) for x,fn in (
        (prefix_t+tlo,floor),(prefix_t+thi,ceil),
        (prefix_m+mlo,floor),(prefix_m+mhi,ceil)))


@cache
def tail_coefficients(a,col):
    """Affine lower and upper payoffs on a high first promotion, direct reward included."""
    lam=S+a
    low,high=rows.child(a,BETA,col,True)[2,1]
    old_lo=BETA*F(low,rows.UNIT)
    old_hi=BETA*F(high,rows.UNIT)
    direct=rows.g(col,rows.N)
    toll_lo,toll_hi=direct+old_lo,direct+old_hi
    return (direct+old_lo+toll_lo/lam,2*toll_lo/lam,
            direct+old_hi+toll_hi/lam,2*toll_hi/(lam+1))


@cache
def star_row(parent,j,col):
    """Constant lower/upper sources plus identical finite multicore coefficients."""
    a=F(1,parent)
    lam=S+a
    K=max(rows.N,N,j)
    assert K<=END
    tlo,thi,mlo,mhi=moments(lam,K)
    intercept_lo,slope_lo,intercept_hi,slope_hi=tail_coefficients(a,col)
    coeff={}
    lo,hi=F(0),F(0)
    bud=F(j-1,j)
    mass=bud/(lam+j-1+bud)
    for k in range(j,K):
        successor=promoted(parent,k)
        direct=rows.g(col,k)
        if k+1<=N:
            coeff[successor]=mass
            bl,bh=F(0),F(0)
        else:
            bl,bh=boundary(successor,col)
        lo+=mass*(direct+bl)
        hi+=mass*(direct+bh)
        mass*=ratio(lam,k)
    lo+=mass*(intercept_lo*tlo-slope_lo*mhi)
    hi+=mass*(intercept_hi*thi-slope_hi*mlo)
    return F(floor(UNIT*lo),UNIT),F(ceil(UNIT*hi),UNIT),coeff


@cache
def dormant_boundary(tree,col):
    D,n=info(tree)
    c=1+len(tree[2])
    X=sum(mark(z)*weight(z) for z in tree[2])
    lv,uv=rows.child(BETA,A,col,False)[D,c]
    coefficient=rows.KAPPA*rows.alpha(col)/((S+A-rows.Q)*(S+BETA-rows.Q))
    scale=A-BETA
    return scale*F(lv,rows.UNIT),scale*(F(uv,rows.UNIT)+coefficient*X)


def dormant_return(col,trees,index,backward,active):
    """Positive dormant activation, with no further fixed-point iterations."""
    dorm={}
    for i in backward:
        t=trees[i]
        if mark(t)!=BETA:
            continue
        D,n=info(t)
        denom=S+A+n-F(1,D)
        low,high=[(A-BETA)*v for v in active[i]]
        for nxt,rate in physical(t):
            lv,uv=dorm[nxt] if nxt in index else dormant_boundary(nxt,col)
            low+=rate*lv
            high+=rate*uv
        dorm[t]=(F(floor(UNIT*low/denom),UNIT),F(ceil(UNIT*high/denom),UNIT))
    K=rows.N
    lam,v=S+A,S+BETA
    low=high=F(0)
    m_a=F(1)/(3+2*lam)
    m_b=F(1)/(3+2*v)
    for k in range(2,K):
        assert m_b>=m_a
        nxt=promoted(16,k)
        if nxt in index:
            al,au=active[index[nxt]]
            dl,du=dorm[nxt]
        else:
            al,au=boundary(nxt,col)
            dl,du=dormant_boundary(nxt,col)
        direct=rows.g(col,k)
        low+=(m_b-m_a)*(direct+al)+m_a*dl
        high+=(m_b-m_a)*(direct+au)+m_a*du
        m_a*=ratio(lam,k)
        m_b*=ratio(v,k)
    ta_lo,ta_hi,ma_lo,ma_hi=[m_a*x for x in moments(lam,K)]
    tb_lo,tb_hi,mb_lo,mb_hi=[m_b*x for x in moments(v,K)]
    tlo,thi=max(0,tb_lo-ta_hi),tb_hi-ta_lo
    mlo,mhi=max(0,mb_lo-ma_hi),mb_hi-ma_lo
    il,sl,iu,su=tail_coefficients(BETA,col)
    low+=il*tlo-sl*mhi
    high+=iu*thi-su*mlo
    # A late activation has exactly one old virgin child at the high-degree boundary.
    X=BETA*weight(newborn(16))
    cA=rows.KAPPA*rows.alpha(col)/(lam-rows.Q)
    direct=rows.g(col,K)
    dl=(A-BETA)*direct/(lam*v)
    du=(A-BETA)*(((direct+cA*rows.z(K))/lam+cA*rows.z(K)/(v-rows.Q))/v
                 +cA*X/(v-rows.Q))
    low+=dl*(ta_lo-2*ma_hi)
    high+=du*ta_hi
    return low,high


def output_bracket(lo,hi):
    return f'[{floor(UNIT*lo)}/{UNIT},{ceil(UNIT*hi)}/{UNIT}]'


def calculate(col):
    trees=retained()
    index={t:i for i,t in enumerate(trees)}
    weights=[weight(t) for t in trees]
    matrix,intrinsic,copies,lower,upper=[],[],[],[],[]
    for t in trees:
        degree,n=info(t)
        denom=S+mark(t)+n-F(1,degree)
        lo=hi=b(t)*rows.g(col,degree)
        coeff=defaultdict(F)
        for nxt,rate in successors(t).items():
            if not nxt[2]:
                l,h,star_coeff=star_row(nxt[0],info(nxt)[0],col)
                lo+=rate*l
                hi+=rate*h
                for m,k in star_coeff.items():
                    coeff[index[m]]+=rate*k
            elif nxt in index:
                coeff[index[nxt]]+=rate
            else:
                l,h=boundary(nxt,col)
                lo+=rate*l
                hi+=rate*h
        full={j:v/denom for j,v in coeff.items()}
        inside={index[nxt]:r/denom for nxt,r in physical(t) if nxt in index}
        rest={j:v-inside.get(j,0) for j,v in full.items() if v>inside.get(j,0)}
        assert all(full[j]>=v for j,v in inside.items())
        matrix.append(full)
        intrinsic.append(inside)
        copies.append(rest)
        lower.append(F(floor(UNIT*lo/denom),UNIT))
        upper.append(F(ceil(UNIT*hi/denom),UNIT))
    backward=sorted(range(len(trees)),key=lambda i:info(trees[i])[1],reverse=True)

    def intrinsic_sweep(rhs):
        out=rhs[:]
        for i in backward:
            out[i]+=sum(v*out[j] for j,v in intrinsic[i].items())
        return out

    copied=intrinsic_sweep([sum(v*weights[j] for j,v in row.items()) for row in copies])
    contraction=max((x/w for x,w in zip(copied,weights)),default=F(0))
    crude=rows.KAPPA*(A+rows.DAMP)
    universal=crude/(crude+S-rows.Q-rows.DAMP)
    assert contraction<1 and universal<1
    if rows.VT is None:
        assert contraction<=universal
    round_i=[[(j,round(UNIT*v)) for j,v in row.items()] for row in intrinsic]
    round_c=[[(j,round(UNIT*v)) for j,v in row.items()] for row in copies]

    def approximate(source):
        toll=[round(x*UNIT) for x in source]
        values=toll[:]
        for _ in range(SWEEPS):
            nxt=[t+sum(m*values[j] for j,m in row)//UNIT
                 for t,row in zip(toll,round_c)]
            for i in backward:
                nxt[i]+=sum(m*nxt[j] for j,m in round_i[i])//UNIT
            values=nxt
        model=[F(x,UNIT) for x in values]
        residual=intrinsic_sweep([x+sum(v*model[j] for j,v in row.items())-model[i]
                                  for i,(x,row) in enumerate(zip(source,matrix))])
        error=max((abs(x)/w for x,w in zip(residual,weights)),default=F(0))
        return model,error/(1-contraction)

    bottom,err_lo=approximate(lower)
    top,err_hi=approximate(upper)
    print('test',col,'retained multi-core states',len(trees),
          'copy contraction <=',ceil(contraction*10**7),'/ 10000000',
          'fixed point weighted error <=',ceil(max(err_lo,err_hi)*UNIT),'/',UNIT)
    result=[]
    for parent in (2,16):
        lo,hi,coeff=star_row(parent,2,col)
        lo+=sum(v*(bottom[index[t]]-err_lo*weight(t)) for t,v in coeff.items())
        hi+=sum(v*(top[index[t]]+err_hi*weight(t)) for t,v in coeff.items())
        result.append((lo,hi))
        print('true row',F(1,parent),output_bracket(lo,hi))
    active=[(max(F(0),x-err_lo*w),min(y+err_hi*w,rows.alpha(col)*w))
            for x,y,w in zip(bottom,top,weights)]
    diff=dormant_return(col,trees,index,backward,active)
    print('true positive dormant row',output_bracket(*diff))
    return result,diff


def scalar_enclosure():
    first,direct_excess=calculate(A)
    other,high_excess=calculate(BETA)
    kA,kB=first[0],other[0]
    L,U=high_excess
    # These prior-assisted caps each have a separate proof in excursion_constants.py.
    caps=((F('.5665'),F('.807760')),(F('.565809'),F('.8423')),
          (F('.565615'),F('.8522')),(F('.565'),F('.8844')),
          (F('.564'),F('.9395')),(F('.5631'),F('.9920')),(F('.563'),F('.9980')))
    usable=[cap for threshold,cap in caps if S>=threshold]
    cap=min([U]+usable)
    if L>=1:
        print('dormant generations have an infinite sum at this parameter')
        return
    low=kA[0]+kB[0]*(1+direct_excess[0])/(1-L)
    if cap>=1:
        print('scalar lower',floor(UNIT*low),'/',UNIT,'upper needs a smaller dormant upper')
        return
    high=kA[1]+kB[1]*(1+direct_excess[1])/(1-cap)
    print('true scalar',output_bracket(low,high),
          'dormant ceiling from', 'current entries' if U==cap else 'verified prior-assisted cap')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--virtual',action='store_true')
    p.add_argument('--size',type=int,default=N)
    p.add_argument('--s',type=F,default=S)
    p.add_argument('--tail',type=int,default=END)
    p.add_argument('--digits',type=int,default=14)
    p.add_argument('--sweeps',type=int,default=SWEEPS)
    p.add_argument('--test',choices=('A','B','half','all','scalar'),default='A')
    arg=p.parse_args()
    S,N=arg.s,arg.size
    assert arg.digits>=1 and arg.sweeps>=0
    UNIT,SWEEPS=10**arg.digits,arg.sweeps
    rows.UNIT=max(10**12,UNIT//100)
    rows.N=max(32,N+1)
    END=max(arg.tail,rows.N)
    rows.S=S
    if arg.virtual:
        import virtual_weight
        rows.virtual_weight=virtual_weight
        rows.VT,rows.KAPPA=virtual_weight.verify(False)
        rows.Q,rows.DAMP=virtual_weight.Q,virtual_weight.DAMP
    assert 3<=N and S>rows.Q+rows.DAMP
    print('s',S,'multi-core cutoff',N,'two-moment tail starts',rows.N,
          'telescope after',END)
    for a in (A,BETA):
        Tlo,Thi,Mlo,Mhi=moments(S+a,rows.N)
        print('relative moment intervals',a,output_bracket(Tlo,Thi),output_bracket(Mlo,Mhi))
    tests={'A':A,'B':BETA,'half':(F(1,2),F(1)),'all':0}
    if arg.test=='scalar':
        scalar_enclosure()
    else:
        calculate(tests[arg.test])
