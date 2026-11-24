"""Directed integer bounds for the nineteen-mark virgin-trigger comparison."""
from fractions import Fraction as F
from functools import cache
from math import ceil, floor
from argparse import ArgumentParser

UNIT=10**12
D=20
BLOCK=(2,3,4,6,10,20)


@cache
def stencil(lam,port):
    """Two endpoints for the old-child inputs and an optional live parent."""
    j=port or 1
    n,q=lam.numerator,lam.denominator
    lines=[]
    for d in range(D-1,1,-1):
        for c in range(d,0,-1):
            L=d-c
            common=2*j*d*q
            b=2*j*q*L
            g=2*j*d*q*L
            G=g+d*q*(j*(c-1)+2*bool(port))
            kill=2*j*d*n
            lines.append((d,c,common,b,g,G,kill+b+g,kill+b+G))
    return lines


def eta(c):
    return F(1,D*(D+1)) if c==1 else F(0)


def table(lam,port,toll,base,tail,upper=False):
    """When toll is given the source is b*toll[d]; otherwise it is base[d,c]."""
    flat=[[0]*(D+2) for _ in range(D+1)]
    for c in range(1,D+1):
        flat[D][c]=(ceil if upper else floor)(tail(c))
    for d,c,common,b,g,G,den,Den in stencil(lam,port):
        source=b*toll[d] if toll is not None else common*base[d][c]
        source+=b*flat[d][c+1]
        jump=flat[d+1][c]
        v,V=source+g*jump,source+G*jump
        flat[d][c]=(max(-((-v)//den),-((-V)//Den)) if upper else min(v//den,V//Den))
    return flat


def root_table(lam,toll):
    return table(lam,0,toll,None,lambda c:F(D-c,D)*toll[D]/(lam+eta(c)))


def blocks(s,R,count,marks=BLOCK):
    bind=[next((i for i in reversed(marks) if i<=d),marks[0]) for d in range(D+1)]
    g={i:int(i==20)*UNIT for i in marks}
    H={i:[[0]*(D+2) for _ in range(D+1)] for i in marks}
    tails=dict.fromkeys(marks,0)
    for _ in range(count):
        least={}
        for j in marks:
            v=H[j]
            least[j]=[[min(v[d][c],v[min(d+1,D)][c]) for c in range(D+2)]
                      for d in range(D+1)]
        new,after={},{}
        for i in marks:
            lam=s+F(1,i)
            w={}
            for j in marks:
                child_lam=s+F(1,j)
                w[j]=table(lam,j,None,least[j],lambda c:
                           F(D-c,D)*tails[j]/((lam+eta(c))*(child_lam+eta(c))))[2][1]
            toll=[0,0]+[g[bind[d]]+w[bind[d]]//d for d in range(2,D+1)]
            new[i]=root_table(lam,toll)
            after[i]=toll[D]
        H,tails=new,after
        for idx,i in enumerate(marks):
            g[i]=max(g[i],floor(H[i][2][1]/R),g[marks[max(0,idx-1)]])
    assert g[20]==UNIT
    return g


def interp(anchor,critical,start,u,target):
    """Check a positive log-convex bound without floating powers."""
    p=(u-start[0])/(u-start[1])
    assert 0<=p<=1
    assert anchor**p.numerator*critical**(p.denominator-p.numerator)<target**p.denominator


def simple_anchor(s,trial,damp=F(3,20),q=F(9,25),hmin=F(7,10)):
    alpha=1/hmin
    ratios=[]
    for i in range(2,D+1):
        lam=s+F(1,i)
        toll=[0,0]+[ceil(UNIT*((trial[d] if d>2 else 0)
                   +alpha/(1+damp*d)/(lam-q))) for d in range(2,D+1)]
        result=table(lam,0,toll,None,lambda c:toll[D]*(1/lam-F(c,D)/(lam+1)),True)
        ratios.append(F(result[2][1],UNIT)/trial[i])
    return max(ratios)


def decreasing_weight():
    q,damp=F('.272'),F('.19')
    unit=10**9
    z=lambda d:1/(1+damp*d)
    h={(D,c):F(ceil(unit*z(D)/q),unit) for c in range(1,D+1)}
    for d in range(D-1,1,-1):
        for c in range(d,0,-1):
            L=d-c
            b=F(L,d)
            val=(L*h[d+1,c]+b*(h.get((d,c+1),0)+z(d)))/(q+L+b)
            h[d,c]=F(ceil(unit*max(F(1,4),h[d+1,c],val)),unit)
    h[2,1]=F(1)
    for d in range(2,D):
        for c in range(1,d+1):
            L=d-c
            assert h[d+1,c]<=h[d,c]
            assert (L*(h[d+1,c]-h[d,c])+F(L,d)*(h.get((d,c+1),0)-h[d,c]+z(d))
                    <=q*h[d,c])
    assert z(D)<=q*h[D,1] and min(h.values())==F('.765931373')
    return min(h.values())


if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('--iterations',type=int,default=24)
    p.add_argument('--all-labels',action='store_true')
    p.add_argument('--critical-only',action='store_true')
    p.add_argument('--refined-only',action='store_true')
    args=p.parse_args()
    u=F('.5665')
    fine=F('.5658')
    s0=F('.56564')
    for s in (s0,fine,u):
        lam=s+F(1,2)
        toll=[0]*D+[UNIT]
        value=root_table(lam,toll)[2][1]
        uniform=table(lam,0,toll,None,lambda c:F(D-c,D)*UNIT/(lam+F(1,420)))[2][1]
        print('s',s,'original root label-20 lower',value,'/ 1e12',
              'uniform tail',uniform,'/ 1e12',
              'root total direct lower',floor(UNIT/(2*(lam+F(1,6)))),'/ 1e12')
    # This input uses the reviewed physical upper .5657979941 < .5658.
    interp(F(50,7),F(1),(s0,F('.51')),fine,F('1.006'))
    marks=tuple(range(2,21)) if args.all_labels else BLOCK
    cases=((u,F(1)),(s0,F('1.006')))
    assert not (args.critical_only and args.refined_only)
    for s,R in (cases[:1] if args.critical_only else cases[1:] if args.refined_only else cases):
        g=blocks(s,R,args.iterations,marks)
        print('s used in lower tables',s,'upper PF multiplier',R,'blocks',
              ' '.join(f'{i}={g[i]}' for i in marks),'units 1e-12')
        removed=F(g[2],UNIT)/(2*(s+F(1,20)+F(3,2)))
        print('complement Perron-norm contraction <=',ceil(UNIT*(R-removed)),'/ 1e12')
    for k in (0,2,4,6):
        trial={j:(1+F(k,20))/(1+F(k,j)) for j in range(2,21)}
        cap=simple_anchor(F('.51'),trial)
        print('one-residence elementary complement anchor at .51, k',k,'<=',ceil(10**6*cap),'/ 1e6')
    if args.iterations>=24:
        hmin=decreasing_weight()
        assert F(1)/(F('.24')*F('.765931373'))<F('5.45')
        interp(F('5.45'),F('.971864'),(F('.5649'),F('.462')),u,F('.9981'))
        print('old physical bracket: complement contraction at .5649 < .9981')
        for k in (0,4,10):
            trial={j:(1+F(k,20))/(1+F(k,j)) for j in range(2,21)}
            cap=simple_anchor(F('.462'),trial,F('.19'),F('.272'),hmin)
            print('one-residence complement anchor at .462, k',k,'<=',ceil(10**6*cap),'/ 1e6')
            if k==4:
                assert cap<F('3.489')
    m=F('.11793629')
    theta=F('.978128')
    assert 1+F('.006')/(1-theta)<F('1.275')
    assert 1+F(10,7)*5/(m*(1-theta))<2800
    assert F(10,7)*F('1.275')/m<16
    assert F(16)*5/(1-theta)<3660
    assert F('1.275')*theta**900<F('2.895e-9')
    assert (F('.5656581227')+2)/(F('.565615')+2)>F('1.00001680')
    assert (F('.5657979941')+2)/(F('.565809')+2)<F('.99999572')
    print('fine bracket: chronological tail factor < 1.275, weight factor < 2800,',
          'scalar size-exit factor < 3660')
