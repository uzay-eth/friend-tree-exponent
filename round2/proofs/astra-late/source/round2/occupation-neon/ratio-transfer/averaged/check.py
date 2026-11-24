#!/usr/bin/env python3
"""Upper stationary prices for two exact paired D18 exception classes."""
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import comb,prod
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'fringe'))
from charges import Charges,Cert

LO=Q('.565676634491')
HI=Q('.565771874464')
A=Q('.565673')
B=Q('.565775')
PLO=Q(113465364604,4*10**11)
PHI=Q('1.356833021311')
QLO=Q(115462837236,4*10**11)
QHI=Q('1.367169995210')
GRID=10**22


def up(x):return Q(-(-x*GRID//1),GRID)


def pair_linear(E,motif,x,theta):
    (a,f),(b,g)=[e.terms(motif,x) for e in E]
    assert all(s[:3]==t[:3] for s,t in zip(f,g,strict=True))
    return a-theta*b,tuple((*u[:3],u[3]-theta*v[3]) for u,v in zip(f,g))


def spine_entry(m):
    # Only m successive genuine tip promotions can enter this face. Any
    # crown gain before entry excludes the face permanently. The crown's
    # total intrinsic monitored intensity is size-1/2.
    return prod(Q(1,2)/(LO+size-Q(1,2)) for size in range(2,m+2))


def spine(E,theta,cap=16,hinge=True):
    a,f=pair_linear(E,'P4',((0,1,1),)*4,theta)
    assert tuple(x[:3] for x in f)==((0,True,1),(3,False,1))
    _,_,_,p=f[0]; _,_,_,w=f[1]
    out=[]
    for sign in (1,-1):
        # Discount LO suffices for each nonnegative hinge. At an overflow
        # degree, arm moves leave the terminal price unchanged.
        tail=up(max(max(Q(0),sign*(a+p*r+w/d)) for r,d in product(
            (Q(0),Q(1,2)),(Q(cap+1),))) /(LO+3))
        # The reciprocal interval contains zero (not physically a mark).
        tail=max(tail,up(max(max(Q(0),sign*(a+p*r)) for r in (Q(0),Q(1,2)))/(LO+3)))
        V={}
        for d in range(cap,1,-1):
            for k in range(d-1,-1,-1):
                l=d-k-1; b=Q(l,d)
                gval=V[d+1,k] if d<cap else tail
                bval=V[d,k+1] if l else 0
                points=[(Q(0),Q(0)),(Q(1,2),Q(0))]
                # reward has one possible kink; endpoint maximum already
                # suffices for an affine numerator on each branch, but
                # check the kink also.
                if p and 0<-(a+w/d)/p<Q(1,2): points.append((-(a+w/d)/p,Q(0)))
                tests=[]
                for (r,_),g,s in product(points,(l+Q(1,2),l+Q(k+1,2)),(LO,) if hinge else (LO,HI)):
                    den=s+3+1/Q(d)+r+g+b
                    flow=sign*(a+p*r+w/d)
                    rhs=(max(Q(0),flow) if hinge else flow)+g*gval+b*bval
                    tests.append((den,rhs))
                V[d,k]=up(max(rhs/den for den,rhs in tests))
                assert all(den*V[d,k]>=rhs for den,rhs in tests)
        out.append(HI*spine_entry(4)*max(Q(0),V[2,0]))
    return out


def binomial_joint(K,j):
    """Dominating rate-one births / independent rate-3/2 lifetimes."""
    # Select j live arrivals from the *same* Poisson process. Campbell's
    # factorial identity leaves a rate-one process for the other arrivals.
    # E[binom(V,j)1(N>=K)|age x] <= (1-exp(-c*x))**j/(j!*c**j)
    # times P(Pois(x)>=K-j), with equality in the dominating model.
    c=Q(3,2)
    M=max(0,K-j)
    from math import factorial
    return sum((-1)**i*comb(j,i)*LO/(LO+c*i)/(1+LO+c*i)**M for i in range(j+1))/(factorial(j)*c**j)


def exact_hub_mass(K):
    # At K arrivals the probability of having lost no virgin is bounded
    # by racing birth (<=1) against cumulative killing (>=3j/2).
    # The actual K-kid all-virgin state is visited at most once.
    return HI*prod(1/(LO+1+Q(3*j,2)) for j in range(K))/(LO+Q(3*K,2))


def vertex_balance():
    # Whole physical stationary balance: pi(k)=1 and E exp(eta*k) grows
    # under births only. Combining the exact mean with arbitrary rate<=1
    # gives a universal improvement over a pure Poisson tail: bound
    # pi(k>=K) by optimizing the first-hit mass at any earlier j.
    out={}
    for K in (8,12,16):
        # Once k>=j, the independent stable-age residual age is exponential.
        # Thus pi(k>=K) <= q**(K-j) pi(k>=j). Sum backward *lower*
        # bounds on the earlier survival probabilities; E[k]=1.
        v=1/sum((1+LO)**(K-j) for j in range(1,K+1))
        out[K]=min(v,1/(1+LO)**K)
    return out


def bare_host(K):
    # Pay every bare D=K+1,k=K host, whether its arms are young or old.
    # Degree overflow never returns. Core gain and promotion use common
    # rho; no independent free probability of arity and bare stock.
    V={}
    n=K+1
    for d in range(n,1,-1):
        for k in range(d-1,-1,-1):
            l=d-k-1; b=Q(l,d)
            r=Q((d,k)==(n,K))
            nxt=V.get((d+1,k),0); same=V.get((d,k+1),0)
            cand=[(r+g*nxt+b*same)/(LO+g+b) for g in (Q(l),l+Q(k+1,2))]
            V[d,k]=up(max(cand))
            assert all(z<=V[d,k] for z in cand)
    return HI*V[2,0]


def hub(E,theta):
    nb=(0,1,())
    z=(0,0,(nb,)*16)
    out=[e.rows(z) for e in E]
    assert set(out[0])==set(out[1])
    corners=[]
    for p in (Q(0),Q(1,2)):
        corners.append(sum(out[0][f][0]+p*out[0][f][1]-theta*(out[1][f][0]+p*out[1][f][1]) for f in out[0]))
    return corners,[exact_hub_mass(16)*max(Q(0),*(s*q for q in corners)) for s in (1,-1)]


def mixing_needed(theta,exception):
    # If both signed inequalities are inferred from one symmetric cone,
    # these are allowed q-widths about the calibration for a bracket gain
    # with only the currently certified denominator extrema. Also price a
    # one-sided test at the nearest endpoint with a prospective 1e-5 move.
    l,r=LO+Q('1e-5'),HI-Q('1e-5')
    def odds(s):return (s-A)/(B-s)
    t_lower=odds(l)*PHI/QLO
    t_upper=odds(r)*PLO/QHI
    print('for 1e-5 per endpoint using published denominator extrema: required lower slack ratio >',float(t_lower),
          'or upper slack ratio <',float(t_upper))
    deficit=(l-A)*PHI-((B-l)*theta*QLO)
    surplus=(r-A)*PLO-((B-r)*theta*QHI)
    print('with unknown denominator ratio, sufficient lower/upper global contrast',float(deficit),float(surplus),
          'and known two-exception errors',*(float(x) for x in exception))
    # This is an illustrative design target, not a proved denominator ratio.
    lo_ratio=Q('.98'); hi_ratio=Q('1.06')
    needl=odds(l)/lo_ratio
    needu=odds(r)/hi_ratio
    # Reuse least possible δ+, on the existing (unmoved) upper endpoint.
    dplus=(B-HI)*QLO
    dplus=max(dplus,Q('3.921982939e-6'))  # independent fringe-ledger floor
    if t_lower<theta<t_upper:
        print('allowed centred hinge using published denominator extrema and existing upper-slack floor',
              float((theta-t_lower)*dplus),float((t_upper-theta)*dplus))
    if needl<theta<needu:
        print('conditional on denominator ratio [.98,1.06]: allowed centred hinge budget',
              float((theta-needl)*dplus),float((needu-theta)*dplus))


def relative_hinge(theta,block,huberr):
    # This targets a quasi-relative bound, rather than requiring the
    # full contrast itself to have either sign. H is per core; for each
    # endpoint evaluate the worst corner of the three valid denominators.
    Rlo=Q('4.410155511e-6')
    Rhi=Q('3.921982939e-6')
    for l,r in [(LO+Q('1e-5'),HI-Q('1e-5'))]:
        a1=(l-A); b1=(B-l)*theta
        a2=(r-A); b2=(B-r)*theta
        assert min(a1,b1,a2,b2)>0
        # At a proposed lower crossing, E = (mu-A)P-theta(B-mu)Q.
        # πδ+>=Rhi forces Q>=Rhi/(B-mu); πδ− similarly forces P.
        rejected_lower=a1*PHI-max(b1*QLO,theta*Rhi)
        rejected_upper=max(a2*PLO,Rlo)-b2*QHI
        # Conditional tolerances for *rest* of all-motif paired hinges.
        lower_tol=-rejected_lower-block[1]-huberr[1]
        upper_tol=rejected_upper-block[0]-huberr[0]
        print('for 1e-5 moves, rest-cone budget (negative, positive) after block and independent slack floors',
              float(lower_tol),float(upper_tol))
    for x in (LO,HI):
        # Concrete feasible accounts in the coarse residual relaxation;
        # the two exact exception signs allow zero there, not a positive
        # mandatory amount. Remaining charges need actual upper control.
        dm,dp=(x-A)*PHI,(B-x)*QHI
        assert dm>=Rlo and dp>=Rhi
        print('bare-cone feasible endpoint contrast at',float(x),float(dm-theta*dp))


def main():
    E=[Charges(Cert(Path(x))) for x in sys.argv[1:]]
    assert [e.c.lam for e in E]==[A,B]
    print('spine exact-entry upper discounted mass',float(spine_entry(4)))
    for K,v in vertex_balance().items():
        b=[binomial_joint(K,j) for j in (0,1,2,3,K)]
        assert min(b)>0 and b[-1]>=exact_hub_mass(K)
        print('arity',K,'mass, first/second/third binomial virgin moment, all-virgin moment',*(float(q) for q in b))
        print('  improved arity mass from exact core-edge balance',float(v))
        print('  canonical bare host, any arms, upper',float(bare_host(K)))
    for theta in (Q(1),Q('1.2415')):
        S=spine(E,theta); corners,H=hub(E,theta)
        block=spine(E,theta,hinge=False)
        assert block[0]==0 and block[1]<S[1]
        print('common multiplier',theta,'spine hinge pi positive/negative upper',*(float(q) for q in S))
        print('  after on-path Poisson correction, spine positive/negative upper',*(float(q) for q in block))
        print('hub K16 exact event pi upper',float(exact_hub_mass(16)),'paired apex corners',*(float(q) for q in corners),
              'hinge pi positive/negative upper',*(float(q) for q in H))
        gap=B-A
        # κ=Φ-+θΦ+, old bracket safe δ total bounds. At a trial m
        # use inequalities on both Φ averages to compare the target sign.
        m=(A+B)/2
        fldlo=(m-A)*PLO-theta*(B-m)*QHI
        fldhi=(m-A)*PHI-theta*(B-m)*QLO
        print('midpoint global contrast threshold from denominator only',float(fldlo),float(fldhi))
        budget=max((s-A)*PHI+(B-s)*QHI for s in (LO,HI))
        print('remaining worst-case sum budget before allocating any charge',float(budget),
              'negative exception fraction of gap',float((S[1]+H[1])/gap))
        mixing_needed(theta,[x+y for x,y in zip(S,H)])
        relative_hinge(theta,block,H)


if __name__=='__main__':main()
