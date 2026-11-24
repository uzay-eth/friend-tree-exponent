"""Exact excess-on-directed-edges certificate (no tree enumeration)."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path

M=run_path(str(Path(__file__).resolve().parents[1]/'mixed-stock/check.py'))
G,U,L,P=(M[k] for k in ('G','U','L','profile'))
A,B=Q('.00634'),Q('.017')
F=Q('.00356')
D2=Q('.00625')
H=Q('.012')


def upper_out(d,e,s,u):
    # The opposite site receives at least the known 1/d across this edge.
    # Symbol 0 only for a hypothetical neighbour of arbitrarily high degree.
    t=Q(1,d*e*(e+1)) if e else 0
    return max(0,s-t)*max(0,u-(L(e+1) if e else 0))


def row(d,c,s,v,u):
    return v+A*(c-2)-B-c*max(0,s*u-B)


def margin(d,ds,high):
    c=len(ds);s,v,u=P(d,tuple(e or high for e in ds))
    return v+A*(c-2)-B-sum(max(0,upper_out(d,e,s,u)-B) for e in ds)


def transfer():
    # Unpaid, degree-3/4 and chain ports have reciprocal minima 0,1/4,1/2.
    # Check the vertices of the port-count simplex at both receive ends.
    z=min((row(d,c,*P(d,(e,)*c))-c*cost,d,c,cost,e)
          for d in range(10,113) for c in range(2,d+1)
          for cost,e in ((0,0),(0,2),(F,4),(F,2),(D2,2)))
    assert z[0]>0
    bracket=M['g']*Q(113,114)+G(2)/2-M['p']/2
    assert bracket<0 and bracket/2+(A-F)*Q(113,2)-2*A-B>0
    assert M['g']*Q(113,228)+2*(A-F)-2*A-B>0
    # At full reciprocal receive, s_l-s<=1/[2d(d+1)]; the old-leaf
    # credit pays this, and pays the opposite debit for l<=3 by the
    # mixed-stock estimate s-s_l<=(l-1)/[2d(d+1)]. For l>=4 it is zero.
    assert G(2)>Q(2,114)
    assert M['g']*Q(113,228)+2*(A-D2)-2*A-B>0
    # A terminal high endpoint can pay one port with ample surplus.
    assert Q(9,11)*M['g']+Q(9,20)*G(2)-max(B,Q(1,22))-A>D2
    print('PASS high donors, finite minimum',z,'and unbounded tail')

    # Count only high slots: low-low edges retain exactly the old rule.
    for d,c in ((3,2),(3,3),(4,4)):
        z=min((margin(d,ds,r)+F*ds.count(0),ds,r)
              for ds in combinations_with_replacement((0,*range(2,10)),c) for r in (0,10))
        assert z[0]>0
        print('PASS transfer recipient',d,c,'minimum',z)


def chains():
    for d in range(10,42):
        for l in range(d):
            assert Q(d+l,2*d*(d+1))*U(d+1,l)<H
    assert Q(1,86)<H
    u=U(3,0);f=max(0,u-L(3));g=G(3)
    def port(s,l):return max(H,s*max(0,u-l))
    def debit(ss,l,r):
        return port(ss[0],l)+port(ss[-1],r)+f*sum(abs(x-y) for x,y in zip(ss,ss[1:]))-2*g*sum(ss)
    def low_terminal(e,s):
        t,v,w=P(e,(2,))
        return v-max(0,(s-t)*max(0,u-L(e+1)),(t-s)*max(0,w-L(3)))
    from itertools import product
    m=Q(1,6)
    ends=[(Q(1,d), L(d+1)) for d in range(10,31)] + [(Q(0),Q(0)),(Q(1,31),Q(0))]
    for k in (1,2,3):
        vals=[]
        for (x,l),(y,r) in product(ends,repeat=2):
            ss= ((x+y)/6,) if k==1 else ((Q(1,2)+x)/6,)+(m,)*(k-2)+((Q(1,2)+y)/6,)
            vals.append(debit(ss,l,r))
        z=max(vals)
        assert z<B+2*D2
        print('PASS suppressed degree-two chain length',k,'max debit',z,'budget',B+2*D2)
        vals=[]
        for x,l in ends:
            for e in range(2,10):
                y=Q(1,e)
                ss=((x+y)/6,) if k==1 else ((Q(1,2)+x)/6,)+(m,)*(k-2)+((Q(1,2)+y)/6,)
                vals.append(low_terminal(e,ss[-1])+2*g*sum(ss)-port(ss[0],l)-f*sum(abs(a-b) for a,b in zip(ss,ss[1:])))
        z=min(vals)
        assert z>A+B
        print('PASS high-to-low-terminal degree-two path length',k,'minimum joint credit',z)
        vals=[]
        for d,e in product(range(2,10),repeat=2):
            x,y=Q(1,d),Q(1,e)
            ss=((x+y)/6,) if k==1 else ((Q(1,2)+x)/6,)+(m,)*(k-2)+((Q(1,2)+y)/6,)
            vals.append(low_terminal(d,ss[0])+low_terminal(e,ss[-1])+2*g*sum(ss)-f*sum(abs(a-b) for a,b in zip(ss,ss[1:])))
        z=min(vals)
        assert z>B+2*A
        print('PASS low-terminal to low-terminal path length',k,'minimum joint credit',z)
    assert g>0  # every further interior site has s=1/6 and zero extra pair debit.


def check():
    # Import only the finite endpoint conclusions and analytic constants of
    # mixed-stock, not its graph fixtures. Its large-d inequalities give
    # C+a(c-2)>1/22 once d>=113; sU<=1/[2(d+1)]<B there.
    M['scalar']()
    assert Q(1,228)<B<Q(1,22)
    dmin=min((row(d,c,*P(d,(r,)*c)),d,c,r)
             for d in range(5,113) for c in range(2,d+1) for r in (0,2))
    assert dmin==(Q(101,50000),5,5,0)
    extra=[min(row(4,c,*P(4,(r,)*c)) for r in (0,2)) for c in (2,3)]
    assert min(extra)>0
    print('PASS 12424 degree/stock corners; minimum',dmin,'stocked d4',*extra)

    tmin=min((v-max(B,s*u)-A,d,r) for d in range(2,10) for r in (0,2) for s,v,u in (P(d,(r,)),))
    assert tmin==(Q(16436047,3000000000000),2,0)
    # For d>=10: terminal-stock analytic bound from mixed-stock and outgoing
    # exposure <=1/[2(d+1)]. This also bounds max(B, own exposure).
    assert Q(9,11)*M['g']+Q(9,20)*G(2)-max(B,Q(1,22))>A
    print('PASS terminal vouchers at all degrees; finite minimum',tmin)

    # Keep all low slots exact; remaining slots have degree >=10 and can
    # be bounded as a single reciprocal interval. Only local multisets.
    for d,c,k,want in ((3,2,1,Q(1173389118653,270000000000000)),
                       (3,3,3,Q(15612335509,4000000000000)),
                       (4,4,2,Q(2503587899,750000000000)),
                       (2,2,1,None)):
        z=min((margin(d,ds,r),ds,r) for ds in combinations_with_replacement((0,*range(2,10)),c)
              if sum(bool(e) for e in ds)>=k for r in (0,10))
        if want is not None:assert z[0]==want>0
        else:assert z[0]<0
        print('local degree multiset',d,c,'at least',k,'low slots:',z)

    # Uniform adverse budgets, usable on every old tree and hence at each
    # deterministic size under the actual law (no stationarity assumed).
    for d,c,bill in ((2,2,B),(3,2,Q('.00414')),(3,3,Q('.01066')),(4,4,Q('.00432'))):
        z=min(margin(d,ds,r) for ds in combinations_with_replacement((0,*range(2,10)),c) for r in (0,10))
        assert z>=-bill
        print('all-neighbour lower ledger',d,c,z,'budget',bill)

    # These are actual local profile deficits for the corrected debit rule,
    # not physical signs. Keep the actual finite opposite degree on each port.
    bad=((2,(2,20)),(3,(111,111)),(3,(111,111,111)))
    for d,ds in bad:
        r=margin(d,ds,0)
        assert r<0
        print('FAIL local ledger only',d,ds,r)
    transfer()
    chains()
    run_path(str(Path(__file__).resolve().parents[1]/'interfaces/low-end.py'))['verify']()
    print('PASS all degree-two chain endpoint types and simultaneous low endpoint budgets')


if __name__=='__main__': check()
