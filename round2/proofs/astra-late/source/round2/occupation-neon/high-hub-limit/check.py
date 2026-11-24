#!/usr/bin/env python3
"""Small exact D18 arm pilot and quantitative high-hub error account."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from math import comb, factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
OCC = HERE.parent
sys.path.insert(0, str(OCC / 'fringe'))
from charges import Cert, Charges, NB, embedding, shape  # noqa: E402
sys.path.insert(0, str(OCC.parent / 'sol-kernel-verifier'))
import verify as v  # noqa: E402

LO, HI = Q('.565673570322'), Q('.565775')
J = (2, 1, ())
FAMILIES = ('hV', 'P4', 'Yt')


def ztuple(n=3):
    return (Q(0),) * n


@cache
def plain(x):
    return (2, x[1], tuple(sorted(map(plain, x[2]))))


def arms(cap):
    todo = [J]
    seen = {J}
    arc = {}
    for x in todo:
        out = defaultdict(Q)
        for y, r in v.row(x)[3].items():
            y = plain(y)
            if v.size(y) <= cap:
                out[y] += r
                if y not in seen:
                    todo.append(y)
                    seen.add(y)
        arc[x] = out
    return sorted(todo, key=lambda x: (-v.size(x), x)), arc


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def scale(a,x):
    return tuple(y*x for y in a)


def mul(a,b):
    out = [Q(0)] * (len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return tuple(out)


@cache
def basis(nodes):
    out = []
    for x in nodes:
        p = (Q(1),)
        for y in nodes:
            if y!=x:
                p = mul(p,(-y/(x-y),1/(x-y)))
        out.append(p)
    return out


def interp(nodes,ys):
    out = ztuple(len(nodes))
    for a,y in zip(basis(nodes),ys):
        out = add(out,scale(a,y))
    return out


def at(p,x):
    a = Q(0)
    for y in reversed(p):
        a = a*x+y
    return a


def bernstein(p, length):
    n = len(p)-1
    q = tuple(x*length**k for k,x in enumerate(p))
    return tuple(sum(q[j]*Q(comb(i,j),comb(n,j)) for j in range(i+1)) for i in range(n+1))


class Rewards:
    def __init__(self, blob, states):
        self.c = Charges(Cert(blob))
        self.stencils = {}
        for x in states:
            p,old = shape(x,1)
            host = (-1,0,*(i+1 for i in p[1:]))
            out = []
            for f,ids in embedding(host):
                if f not in FAMILIES or ids[:2]!=(0,1):
                    continue
                xx = tuple(old[i-1] for i in ids[1:])
                assert all(sum(y)<18 for y in xx)
                ext = []
                for i in ids[1:]:
                    external_children = [j for j in range(1,len(host)) if host[j]==i and j not in ids]
                    ext.append(sum((Q(1,sum(old[j-1])) for j in external_children), Q(0)))
                out.append((f,xx,ext))
            self.stencils[x] = out

    @cache
    def row(self,x,d,k):
        """Actual D18 signed drift in one arm, affine in shared parent+sibling receive."""
        assert d>=32 and 1<=k<d
        a,b = [Q(0)]*3,[Q(0)]*3
        for f,xx,ex in self.stencils[x]:
            tag = ((d-k-1,k,1),*xx)
            if not self.c.supported(f,tag):
                continue
            j = FAMILIES.index(f)
            base, free = self.c.terms(f,tag)
            slope = None
            for i,parent,number,w in free:
                if i==0:
                    assert number==(1 if parent else k-1)
                    if parent:
                        slope = w
                        b[j] += w
                    else:
                        assert w==slope
                else:
                    assert not parent
                    base += w*ex[i-1]
            a[j] += base
        assert all(x+min(Q(0),y*k/2)>=0 for x,y in zip(a,b))
        return tuple(a),tuple(b)

    def polynomials(self,x):
        """Numerators after clearing (1+1/D)^3, k≤18 then affine overflow."""
        if not self.stencils[x]:
            return None
        # The family promotions at this arm cap are entirely inner below the
        # host. At the apex only the gain uses inverse 1/(d+1); numerator
        # degrees at most seven (the apparent 1/x from L cancels).
        ds = tuple(range(33,41))
        nodes = tuple(Q(1,d) for d in ds)
        out = {}
        for k in range(1,self.c.c.K+3):
            rows = [self.row(x,d,k) for d in ds]
            p = tuple(tuple(interp(nodes,[(1+t)**3*ab[i][f] for t,ab in zip(nodes,rows)])
                            for f in range(3)) for i in range(2))
            a,b = self.row(x,83,k)
            assert all(at(p[i][f],Q(1,83)) == (1+Q(1,83))**3*ab[f]
                       for i,ab in enumerate((a,b)) for f in range(3))
            out[k] = p
        arc = tuple(p[0] for p in out[1][0])
        assert all(tuple(p[0] for p in a)==arc and all(p[0]==p[1]==0 for p in b)
                   for a,b in out.values())
        for d,k in ((32,18),(83,29),(101,81)):
            a,b = out[min(18,k)]
            if k>18:
                a = tuple(add(z,scale(add(w,scale(z,-1)),k-18))
                          for z,w in zip(a,out[19][0]))
            p = Q(1,d)
            assert tuple(tuple(at(f,p)/(1+p)**3 for f in part) for part in (a,b))==self.row(x,d,k)
        return arc,out


def stage_bounds(p,inf,M):
    """Envelope for the FIRST-order discounted residual, per actual arm."""
    # q_real − q_inf = D^-1 * numerator / (1+1/D)^3.
    # Keep actual donor sum u in [0,k/2] via both sides of the affine row.
    arc, rows = p if p is not None else (ztuple(),{})
    assert arc==inf
    if not rows:
        return [(ztuple(),ztuple())]*19
    out = []
    for k,(a,b) in rows.items():
        polys = [add(z,scale((Q(1),Q(3),Q(3),Q(1),Q(0),Q(0),Q(0),Q(0)),-f))[1:]
                 for z,f in zip(a,arc)]
        upper = [add(a,scale(b[1:],Q(k,2))) for a,b in zip(polys,b)]
        lo,hi = [],[]
        bound = 1/(1+Q(1,M))**3
        for a,b in zip(polys,upper):
            c = (*bernstein(a,Q(1,M)),*bernstein(b,Q(1,M)))
            l,u = min(c),max(c)
            lo.append(min(l,l*bound)); hi.append(max(u,u*bound))
        out.append((tuple(lo),tuple(hi)))
    return out


def resolvent(states, arcs, reward, s, hazard=Q(0)):
    out = {}
    for x in states:
        a = reward[x]
        for y,r in arcs[x].items():
            a = add(a,scale(out[y],r))
        out[x] = scale(a,1/(s+hazard+v.size(x)-Q(1,v.degree(x))))
    return out


def fixed(a):
    return f'{float(a):.10g}'


def upward(a, n=14):
    f = -(-a*10**n//1)
    return f'{f//10**n}.{f%10**n:0{n}d}'


def downward(a, n=14):
    return '-'+upward(-a,n) if a<0 else format_scaled(a*10**n//1,n)


def format_scaled(f,n):
    return f'{f//10**n}.{f%10**n:0{n}d}'


def clock(s,M,c,depth=4):
    """Iterate inverse-degree moment transport against the same Poisson age."""
    # A_(q,m)=∫ e^-st E[(N)_m (D-1)^(-q)]dt. Every gain decreases
    # f_q by f_q*(1-(1-1/D)^q); N and D never jump simultaneously.
    # gamma is a lower bound on D*[1-(1-1/D)^q]. Then compensate
    # c_t/D f_q <= (c+N) f_(q+1) at the same age.
    out = [Q(1,(M-1)**(depth+1)) * factorial(m)/s**(m+1)
           for m in range(depth+4)]
    for q in range(depth,0,-1):
        gamma = M*(1-(1-Q(1,M))**q)
        norm = Q(1,(M-1)**q)
        cur=[]
        for m in range(len(out)-1):
            drift = (c+m)*out[m]+out[m+1]
            a = ((norm if m==0 else m*cur[m-1])+gamma*drift)/(s+gamma)
            cur.append(min(a,norm*factorial(m)/s**(m+1)))
        out=cur
    return tuple(out)


def common_age(s,M,order=2):
    """E[(N(A))_j 1{Y(A)≥M}], A~Exp(s), using the same stable age."""
    p=Q(1)
    h = [1/s**m for m in range(order+1)]
    for d in range(2,M):
        g=Q(d)-Q(1,2)
        p *= g/(g+s)
        for m in range(1,order+1):
            h[m] += h[m-1]/(g+s)
    return tuple(p*h[m]*factorial(m) for m in range(order+1))


def band_price(states, arcs, inf, M, s, K, pres, upper=False):
    """Two-sided pricing on the young-arm/arity mask, at physical rates."""
    # For any k, discard sign-asymmetric unreliable reciprocals only after
    # grouping the actual charge, and kill at one exterior input proposal.
    out = []
    last = {x:ztuple() for x in states}
    grid=10**23
    for k in range(K,0,-1):
        b0=Q(M-k-1,M)
        assert b0>0
        nxt = {}
        for x in states:
            lo,hi = pres[x][k-1]
            if upper:
                a = tuple(f+max(Q(0),b)/M for f,b in zip(inf[x],hi))
            else:
                a = tuple(max(Q(0),f+min(Q(0),b)/M) for f,b in zip(inf[x],lo))
            for y,r in arcs[x].items():
                a=add(a,scale(nxt[y],r))
            rate=s+v.size(x)-Q(1,v.degree(x))
            obs = nxt.get(v.plus(x),ztuple()) if upper else ztuple()
            def vals(a,new,plus):
                inputs=(Q(0),Q(1,M)) if upper else (Q(1,M),)
                return [(a+b*new+r*plus)/(rate+b+r) for b in (b0,Q(1)) for r in inputs]
            data = [vals(*z) for z in zip(a,last[x],obs)]
            if upper:
                nxt[x]=tuple(Q(-(-max(b)*grid//1),grid) for b in data)
            else:
                nxt[x]=tuple(Q(min(b)*grid//1,grid) for b in data)
        last=nxt
        out.append(nxt)
    return out


def birth_price(M,K,prices,lower_k,upper=False):
    """Discounted compensation over actual host births in degree band [M,∞)."""
    s = LO if upper else HI
    nxt = [ztuple() for _ in range(K+1)]
    for d in range(M,1,-1):
        cur = [ztuple() for _ in range(K+1)]
        for k in range(min(K,d-1),-1,-1):
            L,c = d-k-1,k+1
            b=Q(L,d)
            grow = nxt[k] if d<M else None
            # At entry and thereafter remove all exits at which k exceeds K.
            # For d≥M the unconditional floor b≥1-(k+1)/M suffices;
            # increased rate can help or hurt as k advances, so keep both.
            enter = prices[K-k-1][J] if d==M and lower_k<=k+1<=K else ztuple()
            source = scale(add(enter,cur[k+1]),b) if k<K else ztuple()
            if d==M:
                opt = [scale(add(enter,cur[k+1]) if k<K else ztuple(), beta/(s+beta))
                       for beta in (Q(M-c,M),Q(1))]
            else:
                opt = [scale(add(source,scale(grow,r)),1/(s+b+r))
                       for r in (Q(L),L+Q(c,2))]
            if upper:
                cur[k] = tuple(Q(-(-max(z)*10**23//1),10**23) for z in zip(*opt))
            else:
                cur[k] = tuple(Q(min(z)*10**23//1,10**23) for z in zip(*opt))
        nxt = cur
    return scale(nxt[0],HI if upper else LO)


def polynomial_checks(R, states):
    # Newborn and the bare root are the first two physical strata.
    bare=plain((2,0,(J,)))
    assert bare in states and not R.stencils[J]
    a,b = R.row(bare,32,5)
    assert a[1:]==b[1:]==ztuple()[1:]
    # Guard against accidentally centering at a reset birth or paying an
    # entire candidate class twice: first hV/P4/Yt each start in one arm.
    seen=set()
    for x in states:
        for f,xx,ex in R.stencils[x]:
            assert len(xx)==(2 if f=='hV' else 3) and len(xx)==len(ex)
            seen.add(f)
    assert seen==set(FAMILIES)


def first_pass(blob, cap, Ms, cutoff):
    states, arcs = arms(cap)
    R = Rewards(blob,states)
    polynomial_checks(R,states)
    print('witness',R.c.sign,'states',len(states),'nonzero-shape arms',sum(map(bool,R.stencils.values())),flush=True)
    polynomials = {x:R.polynomials(x) for x in states}
    inf = {x:p[0] if p else ztuple() for x,p in polynomials.items()}
    assert all(min(a)>=0 for a in inf.values())
    assert R.c.c.D==18 and R.c.c.K==17 and R.c.c.lam==(Q('.565673') if R.c.sign==1 else Q('.565775'))
    for s in (HI,):
        u = resolvent(states,arcs,inf,s)
        print('discount',fixed(s),'inert newborn arm',' '.join(map(fixed,u[J])),flush=True)
        for M in Ms:
            c = 5  # explicit conditional fixture (four entering children)
            A = clock(s,M,c)
            Q0,Q1 = [],[]
            pres = {x:stage_bounds(p,inf[x],M) for x,p in polynomials.items()}
            for x,p in polynomials.items():
                stage = pres[x]
                # Beyond the first overflow bucket, envelopes must retain
                # actual k. Get the linear slope from the two overflow rows.
                lu = tuple(zip(*stage[:-1]))
                slope = overflow_slopes(p,M) if p else (ztuple(),ztuple())
                du = add(u.get(v.plus(x),ztuple()),scale(u[x],-1))
                Q0.append((add(tuple(min(a) for a in zip(*lu[0])),du),
                           add(tuple(max(a) for a in zip(*lu[1])),du)))
                Q1.append(slope)
            a = tuple(min(Q(0),min(z[0][i] for z in Q0)) for i in range(3))
            b = tuple(max(Q(0),max(z[1][i] for z in Q0)) for i in range(3))
            aa = tuple(min(Q(0),min(z[0][i] for z in Q1)) for i in range(3))
            bb = tuple(max(Q(0),max(z[1][i] for z in Q1)) for i in range(3))
            k, n = c-1,c-1
            E0 = k*A[0]+A[1]
            # Sup (k−18)_+/k <= 1. For this c=5 fixture, additionally
            # (k−18)_+ <= (N)_2/28 when k≤4+N: N≥15 on the overflow.
            E1 = min(n*n*A[0]+(2*n+1)*A[1]+A[2],
                     ((n+2)*A[2]+A[3])/28)
            rej = scale(u[J],c*A[0]+A[1]); birth = scale(u[J],1/s)
            low = add(add(scale(a,E0),scale(aa,E1)),scale(rej,-1))
            high = add(scale(b,E0),scale(bb,E1))
            # A nonzero material lower probe on young arms, stopped at real
            # apex arity K+1. This is the same D18 objective, not extra credit.
            K=cutoff
            guard_min = band_price(states,arcs,inf,M,s,K,pres)
            fresh = birth_price(M,K,guard_min,1)
            at_least_five = birth_price(M,K,guard_min,5)
            guard_max = band_price(states,arcs,inf,M,LO,K,pres,True)
            upper = birth_price(M,K,guard_max,1,True)
            upper5 = birth_price(M,K,guard_max,5,True)
            assert all(a<=b for a,b in zip(fresh,upper)) and all(a<=b for a,b in zip(at_least_five,upper5))
            print(' M',M,'c',c,'new-arm ideal',' '.join(map(fixed,birth)),
                  'conditional error',tuple(map(downward,low)),tuple(map(upward,high)),
                  'stationary young-arm box',tuple(map(downward,fresh)),tuple(map(upward,upper)),
                  'post-birth k≥5 only',tuple(map(downward,at_least_five)),tuple(map(upward,upper5)),flush=True)
        # The uniform stationary comparison uses µ in the closed prior.
        # u_s and du_s are interval-enclosed by positivity of both values,
        # separately with the sign appropriate to their difference.
        u0 = resolvent(states,arcs,inf,LO)
        for M in Ms:
            X0,X1,X2,X3 = common_age(LO,M,3)
            Xover = min(X1+X2, (X3+2*X2)/17)
            all_rows=[]
            overflow=[]
            for x,p in polynomials.items():
                stage=stage_bounds(p,inf[x],M)[:18]
                du0=add(u.get(v.plus(x),ztuple()),scale(u0[x],-1))
                du1=add(u0.get(v.plus(x),ztuple()),scale(u[x],-1))
                all_rows.append((add(tuple(min(a[i] for a,b in stage) for i in range(3)),du0),
                                 add(tuple(max(b[i] for a,b in stage) for i in range(3)),du1)))
                overflow.append(overflow_slopes(p,M) if p else (ztuple(),ztuple()))
            def sides(rows):
                return (tuple(min(Q(0),min(z[0][i] for z in rows)) for i in range(3)),
                        tuple(max(Q(0),max(z[1][i] for z in rows)) for i in range(3)))
            a,b=sides(all_rows); aa,bb=sides(overflow)
            err_lo=add(add(scale(a,X1/M),scale(aa,Xover/M)),scale(u0[J],-(X0+X1)/M))
            err_hi=add(scale(b,X1/M),scale(bb,Xover/M))
            print(' M',M,'stationary comparison error uniformly in µ',
                  tuple(map(downward,err_lo)),tuple(map(upward,err_hi)),flush=True)


def overflow_slopes(p,M):
    arc,rows = p
    a,b = rows[18]
    c,d = rows[19]
    assert b==d
    low,high = [],[]
    for a,c,b in zip(a,c,b):
        q = add(c,scale(a,-1))[1:]
        assert c[0]==a[0]
        q1 = add(q,scale(b[1:],Q(1,2)))
        vals = (*bernstein(q,Q(1,M)),*bernstein(q1,Q(1,M)))
        low.append(min(Q(0),min(vals)))
        high.append(max(Q(0),max(vals)))
    return tuple(low),tuple(high)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('blob', nargs='+', type=Path)
    p.add_argument('--cap', type=int, default=7, choices=range(6,9))
    p.add_argument('--degrees', type=int, nargs='+', default=[32,64,128])
    p.add_argument('--arity', type=int, default=17, choices=range(5,18))
    a = p.parse_args()
    assert min(a.degrees)>=32
    for blob in a.blob:
        first_pass(blob,a.cap,a.degrees,a.arity)


if __name__=='__main__':
    main()
