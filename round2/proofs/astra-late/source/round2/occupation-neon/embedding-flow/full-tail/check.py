"""Effective uniform-core budgets for every uncharged D18 generator term."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from pathlib import Path
import sys


LO, HI = Q(565673, 10**6), Q(565775, 10**6)
D0, LF, KF = 18, 14000, 180
SCALE = 10**11
NAMES = ('f', 'g', 'hvtab', 'hctab', 'p4tab', 'yttab', 'ybtab', 's3tab')
ORDERS = (1, 2, 3, 3, 4, 4, 4, 4)


def scan(paths):
    """One read-only reduction; no payload copies, fixed-size chunks."""
    import numpy as np
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'higher'))
    from extract import Cert

    for mode, path in zip((2, 1), paths, strict=True):
        c = Cert(path)
        assert (c.mode,c.D,c.LF,c.KF,c.S) == (mode,D0,LF,KF,SCALE)
        maxima = []
        for name in NAMES:
            flat = getattr(c,name).reshape(-1)
            top = 0
            for i in range(0,flat.size,2**20):
                a = flat[i:i+2**20]
                top = max(top, abs(int(a.min())), abs(int(a.max())))
            maxima.append(top)
        stock = c.f[:-2].reshape(2,KF+1,LF+1)
        assert 2*maxima[0]*LF < 2**63
        leaf = max(int(np.abs(np.diff(row)*np.arange(LF,dtype=np.int64)).max())
                   for row in stock.reshape(-1,LF+1))
        assert all(x<=y for x,y in zip(maxima,MAXIMA,strict=True)) and leaf<=LEAF_F
        print(mode, tuple(maxima), leaf, flush=True)


# Verified outward integer bounds in units of 1/SCALE, uniformly for each
# original blob. Replay --scan gives the raw maxima; no payload input needed
# for the mathematical budget replay.
MAXIMA = (167836722289, 26340094623, 10079031623, 13860049440,
          2377223860, 4095258789, 2871334141, 6263404753)
LEAF_F = 146044806770


def geom_moment(k, exponent):
    """E[(N+1)^exponent; N>k] for the Exp(LO)/Poisson mixture."""
    q = 1/(1+LO)
    # Stirling numbers converting powers to factorial powers.
    st = [[1]]
    for m in range(1,exponent+1):
        row = st[-1] + [0]
        st.append([i*row[i]+(row[i-1] if i else 0) for i in range(m+1)])
    moments = [sum(st[m][i]*factorial(i)/LO**i for i in range(m+1))
               for m in range(exponent+1)]
    return q**(k+1)*sum(comb(exponent,i)*(k+2)**(exponent-i)*moments[i]
                        for i in range(exponent+1))


def degree_tail(dyadic):
    # At J=2^j, product P_J <= (3/(J+2))^LO <= 2^-floor(9(j-2)/16).
    return Q(1,2**(9*(dyadic-2)//16))


def terms(t, discount, k, anchored=False):
    """(family, name, pattern order, multiplicity, absolute norm, frozen error*J)."""
    gain_f = Q(LEAF_F,SCALE)*(1+abs(t))
    free_slots = (1,2,3,3,2,3,2,3)
    for name,r,ticks,mfree in zip(NAMES,ORDERS,MAXIMA,free_slots,strict=True):
        m = Q(ticks,SCALE)*(1+abs(t))
        if anchored:
            m *= 2**mfree
        h = 108*m  # 18 times the Bernstein derivative bound.
        if r == 1:
            yield name,'discount',1,1,discount*m,0
            yield name,'leaf gain',1,1,gain_f,0
            yield name,'vertical gain',2,2,m,2*m
            yield name,'promotion',1,1,2*m,0 if anchored else 2*m*(k+1)
            if not anchored:
                yield name,'newborn f',1,1,m,m*(k+1)
            continue
        yield name,'discount',r,1,discount*m,discount*r*h
        yield name,'leaf gain',r,r,34*m, max(h,34*(r-1)*h)
        yield name,'internal gain',r,2*(r-1),m,r*h+(0 if anchored else 2*m)
        yield name,'external gain',r+1,r+1,m,r*h+2*m
        yield name,'promotion',r,r,2*m,2*r*h+2*m*(k+1)
        # Select the newborn role; its parent is in the old pattern.
        yield name,'newborn motif',r-1,r-1,m,(r-1)*h+m*(k+1)


def budget(t, discount, k, j, anchored=False, stationary=False):
    a = [geom_moment(k,r) for r in range(5)]
    p, d = degree_tail(j), 2**j
    arity = omit = frozen = Q(0)
    for _,name,r,count,b,c in terms(t,discount,k,anchored):
        if stationary and name!='discount':
            continue
        multiplicity = r*(k+1)**(r-1)
        arity += count*b*r*a[r-1]
        omit += count*b*multiplicity*p
        frozen += count*c*multiplicity*p/d
    return arity+omit, arity+frozen, arity


def select(t, discount, goal, scheme):
    # Allocate half the tolerance to the arity tail, then minimize the
    # degree exponent at that K. This is a reproducible choice, not an LP optimum.
    anchored = scheme == 'anchor'
    stationary = scheme.startswith('stationary')
    for k in range(1,300):
        arity = budget(t,discount,k,14,anchored,stationary)[2]
        if arity < goal/2:
            break
    for j in range(14,256):
        # The high tail for the anchored objective is dropped, not lumped;
        # the identical triangle-inequality constants also bound that error.
        ans = budget(t,discount,k,j,anchored,stationary)[int(not scheme.endswith('omit'))]
        if ans < goal:
            return k,j,ans,arity
    raise ValueError('increase search bounds')


def controls():
    # Compute geometric moments a second way using generating-function
    # differentiation through the factorial basis, with a different shift.
    def factorial_tail(k,r):
        q=1/(1+LO)
        a=[Q(1)]
        for _ in range(r):
            a=grow(a,k+2)
        return q**(k+1)*sum(c*factorial(i)/LO**i for i,c in enumerate(a))

    def grow(a,shift):
        b=[Q(0)]*(len(a)+1)
        for i,v in enumerate(a):
            b[i]+=(shift+i)*v
            b[i+1]+=v
        return b

    for k,r in product((0,1,17,80),range(5)):
        assert geom_moment(k,r)==factorial_tail(k,r)

    # Verify loss and limit in a cubic shell at two physical censuses;
    # separate checks for target vs donor high, including a gain seam.
    cs=[Q(x) for x in (2,-3,5,-1)]
    m=max(map(abs,cs)); h=108*m
    def weight(d):
        t=Q(18,d)
        return sum(comb(3,i)*t**i*(1-t)**(3-i)*v for i,v in enumerate(cs))
    for d in (18,19,127,14182,2**43+1):
        assert abs(weight(d)-cs[0])<=h/d
        assert abs(weight(d+1)-weight(d))<=h/(d*(d+1))
        assert (d-1)*abs(weight(d+1)-weight(d))<=34*m

    # Old motifs plus newborn ownership must equal a direct jump, also
    # when a donor is outside the motif and a prior leaf becomes a core.
    spec = ((-1,),(-1,0),(-1,0,1),(-1,0,0),(-1,0,1,2),
            (-1,0,1,1),(-1,0,0,2),(-1,0,0,0))
    from itertools import permutations

    def census(p):
        D=[0]*len(p)
        for i,v in enumerate(p[1:],1):
            D[i]+=1;D[v]+=1
        return D,[sum(v==i and D[j]>1 for j,v in enumerate(p)) for i in range(len(p))]

    def weights(p):
        d,k=census(p)
        cores=[i for i in range(len(p)) if d[i]>1]
        out={}
        for f,top in enumerate(spec):
            for ids in permutations(cores,len(top)):
                if not all(p[ids[j]]==ids[v] for j,v in enumerate(top) if j):
                    continue
                if (f in (3,5) and ids[-2]>ids[-1] or
                    f==7 and tuple(sorted(ids[1:]))!=ids[1:]):
                    continue
                out[f,ids]=sum((Q(k[i],d[i]+1) for i in ids),Q(f,7))
        return out

    checked=0
    for p in ((-1,0,0,1,1,3,0,6,2,6,7),(-1,0,0,1,3,4,5,5,4,2)):
        d,k=census(p)
        old=weights(p)
        for recv in range(len(p)):
            new=weights(p+(recv,))
            if d[recv]==1:
                promoter=p[recv]
                assert promoter>=0
                born={key for key in new if recv in key[1]}
                assert all(len(ids)==1 or promoter in ids for _,ids in born)
                touched={key for key in old if promoter in key[1]}
            else:
                born=set()
                touched={key for key in old if recv in key[1]}
            assert set(new)-set(old)==born and set(old)<=set(new)
            assert (sum(new.values())-sum(old.values()) ==
                    sum(new[x]-old[x] for x in touched)+sum(new[x] for x in born))
            for donor in (i for i in range(len(p)) if p[i]==recv or p[recv]==i):
                assert (sum(new.values())-sum(old.values()))/d[donor] == (
                    sum((new[x]-old[x])/d[donor] for x in touched)
                    +sum(new[x]/d[donor] for x in born))
                checked+=1
    print('PASS: 20 arity moments, 5 cubic degrees,',checked,'directed jump/ownership fixtures')


def sparse_payload_controls(paths):
    """Read only controls touched by a mixed-tail physical tree and successors."""
    from collections import Counter
    occ = Path(__file__).resolve().parents[2]
    sys.path[:0] = [str(occ/'higher'/'older'),str(occ/'higher')]
    from extract import Cert, SPEC, embedding, label, up, bud, NB
    from replay import Local

    # Six linked cores, then private stocks. Two tags straddle the D18 seam.
    p = [-1,0,1,2,3,1,2]
    stock = [19,1,1,1,14,1,1]
    kids = [[] for _ in p]
    for j,i in enumerate(p[1:],1):
        kids[i].append(j)
    x = tuple((l,len(k),int(i>0)) for i,(l,k) in enumerate(zip(stock,kids)))
    motifs = embedding(tuple(p))
    seen = Counter()
    for path in paths:
        c = Cert(path)
        e = Local(c)
        old = c.potential(kids,x)
        for i,u in enumerate(x):
            y = label(x,i,up(u))
            v = label(x,i,bud(u))+(NB,)
            nk = [k[:] for k in kids]+[[]]
            nk[i].append(len(x))
            gain = newborn = promotion = Q(0)
            for fam,ids in motifs:
                if i not in ids:
                    continue
                before = e.value(fam,tuple(x[j] for j in ids))
                gain += e.value(fam,tuple(y[j] for j in ids))-before
                promotion += e.value(fam,tuple(v[j] for j in ids))-before
            for fam,ids in embedding(tuple(p)+(i,)):
                if len(x) in ids:
                    assert len(ids)==1 or i in ids
                    value = e.value(fam,tuple(v[j] for j in ids))
                    newborn += value
                    if value:
                        seen[fam]+=1
            assert c.potential(kids,y)-old==gain
            assert c.potential(nk,v)-old==promotion+newborn
    print('PASS: original blob gain and promotion comparisons:',2*len(x)*len(paths),
          'nonzero newborn families:',','.join(f for f in SPEC if seen[f]))


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--scan',nargs=2,type=Path,metavar=('LOWER','UPPER'))
    ap.add_argument('--sparse',nargs=2,type=Path,metavar=('LOWER','UPPER'))
    ap.add_argument('--t',type=Q,default=Q(1))
    ap.add_argument('--s-bound',type=Q,default=HI)
    ap.add_argument('--stationary-gap',type=Q,default=HI-LO,
                    help='upper bound on |mu-s| for stationary-only replay')
    args = ap.parse_args()
    if args.scan:
        scan(args.scan)
        return
    if args.sparse:
        sparse_payload_controls(args.sparse)
        return
    assert min(MAXIMA)>0 and LEAF_F>0 and args.s_bound>=0 and args.stationary_gap>=0
    controls()
    for goal in (Q(1,10**6),Q(1,10**8)):
        for scheme in ('omit','limit','anchor','stationary-omit','stationary-limit'):
            k,j,ans,arity = select(args.t,args.stationary_gap if scheme.startswith('stationary')
                                  else args.s_bound,goal,scheme)
            print(f'goal={float(goal):.0e} scheme={scheme}'
                  f' K={k} J=2^{j} ({2**j}) error={float(ans):.9e}'
                  f' arity={float(arity):.9e}')


if __name__ == '__main__':
    main()
