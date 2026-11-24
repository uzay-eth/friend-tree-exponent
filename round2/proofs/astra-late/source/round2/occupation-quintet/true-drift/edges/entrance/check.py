"""Upper price of a live marked child at the actual host crossing."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
from subprocess import check_output

edge = run_path(str(Path(__file__).resolve().parents[1]/'check.py'))
m = edge['m']
M, K, S, CAP, LO, HI = (edge[k] for k in ('M','K','S','CAP','LO','HI'))
GRID = 10**13
X = 10**6
LOW, HIGH = 565702, 565749
E = 24  # marked child degree cutoff; beyond it retain a constant terminal price
C = 100 # host arity cutoff for an arrival bound, paid separately


def ceil(a,b):
    return -(-a//b)


def terminal(pair):
    e = pair.e
    slots = [e.ids[e.d,c] for c in range(1,e.k+2)]
    assert all(pair.word(i,j)==0 for i in slots+[e.inf] for j in slots+[e.inf])
    price, floor = {}, {}
    for i,(d,c) in enumerate(e.labels):
        if d>=e.d:
            continue
        t = [-(M-e.d)*pair.word(e.inf,i)-e.d*pair.word(j,i) for j in slots]
        assert min(t)>=0 and max(t)<=M*S*CAP
        price[d,c] = ceil(max(t)*GRID,M*S)
        floor[d,c] = min(t)*GRID//(M*S)
    tail = max(v for (d,c),v in price.items() if d>E)
    return price, floor, tail


def child_values(price, tail, upper):
    """Envelope V_d at each actual registration degree d, with c_host<=K.

    Host mark is relaxed, not frozen. At focal gain the two marked cores
    send at their exact current reciprocals. Each gain changes one mark.
    """
    # Enlarged controller: all host arities [2,min(d,K)] and all exteriors.
    # Stop on host overflow in the actual problem. The upper comparison may
    # continue using any allowed arity because all terminal wages are >=0.
    inds = [(q,c) for q in range(E,1,-1) for c in range(q,0,-1)]
    nxt = {(q,c):v for (q,c),v in price.items() if q<=E}
    # The tail price is carried as a fictitious constant, discounted through
    # the remaining host steps at a host rate <= D. It is not paid early.
    tail_v = tail
    born = [0]*(M+1)
    for d in range(M-1,1,-1):
        tail_cur = ceil(tail_v*d*X,d*X+LOW)
        cur = {}
        dh = 2*d
        hc = min(d,K)
        for q,c in inds:
            # Common denominator 2*d*q for hazards. Host gain includes child
            # 1/q; child gain includes host 1/d; the other inputs may be adaptive.
            scale = dh*q
            # If c_host>=2, the prescribed edge is to a child and the host
            # has its own mother. Use min=0+1/q, max from c_host>=2.
            hg0 = (d-hc)*scale+dh
            hg1 = (d-2)*scale+dh+d*q
            cg0 = (q-c)*scale+2*q
            cg1 = cg0+(c-1)*d*q
            bp = (q-c)*dh
            succ = cur[q+1,c] if q<E else tail_cur
            pro = cur[q,c+1] if c<q else 0
            hsucc = nxt[q,c]
            # For the lower bound, kill on any subsequent host promotion
            # (hazard <=1); this preserves the registration's actual arity.
            n, den = corner(hg0,hg1,hsucc,cg0,cg1,succ,bp,pro,scale,upper,
                            0 if upper else scale)
            v = ceil(n,den) if upper else n//den
            cur[q,c] = v
        born[d] = cur[2,1]
        nxt, tail_v = cur, tail_cur
    return born


def corner(g0,g1,x,h0,h1,y,b,z,scale,upper,kill=0):
    corners = [(X*(g*x+h*y+b*z), X*(g+h+b+kill)+(LOW if upper else HIGH)*scale)
               for g in (g0,g1) for h in (h0,h1)]
    n,den = corners[0]
    for a,d in corners[1:]:
        if (a*den>n*d if upper else a*den<n*d):
            n,den = a,d
    v = ceil(n,den) if upper else n//den
    assert all(v*d>=a if upper else v*d<=a for a,d in corners)
    return n,den


def host_sum(born, upper):
    """Price the sum over actual pre-crossing promotion registrations."""
    nxt = [0]*(C+2)
    for d in range(M-1,1,-1):
        cur = [0]*(C+2)
        for c in range(min(C,d),0,-1):
            b = 2*(d-c)
            # At c=C exit, still pay this promotion's marked child, and
            # bound later registrations by arity overflow outside the sweep.
            z = born[d]+cur[c+1]
            n, den = corner(2*d*(d-c),d*(2*d-c),nxt[c],0,0,0,b,z,2*d,upper)
            cur[c] = ceil(n,den) if upper else n//den
        nxt = cur
    return Q(nxt[1],GRID)


def fixtures(pair, price, floor):
    # Host, marked child, mother, and a second child of the host. Add a
    # child of the mark too, so neither gain rate is purely one-way.
    for D,q in ((1998,2),(1999,3),(1999,E),(1999,E+1),(1999,180)):
        graph = m['physical_tree']([(0,1),(0,2),(0,3),(1,4)], [D,q,5,6,7])
        core, marks = m['physical_core'](graph)
        vals = []
        for receiver, donor in ((0,1),(1,0)):
            d,c = marks[receiver]
            exact = sum(Q(1,len(graph[u])) for u in graph[receiver])
            assert d-c+Q(1,len(graph[donor])) <= exact <= d-c+Q(1,len(graph[donor]))+Q(c-1,2)
            vals.append(exact)
        weight = -pair.g((M,marks[0][1]),marks[1])
        if q<=E:
            assert Q(floor[marks[1]],GRID)<=weight<=Q(price[marks[1]],GRID)
        # At the seam the crossing flux is from all donor slots into the
        # host, including the tagged child's own arrow once.
        if D==M-1:
            direct = sum(weight/len(graph[u]) for u in graph[0])
            assert direct==vals[0]*weight
        if q==180:
            assert not weight
    print('PASS 5 marked-edge rate and crossing fixtures')


def main():
    assert LOW*HI.denominator <= HI.numerator*X <= HIGH*HI.denominator
    assert Q(LOW,X)<=LO<=HI<=Q(HIGH,X)
    assert check_output(['git','-C',str(m['WORK']/'source'),'rev-parse','HEAD'],text=True).strip()==m['PIN']
    pair = edge['Pair']()
    price, floor, tail = terminal(pair)
    fixtures(pair,price,floor)
    born = child_values(price, tail, True)
    born_lo = child_values(floor, 0, False)
    raw = host_sum(born, True)
    gap = CAP*m['geom_weight'](C,1)
    p_lo, p_hi, *_ = edge['newborn_bounds'](pair)
    mother = CAP*p_hi
    phys_child = HI*(raw+gap)
    phys_mother = HI*mother
    mass = max(Q(0),p_lo-m['geom_weight'](33,0))
    old_lower = -edge['seam_floor'](pair.verify_faces(),mass,34)
    lower = max(old_lower, LO*host_sum(born_lo, False))
    upper = phys_child+phys_mother
    old_upper = CAP*HI*m['common_age'](M)[0]+CAP*HI*m['common_age'](M)[1]
    for name,v in (('child terminal tail ceiling',Q(tail,GRID)),('child price at registration D=1999',Q(born[1999],GRID)),
                   ('child seam upper',phys_child),('unmapped mother seam upper',phys_mother),
                   ('arity gap',HI*gap),('seam magnitude lower',lower),('seam magnitude upper',upper),
                   ('width improvement',(old_upper-old_lower)/(upper-lower))):
        print(name,edge['fmt'](v))
    assert upper-lower < (old_upper-old_lower)/3
    assert Q('3.6674464e-4') <= lower <= upper <= Q('8.5300112e-3')
    print('payload words:',*(len(p.reads) for p in pair.ps))


if __name__=='__main__':
    main()
