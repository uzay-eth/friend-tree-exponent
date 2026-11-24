"""Selected first promotion and virgin holding: conservative mixed-g pilot."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'higher/older'))
from replay import Cert, Local, NB, HI, LO  # noqa: E402


def price(e):
    # Newborn starts at D=2, k=0. Intrinsic gains from leaves have rate D-1.
    # Ignore beneficial exterior gains when lower-bounding the chosen path;
    # count their entire possible rate in the denominator.
    reach=Q(1); lower=Q(); floors=[]
    for d in range(2,32):
        l=d-1
        hazard=HI+l+Q(l,d)+Q(1,2)
        if d>=18:
            src=(d-2,1,1)
            # The only free role is grandparent send; the virgin has no
            # children, and its parent is included on the edge.
            raw=min(e.charge('g',(src,NB),((rho,0),(0,0))) for rho in (0,Q(1,2)))
            assert raw>=0
            floors.append(raw)
            # Kill on ANY jump of either endpoint. The parent receives
            # <= d-2 + 1/2 + 1/2, promotes at (d-2)/d; the child receives
            # 1+1/d and promotes at 1/2. Holding reward is uniform in rho.
            hold=HI+(d-1)+Q(d-2,d)+1+Q(1,d)+Q(1,2)
            lower+=reach*Q(l,d)/hazard*raw/hold
        reach*=l/hazard
    return min(floors),max(floors),LO*lower


def all_promotions(e):
    # Follow the same physical core until D=32. At each promotion, buy only
    # the first residence of the newly emitted edge. Distinct births give
    # distinct edges, even if their holding intervals overlap.
    grid=10**24
    value={}
    nonvirgin={}
    floors=[]
    for d in range(31,1,-1):
        for k in range(d-1,-1,-1):
            l=d-1-k
            b=Q(l,d)
            reward=Q()
            grown=Q()
            if l and d>=18:
                src=(l-1,k+1,1)
                raw=min(e.charge('g',(src,NB),((p,w),(0,0)))
                        for p,w in product((0,Q(1,2)),(0,Q(k,2))))
                assert raw>=0
                hold=HI+(l-1)+Q(k+2,2)+Q(l-1,d)+1+Q(1,d)+Q(1,2)
                reward=raw/hold
                floors.append((raw,d,k+1))
                # One wanted gain or promotion at the child, with parent
                # fixed, then stop at the next change at either endpoint.
                for child,rate in (((2,0,1),1+Q(1,d)),((0,1,1),Q(1,2))):
                    ell,j,_=child; D=sum(child)
                    raw=min(e.charge('g',(src,child),((p,w),(0,z)))
                            for p,w,z in product((0,Q(1,2)),(0,Q(k,2)),(0,Q(j,2))))
                    assert raw>=0
                    exit_rate=(HI+(l-1)+Q(k+1,2)+Q(1,D)+Q(l-1,d)
                               +ell+Q(1,d)+Q(j,2)+Q(ell,D))
                    grown+=rate/hold*raw/exit_rate
            continuation=b*(value.get((d,k+1),0)+reward)
            up=value.get((d+1,k),0)
            rates=(l,l+Q(k+1,2))
            v=min((g*up+continuation)/(HI+g+b) for g in rates)
            value[d,k]=Q(v*grid//1,grid)
            assert all((HI+g+b)*value[d,k]<=g*up+continuation for g in rates)
            continuation=b*(nonvirgin.get((d,k+1),0)+grown)
            up=nonvirgin.get((d+1,k),0)
            v=min((g*up+continuation)/(HI+g+b) for g in rates)
            nonvirgin[d,k]=Q(v*grid//1,grid)
            assert all((HI+g+b)*nonvirgin[d,k]<=g*up+continuation for g in rates)
    return sum(a>0 for a,_,_ in floors),len(floors),max(floors),LO*value[2,0],LO*nonvirgin[2,0]


if __name__=='__main__':
    for side in ('lo','hi'):
        e=Local(Cert(Path(f'b3f18_{side}.bin')))
        lo,hi,credit=price(e)
        print(side,'exact row floor span',lo,hi,'as floats',float(lo),float(hi),
              'selected first-promotion price per core',float(credit))
        assert credit>0
        n,total,(peak,d,k),credit,grown=all_promotions(e)
        print(side,'all-promotions cut price per core',float(credit),'exact',credit,
              'positive corners',n,'/',total,'peak floor',float(peak),'at D,k',d,k)
        print(side,'nonvirgin one-step price per core',float(grown),'exact',grown)
