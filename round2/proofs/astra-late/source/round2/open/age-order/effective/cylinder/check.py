"""Compact integer checks for the quantitative cylinder; scalar diagnostic fixtures."""

from fractions import Fraction as Q
from functools import cache
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def euler_hit(adj, slots, tag, degree, size, time, steps):
    """Reference rational evaluation of the prescribed finite transition formula.

    `adj` is the placed ordinary tree; `slots` counts its inert exterior
    slots. This deliberately small reference iterates the Euler kernel.
    """
    dt = Q(time, steps)

    @cache
    def value(tree, left):
        d = [len(ns)+(slots[i] if i < len(slots) else 0) for i,ns in enumerate(tree)]
        if d[tag]>=degree:
            return Q(1)
        if len(tree)>=size or not left:
            return Q(0)
        rates = [sum(Q(1,d[w]) for w in ns) for ns in tree]
        stay = 1-dt*sum(rates)
        assert stay>=0
        ans = stay*value(tree,left-1)
        n = len(tree)
        for i,r in enumerate(rates):
            nxt = tree[:i]+(tree[i]+(n,),)+tree[i+1:]+((i,),)
            ans += dt*r*value(nxt,left-1)
        return ans

    return value(tuple(map(tuple,adj)),steps)


def kernel_fixture():
    adj = ((1,), (0,2), (1,))
    # The middle gains at rate 2, or 1+1/3 with two inert slots on leaf
    # 0 (an ordinary representation of a boundary core with two slots).
    # Independent no-hit block at sizes 3 and 4: stay probabilities
    # 2/5 and 1/5 and transfer 1/5; the only further surviving mass
    # is already killed at size 5. For this row each promoted end
    # changes the middle gain intensity from 2 to 3/2.
    hit3 = [Q(0)]
    hit4 = [Q(0)]
    for _ in range(5):
        hit3.append(Q(2,5) + Q(2,5)*hit3[-1] + Q(1,5)*hit4[-1])
        hit4.append(Q(3,10)+Q(1,5)*hit4[-1])
    assert euler_hit(adj,(0,0,0),1,3,5,1,5) == hit3[-1] == Q(2437,3125)
    assert euler_hit(adj,(2,0,0),1,3,5,1,5) < euler_hit(adj,(0,0,0),1,3,5,1,5)
    # One step: starting from size three, the hit/size-four tie must count.
    assert euler_hit(adj,(0,0,0),1,3,4,Q(1,4),1)==Q(1,2)
    assert euler_hit(adj,(2,0,0),1,3,4,Q(1,4),1)==Q(1,3)
    print('finite Euler kernel: hit flag, inert denominators and cutoff tie checked')


def compact():
    beta = Q(97971,800000)
    d, m, t = 2**12, 2**56, 2**24
    T, R = d*t, 2**40
    assert (R-5) >= 8*T  # e < 4 gives factorial decay

    # A block either ends the gap or reaches the prescribed other-leaf stock.
    block = 2*d+128
    blocks = t//block
    assert blocks >= 60*32 and Q(32,31)**32 > 2
    assert 2**64 >= 2*m
    assert Q(20,17)**8 > 2
    lag = d*d*(Q(1,2**60)+Q(2,m+1))
    assert lag < Q(1,2**28)

    age_gap = Q(113,200)-Q(1,4)
    core_moment, degree_moment = Q(5,4)/age_gap, 2/age_gap
    assert 40*core_moment*max(core_moment,degree_moment) < 32**2
    assert Q(5,4)**4 > 2

    # A radius-six prefilter: core-degree cap 2^10, total degree cap 2^640.
    # Root plus six reversed-path terms, bounded each by the worst branching.
    pref = 7*2**50*(Q(512,2**128)+Q(32,2**80))
    assert pref < Q(1,2**22)

    # Outer cap C=2^60, D=2^E. Compare exponents only.
    E = 8*(60*R + 41 + 40)
    assert 2**57 >= 60*R+41+9+30  # C/8 dominates branching and 512.
    assert E//8 >= 60*R+41+5+30
    outer = Q(1,2**30)+Q(1,2**30)
    assert R-5 > 800+30  # localization: V<=2^800, factorial bound.
    assert 640+6*10+12 < 800
    a = E+60*(R+1)+10
    b = a+2*T+32
    # B=2^b; U=2^u is the Euler step count, giving (BT)^2/U=2^-32.
    u = 32+2*(b+36)
    assert T==2**36 and b-a-2*T==32 and u-2*(b+36)==32
    per_tag = Q(1,d)+lag+Q(2,2**30)+Q(4,2**32)
    error = 6*(pref+outer+per_tag)
    assert error < beta/2

    # Spatial k includes one level for reading boundary degrees.
    radius = R+1
    log2_j = 8+42+(radius+1)*(E+1)
    assert log2_j < 2**92
    s = 2**396
    assert s > 515*2**348
    assert s > 20*2**348*(982+2*radius)
    assert s > 25*(71+log2_j)
    assert error < Q('0.001467') and beta-error > Q('0.120996')
    print('weighted stationary error < 0.001467; required < 97971/1600000')
    print('pi L > 0.120996; certified N0 = 2^(2^396)')


def fixtures():
    base = load('physical_scalar', HERE.parents[1]/'preplacement-transfer/check.py')
    q = base.q(2,1,16)
    U = base.p(2,0,112)
    gap = q-U
    assert gap > Q(29,10000)

    def tip(d, l):
        return q/(1+d*d), base.p(d,l,112) if d < 112 else Q(1,l+1)

    # Root (L,c)=(L,2), with two adjacent old cores of specified (D,L).
    # Merge both orientations on root--old-leaf edges. Core--core edges
    # retain original per-root charge, as in G.
    def row(L, nbs):
        d=L+2
        w=(L+sum(Q(1,p) for p,_ in nbs))/(d+1)
        s=w/d
        ans=(w+Q(L,2*d))*gap
        for pn,ln in nbs:
            ans-=s*(tip(d+1,L)[1]-tip(pn+1,ln)[0])
        if L:
            coef=Q(1,2*d)-s
            a,b=tip(d+1,L-1)
            c,e=tip(2,0)
            ans+=L*coef*((a-e) if coef>=0 else (b-c))
        return ans

    for name, L, ns in [('bare/live',0,((2,0),(2,1))),
                        ('one-leaf/live',1,((2,0),(2,1))),
                        ('eight-leaf/live',8,((2,0),(2,1)))]:
        b=row(L,ns)
        # Floor to nanounits, outward.
        lo = (b*10**9)//1
        print(f'{name}: scalar lower {lo/10**9:.9f}; '
              f'additional lower score to reach zero <= {max(0,-lo)/10**9:.9f}')
        assert (b>0)==(L==8)


if __name__=='__main__':
    compact()
    kernel_fixture()
    fixtures()
