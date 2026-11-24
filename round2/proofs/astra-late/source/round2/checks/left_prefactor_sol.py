#!/usr/bin/env python3
"""Exact negative controls for the large-star left-tail prefactor."""
from fractions import Fraction as Q
from math import isqrt


def p_depart(k):
    return Q(2, (k-2)*(k-1))


def entry_clock():
    # Until the centre rings, leaves grow by Yule births, independent of
    # that centre's Exp(1) clock. Conditional weights at deterministic
    # departure time t, where z=e^-t, are NB from two at size K-2.
    for z in (Q(1, 2), Q(1, 3), Q(3, 4)):
        assert (z*z, 2*z*z*(1-z)) != (p_depart(4), p_depart(5))
    # First two conditional Laplace transforms of the departure time,
    # using the actual free holding rates 3 through k-1, not k.
    def discounted(k, r):
        out = p_depart(k)
        for j in range(3, k):
            out *= Q(j, j+r)
        return out
    assert discounted(4, 1)/p_depart(4) == Q(3, 4)
    assert discounted(5, 1)/p_depart(5) == Q(3, 5)
    assert sum(map(p_depart, range(4, 129))) + Q(2, 127) == 1
    # Earliest exit S4: a path. Two end cores each have one leaf.
    # While waiting to reach size 5, *any* jump has rate 4, but only
    # promotion has initial rate 1. After each nonpromotion jump at
    # size 5 without changing C=2, the hazard is at least 7/6.
    hazards = (Q(1, 2)+Q(1, 2), Q(2, 3)+Q(1, 2))
    r = Q(13, 12)
    assert hazards[0] < r < hazards[1] < 4
    assert 4/(4-r) == Q(48, 35) and hazards[1]/(hazards[1]-r) == 14
    print('PASS: dependent departure time/size; Laplace marks 3/4 and 3/5; gated MGF factors 48/35 and 14')


def mixing():
    # Local embedded clock experiment, not a conditional physical law at
    # deterministic time: toss an independent fair choice on a degree-2
    # virgin tick (leaf or outward slot), then propose the next hub ring.
    # Between the two proposals there are no other changes in this toy.
    for d in range(2, 129):
        # Internal promotion yields two branch cores; a slot return gives
        # the physical hub one extra leaf but is null in the ideal branch.
        ok_return = 1-Q(1, d+1)
        ok_promote = 1-Q(1, d)
        posterior = ok_return/(ok_return+ok_promote)
        assert posterior-Q(1, 2) == Q(1, 2*(2*d*d-1))
        assert posterior > Q(1, 2)
        # On failure the ideal mark still exists; including both pieces
        # preserves exactly the original Bernoulli law.
        assert Q(1, 2)*ok_return + Q(1, 2)*(1-ok_return) == Q(1, 2)
    assert ok_return < 1
    # Conditioning a marked arrival sequence on accepted proposals can even
    # distort the count law: with no inter-proposal events, each successful
    # hub promotion preserves d and adds a core neighbour.
    d = 4
    survive = Q(1)
    for rings in range(1, 5):
        survive *= 1-Q(rings, d)
        if rings == 2:
            assert survive == Q(3, 8) != (1-Q(1, d))**2
    assert survive == 0
    print('PASS: success-only mark bias 1/14 at d=2; recurrent hub rejections, 3/8 != 9/16')


def physical_arrows():
    # Reconstruct physical directed source clocks from the labelled tree;
    # remove only the hub send and the outside effect of a return when
    # projecting to one ideal branch. Reach later states by actual arrows.
    def rates(edges):
        arcs = edges + [(y,x) for x,y in edges]
        return {(x,y): Q(1,sum(z==x for z,_ in arcs)) for x,y in arcs}
    for d in (2,3,8,32):
        tree = [(0,i) for i in range(1,d+1)]+[(1,d+1)]
        new = d+2
        a = rates(tree)
        assert (a[0,1],a[1,0],a[1,d+1],a[d+1,1]) == (Q(1,d),Q(1,2),Q(1,2),1)
        # An actual return 1->0 creates a new hub leaf; it leaves the
        # branch generator invariant and changes the next rejection rate.
        b = rates(tree+[(0,new)])
        assert b[0,1] == Q(1,d+1)
        assert (b[1,0],b[1,d+1],b[d+1,1]) == (a[1,0],a[1,d+1],a[d+1,1])
        # An actual hub promotion from the same departure seed: same hub
        # degree, one more core neighbour. The newborn is again a virgin.
        c = rates(tree+[(2,new)])
        assert c[0,1]+c[0,2] == Q(2,d)
        assert (c[2,0],c[2,new],c[new,2]) == (Q(1,2),Q(1,2),1)
    print('PASS: four exact departure-tree arrow projections, live return and new virgin immigration')


def normalization():
    # For r=2 the exact free Gamma factor is k²/((k-1)(k-2)).
    for k in range(4, 257):
        factor = Q(k*k, (k-1)*(k-2))
        assert 1 < factor <= Q(8, 3)
    # Independence is required of W and fixed-size D, not of W and Y.
    # Toy with mu=1/2: W/k ∈{1/4,7/4}, E(W/k)=1 and D≡1.
    w = (Q(1, 4), Q(7, 4))
    e = lambda xs: sum(xs)/2
    assert e(w) == 1 and e(1/z for z in w) == Q(16, 7)
    # If one redraws an independent W after observing Y, the inferred
    # E[D^-2] becomes E[Y^-2] E[W/k], visibly wrong.
    assert e(1/z for z in w)*e(w) != 1
    # Correlate the toy W with D^-2 without changing either marginal.
    inv = (Q(1),Q(3))
    actual = e(x/z for x,z in zip(inv,w))
    assert actual == Q(20, 7) != e(inv)*e(1/z for z in w)
    print('PASS: Gamma inverse factor endpoint 8/3; independent-redraw mutant 16/7 != 1; no marginal factorization')


def inverse_and_mixture():
    # Toy at mu=1/2, alpha=2. D_k=1 except for an atom k^-1/2
    # of probability k^-1. D_k⇒1, but critical inverse moments are
    # bounded and not UI. The rare atom starts contributing at k~sqrt(M)
    # when h=1/M, so it adds a constant to the departure mixture.
    for j in (4, 16, 32):
        M = j*j
        assert isqrt(M)==j and 1-Q(1, j)+Q(1, j)*j == 2-Q(1, j)
        atom = M*sum(p_depart(k)/k for k in range(j, 8*j+1))
        # 2/[k(k-1)(k-2)] telescopes, so this also checks the infinite sum.
        rest = Q(M, (8*j)*(8*j-1))
        assert atom+rest == Q(M, (j-1)*(j-2)) > 1
        # Normal atom appears at k>=M; subtract the mass switched to the
        # rare atom. Limit is 2+1=3 rather than the weak-limit prediction 2.
        total = Q(2*M, M-2)-Q(M, (M-1)*(M-2))+atom+rest
        assert 3 < total < 3+Q(9, j)
        # At k=j the negative moment of order 2r pays j^(r-1),
        # unbounded for every r>1; use the exact r=2 case.
        assert Q(1,j)*j*j == j
    print('PASS: critical-only inverse moment mutant changes scaled left-tail constant from 2 to 3')


if __name__ == '__main__':
    entry_clock()
    mixing()
    physical_arrows()
    normalization()
    inverse_and_mixture()
