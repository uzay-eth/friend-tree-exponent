"""Bounded independent checks of initiator colors and conditional preplacement."""
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial


def add(tree, dest):
    n = len(tree)
    return tuple(row + (n,) if i == dest else row for i, row in enumerate(tree)) + ((dest,),)


def expand(tree, reds, tag, word, mask, blue_index=0):
    if not word:
        return Q(1), 1
    rd = bool(word[0])
    senders = [j for j in range(len(tree)) if (j in reds) == rd]
    good = lambda g, rr: sum(g[j] == (tag,) for j in rr)
    old = 1-Q(good(tree, reds), len(reds))
    total = drift = mass = Q(0)
    count = 0
    for src in senders:
        for dst in tree[src]:
            weight = Q(1, len(senders)*len(tree[src]))
            new = add(tree, dst)
            next_red = reds | {len(tree)} if rd else reds
            diff = 1-Q(good(new, next_red), len(next_red))-old
            drift += weight*diff
            total += weight
            if mask and not rd and blue_index < 2 and src == tag:
                continue
            p, rows = expand(new, next_red, tag, word[1:], mask, blue_index+int(not rd))
            mass += weight*p
            count += rows
    assert total == 1
    assert drift <= (0 if rd else Q(1, len(senders)*len(reds)))
    return mass, count


def colors():
    count = 0
    for base in (((1,), (0,)), ((1,), (0,2), (1,3), (2,)),
                 ((1,2,3,4), (0,), (0,), (0,), (0,))):
        for u in (0, len(base)-1):
            tag = len(base)
            tree = add(add(base, u), tag)
            b = len(tree)-1
            for word in product((0,1), repeat=4):
                p, n = expand(tree, {tag+1}, tag, word, True)
                h = min(word.count(0), 2)
                assert p == Q(b-1, b+h-1)
                count += n
    print(f'PASS {count} masked colored histories with reciprocal drift and exact prefix mass')


def reciprocal_and_zero_stock():
    for theta, f in product((Q(1,2), Q(1,3), Q(2,3)), range(1,5)):
        u = 1-theta
        p = u**(f+1)
        partial = prob = Q(0)
        for s in range(60):
            assert p == comb(s+f,f)*u**(f+1)*theta**s
            prob += p
            partial += p/(s+1)
            p *= theta*Q(s+f+1,s+1)
        exact = u*(1-u**f)/(theta*f)
        assert partial < exact <= partial+(1-prob)/61
    old = ((1,), (0,2), (1,3), (2,))
    u, v = 0, 1
    charges = [sum(Q(1,len(old[w])) for w in nbrs)/len(nbrs)/(len(nbrs)+1) for nbrs in old]
    assert charges[v] == charges[u]  # a path endpoint tie is not an adverse fixture
    old = ((1,), (0,2), (1,3,4), (2,), (2,))
    charges = [sum(Q(1,len(old[w])) for w in nbrs)/len(nbrs)/(len(nbrs)+1) for nbrs in old]
    assert charges[u] > charges[v]
    parallel = add(add(old, u), v)
    assert parallel[u] == (v,len(old)) and len(parallel[v]) == 3
    print('PASS 12 geometric reciprocal brackets and nonzero adverse zero-sibling endpoint')


def budget():
    log3 = sum(((-1)**(j+1)*Q(1,2)**j/j for j in range(1,8)), Q(0))
    c = Q(83,250)*Q(297,500)*Q(9,25)
    assert 1-log3 >= Q(297,500) and Q(4,11)*Q(254,255)>Q(9,25)
    assert 800*c*c > 4
    t = Q(5,6)
    log16 = sum(((-1)**(j+1)*Q(5,11)**j/j for j in range(1,12)), Q(0))
    e = sum(((-t)**j/factorial(j) for j in range(12)), Q(0))
    b = 232
    c = (b-1)/(Q(16,5)*b+1)*(1-log16/t)*e*(1-t*t/(2*(b-t)))
    assert Q(157,50)*232*c*c > 4
    print('PASS original 256 and tightened 232 cutoffs with rational bounds')


if __name__ == '__main__':
    colors()
    reciprocal_and_zero_stock()
    budget()
