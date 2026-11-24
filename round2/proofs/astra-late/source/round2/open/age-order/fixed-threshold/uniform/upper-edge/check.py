"""Bounded exact fixtures for the outside-clock change of measure."""
from fractions import Fraction as Q
from itertools import product


def birth(tree, dest):
    n=len(tree)
    return tree[:dest]+(tree[dest]+(n,),)+tree[dest+1:]+((dest,),)


def weights(tree):
    return [sum(Q(1,len(tree[y])) for y in x) for x in tree]


def walk(tree, waits, f=0):
    if not waits:
        return 1
    R=Q(len(tree[0]))-Q(1,2)
    B=Q(len(tree))-R
    delta=R-weights(tree)[0]
    ref=actual=block=Q(1)
    for j in range(waits[0]):
        r=weights(tree)[0]
        assert R+j-r==delta
        ref*=(R+j)/len(tree)
        actual*=r/len(tree)
        block*=1-delta/(R+j)
        tree=birth(tree,0)
    assert B==len(tree)-(R+waits[0])
    w=weights(tree)
    assert sum(w[1:])==B+delta
    cells=0
    for dest in range(1,len(tree)):
        q=ref*(B/len(tree))*w[dest]/(B+delta)
        p=actual*w[dest]/len(tree)
        assert q*block*(1+delta/B)==p
        new=birth(tree,dest)
        assert R+waits[0]-weights(new)[0]>=delta
        cells+=walk(new,waits[1:],f+1)
    return cells


def harmonic_controls():
    count=0
    for r,g in product((Q(1,2),Q(5,2),Q(101,2),Q(10001,2)),(0,1,2,25,1000)):
        h=sum((1/(r+j) for j in range(g)),Q(0))
        lo=sum((1/(r+j+1) for j in range(g)),Q(0))
        sq=sum((1/(r+j)**2 for j in range(g)),Q(0))
        diff=1/r-1/(r+g)
        assert h-lo==diff
        assert sq<=(1+1/r)*diff
        count+=1
    return count


def geom_controls():
    for u in (Q(1,2),Q(1,10),Q(1,10**24)):
        theta=1-u
        mean=theta/u
        factorial_second=2*theta**2/u**2
        assert u*mean==theta
        assert u**2*(mean+factorial_second-mean**2)==theta
        # Prefix conditioning adds k to G_0 and hence exactly uk to uR_f.
        k=1//u
        assert 0<u*k<=1 and u*k>1-u


if __name__=='__main__':
    seed=((1,),(0,2),(1,3),(2,))
    cases=sum(walk(seed, waits) for waits in product(range(3),repeat=3))
    print(f'PASS {cases} exact stopped prefixes: block likelihood and reciprocal deficit')
    print(f'PASS {harmonic_controls()} reciprocal-telescoping controls, including zero and long waits')
    geom_controls()
    print('PASS scaled geometric moments and prefix shifts at u=1/2, 1/10, 10^-24')
