"""Whole-tree first-event and exhaustive finite-box checks for the joint backup."""
from fractions import Fraction as F
from itertools import product

from joint_bellman import S,UNIT,extremum,row_terms
from sector import BETA,hy,info,mark,newborn,physical,plus


def score(tree):
    degree,n=info(tree)
    return UNIT*F(n+tree[1]+7*(tree[0]==16),degree+1)


if __name__=='__main__':
    number=0
    boxes=0
    for tag in (2,16):
        for d in (7,15,16,24):
            j=newborn()
            h=(2,d-1,())
            trees=((tag,d-3,(j,plus(j))), (tag,0,((2,d-2,(j,)),)),
                   (tag,d-2,((2,0,(j,)),)), (tag,2,tuple(sorted((j,h)))))
            for t in trees:
                N,den,terms=row_terms(t,lambda y:(score(y),score(y)+UNIT))
                degree,n=info(t)
                p=F(1,degree)
                real=physical(t)
                for side in (0,1):
                    f=lambda y:score(y)+side*UNIT
                    value=UNIT*F(t[1],degree)*hy.rows.g(BETA,degree)
                    value+=sum(rate*f(y) for y,rate in real)
                    value+=sum((mark(y)-p)*f(y)+p*f(plus(y)) for y in t[2])
                    independent=S+mark(t)+sum(rate for _,rate in real)
                    assert value==N[side]+sum(x[side+1] for x in terms)
                    assert independent==den+sum(x[3] for x in terms)
                    number+=1
                    choices=[(x[side+1],x[3]) for x in terms]
                    candidates=[(N[side]+sum(flag*a for flag,(a,b) in zip(bits,choices)))/
                                (den+sum(flag*b for flag,(a,b) in zip(bits,choices)))
                                for bits in product((0,1),repeat=len(choices))]
                    assert min(candidates)==extremum(N[side],den,choices,False)
                    assert max(candidates)==extremum(N[side],den,choices,True)
                    boxes+=2
    print('whole-tree source and denominator identities',number,'box extrema against all vertices',boxes)
