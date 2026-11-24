"""Small exact controls for the fixed-interior change of measure."""
from fractions import Fraction as Q
from itertools import product


def grow(tree, target):
    n=len(tree)
    return tree[:target]+(tree[target]+(n,),)+tree[target+1:]+((target,),)


def next_rows(tree, a):
    n=len(tree); mark=0
    R=len(tree[mark])-1+a; B=n-R
    weights=[sum(Q(1,len(tree[src])) for src in nbrs) for nbrs in tree]
    r=weights[mark]; delta=R-r
    assert delta>=0 and sum(weights)==n
    for j,w in enumerate(weights):
        q=R/n if j==mark else B/n*w/(B+delta)
        ell=1-delta/R if j==mark else 1+delta/B
        assert q*ell==w/n
        yield grow(tree,j),q,ell,w/n,j==mark


def history_control(seed, a, depth):
    stack=[(seed,Q(1),Q(1),Q(1),())]
    cases=0
    for _ in range(depth):
        nxt=[]
        for tree,ref,lik,physical,word in stack:
            rows=list(next_rows(tree,a))
            assert sum(q for _,q,_,_,_ in rows)==1
            assert sum(q*l for _,q,l,_,_ in rows)==1
            for t,q,l,p,h in rows:
                assert ref*q*lik*l==physical*p
                nxt.append((t,ref*q,lik*l,physical*p,word+(h,)))
            cases+=len(rows)
        stack=nxt
    b=len(seed)-a
    for word in product((False,True),repeat=depth):
        R=a; B=b; expected=Q(1)
        for step in word:
            expected*=(R if step else B)/(R+B)
            R+=step; B+=not step
        assert sum(q for _,q,_,_,w in stack if w==word)==expected
    assert sum(q*l for _,q,l,_,_ in stack)==1
    return cases


def drift_control():
    cases=0
    for v, i, s, theta in product((4,17,100), (3,89), (0,1,2,3), (Q(1,5),Q(1,2),Q(4,5))):
        a=Q(1,2); R=a+s; B=v+i-R; z=R/(v+i)
        delta=Q(1,7)
        drift=theta*(-delta/R)+(1-theta)*delta/B
        assert drift==delta*(z-theta)/((v+i)*z*(1-z))
        assert drift<0 if z<theta else drift>=0
        # Clock derivative: t=log(1+(1-theta)y), exp(t)-1=(1-theta)y
        y=Q(i,v); et=1+(1-theta)*y
        assert (1-theta)/et-1/y==-(1-theta)/et/(et-1)
        cases+=1
    return cases


if __name__=='__main__':
    leaf_parent=((1,),(0,2),(1,3),(2,))
    core_parent=((1,),(0,2,3),(1,),(1,))
    n=history_control(leaf_parent,Q(1,2),4)+history_control(core_parent,Q(1,3),4)
    print(f'PASS {n} physical/ref arrow descendants and all depth-4 Pólya words')
    print(f'PASS {drift_control()} signed drift and clock identities')
