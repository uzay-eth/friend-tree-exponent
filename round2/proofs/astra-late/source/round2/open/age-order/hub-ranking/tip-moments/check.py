"""Rational replay of physical Gamma caps and the discounted promotion race."""
from fractions import Fraction as Q
from math import factorial, prod


def ap(a, p):
    return prod((a+j for j in range(p)), start=Q(1))/factorial(p)


def cap(k, p, sharp):
    a = 1+Q(1,k)
    plain = ap(a,p)
    assert plain == (a*ap(a+1,p)+plain)/(p+1+a)
    if not sharp:
        return plain
    race = (a*ap(a+1,p)+(ap(a-Q(1,2),p)+ap(1+Q(1,k+1),p))/2)/(p+1+a)
    assert 0 < race < plain
    return race


def bounds(p, sharp, stop=64):
    upper = cap(stop,p,sharp)/(1+p*stop)
    lower = Q(0)
    for k in range(stop-1,1,-1):
        b = cap(k,p,sharp)
        assert b >= cap(k+1,p,sharp)
        rate = Q(1,k)
        previous = upper
        upper = (rate*b+upper)/(p+1+rate)
        lower = (rate*b+lower)/(p+1+rate)
        assert previous <= upper <= b/(1+p*k)
        assert upper-lower <= Q(1,2)**(stop-k)
        if k in (2,3):
            yield k,upper


def physical(k):
    # Parent has k neighbours including tag. Before promotion the tag's own
    # arrow grows parent and only the parent's specific slot can promote tag.
    adj = [[1,*range(2,k+1)],[0]]+[[0] for _ in range(k-1)]
    incoming = lambda x: sum(Q(1,len(adj[y])) for y in adj[x])
    assert incoming(1)==Q(1,k)
    adj[0].append(len(adj)); adj.append([0])  # actual return arrow
    assert incoming(1)==Q(1,k+1)
    adj[1].append(len(adj)); fresh=len(adj); adj.append([1])
    assert incoming(1)==1+Q(1,k+1)
    adj[fresh].append(len(adj)); adj.append([fresh])  # tag promotes fresh leaf
    assert incoming(1)==Q(1,2)+Q(1,k+1)


if __name__=='__main__':
    budgets={(1,True):{2:387,3:265},(2,True):{2:302,3:185},
             (1,False):{2:414,3:284},(2,False):{2:325,3:199}}
    for (p,sharp),limits in budgets.items():
        for k,value in bounds(p,sharp):
            assert value < Q(limits[k],1000)
            print(f'PASS p={p} parent={k} sharp={sharp}: bound < {limits[k]/1000:.3f}')
    for k in range(2,66):
        physical(k)
    print('PASS 64 physical arrow fixtures and discounted Gamma-race identities')
