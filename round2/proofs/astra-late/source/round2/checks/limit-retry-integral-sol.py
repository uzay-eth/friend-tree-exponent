"""Exact checks for defective clock, age integration and fold obstructions."""
from fractions import Fraction as Q


def clock():
    r = Q(747, 10000)
    # q(w) and q(w)^2 from independent geometric power series.
    def coef(j):
        return Q(1, 2**(j+1)) + sum(Q(1, 2**(k+1))*Q(1, 2**(j-k+1)) for k in range(1,j))
    weights=[coef(j) for j in range(1,9)]
    remain=Q(3,4)-sum(weights)
    assert remain==Q(11,1024)
    # Integrate the defective density represented by exponential components.
    f=sum(w*j/(j-2*r) for j,w in enumerate(weights,1))+remain*9/(9-2*r)
    rho=(1+f/(1-Q(3,2)*r))/(2-r)
    assert rho < Q(99991,100000)
    sig=r/(1+2*r)
    assert sig==Q(747,11494) and sig-Q(13,278)==Q(14561,798833)
    assert 1/(1-sig)==Q(11494,10747)
    assert Q(13,200)/(1+2*Q(13,200))==Q(13,226)


def age():
    rows=0
    for b in (1,2,4,8,16,32):
        q=Q(1,b)
        u=1-q
        # (1+Z)^-1 series versus the integral-series for -log q.
        for k in (2,4,10,30):
            partial=sum((q*u**(j-1)/(j+1) for j in range(1,k+1)),Q(0))
            assert partial*u*u==q*sum((u**j/j for j in range(2,k+2)),Q(0))
            assert partial+u**k/(k+2)<=1
            rows+=1
    # A finite weighted reached-setup measure. Work with y=e^r and integer ages.
    y=Q(5,4)
    setups={0:Q(1,2),2:Q(1,3),6:Q(1,9),9:Q(1,27)}
    for t in range(15):
        late=sum(w for s,w in setups.items() if s>t)
        age_cost=sum(w/y**(t-s) for s,w in setups.items() if s<=t)
        assert late+age_cost <= sum(w*y**s for s,w in setups.items())/y**t
    return rows


def fold():
    # Uniform sqrt-location coordinate; every dyadic shell pays a fixed amount.
    for k in range(40):
        upper=Q(1,2**k)
        assert (upper-upper/2)/upper == Q(1,2)
    # Planted nested-mask error: in a uniform two-bit space the later event
    # is a subset of the earlier, not an independent new failure.
    first={0,1}; later={0}
    assert Q(len(first & later),4) > Q(len(first),4)*Q(len(later),4)


if __name__=='__main__':
    clock()
    print('PASS exact defective-clock and full-law exponent arithmetic')
    print(f'PASS {age()} geometric reciprocal partial sums and 15 integrated-age controls')
    fold()
    print('PASS 40 critical fold shells and planted nested-mask error')
