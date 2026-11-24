"""All-size serial and favourable sharp tails using a uniform monitored clock."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from math import comb, factorial, lcm
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
BASE = run_path(str(HERE.parent/'check.py'))
S = 10**15


def ceil(x):
    return -((-x)//1)


def poisson_cdf(times, mass_lo, mass_hi, xlo, xhi, unit):
    """Directed Poisson CDFs through times-1 at a single clock threshold."""
    lo_total = hi_total = 0
    result = []
    for power in range(times):
        if power:
            mass_lo = mass_lo*xlo.numerator//(power*xlo.denominator)
            mass_hi = -((-mass_hi*xhi.numerator)//(power*xhi.denominator))
        lo_total += mass_lo
        hi_total += mass_hi
        result.append((min(lo_total,unit),min(hi_total,unit)))
    return result


def kernel(left, right, shapes, rate, times, mesh, minimum, maximum=None):
    """Lower enclosed normalized binomial masses at each virtual event time."""
    unit = 10**80
    ratio = Q(mesh+1,mesh)
    tlo = left if maximum is None else left/(1-left/maximum)
    thi = right/(1-right/minimum)
    xi = Q(1,2*mesh+1)
    hlow = 2 * (xi + xi**3/3 + xi**5/5)
    hhigh = hlow + 2*xi**7/(7*(1-xi*xi))
    step_lo = unit//ratio**rate
    step_hi = ceil(unit/ratio**rate)
    min_pow = max_pow = exp_lo = exp_hi = unit
    previous = [unit]*times
    result = [[0]*shapes for _ in range(times)]
    i=0
    while thi*(Q(max_pow,unit)-1) < minimum:
        i+=1
        next_min = min_pow*(mesh+1)//mesh
        next_max = -((-max_pow*(mesh+1))//mesh)
        z = thi*(Q(next_max,unit)-1)
        base = max(0,1-z/minimum)
        limit = BASE['exponential_lower'](z,unit)
        basis = [S*comb(minimum-1,v)*(tlo*min_pow/(minimum*unit))**v
                 *min(base**(minimum-1-v),limit)//1 for v in range(shapes)]
        exp_lo = exp_lo*step_lo//unit
        exp_hi = -((-exp_hi*step_hi)//unit)
        endpoints = poisson_cdf(times,exp_lo,exp_hi,i*rate*hlow,i*rate*hhigh,unit)
        for n,(lo,hi) in enumerate(endpoints):
            mass=max(0,previous[n]-hi)
            previous[n]=lo
            if mass:
                dest=result[n]
                for v,number in enumerate(basis):
                    dest[v]+=mass*number
        min_pow,max_pow=next_min,next_max
    return [[x//unit for x in row] for row in result]


def model(degree, neighbour, rate):
    stock=[('s',d,j) for d in range(2,degree) for j in range(1,d)]
    waiting=[('r',d,k) for d in range(1,degree) for k in range(2,neighbour)]
    host=[('h',2,k) for k in range(2,neighbour)]
    states=stock+waiting+host
    index={s:i for i,s in enumerate(states)}
    rows=[]
    for state in states:
        kind,d,j=state
        if kind=='s':
            death=Q(j,d)
            branches=[(Q(j),('s',d+1,j+1)),
                      (death,('s',d,j-1) if j>1 else ('r',d,2))]
            sources=[(1-death,j)]
        elif kind=='r':
            ret=Q(1,j)
            branches=[(ret,('s',d+1,1)),(j-1+Q(1,d),('r',d,j+1))]
            sources=[(ret,1)]
        else:
            ret=Q(1,j)
            next_host=('h',2,j+1) if j+1<neighbour else ('s',2,1)
            branches=[(1+ret,('s',3,2)),(Q(1,2),('r',2,2)),(j-Q(1,2),next_host)]
            sources=[(ret,2),(Q(1,2),1)]
        hold=rate-sum(r for r,_ in branches)
        assert hold>=0
        branches.append((hold,state))
        den=lcm(*(x.denominator for x,_ in branches+sources))
        moves=[(int(r*den),index[q]) for r,q in branches if q in index]
        extras=[(int(r*den),j) for r,j in sources]
        rows.append((moves,extras,den*rate))
    return index,rows


def prices(left,right,minimum=21,maximum=None,degree=14,neighbour=16,rate=18,times=110,mesh=128):
    assert rate<=minimum and max(degree,neighbour+1)<=rate and right<minimum
    index,rows=model(degree,neighbour,rate)
    values=[0]*len(rows)
    for basis in reversed(kernel(left,right,degree,rate,times,mesh,minimum,maximum)):
        jumps=[0]+[sum(basis[v]//factorial(j-1-v) for v in range(j))
                   for j in range(1,degree+1)]
        values=[(sum(k*values[a] for k,a in moves)+sum(k*jumps[a] for k,a in extras))//den
                for moves,extras,den in rows]
    e=BASE['POISSON_AT_ZERO']
    def tail(state,baseline):
        return Q((e*(baseline*S+values[index[state]]))//1,S)
    middle=tail(('s',2,1),1)
    fav={k:tail(('r',1,k),0) for k in range(2,neighbour)}
    parent={k:max(middle,tail(('h',2,k),1)) for k in range(2,neighbour)}
    return middle,fav,parent


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('left',type=Q)
    p.add_argument('right',type=Q)
    p.add_argument('--minimum',type=int,default=21)
    p.add_argument('--maximum',type=int)
    a=p.parse_args()
    mid,fav,host=prices(a.left,a.right,a.minimum,a.maximum)
    print('physical prices',a.left,a.right,'middle',mid,'parent 2',host[2],
          'parent 3',host[3],'favourable 2',fav[2],'favourable 3',fav[3])
