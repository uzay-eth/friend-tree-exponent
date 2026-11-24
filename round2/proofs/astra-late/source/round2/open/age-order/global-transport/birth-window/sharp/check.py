"""Rational deadline budgets for the physical birth-scale comparison."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path


def exp_lo(x, terms=160):
    value = term = Q(1)
    for j in range(1, terms+1):
        term *= x/j
        value += term
    return value


def return_mgf(d, theta, start=2):
    v = Q(1)
    for k in range(31,start-1,-1):
        rate = k-1+Q(1,d)
        value = (Q(1,k)+rate*v)/(Q(1,k)+rate-theta)
        assert value >= v >= 1
        v = value
    return v


def terminal_mgf(theta,d,K):
    # Pessimistic exponential-time potential at both exits 0,K.
    a,b = [Q(0)],[Q(1)]
    for j in range(1,K):
        q = Q(j,d)
        den = j+q-theta-q*a[-1]
        assert den > 0
        a.append(j/den)
        b.append(q*b[-1]/den)
    row = [Q(0)]*K+[Q(1)]
    for j in range(K-1,0,-1):
        row[j] = a[j]*row[j+1]+b[j]
    row[0] = Q(1)
    for j in range(1,K):
        assert row[j]>=1
        assert j*row[j+1]+Q(j,d)*row[j-1] == (j+Q(j,d)-theta)*row[j]
    # Round up once to keep the outer rational recursions small.
    return [Q(-((-v*10**12)//1),10**12) for v in row]


def time_mgf(theta,K):
    """Charge real positive-stock time and the capped virtual return races."""
    nxt = terminal_mgf(theta,20,K)
    for d in range(19,1,-1):
        row = [return_mgf(d,theta)*nxt[1]]
        for j in range(1,d):
            gain = nxt[j+1]
            row.append(max(gain,(j*gain+Q(j,d)*row[-1])/(j+Q(j,d)-theta)))
        if d == 3:
            three = row
        nxt = row
    # Charge either the initial parallel return or the initial serial host race.
    worst = max(nxt[0],nxt[1],three[1],three[2])
    assert all(return_mgf(1,theta,k)<=return_mgf(1,theta) for k in range(2,32))
    return max(return_mgf(1,theta)*nxt[1],worst/(1-theta))


def main():
    frozen = run_path(str(Path(__file__).resolve().parents[1]/'check.py'))
    frozen['check']()
    theta,K,T = Q(1,2),96,46
    time_error = time_mgf(theta,K)/exp_lo(theta*T)
    # Gamma(K,1) > 48 and Exp(1) < 16 give remaining mass > 32.
    mass_error = (1-Q(48,10000))**(-10000)/2**K + 1/exp_lo(Q(16))
    size_error = Q(2)**3/exp_lo(Q(24))
    eta = time_error+mass_error+size_error
    assert 6*eta < frozen['delta']
    assert eta < min(frozen['G'](2), frozen['G'](3), Q(1,400))
    founder_theta = Q(1,4)
    founder = terminal_mgf(founder_theta,2,K)[2]/exp_lo(founder_theta*T)+mass_error+size_error
    assert founder < Q(127,1220)
    # exp(46)<2^67 and 16*2^-66 = 32*2^-67.
    assert (1-Q(T,67000))**(-1000) < 2
    print('PASS accumulated time, terminal mass, size:',
          *(f'{float(x):.4g}' for x in (time_error,mass_error,size_error)))
    print(f'PASS signed endpoint error < {float(eta):.6g}, time cutoff {T}, stock {K}')
    print(f'PASS y <= 2^-66, including founders (founder error < {float(founder):.6g})')


if __name__ == '__main__':
    main()
