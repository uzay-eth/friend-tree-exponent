"""Successful stock payoff, no common deadline, with exact rational budgets."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
from runpy import run_path


ROOT = Path(__file__).resolve().parents[1]
SHARP = run_path(str(ROOT/'sharp/check.py'))
FROZEN = run_path(str(ROOT/'check.py'))


def stock_cost(K):
    """Tail atom and upper inverse half moment for the Poisson-Gamma payoff."""
    weights = [Q(1,factorial(k)) for k in range(K)]
    e_lo = sum(Q((-1)**k,factorial(k)) for k in range(18))
    e_hi = e_lo+Q(1,factorial(18))
    # At shape one use sqrt(pi)<9/5; at shape r>=2 use Cauchy–Schwarz
    # and an upward integer grid root of 1/(r−1).
    def inv(r):
        if r == 1:
            return Q(9,5)
        return next(Q(i,1000) for i in range(1,1001) if (r-1)*i*i>=1000**2)
    return 1-e_lo*sum(weights), e_hi*sum(w*inv(K-k) for k,w in enumerate(weights))


def successful_terminal(theta,d,K):
    a = [Q(0)]
    for j in range(1,K):
        q=Q(j,d)
        den=j+q-theta-q*a[-1]
        assert den>0
        a.append(j/den)
    row=[Q(0)]*K+[Q(1)]
    for j in range(K-1,0,-1):
        row[j]=a[j]*row[j+1]
        q=Q(j,d)
        # Check the zero-boundary tilted equation explicitly.
        assert (j+q-theta)*row[j]==j*row[j+1]+q*a[j-1]*row[j]
    return [Q(-((-v*10**12)//1),10**12) for v in row]


def source_times(K,theta=Q(1,2)):
    def ret(d,theta,start=2):
        v=Q(0)
        for k in range(31,start-1,-1):
            b=k-1+Q(1,d)
            v=max(Q(1),v,(Q(1,k)+b*v)/(Q(1,k)+b-theta))
        return v
    row = successful_terminal(theta,20,K)+[Q(1)]*20
    rows={20:row}
    for d in range(19,1,-1):
        nxt = row
        row = [ret(d,theta)*nxt[1]]
        for j in range(1,d):
            gain = nxt[j+1]
            row.append(max(gain,(j*gain+Q(j,d)*row[-1])/(j+Q(j,d)-theta)))
        rows[d]=row

    def serial(k):
        # W/ host contribution: initial clock of rate at least 3/2+1/k.
        # Larger host-growth rate preempts the stage and pays the base state.
        a = Q(1,k) if k is not None else Q(0)
        rate = Q(3,2)+a
        payoff = rows[3][2]*(1+a)+rows[2][0]/2
        host = rows[2][1]
        trial = max(payoff/(rate-theta),host)
        return trial

    def parallel(k):
        return ret(1,theta,k)*rows[2][1] if k<32 else Q(0)

    return serial,parallel


def main():
    FROZEN['check']()
    K=16
    atom,inv=stock_cost(K)
    serial,parallel=source_times(K)
    # The host transform is a linear-fractional function of 1/k;
    # include the exact endpoint at infinite host degree.
    E_s=max(serial(2),serial(None))
    E_e=parallel(2)
    root=Q(1,650_000_000)
    val=root*root
    size=Q(8)/SHARP['exp_lo'](Q(24))
    eta_s=atom+4*inv*E_s*root+size
    eta_e=atom+4*inv*E_e*root+size
    core=Q(3,2)*(eta_s+eta_e)
    chain=eta_s+Q(4,3)*eta_e
    assert max(core,chain)<FROZEN['delta']
    assert eta_s<min(FROZEN['G'](2),FROZEN['G'](3))
    assert root*root < Q(1,2**58)
    assert 16*val < Q(1,2**52)
    founder=size+atom+SHARP['terminal_mgf'](Q(1,4),2,K)[2]*Q(1,2**13)
    assert founder<Q(127,1220)
    print('K, atom, inverse moment, serial and parallel transforms:', K,
          *(f'{float(x):.9g}' for x in (atom,inv,E_s,E_e)))
    print(f'PASS y <= {val}; retained loss < {float(core):.9g}, chain loss < {float(chain):.9g}')
    print(f'PASS founder loss < {float(founder):.9g}')

    # Hostile bare/direct obstruction to using the same uniform loss budget
    # even with no retry time (j=1 gives a stock-one exponential payoff).
    P,A,B = (FROZEN[k] for k in ('P','A','B'))
    s,C,u = P(2,(0,))
    assert s==Q(1,6)
    slack = C-A-max(B,s*u)
    assert 0<slack<Q(1,100000)
    print(f'CONTROL degree-two terminal, opposite reciprocal zero: slack {float(slack):.9g}')
    # In the credited degree-20 stock-one exit, replacing survival by H
    # loses at least 19/20 − e^-1, even at zero threshold. Do not cash there.
    assert Q(19,20)-Q(3,8)>Q(1,2)
    # At a stock-one entrance on the limiting birth scale, even the no-wait
    # Poisson-beta payoff loses e^-1(1-e^-y). This is not a physical upper.
    y=Q(1,100)
    immediate=Q(36787,100000)*y/(1+y)
    assert Q(1,2)*immediate>300*slack
    print(f'CONTROL y=.01: immediate stock-one half-weight loss > {float(immediate/2):.9g}')


if __name__ == '__main__':
    main()
