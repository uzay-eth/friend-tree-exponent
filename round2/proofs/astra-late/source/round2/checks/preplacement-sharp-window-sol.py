"""Independent rational exit-clock and birth-scale budgets."""

from fractions import Fraction as Q


def exits(d, tilt, K):
    # Compute the two Dirichlet solutions by forward shooting, not elimination.
    A,B = [Q(0),Q(1)],[Q(1),Q(0)]
    for j in range(1,K):
        a = 1+Q(1,d)-tilt/j
        A.append(a*A[-1]-A[-2]/d)
        B.append(a*B[-1]-B[-2]/d)
    assert A[-1]>0
    x=(1-B[-1])/A[-1]
    values=[a*x+b for a,b in zip(A,B)]
    assert values[0]==values[-1]==1 and min(values)==1
    for j in range(1,K):
        assert (j+Q(j,d)-tilt)*values[j]==j*values[j+1]+Q(j,d)*values[j-1]
        h=lambda j: (1-Q(1,d)**j)/(1-Q(1,d)**K)
        assert j*(h(j+1)-h(j))+Q(j,d)*(h(j-1)-h(j))==0
    return values


def whole_return(d, tilt):
    # Forward measure through the nominal race, stopping on either exit.
    out=Q(0)
    alive=Q(1)
    for k in range(2,32):
        send=Q(1,k)
        gain=k-1+Q(1,d)
        den=send+gain-tilt
        assert den>0
        out+=alive*send/den
        alive*=gain/den
    return out+alive


def full(tilt, terminal):
    V={(20,j):v for j,v in enumerate(terminal)}
    for d in range(19,1,-1):
        V[d,0]=whole_return(d,tilt)*V[d+1,1]
        for j in range(1,d):
            g,l=V[d+1,j+1],V[d,j-1]
            V[d,j]=max(g,(d*g+l)/(d+1-d*tilt/j))
            assert j*(g-V[d,j])+Q(j,d)*(l-V[d,j])+tilt*V[d,j]<=0
    return max(whole_return(1,tilt)*V[2,1],max(V[k] for k in ((2,0),(2,1),(3,1),(3,2)))/(1-tilt))


def exp_bounds(x):
    t=out=Q(1)
    for k in range(1,161):
        t*=x/k
        out+=t
    ratio=x/161
    assert 0<ratio<1
    return out,out+t*ratio/(1-ratio)


def main():
    T,K=46,96
    tilt=Q(1,2)
    retry=full(tilt,exits(20,tilt,K))/exp_bounds(tilt*T)[0]
    founder=exits(2,Q(1,4),K)[2]/exp_bounds(Q(T,4))[0]
    mass=exp_bounds(Q(48))[1]/2**96+1/exp_bounds(Q(16))[0]
    size=8/exp_bounds(Q(24))[0]
    assert exp_bounds(Q(T))[1]<2**67
    assert retry<Q('4.492e-7')
    eta=retry+mass+size
    assert eta<Q('5.72e-7') and 6*eta<Q(1,200000)
    assert founder+mass+size<Q('.000083')<Q(127,1220)
    print('PASS 190 exact transform equations and 190 hit-probability equations')
    print(f'PASS independent bounds: retry {float(retry):.9g}; loss {float(eta):.9g}; founder loss {float(founder+mass+size):.9g}')
    print('PASS strict mass threshold 2^-66; both frozen owner budgets fit')


if __name__=='__main__':
    main()
