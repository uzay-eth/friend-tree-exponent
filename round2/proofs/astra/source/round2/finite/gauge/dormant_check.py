"""Exact finite killed replay for the binary active/dormant generator."""
from fractions import Fraction as F
from transfer_check import A, BETA, b, info, mark, newborn, physical, plus

S = F(2)
SHAPES = {(L, ()) for L in range(1,16)}
for k in (2,3,8,16):
    y = newborn(2 if k<16 else 16)
    for child in (y,plus(y),(y[0],0,(newborn(),))):
        SHAPES.add((k-2,(child,)))
ACTIVE = sorted((a,L,ch) for L,ch in SHAPES for a in (2,16))
DORMANT = [t for t in ACTIVE if t[0]==16]


def solve(rows, rhs):
    """Exact multiple-right-hand-side elimination, with sparse zero skipping."""
    n = len(rows)
    width = len(rhs[0])
    matrix = [row[:] + pay[:] for row,pay in zip(rows,rhs)]
    for j in range(n):
        pivot = next(i for i in range(j,n) if matrix[i][j])
        matrix[j],matrix[pivot] = matrix[pivot],matrix[j]
        divisor = matrix[j][j]
        matrix[j][j:] = [x/divisor for x in matrix[j][j:]]
        for i in range(n):
            if i==j or not matrix[i][j]:
                continue
            factor = matrix[i][j]
            for k in range(j,n+width):
                matrix[i][k] -= factor*matrix[j][k]
    return [row[n:] for row in matrix]


def matrix(dormant, full):
    states = [(False,t) for t in ACTIVE]
    if dormant:
        states += [(True,t) for t in DORMANT]
    index = {st:i for i,st in enumerate(states)}
    rows,rhs = [],[]
    for asleep,t in states:
        row = [F(0)]*len(states)
        row[index[asleep,t]] = S+(A if asleep else mark(t))
        for nxt,rate in physical(t):
            row[index[asleep,t]] += rate
            if (asleep,nxt) in index:
                row[index[asleep,nxt]] -= rate
        if asleep:
            row[index[False,t]] -= A-BETA
        else:
            p = 1/F(info(t)[0])
            for y in t[2]:
                for nxt,rate in ((y,mark(y)-p),(plus(y),p)):
                    if (False,nxt) in index:
                        row[index[False,nxt]] -= rate
            high = info(t)[0]>=16
            if dormant and high:
                row[index[True,newborn(16)]] -= b(t)
            if full:
                row[index[False,newborn(2 if dormant or not high else 16)]] -= b(t)
        high = info(t)[0]>=16
        rhs.append([F(0),F(0),F(0)] if asleep
                   else [b(t)*int(not high),b(t)*int(high),F(1)])
        rows.append(row)
    return dict(zip(states,solve(rows,rhs)))


if __name__=='__main__':
    plain = matrix(False,False)
    enlarged = matrix(True,False)
    kA,kB,_ = plain[False,newborn()]
    highA,highB,_ = plain[False,newborn(16)]
    dA,dB = highA-kA,highB-kB
    assert 0<dB<1 and dA>0
    assert enlarged[True,newborn(16)][1] == dB/(1-dB)
    scalar = kA+kB*(1+dA)/(1-dB)
    assert scalar == sum(enlarged[False,newborn()][:2])
    original_full = matrix(False,True)
    converted_full = matrix(True,True)
    for t in ACTIVE:
        assert original_full[False,t][2] == converted_full[False,t][2]
    print('spectral parameter',S,'active states',len(ACTIVE),'dormant states',len(DORMANT))
    print('positive excess, scalar return and full occupation: exact equalities',2+len(ACTIVE))
    print('finite killed scalar return',scalar)
