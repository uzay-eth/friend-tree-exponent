"""Check elimination against an independent finite killed first-event solve."""
from fractions import Fraction as F
from dormant_check import solve
import excursion_rows as ex
from transfer_check import A,BETA,b,info,mark

ex.N=5
s=F('0.565615')
last=18
multi=ex.retained()
stars=[ex.star(a,k) for a in (2,16) for k in range(2,last+1)]
states=multi+stars
idx={t:i for i,t in enumerate(states)}
mi={t:i for i,t in enumerate(multi)}
cols=(A,BETA)


def charge(t,g):
    return int((info(t)[0]<16)==(g==A))


raw,pays=[],[]
for t in states:
    D,n=info(t)
    row=[F(0)]*len(states)
    row[idx[t]]=s+mark(t)+n-F(1,D)
    for nxt,r in ex.successors(t).items():
        if nxt in idx:
            row[idx[nxt]]-=r
    raw.append(row)
    pays.append([b(t)*charge(t,g) for g in cols])
original=solve(raw,pays)

# Each star is eliminated from its own gain/bud rates, independently of the solved rows.
star={}
for t in stars:
    a,L,_=t
    j=1+L
    lam=s+mark(t)
    prob=F(j-1,j)/(lam+j-1+F(j-1,j))
    coeff={}
    reward=[F(0),F(0)]
    for k in range(j,last+1):
        nxt=ex.promoted(a,k)
        reward[int(k>=16)]+=prob
        if nxt in mi:
            coeff[mi[nxt]]=prob
        prob*=ex.ratio(lam,k)
    star[t]=(coeff,reward)

matrix,rhs=[],[]
for t in multi:
    D,n=info(t)
    row=[F(0)]*len(multi)
    row[mi[t]]=s+mark(t)+n-F(1,D)
    reward=[b(t)*charge(t,g) for g in cols]
    for nxt,r in ex.successors(t).items():
        if nxt in mi:
            row[mi[nxt]]-=r
        elif nxt in star:
            cc,h=star[nxt]
            for j,v in cc.items():
                row[j]-=r*v
            for j in range(2):
                reward[j]+=r*h[j]
    matrix.append(row)
    rhs.append(reward)
small=solve(matrix,rhs)
for t in states:
    if t in star:
        coeff,x=star[t]
        result=[x[i]+sum(v*small[j][i] for j,v in coeff.items()) for i in range(2)]
    else:
        result=small[mi[t]]
    assert result==original[idx[t]]
j=ex.star(2,2)
assert original[idx[j]][0]>0 and original[idx[j]][1]>0
print('finite reference star degree limit',last,'multicore size limit',ex.N)
print('raw states',len(states),'Schur states',len(multi),'exact row-column equalities',2*len(states))
print('raw positive high-label entry',original[idx[j]][1])
