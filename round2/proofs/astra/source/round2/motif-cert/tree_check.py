"""End-to-end rational test of the signed generator, including FY and both physical reserves."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from pathlib import Path
from low_oracle import Local, INF, NB, bud, degree, gain, inv, probes
from oracle import ARITY, Tree, newborn_extensions, replace
from physical_check import trees
from tail import verify as verify_tail
import argparse
import subprocess


def cases(c,nmax):
    yield from trees(c,nmax)
    # The two FY support exits, degree seams, corner T1 and a physical double tail.
    for dp in (21,22,23,24,25,71,72,73,1399,1400,1401,2000003):
        for dw in (3,72,79):
            # Edges both above and below a donor; a cherry, a fork, a path and a Y.
            yield Tree((-1,0,1,2,2,4,4,5),[0,3,dp-3,4,dw-3,2,4,1])
    for k in (10,12,13,24,25,26):
        parents=[-1,0,1,2]+[3]*k
        for D in (max(28,k+1),72,81):
            yield Tree(parents,[0,4,29,D-k-1]+[1+(j%4) for j in range(k)])


def phi(t,c):
    return t.potential(c)+sum(map(c.f,t.s[1:]))


def gap(jump,inverse):
    return max(Q(0),jump)/2-inverse*jump


def check_one(t,c,target,collect,tally,e,lower,upper,plant):
    old=t.active(c)
    signed=c.sign
    base=phi(t,c)
    number=len(t.s)-1
    if not lower*number<=base<=upper*number: raise ArithmeticError("tree comparison failed")
    direct=Q(0)
    for x in range(1,len(t.s)):
        direct+=t.r[x]*(phi(t.after(x,False),c)-base)
        if t.b[x]: direct+=t.b[x]*(phi(t.after(x,True),c)-base)
    adjusted=defaultdict(Q)
    reserve=defaultdict(Q)
    edges={frozenset(v):v for f,v,_,_ in old.values() if f=="E"}
    # Every vertex-to-edge transfer has an actual incident edge, including the dummy.
    for x,s in enumerate(t.s[1:],1):
        key=frozenset([x])
        change=c.f(gain(s))-c.f(s)
        born=(c.f(bud(s))-c.f(s)+c.f(NB)+c.weight("E",(bud(s),NB))) if t.b[x] else 0
        adjusted[key]=signed*(t.r[x]*change+t.b[x]*born-target*c.f(s))
        for y in [t.p[x]]+t.ch[x]:
            if degree(s)<c.l2.d0:
                b=signed*(inv(t.s[y])-c.center(0 if y==t.p[x] else 1,s))*change
                adjusted[key]-=b;adjusted[frozenset([x,y])]+=b
            elif degree(s)<c.par["DF"]:
                reserve[key]+=gap(signed*change,inv(t.s[y]))
    for support,(f,vs,st,w) in old.items():
        adjusted[support]-=signed*target*w
        for pos,x in enumerate(vs):
            if x==0: continue
            state=st[pos]
            delta=signed*(c.weight(f,replace(st,pos,gain(state)))-w)
            adjusted[support]+=t.r[x]*delta
            if t.b[x]:
                post=replace(st,pos,bud(state))
                bdiff=c.weight(f,post)-w
                if ARITY[f]<4:
                    newborn,labels=newborn_extensions(f,post,pos)
                    bdiff+=c.weight(newborn,labels)
                adjusted[support]+=signed*t.b[x]*bdiff
            for y in [t.p[x]]+t.ch[x]:
                if y in support: continue
                extended=support|{y}
                iy=inv(t.s[y]); from_parent=y==t.p[x]
                if ARITY[f]==4:
                    if degree(state)<c.par["D4"]:
                        reserve[support]+=gap(delta,iy)
                    continue
                if extended in old:
                    # A true centering transfer is one item on both motif rows.
                    if f=="E": cent=4 if old[extended][0]=="C" else (2 if from_parent else 3)
                    else: cent=0 if from_parent else 1
                    move=(iy-c.center(cent,state))*delta
                    adjusted[support]-=move;adjusted[extended]+=move
                elif f=="E" and degree(state)<c.l2.d0:
                    ty=1 if pos==1 else (0 if from_parent else 2)
                    amt=signed*c.fy(ty,c.l2.id[state],iy)
                    donor=frozenset((x,y)); assert donor in edges
                    adjusted[support]-=amt
                    if plant!="fy" or ty!=2: adjusted[donor]+=amt
                    tally["FY incidences"]+=1
                    tally["nonzero FY incidences"]+=amt!=0
    actual=signed*(direct-target*base)
    if sum(adjusted.values())!=actual: raise ArithmeticError("complete physical drift transfer identity failed")
    total_gap=sum(reserve.values())
    if min(reserve.values(),default=0)<0: raise ArithmeticError("negative physical gap")
    if actual+total_gap>number*e: raise ArithmeticError("whole-tree additive bound with both physical credits failed")
    for key,v in adjusted.items():
        f,st=(("vertex",(t.s[next(iter(key))],)) if len(key)==1 else (old[key][0],old[key][2]))
        g=reserve[key]
        if g:
            tally["vertex credit" if f=="vertex" else "four-motif credit"]+=1
        reward=v+g*(2 if plant=="double-gap" else 1)
        record=(f,st)
        if record not in collect or reward>collect[record]: collect[record]=reward
    tally["trees"]+=1;tally["actual events"]+=sum(1+bool(b) for b in t.b[1:]);tally["local rows"]+=len(adjusted)


def run(path,low,high,target,e,lo,hi,nmax,plant):
    c=Local(path);records={};tally=Counter()
    for t in cases(c,nmax): check_one(t,c,target,records,tally,e,lo,hi,plant)
    corner_set=set();low_queries=set()
    for f,st in records:
        if f=="vertex":
            if degree(st[0])<c.par["DF"]: low_queries.add(("v",*st[0]))
        elif f=="E": low_queries.update(("g",*x) for x,w in c.corners(f,st))
        else: corner_set.update((f,x) for x,w in c.corners(f,st))
    lookup=probes(c,low,path,target,sorted(low_queries))
    keys=sorted(corner_set)
    inp="".join(f+" "+" ".join(map(str,x))+"\n" for f,x in keys)
    proc=subprocess.run([high,path,"--query"],input=inp,text=True,capture_output=True,check=True)
    rational=list(map(Q,proc.stdout.splitlines()))
    if len(rational)!=len(keys): raise ArithmeticError("short high corner response")
    high_lookup=dict(zip(keys,rational))
    tail=verify_tail(path,target)["cap"]
    for (f,st),r in records.items():
        if f=="vertex": cap=lookup["v",*st[0]] if degree(st[0])<c.par["DF"] else tail
        elif f=="E": cap=sum(w*lookup["g",*x] for x,w in c.corners(f,st))
        else:
            # The exact query returns a rational at the header exponent.
            cap=sum(w*high_lookup[f,x] for x,w in c.corners(f,st))+(c.lam-target)*c.sign*c.weight(f,st)
        if cap<r: raise ArithmeticError(f"credited physical row exceeds corners: {f} {st} gap={cap-r}")
    assert tally["nonzero FY incidences"] and tally["vertex credit"] and tally["four-motif credit"]
    print("complete physical tree audit",dict(tally),"physical types",len(records),"low probes",len(lookup),
          "GMP high corners",len(keys),"PASS")


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("snapshot");p.add_argument("--target",type=Q,required=True)
    for arg in ("cap","lo","hi"):p.add_argument("--"+arg,required=True,type=int)
    p.add_argument("--low",default="round2/outputs/mc-low");p.add_argument("--high",default="round2/outputs/mh-check")
    p.add_argument("--nmax",type=int,default=5);p.add_argument("--plant",choices=["fy","double-gap"])
    a=p.parse_args()
    run(a.snapshot,a.low,a.high,a.target,a.cap,a.lo,a.hi,a.nmax,a.plant)
