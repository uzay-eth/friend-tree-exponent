#!/usr/bin/env python3
"""Immutable seven-core regression of the withdrawn first-body witness."""
from pathlib import Path
import sys

from replay import Census, D, Charges, Q, NB, shape, label, up, bud
from check import drift

PARENTS=(-1,0,1,2,1,0,5)
LEAVES=(62,63,0,1,4,65,1)
TRIAL=Q('.565770220167')


def main(path):
    c=D.Cert(path);assert c.mode==1
    cs=Census(c)
    floors=list(cs.floors);floors[1]*=1000
    margin,_,z,tail=cs.solve(TRIAL,{},floors=floors,legacy=True)
    # Match original rounding, zero cutoff and founder treatment exactly.
    children=[[j for j,p in enumerate(PARENTS) if p==i] for i in range(len(PARENTS))]
    def tree(i):return 0,LEAVES[i],tuple(map(tree,children[i]))
    p,x=shape(tree(0),0);assert p==PARENTS
    delta=c.lam-TRIAL
    def correction(x):
        out=-delta*cs.M/TRIAL
        for y in x[1:]:
            d,k=sum(y),y[1]
            out+=z[d,k] if (d,k) in z else tail[min(k,c.KF+1)]
            out-=z[2,0]
        return out
    old=correction(x);dr=-TRIAL*old
    for i,y in enumerate(x):
        ds=children[i]+([p[i]] if i else [])
        g=y[0]+sum(Q(1,sum(x[j])) for j in ds)
        dr+=g*(correction(label(x,i,up(y)))-old)
        if y[0]:dr+=Q(y[0],sum(y))*(correction(label(x,i,bud(y))+(NB,))-old)
    base=-drift(c,tree(0))-delta*c.potential(children,x)/c.S
    bad=base+dr
    assert margin>0 and bad<Q('-0.0008')
    print('parents',PARENTS,'leaves',LEAVES,'s',TRIAL)
    print('claimed local margin',margin,'base signed drift',base,'correction drift',dr)
    print('actual signed upper residual',bad, 'decimal',float(bad))


if __name__=='__main__':main(Path(sys.argv[1]))
