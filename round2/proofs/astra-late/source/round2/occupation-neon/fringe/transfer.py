#!/usr/bin/env python3
"""Positive 19-mark crown/ancestry transfer for physical D18 charges."""
import argparse
from collections import Counter
from functools import cache
from hashlib import file_digest
from pathlib import Path
import subprocess
import sys

from charges import Cert, Charges, Q, SPEC
from arms import Young

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'sol-kernel-verifier'))
import verify as v

HI=Q('.565772677869'); LO=Q('.565675597658')
B=10**12
DEN=(Q('1.356833021311'),Q('1.367169995210'))
LAMS=(Q('.565673'),Q('.565775'))
HASHES=('9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
        '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41')
FBASE=(Q('.000002512846393'),Q('.000002025752646'))
HBASE=(Q('.000001276902849680'),Q('.000001501029458674'))


@cache
def geometry(s):return (0,s[1],tuple(sorted(map(geometry,s[2]))))


@cache
def geometrical(s):return 1+sum(map(geometrical,s[2]))


def domain(n,degree,arm):
    seen=set(v.J); todo=list(seen)
    while todo:
        x=todo.pop()
        for y,_ in v.transitions(x):
            if y in seen:continue
            if v.size(y)>n and not (v.degree(y)<=degree and sum(map(v.size,y[2]))<=arm):continue
            seen.add(y);todo.append(y)
    return sorted(seen,key=lambda x:(-v.size(x),-geometrical(x),x))


class Reward:
    def __init__(self,c,gate=0,low_host=False):self.c=c;self.gate=gate;self.low_host=low_host

    @cache
    def own(self,x):
        # Choose a joint endpoint for all embeddings of the same family.
        out={} if self.gate and v.size(x)<=self.gate else self.c.rows(x)
        if self.low_host and v.degree(x)>=18:
            out={k:r for k,r in out.items() if k not in ('g','hC','S3')}
        return tuple(tuple(out.get(k,(Q(0),Q(0)))) for k in SPEC)

    def values(self,x):
        u=self.own(x)
        p=Q(1,v.degree(x)); n=v.size(x)
        child=[self.own(t) for t in x[2]]
        def get(first):
            a,b=map(sum,zip(*(u[:1] if first else u[1:])))
            val=a+min(0,b/2)
            for row in child:
                a,b=map(sum,zip(*(row[:1] if first else row[1:])))
                val+=b*p-min(0,b/2)
            return val
        # First body and all seven other motifs are disjoint charges.
        # For a gated channel the *charged apex's* size is tested, not
        # the parent crown size where the correction is paid.
        return get(True),get(False)


def compile_rows(binary,xs,certs,gate=0):
    ix={x:i for i,x in enumerate(xs)}
    channels=19+4*len(certs)
    lines=[f'{len(xs)} {channels} 36',' '.join(str(ix[x]) for x in v.J)]
    rewards=tuple((Reward(c),Reward(c,gate,True)) for c in certs)
    tables={}
    for i,x in enumerate(xs):
        den,j,b,intr,copy=v.row(x);den+=HI
        arcs=Counter()
        for block in (intr,copy):
            for y,r in block.items():
                if y in ix:arcs[ix[y]]+=r
        z=geometry(x)
        if z not in tables:tables[z]=tuple(t for rr in rewards for r in rr for t in r.values(z))
        forcing=[Q(0)]*19+list(tables[z]);forcing[j-2]=b
        # Enclose normalized first-step coefficients, avoiding an integer
        # LCM of the denominator at the exact 12-place discount and all
        # live degrees. Downwards on every positive path factor.
        row=[B,len(arcs)]
        for y,a in arcs.items():
            t=a*B//den
            assert t>0 and y!=i
            row.extend((y,t))
        for a in forcing:
            assert a>=0
            row.append(a*B**2//den)
        lines.append(' '.join(map(str,row)))
    p=subprocess.run([binary],input='\n'.join(lines)+'\n',text=True,capture_output=True)
    print(p.stderr.strip());p.check_returncode()
    ds=[[Q(int(k),B) for k in line.split()] for line in p.stdout.splitlines()]
    assert len(ds)==19 and all(len(r)==channels for r in ds)
    print('physical rows',len(xs),'charge geometries',len(tables),flush=True)
    return ds


def mixture(ds):
    # For the true stationary vector ell K_mu=ell and lower matrix M,
    # ell>=ell M, ell 1=1. Since ell(t-Mt)>=0 for t>=0, use min(a-t+Mt).
    import numpy as np
    from scipy.optimize import linprog
    M=[r[:19] for r in ds]
    I=[[Q(i==j)-M[i][j] for j in range(19)] for i in range(19)]
    ans=[]
    for col in range(19,len(ds[0])):
        a=[r[col] for r in ds]
        opt=linprog([*([0]*19),-1],A_ub=np.array([[*map(float,r),1.] for r in I]),
                    b_ub=[float(x*10**6) for x in a],bounds=[(0,None)]*19+[(None,None)],method='highs')
        assert opt.success
        t=[max(Q(0),Q(round(x*10**14),10**20)) for x in opt.x[:19]]
        z=max(min(a),min(ai-sum(b*x for b,x in zip(r,t)) for ai,r in zip(a,I)))
        assert min(a)<=z<=max(a)
        ans.append(LO*z)
        print('physical occupation column',col,'min/mix/max',v.bracket(LO*min(a),LO*min(a),15),
              v.bracket(LO*z,LO*z,15),v.bracket(LO*max(a),LO*max(a),15))
    print('lower matrix row sums',v.bracket(min(sum(r) for r in M),max(sum(r) for r in M)))
    return ans


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary');p.add_argument('low',type=Path);p.add_argument('high',type=Path)
    p.add_argument('--size',type=int,default=8);p.add_argument('--degree',type=int,default=22)
    p.add_argument('--arm',type=int,default=4)
    p.add_argument('--verify-sha',action='store_true')
    p.add_argument('--young',type=int,nargs=3,default=(120,14,8),metavar=('D','K','M'))
    a=p.parse_args()
    blobs=(a.low,a.high); stamps=tuple((x.stat().st_ino,x.stat().st_size,x.stat().st_mtime_ns) for x in blobs)
    if a.verify_sha:
        for x,digest in zip(blobs,HASHES):
            with x.open('rb') as file:assert file_digest(file,'sha256').hexdigest()==digest
        print('both pinned D18 SHA-256 checks passed',flush=True)
    certs=tuple(Charges(Cert(x)) for x in blobs)
    assert tuple(c.c.lam for c in certs)==LAMS and tuple(c.sign for c in certs)==(1,-1)
    assert all(c.c.stored==(3495091968-96)//8 and c.c.D==18 and c.c.S==10**11 for c in certs)
    ds=compile_rows(a.binary,domain(a.size,a.degree,a.arm),certs,12)
    vals=mixture(ds)
    end=[]
    for j in range(2):
        f=max(vals[4*j],FBASE[j])
        young=Young(certs[j],*a.young,degree3=True).solve()
        h=max(vals[4*j+1],HBASE[j]+vals[4*j+3]+young)
        s=LAMS[j]+(-1)**j*(f+h)/DEN[j];end.append(s)
        print('lambda side',j,'first and all higher (disjoint)',v.bracket(f,f,15),v.bracket(h,h,15),
              'uncapped crown young-arm',v.bracket(young,young,15),
              'outward endpoint',v.bracket(s,s))
    print('joint bracket',v.bracket(*end),'fraction of original width paid',float(1-(end[1]-end[0])/(LAMS[1]-LAMS[0])))
    assert stamps==tuple((x.stat().st_ino,x.stat().st_size,x.stat().st_mtime_ns) for x in blobs)


if __name__=='__main__':main()
