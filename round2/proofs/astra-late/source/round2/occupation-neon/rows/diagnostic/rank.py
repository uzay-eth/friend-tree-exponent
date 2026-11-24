"""Census of nonnegative charged residual classes on the raw diagnostic sample."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from pathlib import Path
import argparse
import sys

from probe import sparse, evaluator, snapshots, physical_view

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'higher/older'))
from replay import Cert, Local, NB  # noqa: E402


def band(x):
    d=sum(x)
    return next(b for b in (3,4,5,6,9,13,18,32,65,float('inf')) if d<b)


def group(f, ids, st, size, nc):
    a=ids[0]; x=tuple(st[j] for j in ids)
    if a==0: return ('founder',f)
    if f=='f': return ('f',min(x[0][1],5),band(x[0]))
    if f=='g': return ('g',band(x[0]),band(x[1]),x[1]==NB)
    return ('higher',f,'size<=12' if size[a]<=12 else f'{min(nc[a],5)}core',
            sum(y==NB for y in x))


def census(e, seed, runs, n):
    score=evaluator(e.c)
    mass=defaultdict(Q); hits=Counter()
    faces=defaultdict(Q); face_hits=Counter()
    phi=Q(); core=0; stream=[]
    for s in range(seed,seed+runs):
        par=next(snapshots(s,(n,)))
        ch,st,_=physical_view(par)
        p=[-1]*len(st)
        for a,children in enumerate(ch):
            for b in children:p[b]=a
        base,_,_,raw=sparse(e.c,score,par)
        q=sum(base)/e.c.S
        rows,_=e.tree(p,st)
        assert sum(v for _,_,v in rows)==sum(map(sum,raw)),(n,s)
        sz=[1+x[0] for x in st]; nc=[1]*len(st)
        for b in range(len(st)-1,0,-1):sz[p[b]]+=sz[b]; nc[p[b]]+=nc[b]
        assert sz[0]==n
        trial=defaultdict(Q)
        for f,ids,v in rows:
            key=group(f,ids,st,sz,nc)
            trial[key]+=v; hits[key]+=1
            if key[0]=='higher' and key[2:] == ('5core',0):
                degrees = [sum(st[j]) for j in ids]
                j=max(range(len(ids)),key=degrees.__getitem__)
                face=(f,j,band(st[ids[j]]),min(st[ids[j]][1],5))
                faces[face]+=v; face_hits[face]+=1
        for key,v in trial.items(): mass[key]+=v
        stream.append((q,trial)); phi+=q;core+=len(st)
    print(f'n={n} runs={runs} cores={core} charged/raw exact equality on every tree '
          f'total_R/phi={float(sum(mass.values())/phi):.10g}')
    for cell in ('f','g','higher','founder'):
        chosen={k:v for k,v in mass.items() if k[0]==cell}
        print(cell,'R/phi',float(sum(chosen.values())/phi))
        for k,v in sorted(chosen.items(),key=lambda x:x[1],reverse=True)[:10]:
            vals=[float(t.get(k,0)/q) for q,t in stream]
            print(' ',k,'hits',hits[k],f'R/phi {float(v/phi):.9g}',
                  f'mean {float(v/hits[k]):.9g}',f'per-tree [{min(vals):.6g},{max(vals):.6g}]')
    print('disjoint first-body slices R/phi by k and degree 2..17 / 18..64 / 65+:')
    for k in range(6):
        vals=[sum((v for key,v in mass.items() if key[0]=='f' and key[1]==k and lo<key[2]<=hi),Q())
              for lo,hi in ((0,18),(18,65),(65,float('inf')))]
        print(' ',k,*(f'{float(v/phi):.8g}' for v in vals))
    print('higher 5+ cores / no virgin slots: family, first max-degree slot (0=apex), degree upper-exclusive, min(k,5):')
    for k,v in sorted(faces.items(),key=lambda x:x[1],reverse=True)[:12]:
        print(' ',k,'hits',face_hits[k],f'R/phi {float(v/phi):.9g}',f'mean {float(v/face_hits[k]):.9g}')
    sys.stdout.flush()


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--seed',type=int,default=20261124)
    p.add_argument('--runs',type=int,default=8)
    p.add_argument('--sizes',type=int,nargs='+',default=[2048,16384])
    a=p.parse_args()
    print('SIMULATION DIAGNOSTIC, NOT A BOUND; seed0',a.seed,flush=True)
    for side in ('lo','hi'):
        print(side,flush=True)
        e=Local(Cert(Path(f'b3f18_{side}.bin')))
        for n in a.sizes:census(e,a.seed,a.runs,n)
        e.terms.cache_clear(); e.value.cache_clear()


if __name__=='__main__':main()
