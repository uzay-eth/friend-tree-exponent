#!/usr/bin/env python3
"""Project the frozen conditioned f residual onto the common sided census."""
import argparse
import ast
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from json import dumps
from pathlib import Path

HERE = Path(__file__).resolve().parent
FROZEN = 'cca747e9f24479f200d30fdfaf58f15f6022de65158649273ca3b961c034b5f7'
D, C = 360, 18
ROOTS = {'lower': (1322979352773609470, 1020282677750521357),
         'upper': (1755696543328942415, 1457540240862028342)}


def load():
    path=HERE.parent/'replay.py'
    raw=path.read_bytes()
    assert sha256(raw).hexdigest()==FROZEN
    spec=spec_from_file_location('frozen_conditioned',path)
    m=module_from_spec(spec);spec.loader.exec_module(m)
    # Add only an observer after each committed value; no Bellman algebra or
    # successor order is replaced, and no previously published file is edited.
    tree=ast.parse(raw)
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='conditional')
    slab=next(n for n in fn.body if isinstance(n,ast.For))
    loop=slab.body[0].body[-1].body[0].body
    assert ast.unparse(loop[-2])=='lookup[d, c, m, t] = (min if up else max)(v, W[d][c])'
    loop.extend(ast.parse('_capture[d,c,m,t] = lookup[d,c,m,t]').body)
    m._capture={}
    exec(compile(ast.fix_missing_locations(ast.Module([fn],[])),str(path),'exec'),vars(m))
    return m


def corners(m,s,W):
    vals=m._capture; G=m.G
    def val(d,c,p,t):
        return Q(vals[d,c,p,t] if d<=D and c<=C else W[d][c],G)
    out=defaultdict(dict)
    num=0
    for (d,c,p,t),raw in vals.items():
        j,i=t; k=c-1; L=d-c; V=Q(raw,G)
        a,b=s.row_ticks(d,c); a+=m.base.GOLD[s.side]['cap_f']
        gain=val(d+1,c,p,t)-V
        a+=L*gain-m.base.HI*V; b+=gain
        if L: a+=Q(L,d)*(val(d,c+1,p,(2,1))-V)
        if p==2:
            jump=val(d,c,-3,t)-V
            a+=(Q(1,d)+(Q(1,2) if jump<0 else 0))*jump
        if i:
            jump=val(d,c,p,(j+1,i) if j<3 else (4,0))-V
            a+=(j-i+Q(1,d)+(Q(i-1,2) if jump<0 else 0))*jump
            if j>i: a+=Q(j-i,j)*(val(d,c,p,(j,i+1))-V)
        vs=(Q(1,2),) if p==2 else (Q(0),Q(1,3))
        if not j: us=(Q(0),)
        elif i: us=(Q(1,j),Q(1,j)+Q(k-1,2))
        else: us=(Q(0),Q(1,4)+Q(k-1,2))
        row=out[d,k]
        for u,v in product(us,vs):
            x=a+b*(u+v)
            assert x>=0
            row[u,v]=min(x,row.get((u,v),x));num+=1
    return out,num


def fit(points,k):
    """A nonnegative affine minorant on the entire common census box.

    Try the constant, then monotone ramps from all four box corners. Every
    candidate is limited by the exact hidden corners, not by sampled receives.
    """
    best=(min(points.values()),Q(0),Q(0))
    U,V=Q(k,2),Q(1,2)
    center=lambda x: x[0]+U*x[1]/2+V*x[2]/2
    for du,dv in product((-1,1),repeat=2):
        u0=U if du<0 else 0; v0=V if dv<0 else 0
        for wu,wv in ((1,0),(0,1),(1,1),(1,3),(3,1)):
            B=du*wu/U if U else Q(0); C0=Q(dv*wv,V)
            if not (B or C0): continue
            A=-B*u0-C0*v0
            candidates=[r/g for (u,v),r in points.items() if (g:=A+B*u+C0*v)>0]
            if not candidates: continue
            r=min(candidates)
            line=(r*A,r*B,r*C0)
            if center(line)>center(best): best=line
    A,B,C0=best
    assert all(0<=A+B*u+C0*v for u,v in product((0,U),(0,V)))
    assert all(A+B*u+C0*v<=r for (u,v),r in points.items())
    return best


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',action='store_true')
    args=parser.parse_args()
    m=load(); lo,hi=m.base.LO,m.base.HI
    for side in ('lower','upper'):
        with m.base.Snapshot(m.base.BASE/f'quintet_{side}.bin',side) as s:
            W=m.broad(s,12000,30,D+1,False)
            roots,n,_=m.conditional(s,W,D,C,2,3,False)
            assert n==len(m._capture)==72322 and tuple(v*m.G for v in roots)==ROOTS[side]
            mixed=m.birth_mix(roots,False)
            expected={'lower':'610545.766447640','upper':'857402.786257122'}
            assert m.decimal(lo*mixed,False)==expected[side]
            points,checks=corners(m,s,W)
            rows=[]; groups=defaultdict(int)
            for (d,k), ps in sorted(points.items()):
                ar=fit(ps,k)
                if any(ar):
                    rows.append([d,k,*map(str,ar)])
                    groups['positive']+=1
                    groups['through64']+=d<=64
                    groups['receive_sensitive']+=bool(ar[1] or ar[2])
            s.check_unchanged()
            ans=dict(schema='quintet-conditioned-f-reduced-v1', side=side,
                conditioned_commit='fd63517', conditioned_replay_sha256=FROZEN,
                source_commit='b04cc50c1e391e0420eb7b006ac781d484c68543',
                payload_sha256=m.base.GOLD[side]['verify']['sha256'], scale=s.scale,
                cap_f_ticks=m.base.GOLD[side]['cap_f'], LO=str(lo),HI=str(hi),
                degree=D, arity=C, rows=rows,
                birth_V_ticks={'D=2':str(roots[0]),'D>=3':str(roots[1])},
                beta_D2_bounds=[str(Q(1,2)/(hi+2)),str(Q(1,2)/(lo+Q(3,2)))],
                birth_floor_ticks=str(lo*mixed), hidden_corners=checks)
        m._capture.clear()
        dest=HERE/f'{side}.json'; raw=(dumps(ans,separators=(',',':'))+'\n').encode()
        if args.replay: assert dest.read_bytes()==raw
        else: dest.write_bytes(raw)
        print(side,'hidden corners',checks,'census rows',len(points),dict(groups),
              'source floor',lo*mixed,'hash',sha256(raw).hexdigest(),flush=True)


if __name__=='__main__': main()
