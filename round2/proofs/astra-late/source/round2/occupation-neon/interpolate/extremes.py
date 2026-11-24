#!/usr/bin/env python3
"""Exact physical counterexamples for frozen-tail interpolation at target s."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ratio-transfer'))
from replay import Pair, Cert, J, tree, train

TARGET=Q('.56572365')


def show(x):return f'{float(x):.13g}'


def metric(pair,z,s):
    (dl,du),(vl,vu)=pair.full(z,True)
    return (vl,vu,dl-(s-pair.c[0].lam)*vl,(pair.c[1].lam-s)*vu-du)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lower',type=Path);p.add_argument('upper',type=Path)
    p.add_argument('--s',type=Q,default=TARGET)
    a=p.parse_args()
    pair=Pair(Cert(a.lower),Cert(a.upper))
    fixtures={'seed':tree(2),'linear-train':train(6,0),'hub16-train':train(6,16),
              'tail-star':tree(1000000),'tail-chain':tree(1000000,tree(1000000,tree(1000000))),
              'tail-cherry':tree(1000000,*([tree(1000000)]*2)),
              'tail-burst':tree(1000000,*([tree(1000000)]*18)),
              'all8':tree(1,tree(1,tree(1,J,J),tree(1,J),J),tree(1,J),J)}
    for name,z in fixtures.items():
        vl,vu,rl,ru=metric(pair,z,a.s)
        print(name,'Phi-',show(vl),'Phi+',show(vu),'R-/Phi-',show(rl/vl),'R+/Phi+',show(ru/vu),
              'R-+R+',show(rl+ru),'scale K breakpoint',show(-rl/ru) if ru else 'none',flush=True)

    # These probe the exact degree tail without extrapolating a decimal fit.
    # Every existing center/source annotation is irrelevant to this check.
    for l in (14001,1000000,10**12):
        for name,z in (('tail-star',tree(l)),('tail-chain',tree(l,tree(l,tree(l)))),
                       ('tail-burst',tree(l,*([tree(l)]*18)))):
            v,w,r,t=metric(pair,z,a.s)
            print(name,l,'signs',int(r>0)-int(r<0),int(t>0)-int(t<0),
                  'witness shift',show(r/v),show(t/w),flush=True)
    # Corners are computed by two exact all8 implementations above; no
    # interpolation in D or assumption about the tail limit is required.
    if a.s==TARGET:
        signs=[]
        for l in (14001,1000000):
            chain=tree(l,tree(l,tree(l)))
            v,w,r,t=metric(pair,chain,a.s);signs.append((r/v,t/w))
        assert max(signs[0])<Q('-4.5e-6') and min(signs[1])>Q('5.1e-6')
        _,_,r,t=metric(pair,fixtures['all8'],a.s)
        assert r<0 and t<0
        print('exact obstruction: every positive mix normalized by its own Phi has',
              'residual/Phi < -4.5e-6 and > +5.1e-6 on the two frozen chains',flush=True)
        def bare_chain(stock,m):
            z=tree(stock)
            for _ in range(m-1):z=tree(stock,z)
            return z
        # Finite interaction range makes each homogeneous train exactly
        # affine after its two boundary neighborhoods cease to overlap.
        for stock,compare in ((14001,lambda x: x<Q('-4.5e-6')),
                              (1000000,lambda x:x>Q('5.1e-6'))):
            checks=[metric(pair,bare_chain(stock,m),a.s) for m in (8,9,10)]
            delta=tuple(b-a for a,b in zip(checks[0],checks[1]))
            assert delta==tuple(b-a for a,b in zip(checks[1],checks[2]))
            v,w,r,t=delta
            assert min(v,w)>0 and compare(r/v) and compare(t/w)
            print('homogeneous bulk',stock,'exact per-core R/Phi',r/v,t/w,flush=True)


if __name__=='__main__':main()
