#!/usr/bin/env python3
"""Short exact checks of free attachment and high-core growth fixtures."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import importlib.util


def graph(parents):
    out = [set() for _ in range(len(parents)+1)]
    for x, parent in enumerate(parents, 1):
        out[x].add(parent)
        out[parent].add(x)
    return out


def cap_check():
    layer = {(0,0)}
    total = 0
    for n in range(3,9):
        nxt = set()
        for parents in layer:
            g = graph(parents)
            deg = list(map(len,g))
            c = Counter(deg)
            assert max(deg) <= 3 and c[3] == c[1]-2
            assert 2*(n-c[1]) >= n-2
            d = 1 if 4*c[1] >= n else 2
            # Reconstruct target probabilities from every individual source.
            prob = sum((F(1,n*deg[x]) for x,h in enumerate(g) for y in h if deg[y]==d),F())
            assert prob >= F(1,12)
            nxt.update(parents+(x,) for x,v in enumerate(deg) if v <= 2)
            total += 1
        layer = nxt
    # A delayed star can cross arbitrarily large core counts by elongating
    # a single arm. Explicitly check its count and each channel's probability.
    n=12
    p=(0,)*(n-1)
    end=1
    for step in range(18):
        g=graph(p)
        assert sum(len(x)>1 for x in g)==1+step
        assert len(g[end])==1
        src,=g[end]
        assert F(1,len(g)*len(g[src]))>0
        p+=(end,)
        end=len(p)
    print(f'PASS {total} degree-capped labelled histories and 18 delayed-star arm steps')


def leaf_check():
    # Reuse the audited arrow runner, testing the special edge seed whose
    # attachment site is initially a genuine degree-one leaf.
    path = Path(__file__).resolve().parents[2]/'influence'/'check'/'coupling.py'
    spec=importlib.util.spec_from_file_location('parity_audit',path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    B=frozenset({mod.edge(0,1)})
    A=B|{mod.edge(0,2)}
    assert mod.invariant(A,B,None)=={0}
    assert (len(mod.neighbors(B)[0]),len(mod.neighbors(A)[0]))==(1,2)
    layer={(A,B)}
    count=0
    for _ in range(4):
        for A,B in layer:
            mod.audit(A,B,None)
            count+=1
        layer={mod.step(A,B,(sg,x,y)) for A,B in layer for sg,x,y,_ in mod.rates(A,B,None)}
    print(f'PASS {count} exact paired free states from initial genuine-leaf promotion')


if __name__=='__main__':
    leaf_check()
    cap_check()
