"""Independently sum exact physical quartet J and drift on reachable and decorated trees."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from pathlib import Path
from subprocess import run
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "occupation"))
from check_five import Certificate, geometry, quartet, subsets, transitions

S = 2**40
U = F(1133, 2000)


@cache
def keys(state):
    parents, _ = state
    _, degree, core = geometry(state)
    rows = {}
    for n in range(1,6):
        for members, shape in subsets(parents,n):
            sh = 0
            if n > 1:
                for i, ch in enumerate(shape, 1):
                    sh = i*sh + int(ch)
            rows[members] = (n, sh, *(z for j in members for z in (degree[j],core[j])))
    return rows


def fixtures():
    layer = {((-1,), (2,))}
    base = set()
    for _ in range(3, 11):
        base.update(layer)
        layer = {nx for st in layer for nx, rate, b in transitions(st)}
    for par in [(-1,0,1,2,3,4),(-1,0,1,1,2,3),(-1,0,1,2,2,3),(-1,0,1,1,1,2),(-1,0,1,1,2,2,3,3)]:
        c = [par.count(i)+(i!=0) for i in range(len(par))]
        for i in range(1,len(par)):
            for d in (5,8,10,13,14,15,33,191,1301):
                leaves = [max(int(k==1),4-k) for k in c]
                leaves[i] = max(0,d-c[i])
                base.add((par,tuple(leaves)))
                for j in range(i+1,len(par)):
                    pair = list(leaves)
                    pair[j] = max(0,19-c[j])
                    base.add((par,tuple(pair)))
    return base


def check(side, base, queries, request, build, mode):
    snapshot = HERE.parent / "outputs" / f"mc-{side}.bin"
    output = build / f"scores-{side}.tsv"
    with (build/f"fixture-{side}.log").open("w") as log:
        run([str(build/"reserve"),str(snapshot),mode,"3",str(request),str(output)],check=True,stdout=log)
    values = [int(s) for s in output.read_text().splitlines()]
    born = values.pop(0)
    table = dict(zip(queries, values, strict=True))
    assert born == table[1,0,2,1]
    cert = Certificate(snapshot)
    counts = Counter()
    @cache
    def score(state):
        return sum(table[key] for key in keys(state).values())
    @cache
    def gradient(family, state, i):
        after = list(state)
        l,k = after[i]
        after[i] = (l+1,k)
        return cert.sign*(cert.weight(family,tuple(after))-cert.weight(family,state))
    for st in base:
        par, leaves = st
        near, degree, core = geometry(st)
        actual = F(0)
        for members, _ in subsets(par,4):
            family, ordered = quartet(par,members)
            state = tuple((leaves[j],core[j]-1) for j in ordered)
            for role,x in enumerate(ordered):
                g = gradient(family,state,role)
                for z in near[x]:
                    if z==0 or z in ordered:
                        continue
                    j = max(0,g)/F(2)-g/degree[z]
                    assert j>=0
                    actual+=j
                    if j:
                        counts[family]+=1
                        counts["outside-current-table"] += max(degree[k] for k in ordered)>=14
                    if max(0,-g)/F(2)+g/degree[z] > j:
                        counts["wrong-sign overcharge"]+=1
        t = score(st)
        gt = births = F(0)
        for nx,r,bud in transitions(st):
            gt+=r*(score(nx)-t)
            if bud:births+=r
            counts["elementary successors"]+=1
        gap = actual-(U*t-gt+born*births)
        assert gap>=0,(side,st,gap)
        counts["inflated reserve rejected"]+=gap<S*births/1000
    assert all(counts[k]>0 for k in ("P","F","Y","T","outside-current-table","wrong-sign overcharge","inflated reserve rejected"))
    print(side,mode,len(base),"trees",len(queries),"query configurations",dict(sorted(counts.items())),flush=True)


if __name__ == "__main__":
    build=Path(sys.argv[1])
    mode=sys.argv[2] if len(sys.argv)>2 else "small"
    assert mode in ("small","small-split")
    base=fixtures()
    queries=sorted({key for st in base for nx in (st,*[n for n,r,b in transitions(st)]) for key in keys(nx).values()})
    request=build/"queries.tsv"
    request.write_text("".join(" ".join(map(str,k))+"\n" for k in queries))
    for side in ("lower","upper"):
        check(side,base,queries,request,build,mode)
