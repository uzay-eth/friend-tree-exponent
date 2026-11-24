"""Test physical seam and pooled gradients and compare the full census on an exact tree."""
from fractions import Fraction as F
from pathlib import Path
from random import Random
from subprocess import run
import sys
from check import HERE, Certificate, geometry, quartet, subsets

FAMILIES=("P","F","Y","T")
KIDS=((1,1,1,0),(2,0,1,0),(1,2,0,0),(3,0,0,0))
ENDS=((0,3),(1,3),(0,2,3),(1,2,3))


def observations():
    rng=Random(53275101)
    for j,shape in enumerate(FAMILIES):
        for index in range(160):
            d=[];c=[]
            for role,kids in enumerate(KIDS[j]):
                k=kids+rng.choice([0,0,0,1,2])
                ds=(2,3,4,6,8,13,14,15,21,41,105,1031) if role in ENDS[j] else (3,4,5,6,7)
                deg=max(k+1,rng.choice(ds));d.append(deg);c.append(k+1)
            actual=tuple((a-b,b-1) for a,b in zip(d,c))
            for site in range(4):
                exact=[(a,b,0,0) for a,b in zip(d,c)]
                yield j,shape,site,tuple(exact),actual,None
                partial=exact.copy()
                changed=[]
                for endpoint in ENDS[j]:
                    if site==endpoint:continue
                    coarse=exact.copy()
                    cut=min(d[endpoint],rng.choice((5,7,11,14,19)))
                    cmin=max(KIDS[j][endpoint]+1,min(c[endpoint],rng.choice((2,3,4,5))))
                    clamped=(index%2==0 or cmin<c[endpoint])
                    coarse[endpoint]=(cut,cmin,1,int(clamped))
                    yield j,shape,site,tuple(coarse),actual,endpoint
                    partial[endpoint]=coarse[endpoint];changed.append(endpoint)
                if len(changed)>=2:
                    yield j,shape,site,tuple(partial),actual,tuple(changed)


def test_slopes(build):
    cases=list(observations())
    request=build/"slope-requests.txt"
    request.write_text("".join(" ".join(map(str,(j+3,site,*(t for a in cells for t in a))))+"\n"
                               for j,f,site,cells,state,ep in cases))
    for side in ("upper","lower"):
        snap=HERE.parent/"outputs"/f"mc-{side}.bin"
        cert=Certificate(snap);output=build/f"slope-{side}.txt"
        run([str(build/"gradient_probe"),str(snap),str(request),str(output),"multiple"],check=True)
        result=[int(v) for v in output.read_text().splitlines()]
        exact=coarse=positive=bad_inf=multiple=multiple_positive=0
        for (_,fam,site,cell,state,ep),v in zip(cases,result,strict=True):
            def gradient(s):
                after=list(s);l,k=after[site];after[site]=(l+1,k)
                return cert.sign*(cert.weight(fam,tuple(after))-cert.weight(fam,tuple(s)))
            a=gradient(state)
            if ep is None:
                assert v==int(a),(fam,site,state,v,a)
                exact+=1
            else:
                assert (not v) or (v*a>0 and abs(v)<=abs(a)),(fam,site,cell,state,v,a)
                coarse+=1;positive+=v!=0
                endpoints=ep if isinstance(ep,tuple) else (ep,)
                if len(endpoints)>1:multiple+=1;multiple_positive+=v!=0
                inf=list(state)
                for p in endpoints:inf[p]=None
                wrong=gradient(inf)
                j=lambda x:max(x,0)/F(2)-x/3
                bad_inf+=j(wrong)>j(a)
        assert bad_inf and positive and multiple_positive
        print(side,exact,"exact quartet slopes,",coarse,"compatible pooled comparisons,",positive,
              "nonzero pooled bounds;",multiple,"multiple endpoint cases,",multiple_positive,
              "nonzero;",bad_inf,"infinity-only donor overcharges",flush=True)


def test_census(build):
    par=(-1,0,1,1,1,2,2,3,4,5,6,6,8,9)
    core=[par.count(i)+(i!=0) for i in range(len(par))]
    leaves=tuple(max(1,deg-c) for deg,c in zip((413,5,6,15,13,4,7,2,4,5,32,7,2,81),core))
    file=build/"exact-census-tree.tsv"
    file.write_text("".join(f"{p} {l}\n" for p,l in zip(par,leaves)))
    out=HERE.parent/"outputs"
    res=run([str(build/"census"),str(out/"mc-lower.bin"),str(out/"mc-upper.bin"),"--tree",str(file)],
            check=True,capture_output=True,text=True).stdout.splitlines()
    fields=res[0].split();records=[dict(zip(fields,line.split(),strict=True)) for line in res[1:]]
    near,ds,cs=geometry((par,leaves))
    for side in ("upper","lower"):
        cert=Certificate(out/f"mc-{side}.bin");values={f:[F(),F()] for f in FAMILIES}
        for group,_ in subsets(par,4):
            f,order=quartet(par,group);x=tuple((leaves[i],cs[i]-1) for i in order)
            for j,i in enumerate(order):
                nx=list(x);l,k=nx[j];nx[j]=(l+1,k)
                a=cert.sign*(cert.weight(f,tuple(nx))-cert.weight(f,x))
                for z in near[i]:
                    if z not in order:values[f][int(z==0)]+=max(0,a)/F(2)-a/ds[z]
        for r in records:
            if r["side"]!=side:continue
            for label,w in zip(("J","founder"),values[r["family"]]):
                assert abs(F(r[label])-w)<=max(1,w)/10**14,(side,r,w)
        print(side,"four physical family totals and four founder-donor totals agree with Fraction",flush=True)


if __name__=="__main__":
    build=Path(sys.argv[1])
    test_slopes(build);test_census(build)
