from collections import defaultdict
from pathlib import Path

H=Path(__file__).resolve().parent
base={}; trees={}
for line in (H.parent/'upper.tsv').read_text().splitlines():
    a=line.split('\t')
    if a[0]=='tree': path=a[1]
    if a[:2]==['family','s4']:base[path]=(int(a[2]),float(a[6]))
for line in (H/'s4-sample.tsv').read_text().splitlines():
    a=line.split('\t')
    if a[0]=='tree':
        _,path,c,n,v=a
        assert base[path][0]==int(n) and abs(base[path][1]-float(v))<2e-12
        trees[path]=(int(c),int(n),float(v),[])
    else:trees[path][3].append((tuple(map(int,a[1:-2])),int(a[-2]),float(a[-1])))
assert len(trees)==9
for held in (False,True):
    tt=[v for path,v in trees.items() if path.endswith('numeric4-root-tree.tsv')==held]
    cores=sum(v[0] for v in tt)
    assert cores==(235709 if held else 896887)
    rows=[r for v in tt for r in v[3]]
    assert sum(n for _,n,_ in rows)==sum(v[1] for v in tt)
    print('held-out' if held else 'pooled',cores,'all S4',sum(v[1] for v in tt),sum(v[2] for v in tt)*1e6/cores)
    def show(name,fn,key,top=None):
        d=defaultdict(lambda:[0,0.])
        for k,n,v in rows:
            if fn(k):
                p=d[key(k)];p[0]+=n;p[1]+=v*1e6/cores
        print(name,'total',round(sum(v[1] for v in d.values()),8),
              [(k,n,round(v,8)) for k,(n,v) in sorted(d.items(),key=lambda r:-r[1][1])[:top]])
    clean=lambda k:not k[0] and k[1]
    oi=('mother max','child max','tied max')
    show('all by founder/mother included/max role',lambda k:True,lambda k:(k[0],k[1],oi[k[2]]))
    show('mother + 3 children, max arm >=19 by role/band',lambda k:clean(k) and k[12],lambda k:(oi[k[2]],k[6]))
    show('mother + 3 children, old+young (max>=19, min child<=5), max and 2nd-max bands',
         lambda k:clean(k) and k[12] and k[8]==5,lambda k:(oi[k[2]],k[6],k[7]),24)
    show('old+young exact hub D/c/ordinary leaves',lambda k:clean(k) and k[12] and k[8]==5,
         lambda k:(oi[k[2]],k[4],k[5],k[4]-k[5]),20)
    show('old+young separately by hub extra child count',lambda k:clean(k) and k[12] and k[8]==5,
         lambda k:(oi[k[2]],k[5]-4))
    show('old child rank (among included children, increasing birth order)',lambda k:clean(k) and k[12] and k[2]==1,
         lambda k:(k[3],k[6],k[11]),15)
    show('mother-max, 1 vs 2 vs 3 young included children',lambda k:clean(k) and k[12] and k[2]==0,
         lambda k:(k[6],k[10],k[9]),20)
