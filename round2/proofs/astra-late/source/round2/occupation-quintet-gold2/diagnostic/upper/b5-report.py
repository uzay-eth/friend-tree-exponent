from collections import defaultdict
from pathlib import Path

H=Path(__file__).resolve().parent
base={}; trees={}
for line in (H.parent/'upper.tsv').read_text().splitlines():
    a=line.split('\t')
    if a[0]=='tree': path=a[1]
    if a[:2]==['family','b5']:base[path]=(int(a[2]),float(a[6]))
with (H/'b5-sample.tsv').open() as inp:
    for line in inp:
        a=line.rstrip().split('\t')
        if a[0]=='tree':
            _,path,c,n,v=a
            assert base[path][0]==int(n) and abs(base[path][1]-float(v))<2e-12,(path,base[path],n,v)
            trees[path]=(int(c),int(n),float(v),[])
        else:trees[path][3].append((tuple(map(int,a[1:-2])),int(a[-2]),float(a[-1])))
assert len(trees)==9
orient=('hub','short tip','interior arm','far tip')
hi=('short child','short parent','far tip','tied')
for held in (False,True):
    tt=[v for path,v in trees.items() if path.endswith('numeric4-root-tree.tsv')==held]
    cores=sum(v[0] for v in tt)
    assert cores==(235709 if held else 896887)
    rows=[r for v in tt for r in v[3]]
    assert sum(n for _,n,_ in rows)==sum(v[1] for v in tt)
    print('held-out' if held else 'pooled',cores,'all B5',sum(v[1] for v in tt),sum(v[2] for v in tt)*1e6/cores)
    def show(name,fn,key,top=None):
        d=defaultdict(lambda:[0,0.])
        for k,n,v in rows:
            if fn(k):
                p=d[key(k)];p[0]+=n;p[1]+=v*1e6/cores
        print(name,'total',round(sum(v[1] for v in d.values()),8))
        for k,(n,v) in sorted(d.items(),key=lambda r:-r[1][1])[:top]:print(' ',k,n,round(v,8))
    clean=lambda k:not k[0]
    mixed=lambda k:clean(k) and k[11] and k[10]>=1
    show('root role (all physical)',lambda k:True,lambda k:(k[0],orient[k[1]]))
    show('max tip role and degree, max>=19',lambda k:clean(k) and k[11],lambda k:(orient[k[1]],hi[k[2]],k[7]),24)
    show('old+young (max tip>=19, min tip<=5), root and highest tip',mixed,lambda k:(orient[k[1]],hi[k[2]]))
    show('old+young internal labels (root, hub degree/arity, arm degree/arity)',mixed,
         lambda k:(orient[k[1]],k[3:5],k[5:7]),25)
    show('old+young internal ordinary leaf counts (root, hub, arm)',mixed,
         lambda k:(orient[k[1]],k[3]-k[4],k[5]-k[6]),18)
    show('old+two young and exact low internal cut, highest tip role and band',
         lambda k:mixed(k) and k[10]>=2 and k[3]<=5 and k[5]<=5,lambda k:(orient[k[1]],hi[k[2]],k[7]),24)
