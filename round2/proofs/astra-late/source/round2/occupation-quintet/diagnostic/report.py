#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path
import hashlib

HERE=Path(__file__).resolve().parent
F=('f','g','h3','p4','s3','p5','b5','s4')
DEG=('2–5','6–11','12–18','19–37','38–179','180–1998','1999','2000–2011','2012–9999','10000+')
AR=('0','1','2–4','5–8','9–15','16–26','27+')
DIRECTION=('mother pays','mother donates','older child pays','younger child pays')


def report(side):
    rows=[]
    for line in (HERE/f'{side}.tsv').read_text().splitlines():
        tag,*fields=line.split('\t')
        if tag=='tree':
            path,mode,c,n,p,charge,raw,direct=fields
            assert int(mode)==(2 if side=='lower' else 1)
            rows.append(dict(path=Path(path),c=int(c),n=int(n),p=float(p),charge=float(charge),
                             err=max(abs(float(charge)-float(raw)),abs(float(charge)-float(direct))),
                             fam={},bins={}))
        elif tag=='family':
            rows[-1]['fam'][fields[0]]=list(map(float,fields[1:]))
        elif tag=='bin':
            rows[-1]['bins'][tuple(map(int,fields[:8]))]=(int(fields[8]),float(fields[9]))
    assert len(rows)==9 and rows[-1]['path'].name=='numeric4-root-tree.tsv'
    print(side,'inputs')
    for r in rows:
        with r['path'].open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
        print(r['path'].name,r['n'],r['c'],digest)
    for label,rs in (('pooled age20',rows[:8]),('held-out age20',rows[-1:])):
        c=sum(r['c'] for r in rs)
        p=sum(r['p'] for r in rs)
        cap=(19186 if side=='lower' else 19783)/2**40
        mmin=(99882299906 if side=='lower' else 155236771435)/2**40
        k=1e6/c
        print(label,c,'Phi/C',p/c,'max generator absolute delta',max(r['err'] for r in rs))
        total=sum(r['charge'] for r in rs)
        print('net physical ppm/core',-total*k,'nominal gap ppm/potential',-1e6*total/p,
              'cap-to-observed reference gap ppm',1e6*cap/mmin-1e6*total/p,
              'cap numerator ppm/core',cap*1e6)
        if len(rs)>1:
            print('individual net ppm/core',[round(-1e6*r['charge']/r['c'],6) for r in rs])
        print('family n/core potential/core raw_ppmc charged_ppmc neg_ppmc pos_ppmc knot_or_box_gap_ppmc edge18_cap_slack_ppmc')
        for f in F:
            sums=[sum(r['fam'][f][i] for r in rs) for i in range(8)]
            print(f,*(f'{v*(1/c if i<2 else k):.9g}' for i,v in enumerate(sums)))
        bins=defaultdict(float)
        for r in rs:
            for key,(_,v) in r['bins'].items():bins[key]+=v*k
        def summary(name, pred, project):
            a=defaultdict(float)
            for t,v in bins.items():
                if pred(t):a[project(t)]+=v
            print(name,'total',f'{sum(a.values()):.8f}', 'top',[(key,round(v,7)) for key,v in sorted(a.items(),key=lambda kv:-kv[1])[:15]])
        summary('FY by mask',lambda t:t[0]==1,lambda t:t[6])
        remaining=lambda t:t[0]==1 and (t[6]==0 or (side=='upper' and t[6]==3))
        summary('FY residual outside named masks by direction and degree',remaining,
                lambda t:(DIRECTION[t[1]],*(DEG[i] for i in t[2:5]),AR[t[5]],AR[t[7]]))
        summary('FY residual outside named masks by coarse payer/donor degree',remaining,
                lambda t:(DIRECTION[t[1]],t[3]>0,t[4]>0))
        summary('uncovered by family/direction/receiver degree/donor degree',lambda t:t[0]==2,
                lambda t:(F[t[1]],'mother donor' if not t[2] else 'child donor',DEG[t[3]],DEG[t[4]],AR[t[5]]))
        summary('uncovered by family/direction',lambda t:t[0]==2,lambda t:(F[t[1]],'mother donor' if not t[2] else 'child donor'))
        summary('negative charged higher by family/max degree/max arity',lambda t:t[0]==3 and t[1]>=2,
                lambda t:(F[t[1]],DEG[t[2]],AR[t[3]]))
        summary('negative charged higher by family',lambda t:t[0]==3 and t[1]>=2,lambda t:F[t[1]])
        summary('H3 mixed support min-tip <=11 max-tip 19–1998 middle c2–8 outside old masks',
                lambda t:t[0]==5 and t[1]==2 and t[2]<=1 and 3<=t[3]<=5 and t[5] in (2,3) and not(t[6] or t[7]),
                lambda t:(DEG[t[3]],DEG[t[4]],AR[t[5]]))
        summary('H3 broad support max-tip 12–1998 all middle through c8 outside old masks',
                lambda t:t[0]==5 and t[1]==2 and 2<=t[3]<=5 and t[5] in (2,3) and not(t[6] or t[7]),
                lambda t:(DEG[t[3]],AR[t[5]]))
        summary('H3 priority support max-tip 19–1998 middle c2–4 outside old masks',
                lambda t:t[0]==5 and t[1]==2 and 3<=t[3]<=5 and t[5]==2 and not(t[6] or t[7]),
                lambda t:DEG[t[3]])
        summary('quintet slices by family and max tip degree',lambda t:t[0]==6,lambda t:(F[t[1]],DEG[t[2]]))
        summary('edge cap slack after FY and plain donor gaps, by unordered endpoint bands',lambda t:t[0]==7,
                lambda t:(DEG[t[2]],DEG[t[3]],AR[t[4]]))
        summary('singleton slack by degree / c',lambda t:t[0]==4,lambda t:(DEG[t[2]],AR[t[3]]))
        summary('founder-containing charged',lambda t:t[0]==0 and t[6]==1,lambda t:F[t[1]])


for side in ('lower','upper'):
    report(side)
