from collections import defaultdict
from pathlib import Path

HERE=Path(__file__).resolve().parent
D=(3,5,8,11,18,24,37,65,100,179,360,1988,1998,1999,2000,2011,9999,10**30)
OR=('mother high','child high, mother low','two children')
PORT=('middle','low tip','high tip','uncovered mother','uncovered child',
      'outgoing mother transfer','outgoing child transfer','incoming low-tip transfer',
      'incoming high-tip transfer','covered mother: actual','covered child: actual')


def read(side,held):
    check={}
    for line in (HERE/f'{side}.tsv').read_text().splitlines():
        a=line.split('\t')
        if a[0]=='tree':path=a[1]
        if a[:2]==['family','h3']:check[path]=tuple(float(a[i]) for i in (2,5,6,8))
    maps=[defaultdict(lambda:[0.]*4) for _ in range(2)]
    rates=defaultdict(lambda:[0.]*4)
    cores=seen=0
    with (HERE/f'h3-{side}.tsv').open() as stream:
        for line in stream:
            a=line.rstrip().split('\t')
            if a[0]=='tree':
                old=check[a[1]]; new=tuple(map(float,a[4:]));assert new[0]==old[0]
                assert max(abs(x-y) for x,y in zip(old,new))<2e-11,(a,old,new)
                keep=(a[1].endswith('/numeric4-root-tree.tsv')==held)
                if keep:cores+=int(a[3]);seen+=1
                continue
            if not keep:continue
            if a[0]=='rate':
                dest=rates[tuple(map(int,a[1:5]))]
                for i,x in enumerate(a[5:]):dest[i]+=float(x)
                continue
            assert a[0] in ('row','port')
            dest=maps[a[0]=='port'][tuple(map(int,a[1:9]))]
            for i,x in enumerate(a[9:]):dest[i]+=float(x)
    assert seen==(1 if held else 8)
    for d in maps:
        for v in d.values():
            for i in (1,2,3):v[i]*=1e6/cores
    return cores,*maps,rates


def run(side):
    for label in ('pooled','heldout'):
        n,row,port,rates=read(side,label=='heldout')
        def group(src,fn,key):
            d=defaultdict(lambda:[0.]*4)
            for k,v in src.items():
                if fn(k):
                    for i,a in enumerate(v):d[key(k)][i]+=a
            return {k:(int(v[0]),*(round(x,8) for x in v[1:])) for k,v in sorted(d.items())}
        def say(name,fn,key=lambda k:'all',src=row):print(name,group(src,fn,key))
        clean=lambda k:k[6]==0 and not k[7]
        priority=lambda k:clean(k) and 18<D[k[1]]<=1998 and k[5]<=4
        print(side,label,n,'(count, signed charged or port, negative row, uncovered gap), per-million/core')
        say('all orientation (strict degree ordering; ties separate)',lambda k:True,lambda k:'tie' if k[7] else OR[k[0]])
        say('founder / strict old chain / high-middle mask flags (bits 1 / 2 / 4)',lambda k:True,lambda k:k[6])
        say('priority by orientation',priority,lambda k:OR[k[0]])
        say('priority high-tip max band by orientation',priority,lambda k:(OR[k[0]],D[k[1]]))
        say('priority low tip D and actual leaf stock L (99 means >=99)',priority,
            lambda k:(OR[k[0]],D[k[2]],'0' if k[3]%100==0 else '1–2' if k[3]%100<=2 else '3–5' if k[3]%100<=5 else '6+'))
        say('priority low tip actual L and c (low D<=11)',lambda k:priority(k) and D[k[2]]<=11,
            lambda k:(OR[k[0]],k[3]%100,k[3]//100%100))
        say('mixed seam, low tip <= 11, no middle cut',lambda k:clean(k) and D[k[2]]<=11 and D[k[1]]>=1998,
            lambda k:(D[k[1]],OR[k[0]]))
        say('mixed seam, low tip <= 11, middle c <= 4',lambda k:clean(k) and D[k[2]]<=11 and D[k[1]]>=1998 and k[5]<=4,
            lambda k:(D[k[1]],OR[k[0]]))
        for hi in (37,65,100,179,360):
            for lo in (5,11,18,24,37,65):
                if lo>=hi:continue
                for md in (5,8,11,18,37):
                    fn=lambda k:priority(k) and D[k[1]]<=hi and D[k[2]]<=lo and k[4]<=md
                    val=group(row,fn,lambda k:'all').get('all',(0,0,0,0))
                    if val[2]>=1.80 or (hi,lo,md) in ((65,11,8),(100,18,8),(179,18,8)):
                        print('cap',hi,lo,md,val,group(row,fn,lambda k:OR[k[0]]))
        small=lambda k:priority(k) and k[0]==0 and k[4]<=5 and D[k[1]]<=65 and D[k[2]]<=18
        def top(name,fn,key):
            grouped=group(row,fn,key)
            vals=sorted(grouped.items(),key=lambda kv:-kv[1][2]);print(name,'top',vals[:20])
        top('small mother-high role arities (mother,middle,child)',small,lambda k:(k[3]//10000,k[5],k[3]//100%100))
        top('small mother-high actual child stock (degree cutoff, L, c)',small,
            lambda k:(D[k[2]],k[3]%100,k[3]//100%100))
        say('small mother-high child degree split',small,lambda k:D[k[2]])
        say('small mother-high by middle degree / arity',small,lambda k:(k[4],k[5]))
        say('small mother-high high degree split',small,lambda k:D[k[1]])
        say('small mother-high exact virgin child',(lambda k:small(k) and k[3]%10000==101),lambda k:D[k[1]])
        def portkey(k):return OR[k[0]],PORT[k[3]],PORT[k[4]]
        say('all actual/corrected ports',lambda k:k[6]<8,portkey,port)
        say('small 65/18/5 mother-high actual/corrected ports',
            lambda k:k[6]==8 and D[k[1]]<=65 and D[k[2]]<=18,portkey,port)
        for c in sorted({k[0] for k in rates}):
            vs=[v for k,v in rates.items() if k[0]==c]
            N,p,kids,rew=map(sum,zip(*vs))
            uplow=sum(v[0] for k,v in rates.items() if k[0]==c and D[k[1]]>=100)
            childhigh=sum(v[0] for k,v in rates.items() if k[0]==c and k[2]>=max(1,c-2))
            joint=sum(v[0] for k,v in rates.items() if k[0]==c and D[k[1]]>=100 and k[2]>=max(1,c-2))
            print('high mother receive',c,int(N),'mean parent/other children',round(p/N,8),round(kids/N,8),
                  'freq parent D>65 / child sum>=(c-2)/4 / joint',round(uplow/N,6),round(childhigh/N,6),round(joint/N,6),
                  'neg reward ppm/core',round(rew*1e6/n,8))
        say('mixed (low <=11, high 19–1998) actual/corrected ports',
            lambda k:clean(k) and D[k[2]]<=11 and 18<D[k[1]]<=1998,portkey,port)
        say('mixed seam actual/corrected ports by high band',
            lambda k:clean(k) and D[k[2]]<=11 and D[k[1]] in (1999,2000,2011),
            lambda k:(D[k[1]],*portkey(k)),port)


for s in ('lower','upper'):run(s)
