#!/usr/bin/env python3
"""Exact independent set/deletion accounting and literal ordinary-tree generator."""
from bisect import bisect_right
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations
from os import close, pread
from pathlib import Path
from struct import unpack
from subprocess import check_output
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'attacks/neon2-quintet/physical'))
from check import Payload, WORK, NAMES, physical_core, subsets

# One branching, physically attainable topology, four distinct degree fixtures.
# Each configuration hits different inner/tail, stock and singleton seams.
PAR=(-1,0,1,2,3,1,1,1,0,8,8,3,0)
DEGREES=(
    (2000,5,2,3,2,2,3,5,65,2,200,179,190),
    (1999,6,3,3,2,11,2,181,38,2,3,2,2001),
    (2001,5,3,5,180,2,2,48,37,2,183,6,180),
    (179,5,2,4,2,2,5,12,8,2,3,2,3),
)


def dump(path,ds,par=PAR):
    c=Counter(par[1:]); c.update(range(1,len(par)))
    path.write_text(''.join(f'{p} {ds[i]-c[i]}\n' for i,p in enumerate(par)))


def evaluate(p,path):
    rows=[tuple(map(int,line.split())) for line in path.read_text().splitlines()]
    tr=[set() for _ in rows]
    for j,(i,_) in enumerate(rows[1:],1):tr[i].add(j);tr[j].add(i)
    for i,(_,l) in enumerate(rows):
        for _ in range(l):tr[i].add(len(tr));tr.append({i})
    sig=1 if p.h['mode']==1 else -1
    lam=Q(p.h['P'],p.h['Q']); sc=p.h['S'];g=p.ends[2]
    knot=tuple(Q(k,p.h['Q_K']) for k in unpack('<24q',pread(p.fd,192,304))[:p.h['nknots']])
    assert knot[0]==0 and knot[-1]==Q(1,2)
    def r1(m):return Q(p.word(g.ids[m]),p.h['S_r'])
    def fy(m,x):
        i=min(bisect_right(knot,x)-1,len(knot)-2);assert i>=0
        j=g.exact_count+g.ids[m]*len(knot)+i
        t=(x-knot[i])/(knot[i+1]-knot[i]);return (1-t)*p.word(j)+t*p.word(j+1)
    def universe(t):
        a,m=physical_core(t)
        # Brute combinations on 13-14 cores; subsets() supplies wire lookup
        # and geometric role recovery but does not drive C++ enumeration.
        val=subsets(p,a,m)
        brute={}
        for sz in range(1,6):
            for key in combinations(a,sz):
                seen={key[0]}
                while True:
                    nxt=seen|{v for u in seen for v in a[u] if v in key}
                    if nxt==seen:break
                    seen=nxt
                if len(seen)==sz:
                    aa={u:a[u]&seen for u in seen}
                    brute.update((k,z) for k,z in subsets(p,aa,m).items() if len(k[1])==sz)
        assert brute==val
        return a,m,{s:(k,Q(w)) for (k,s),(w,_) in val.items()}
    a,m,old=universe(tr)
    def value(t,s):return t.get(s,('',0))[1]
    def birth(j):t=tr.copy()+[{j}];t[j]=t[j]|{len(tr)};return t
    w=defaultdict(Q); raw=defaultdict(Q); direct=defaultdict(Q); charged=defaultdict(Q)
    rates=defaultdict(Q)
    for u,nb in enumerate(tr):
        for v in nb:rates[v]+=Q(1,len(nb))
    assert sum(rates.values())==len(tr)
    gains={u:universe(birth(u))[2] for u in a}
    delta=lambda u,s:value(gains[u],s)-value(old,s)
    buds={u:next((v for v in tr[u] if v not in a),None) for u in a}
    promot={v:universe(birth(v))[2] for v in buds.values() if v is not None}
    for s,(k,z) in old.items():w[k]+=z;raw[s]=-sig*lam*z;direct[k]-=sig*lam*z
    for u in a:
        for s,(k,_) in old.items():
            if u in s:raw[s]+=sig*rates[u]*delta(u,s)
        for s,(k,_) in gains[u].items():assert s in old
        for s in old:direct[old[s][0]]+=sig*rates[u]*delta(u,s)
        j=buds[u]
        if j is None:continue
        ym=promot[j]; rr=rates[j]*(len(tr[u])-len(a[u]))
        for s,(k,_) in old.items():
            v=sig*rr*(value(ym,s)-value(old,s));raw[s]+=v;direct[k]+=v
        for s,(k,z) in ym.items():
            if j not in s:continue
            base=tuple(v for v in s if v!=j) or (u,)
            raw[base]+=sig*rr*z;direct[k]+=sig*rr*z
    charged.update(raw);moves=defaultdict(Q); fy_both=Counter(); fy_gap=Q(0)
    for s,(k,_) in old.items():
        for u in s:
            d,c=m[u]; ex=d<g.d and 1<=c<=g.k
            for v in a[u]-set(s):
                add=tuple(sorted((*s,v)));x=Q(1,m[v][0]);j=sig*delta(u,s)
                if add in old and ex:
                    inc=(x-r1(m[u]))*j;tgt=old[add][0]
                    charged[s]-=inc;charged[add]+=inc;moves[k,tgt]+=inc
                elif k=='g' and ex and c>=2 and m[u] not in p.ins['3m']:
                    inc=sig*fy(m[u],x);alt=tuple(sorted((u,v)))
                    charged[s]-=inc;charged[alt]+=inc
                    gap=max(t*j-sig*fy(m[u],t) for t in knot)-(x*j-inc)
                    assert gap>=0;fy_gap+=gap;fy_both['slots']+=1;fy_both['positive']+=gap>0
    charged_f=defaultdict(Q)
    own_raw=defaultdict(Q); counts=Counter()
    for s,(k,_) in old.items():charged_f[k]+=charged[s];own_raw[k]+=raw[s];counts[k]+=1
    assert sum(charged_f.values())==sum(direct.values())==sum(own_raw.values())
    if path.stem=='fixture3':
        s=(1,9,10)
        assert old[s][0]=='h3' and charged[s]<0
        print('target H3',s,'marks',tuple(m[u] for u in s),'charged ticks exact',charged[s])
    return counts,w,own_raw,charged_f,moves,fy_gap,fy_both,len(tr)


def main():
    for side in ('lower','upper'):
        path=WORK/f'quintet_{side}.bin'; p=Payload(path)
        all_mov=set(); all_fy=Counter()
        for i,ds in enumerate(DEGREES):
            f=HERE/f'fixture{i}.tsv'
            if i==3: dump(f,(2,*ds),(-1,*(j+1 for j in PAR)))
            else: dump(f,ds)
            counts,w,raw,ch,moves,fy,checks,n=evaluate(p,f)
            out=check_output([HERE/'census',path,f],text=True).splitlines()
            seen=set()
            for line in out:
                tag,*a=line.split('\t')
                if tag=='family':
                    k=a[0];assert int(a[1])==counts[k];seen.add(k)
                    for z,t in zip(a[2:5], (w[k],raw[k],ch[k])):
                        assert abs(Q(z)-t/p.h['S'])<Q('3e-11'),(side,i,k,z,t)
                if tag=='move':
                    k,l=map(lambda j:NAMES[int(j)],a[:2]); x=moves[k,l]/p.h['S']
                    assert abs(Q(a[2])-x)<Q('3e-11')
                    if x:all_mov.add((k,l))
            assert seen==set(NAMES)
            cppfy=sum(Q(line.split('\t')[-1]) for line in out if line.startswith('bin\t1\t'))
            assert abs(cppfy-fy/p.h['S'])<Q('3e-11')
            all_fy.update(checks)
            print(side,i,'N',n,'counts',dict(counts),'FY',dict(checks),'total exact',sum(ch.values())/p.h['S'])
        assert len(all_mov)==8 and all_fy['positive']>0,(all_mov,all_fy)
        print(side,'all signed transfer channels',sorted(all_mov),'FY checks',dict(all_fy))
        close(p.fd)


if __name__=='__main__':main()
