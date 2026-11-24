#!/usr/bin/env python3
"""Physical D20 row checks: exact labelled arrows, transfer cancellation and charge."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import json
from pathlib import Path

from extract import (COMMIT, SOURCE, Cert, Charge, HASH, NB, ROOT, Rows, SPEC, audit,
                     bands, bud, children, embeddings, gain, pick, pinned_sources,
                     recharge, replace)


def fee(e,f,y,ids,p,x):
    a,ports=e.affine(f,y)
    for i,parent,_,s in ports:
        j=ids[i]
        send=Q(1,sum(x[p[j]])) if parent else sum((Q(1,sum(x[v])) for v in range(j+1,len(x))
                                                   if p[v]==j and v not in ids),Q(0))
        a+=s*send
    return a


def tree_check(c,e,q):
    checked=Counter(); cases=[]
    frontier=[(-1,0,0)]
    for n in range(3,7):
        cases+=frontier
        if n<6: frontier=[p+(j,) for p in frontier for j in range(n)]
    recipes=[((-1,0,1,2),(1,0,0,1)),((-1,0,0,2,2,0),(1,1,0,1,1,1)),
             ((-1,0,1,1,1),(1,0,1,1,1)),((-1,0,1,1),(1,c.D-4,1,1)),
             ((-1,0,1,1),(1,c.D-3,1,1)),((-1,0,1,1),(1,c.D-2,1,1)),
             ((-1,0,1,2,3),(1,17,16,16,1)),((-1,0,1,2,3),(1,18,18,18,1)),
             ((-1,0,1,2,3),(1,1,100,50,1)),((-1,0,1,2,3),(1,30,1,1,50))]
    for arity in (c.K,c.K+1): recipes.append(((-1,0)+(1,)*arity,(1,1)+(1,)*arity))
    cases.extend(audit.build(p,l) for p,l in recipes)
    max_slack=Q(0); wrong_rate_detected=False
    for par in cases:
        ch,x,core=audit.physical_view(par); ids={v:i for i,v in enumerate(core)}
        p=(-1,*(ids[par[v]] for v in core[1:])); x=tuple(x)
        parts=c.potential(ch,x,parts=True); old=sum(parts)
        assert sum(k-h for l,k,h in x)==0
        assert c.c0*len(x)<=old/c.S<=c.C0*len(x)
        assert old<=sum(q.exact(*s) for s in x)
        assert old==sum(e.value(f,pick(x,v)) for f,v in embeddings(p))
        by_arrow=Q(0)
        arrows=audit.arrows(par); assert sum(arrows.values())==len(par)
        wrong=audit.arrows(par,wrong=True); wrong_arrow=Q(0)
        for j,r in arrows.items():
            nxt=c.potential(*audit.physical_view(par+(j,))[:2])
            by_arrow+=r*(nxt-old)
            wrong_arrow+=wrong[j]*(nxt-old)
        wrong_rate_detected |= by_arrow!=wrong_arrow
        signed=c.sign*(by_arrow-c.lam*old)/c.S
        rows=Q(0)
        for f,v in embeddings(p):
            y=pick(x,v)
            if not e.live(f,y): continue
            reward=fee(e,f,y,v,p,x)
            assert reward>=0,(f,y,reward)
            rows+=reward; checked[f]+=1
        assert signed==rows,(p,x,signed,rows)
        max_slack=max(signed,max_slack)

    # High stock without rebuilding tens of thousands of ordinary-vertex
    # trees. Compute directly from the physical core/leaf events.
    stock=[(1,c.LF,1,1),(1,c.LF+1,1,1),(1,c.LF+19,1,1)]
    for ls in stock:
        p=(-1,0,1,1); ch=children(p)
        x=tuple((l,len(ch[i]),int(i>0)) for i,l in enumerate(ls))
        old=c.potential(ch,x); rhs=Q(0)
        for i,s in enumerate(x):
            ns=ch[i]+([p[i]] if i else []); rate=s[0]+sum((Q(1,sum(x[v])) for v in ns),Q(0))
            rhs+=rate*(c.potential(ch,replace(x,i,gain(s)))-old)
            if s[0]:
                nn=[kids[:] for kids in ch]+[[]];nn[i].append(len(x))
                rhs+=Q(s[0],sum(s))*(c.potential(nn,replace(x,i,bud(s))+(NB,))-old)
        value=c.sign*(rhs-c.lam*old)/c.S
        assert value>=0 and value==sum(fee(e,f,pick(x,v),v,p,x) for f,v in embeddings(p))
        assert old<=sum(q.exact(*s) for s in x)
    assert set(checked)==set(SPEC) and wrong_rate_detected
    print(c.side,len(cases),'labelled arrow trees and',len(stock),'compressed stock trees;',
          'checked occurrences',dict(checked),'wrong-recipient rate rejected;',
          'max signed physical slack',float(max_slack),flush=True)


def row_check(c,e,q,path):
    data=json.loads(path.read_bytes()); assert data['side']==c.side and data['payload_sha256']==HASH[c.side]
    assert data['source_commit']==COMMIT and data['source_blobs']==SOURCE
    assert data['root_ticks']==list(ROOT[c.side])
    assert {(h,k):v for h,k,v in data['charge_by_h_k_ticks']}==q.by_k
    for l,k,h,a,b,z,w in data['f_rows']:
        A,B=e.first_body((l,k,h)); assert (A,B,A,A+B*Q(k+h,2))==tuple(map(Q,(a,b,z,w)))
    assert json.loads(json.dumps(bands(c,len(data['bands'])-1),default=str))==data['bands']
    assert json.loads(json.dumps(recharge(c,q,data['recharge']['degree']),default=str))==data['recharge']
    # Mixed tail and degree seam, plus ordinary inner edge and every higher
    # family. Check true affine corners, including their signs.
    n=0
    for f,p in SPEC.items():
        ch=children(p)
        for ds in ((3,)*len(p), (19,)*len(p), (20,)*len(p), (21,)*len(p), (64,)*len(p)):
            for slot in range(len(p)):
                x=tuple((max(d-len(ch[i])-1,0)+(int(i==slot)),len(ch[i]),1) for i,d in enumerate(ds))
                if not e.live(f,x):continue
                a,ps=e.affine(f,x); assert a+sum(min(0,b*k/2) for _,_,k,b in ps)>=0
                n+=1
    print(c.side,len(data['f_rows']),'replayed signed f rows;',n,'synthetic motif corner rows;',
          'all 404 recharge caps, bands and physical recharge equal to authenticated export',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('path',type=Path);p.add_argument('side',choices=HASH)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--tables',type=Path,required=True)
    a=p.parse_args()
    if not __debug__:p.error('assertions required')
    pinned_sources(a.source); c=Cert(a.path,a.side); e=Rows(c); q=Charge(c)
    row_check(c,e,q,a.tables); tree_check(c,e,q);c.check_unchanged()


if __name__=='__main__':main()
