"""Gold2 host validation and exact stationary objective at mu=99/175."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from json import dumps
from pathlib import Path
from reader import R2, TARGET, pinned, open_host

CONDITIONED_REPORT='cc6c8b3c3c9c5d52917e41b17c18fad2018bb57e'
PAYER_REPORT_SHA='0020b89ff3d30dcf654225d7f37626d62c8f21b60312b7b2b342036faac6904f'
PAYER_TICKS=Q(1046834122639903380009,20000000000000000)


def fixtures(w):
    p=w.physical
    count=0
    for name,tree in p.fixtures(w.h):
        core,marks=p.physical_core(tree)
        terms=p.subsets(w.p,core,marks)
        for root in {min(core),max(core)}:
            parent={root:None}; todo=[root]
            for u in todo:
                for v in core[u]:
                    if v not in parent: parent[v]=u;todo.append(v)
            owned=dict.fromkeys(core,Q(0))
            for (family,vertices),(value,roles) in terms.items():
                if family=='g':
                    a,b=vertices; hosts=(a if parent[a]==b else b,)
                elif family=='p4': hosts=roles
                elif family=='p5': hosts=roles[1:2]
                else: hosts=roles[:1]
                for v in hosts: owned[v]+=Q(value,len(hosts))
            assert sum(owned.values())==sum(x for x,_ in terms.values())
            for u,(d,c) in marks.items():
                lo,hi=w.comparison_ticks(d,c)
                assert lo<=owned[u]+(w.alpha*(c-2) if c else 0) and owned[u]<=hi,(name,u,d,c)
                if 1<=c and d<w.DF:
                    l,h=w.finite(d);assert (lo,hi)==(l[c-1],h[c-1])
                count+=1
    return count


def objective(w,limit,sat,birth,edges,reduced=None):
    """B=independently certified f floor plus disjoint original edge rebates, ticks."""
    cells=w.cells(limit,sat)
    m=pinned('gold2_physical_lp',R2/'occupation-neon/pressure-memory/birth-law/replay.py',
             '7dd40ac41e5d95ae41843da51de436fd45b781f1')
    m.LO,m.HI=w.prior;m.SAT=sat
    lp=m.build(limit,sided=True)
    assert all(v>=0 for v in (birth,edges))
    def price(fn,scale):
        def obj(k,t):
            if k!='occ':return 0
            v=fn(*cells[t[:2]],t)
            assert v>=0
            return Q(v,scale)
        return max(0,w.prior[0]*scale*lp.lower(obj))
    L,M,H=w.global_low,w.global_high,w.higher_cap
    L2=L+price(lambda l,u,h,t:l-L,w.scale)
    M2=M-price(lambda l,u,h,t:M-u,w.scale)
    shift=w.sign*(TARGET-w.lam)
    anchor=M if shift>=0 else L
    def slack(t):
        return 0 if reduced is None else reduced(*t)
    h_saved=price(lambda l,u,h,t:H-h+slack(t),10**7)
    joint=price(lambda l,u,h,t:H-h+slack(t)+abs(shift)*(M-u if shift>=0 else l-L),10**7)
    base=birth+edges-w.golden['err_per_core']
    margin=base+joint-shift*anchor
    net=base+h_saved
    separate=w.lam+w.sign*net/(M2 if net>=0 else L2)
    if shift>=0:
        # Independently valid pricing may be combined by max, not addition.
        margin=max(margin,net-shift*M2)
    print('LP rows/columns',len(lp.rows),len(lp.cols),'cells',len(cells),flush=True)
    print('M lower/upper',L2,M2,'approx',float(L2),float(M2),flush=True)
    print('H savings / joint rebate ticks',h_saved,joint,'approx',float(h_saved),float(joint),flush=True)
    print('signed shift',shift,'external f / edge ticks',birth,edges,'target margin',margin,'approx',float(margin),flush=True)
    print('target margin >0:',margin>0,'additional ticks needed at most',max(0,-margin),'separate endpoint',float(separate),flush=True)
    return dict(side=w.side, limit=limit, saturation=sat, target=str(TARGET), nominal=str(w.lam),
                prior=list(map(str,w.prior)), f_floor=str(birth), edge_floor=str(edges),
                reduced_f=reduced is not None, cap_error=w.golden['err_per_core'],
                denominator=list(map(str,(L2,M2))), saved_plus_reduced=str(h_saved),
                joint_rebate=str(joint), signed_shift=str(shift), margin=str(margin),
                separate_endpoint=str(separate), rows=len(lp.rows), columns=len(lp.cols),
                cells={f'{d},{k}':list(map(str,v)) for (d,k),v in cells.items()})


def reduced_f(w,limit):
    """A Gold2-only lower lifetime test; price its actual receive residual once."""
    G=10**12; disc=w.prior[1]; C=30; tail=[0]*(C+2); kept={}
    # Zero at exit is a valid lower test since original f slack is nonnegative.
    # Stop at D=3999, with ample room beyond every finite LP cell.
    for d in range(w.DF-1,1,-1):
        cur=[0]*(C+2)
        for c in range(min(d,C),0,-1):
            a,b=w.slack_ticks(d,c); L=d-c
            gain=tail[c]; bud=cur[c+1]
            p=Q(L,d); den=2*d*disc.denominator
            b0=p*den; base=(a*G//1)*den+b0*bud
            n=base+L*den*gain; k=disc.numerator*2*d+b0+L*den
            dx=c*d*disc.denominator
            cur[c]=max(0,min(n//k,(n+dx*(gain+b*G//1))//(k+dx)))
            v=cur[c]
            assert all((a+b*r)*G+(L+r)*(gain-v)+p*(bud-v)>=disc*v for r in (0,Q(c,2)))
        if d<=limit:
            for c in range(1,d+1):
                a,b=w.slack_ticks(d,c)
                if c<=C:
                    v=cur[c]; jump=tail[c]-v
                    a+=Q((d-c)*jump,G)+Q(d-c,d)*Q(cur[c+1]-v,G)-disc*Q(v,G)
                    b+=Q(jump,G)
                assert min(a,a+Q(c,2)*b)>=0
                kept[d,c-1]=a,b
        tail=cur
    print('reduced f birth floor',w.prior[0]*Q(tail[1],G),'rows',len(kept),flush=True)
    def value(d,k,b,p,u,v):
        a,s=kept.get((d,k),(0,0));return a+s*(u+v)
    return w.prior[0]*Q(tail[1],G),value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--side',choices=('lower','upper'),action='append');p.add_argument('--limit',type=int,default=24)
    p.add_argument('--saturation',type=int,default=7)
    p.add_argument('--f-floor',type=Q);p.add_argument('--edge-floor',type=Q,default=Q(0))
    p.add_argument('--conditioned-f',action='store_true',help='use published downward-rounded Gold2 conditioned floors for both signs')
    p.add_argument('--completed-lower-edges',action='store_true',help='bind pinned D18, absorb12 and wide-payer receipts')
    p.add_argument('--reduced-f',action='store_true',help='use a derived birth plus actual reduced f; replaces --f-floor')
    p.add_argument('--work',type=Path)
    p.add_argument('--out',type=Path,help='save exact inputs, outputs and census cells under host/')
    a=p.parse_args()
    if sum((a.f_floor is not None,a.conditioned_f,a.reduced_f))>1:
        p.error('choose one singleton floor method')
    if a.edge_floor and a.completed_lower_edges: p.error('choose one edge floor input')
    from reader import HERE, RECEIPT
    if a.out and a.out.resolve().parent!=HERE: p.error('--out must be inside host/')
    results={}
    for side in a.side or ('lower','upper'):
        with open_host(side,**({'work':a.work} if a.work else {})) as w:
            print(side,'authenticated; actual-degree header',w.DF,w.d2,w.K2,flush=True)
            print('physical local comparisons',fixtures(w),'potential words reduced',w.weight_words,flush=True)
            for d,c in ((2,1),(12,4),(251,22),(251,23),(252,2),(3999,2),(4000,2),(10**8,50)):
                x,y=w.slack_ticks(d,c); assert min(x,x+y*Q(c,2))>=0
            print('singleton all-stock c=1',w.stock_floor(1),flush=True)
            if a.reduced_f:
                b,fun=reduced_f(w,a.limit)
            elif a.conditioned_f:
                report=HERE.parent/'lifetime/README.md'
                import re
                text=report.read_text()
                from hashlib import sha1
                raw=text.encode()
                assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==CONDITIONED_REPORT
                floors=re.search(r'Conditioned floor \| \*\*([\d.]+)\*\* \| \*\*([\d.]+)\*\*',text).groups()
                b,fun=Q(floors[side=='upper']),None
            else: b,fun=a.f_floor or Q(0),None
            edge=a.edge_floor; roots={}; sources={}
            if a.completed_lower_edges and side=='lower':
                import re
                for f,h in (('edge18','1061d512219da5a38f264c07b4082f1db391b213a885f4cd955ea691946e2948'),
                            ('absorb12','eded191d098bb2c9296c5d7d2861a051bc303bce122953c518ebf1981d22d521')):
                    name=f'edges/{f}.log'
                    raw=(HERE.parent/name).read_bytes(); assert sha256(raw).hexdigest()==h
                    sources[name]=h
                    assert raw.decode().splitlines()[-1].startswith('PASS')
                    pat=r'credit ticks/core (\S+)' if f=='edge18' else r'^TOTAL (\S+)'
                    roots[f]=Q(re.search(pat,raw.decode(),re.M)[1])*(1 if f=='edge18' else w.scale)
                name='payer/README.md'
                raw=(HERE.parent/name).read_bytes()
                assert sha256(raw).hexdigest()==PAYER_REPORT_SHA
                normalized=Q(re.search(r'Exact normalized value: `([0-9]+/[0-9]+)`',raw.decode())[1])
                roots['old_payer']=w.scale*normalized
                assert roots['old_payer']==PAYER_TICKS
                sources[name]=PAYER_REPORT_SHA
                edge=sum(roots.values())
            results[side]=objective(w,a.limit,a.saturation,b,edge,fun)
            if roots:
                results[side]['edge_components_ticks']={k:str(v) for k,v in roots.items()}
                results[side]['edge_receipt_sha256']=sources
                results[side]['old_payer_normalized']=str(normalized)
            results[side]['payload_sha256']=w.docs[side+'_sha256']
            if a.conditioned_f: results[side]['f_source_git_blob']=CONDITIONED_REPORT
    if a.out:
        data=(dumps(dict(original_receipt_sha256=RECEIPT,results=results),indent=2,sort_keys=True)+'\n').encode()
        a.out.write_bytes(data)
        print('exact export SHA-256',sha256(data).hexdigest(),flush=True)
    print('PASS exact Gold2 host, overflow and joint inputs',flush=True)


if __name__=='__main__':main()
