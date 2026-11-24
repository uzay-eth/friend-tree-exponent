"""Independent exact motif/transfer/event/charge checks on core-compressed trees.

Only the byte decoder is reused; the shipped potential and drift functions are unused.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cert'))
from tree_identity_check import load_cert


def run(C):
    d,kt,lf,kf=(C[x] for x in ('D0','KT','LF','KF'))
    inf=C['INF_ID']; new=C['lut'][1,0,1]; S=C['S']; Sr=C['S_r']
    elm=lambda name,j:C[name][j]
    p=C['p_to_pi'].__getitem__; w=C['w_to_wi'].__getitem__
    n=C['nchs']; nt=kt+2
    pair=lambda a,b: sum(comb(n-i,1) for i in range(a))+b-a
    triple=lambda a,b,c: comb(n+2,3)-comb(n-a+2,3)+comb(n-a+1,2)-comb(n-b+1,2)+c-b
    @lru_cache(maxsize=None)
    def table(tag, s):
        if tag=='g': return elm('g',s[0]*C['NG']+s[1])
        if tag=='hV':
            a,b,c=s; u=C['hu_of_g'][b]
            return elm(tag,(p(a)*C['NHU']+u)*n+w(c)) if u>=0 else 0
        if tag=='hC':
            a,b,c=s; u=C['cu_of_g'][a]
            return elm(tag,u*comb(n+1,2)+pair(*sorted((w(b),w(c))))) if u>=0 else 0
        if tag=='qP4':
            a,b,c,e=s; N=C['NP4']
            return elm(tag,((p(a)*N+b)*N+c)*n+w(e))
        if tag=='qYt':
            a,b,c,e=s
            return elm(tag,(p(a)*C['NYT']+b)*comb(n+1,2)+pair(*sorted((w(c),w(e)))))
        if tag=='qYb':
            a,b,c,e=s
            return elm(tag,((a*C['NYBW']+b)*n+w(c))*n+w(e))
        if tag=='qS3':
            a,b,c,e=s
            return elm(tag,a*comb(n+2,3)+triple(*sorted((w(b),w(c),w(e)))))
        raise AssertionError(tag)

    @lru_cache(maxsize=None)
    def embed(L,k,h):
        D=L+k+h
        if D<d: return ((C['lut'][h,k,L],Q(1)),)
        t=Q(d,D)
        ids=(inf,*(C[f't{j}'](k,h) for j in (1,2,3)))
        return tuple((idx,comb(3,j)*t**j*(1-t)**(3-j)) for j,idx in enumerate(ids))

    def mix(tag,es,func=None):
        return sum(prod(z for _,z in terms)*(table if func is None else func)(tag,tuple(i for i,_ in terms))
                   for terms in product(*es) if all(z for _,z in terms))

    def snapshot(par,leaves):
        ch=[[] for _ in par]
        for i,p0 in enumerate(par[1:],1): ch[p0].append(i)
        k=[len(z) for z in ch]; h=[int(i>0) for i in range(len(par))]
        D=[a+b+c for a,b,c in zip(leaves,k,h)]
        assert min(D)>=2
        em=[embed(*x) for x in zip(leaves,k,h)]
        nbr=[([par[u]] if h[u] else [])+ch[u] for u in range(len(par))]
        assert all(sum(q*(Q(1,D[u]) if D[u]<d else Q(j,3*d)) for j,(_,q) in enumerate(em[u]))==Q(1,D[u]) for u in range(len(par)))
        assert sum(Q(a-b,4) for a,b in zip(k,h))==0
        m={}; charge=[Q(0)]*len(par)
        def add(tag, vs, to, fixed=None):
            x=(tag,tuple(sorted(vs))); assert x not in m
            es=[((fixed[j],Q(1)),) if fixed and j in fixed else em[u] for j,u in enumerate(vs)]
            m[x]=mix(tag,es); charge[to]+=m[x]
        for u in range(len(par)):
            a,b,c=leaves[u],k[u],h[u]; p0=par[u]
            i=(C['VSTAR1'] if c else C['VSTAR0']) if b>kf else ((0 if c else kf+1)+b)*(lf+1)+min(a,lf)
            m['f',(u,)]=elm('f',i); charge[u]+=m['f',(u,)]
            if c: add('g',(p0,u),u)
            if c and 1<=b<=kt:
                for v in ch[u]: add('hV',(p0,u,v),v)
            if 2<=b<=kt:
                for v,z in combinations(ch[u],2): add('hC',(u,v,z),u)
            ui=C['p4_lut'].get((b,a)) if c else None
            if ui is not None and p0>0:
                pi=C['p4_lut'].get((k[p0],leaves[p0]))
                if pi is not None:
                    for v in ch[u]: add('qP4',(par[p0],p0,u,v),v,{1:pi,2:ui})
            yi=C['yt_lut'].get((b,a)) if c else None
            if yi is not None:
                for v,z in combinations(ch[u],2): add('qYt',(p0,u,v,z),u,{1:yi})
            yi=C['ybu_lut'].get((c,b,a))
            if yi is not None:
                for v in ch[u]:
                    zi=C['ybw_lut'].get((k[v],leaves[v]))
                    if zi is not None:
                        for z in ch[u]:
                            if z!=v:
                                for j in ch[v]: add('qYb',(u,v,z,j),j,{0:yi,1:zi})
            si=C['s3_lut'].get((c,b,a))
            if si is not None:
                for v,z,j in combinations(ch[u],3): add('qS3',(u,v,z,j),u,{0:si})
        assert sum(charge)==sum(m.values())
        return m,k,h,D,em,nbr,charge

    def checks(par,L,name):
        m,k,h,D,em,nbr,charge=snapshot(par,L)
        core=len(par); nb=core; phi=sum(m.values()); lphi=Q(0)
        c0,C0,phi3=((Q(113465364604,4),184124847315,146042017510) if C['mode']==2 else
                   (Q(115462837236,4),195820675182,146044806770))
        assert all(c0<=val+Q(S*(a-b),4) and val<=C0 for val,a,b in zip(charge,k,h))
        if core==1 and L==[2]: assert phi==phi3
        assert c0*core<=phi<=C0*core
        assert sum(L)+sum(map(len,nbr))==sum(D)
        moved=Counter(); births=Counter(); levels=defaultdict(Q)
        for u in range(core):
            # Reweight leaf-gain coefficients by moving each external inverse-degree
            # deviation to the unique motif enlarged by its source neighbour.
            plus=L.copy(); plus[u]+=1
            up=snapshot(par,plus)[0]
            assert up.keys()<=m.keys()
            gain=L[u]+sum(Q(1,D[z]) for z in nbr[u]); bud=Q(L[u],D[u]); delta=sum(up.values())-phi
            lphi+=gain*delta
            local={key:gain*(up.get(key,0)-val) for key,val in m.items()}
            destinations=defaultdict(list)
            for key in m: destinations[key[1]].append(key)
            for key,val in m.items():
                tag,vs=key
                if u not in vs or len(vs)>=4: continue
                for z in nbr[u]:
                    if z in vs: continue
                    tags=destinations[tuple(sorted((*vs,z)))]; assert len(tags)<=1
                    if not tags: continue
                    dest=tags[0]; adj=Q(1,D[z]); parent=(z==par[u]); rate='par' if parent else 'ch'
                    tailidx=(0 if h[u] else nt)+min(k[u],kt+1)
                    coeff=Q(0)
                    invz = Q(1,D[z]) if D[z]<d else None
                    inv_controls=[Q(j,3*d) for j in range(4)] if invz is None else [invz]
                    if D[u]<d:
                        rb=Q(elm('r'+rate,em[u][0][0]),Sr)
                        coeff=up.get(key,0)-val
                        term=(adj-rb)*coeff; kind='inner'
                        received=sum(q*(r-rb)*coeff for (_,q),r in zip(em[z],inv_controls))
                    elif tag=='f':
                        if k[u]>kt: continue
                        rb=Q(elm('tr'+rate,tailidx),Sr)
                        a,b=elm('c2_f',tailidx),elm('c3_f',tailidx)
                        coef=Q(a,D[u]**2)+Q(b,D[u]**3)
                        assert coef==sum(q*v for (_,q),v in zip(em[u],(0,0,Q(a,3*d*d),Q(a,d*d)+Q(b,d**3))))
                        term=(adj-rb)*coef; kind='tail-c2/c3'
                        received=sum(q*qz*cv*(r-rb) for ((_,q),cv),((_,qz),r) in product(
                            zip(em[u],(0,0,Q(a,3*d*d),Q(a,d*d)+Q(b,d**3))),zip(em[z],inv_controls)))
                    elif tag=='g':
                        if k[u]>kt: continue
                        rb=Q(elm('tr'+rate,tailidx),Sr)
                        a,b=vs; other=b if u==a else a
                        def edge(q,i): return elm('g',(q if u==a else i)*C['NG']+(i if u==a else q))
                        control=em[u][1][0]
                        diff=sum(weight*Q(edge(inf,i)-edge(control,i),d) for i,weight in em[other])
                        ref=sum(weight*(elm('Yp_ref',tailidx*n+w(i)) if u==a else elm('Yw_ref',p(i)*nt+k[u])) for i,weight in em[other])
                        term=(adj-rb)*(em[u][2][1]*diff+em[u][3][1]*ref); kind='tail-T2/Yref'
                        received=Q(0)
                        for i,weight in em[other]:
                            f2=Q(edge(inf,i)-edge(control,i),d)
                            f3=elm('Yp_ref',tailidx*n+w(i)) if u==a else elm('Yw_ref',p(i)*nt+k[u])
                            for (_,qz),r in zip(em[z],inv_controls):
                                received+=weight*qz*(r-rb)*(em[u][2][1]*f2+em[u][3][1]*f3)
                    else:
                        # The extra vertex must extend an inner slot, not merely
                        # make a shape which has some other inner slot.
                        raise AssertionError((key,dest,u))
                    assert term==received
                    local[key]-=term; local[dest]+=received
                    moved[tag,dest[0],kind]+=1
            assert sum(local.values())==gain*delta
            if bud:
                minus=L.copy(); minus[u]-=1
                mb,kb,hb,Db,emb,_,_=snapshot(par+[u],minus+[1])
                assert Db[u]==D[u] and emb[nb]==((new,Q(1)),)
                if D[u]>=d: assert tuple(x[1] for x in em[u])==tuple(x[1] for x in emb[u])
                local_b={key:mb.get(key,0)-val for key,val in m.items()}
                for key in mb.keys()-m.keys():
                    tag,vs=key
                    assert tag=='f' or nb in vs
                    if tag=='f': dest=('f',(u,)); assert vs==(nb,)
                    else:
                        old=tuple(x for x in vs if x!=nb)
                        ds=destinations[old]; assert len(ds)==1,(key,old,ds)
                        dest=ds[0]
                    local_b[dest]+=mb[key]
                    births[dest[0],tag,'inner' if D[u]<d else 'tail']+=1
                if k[u]==kt:
                    gone=[key for key in m.keys()-mb.keys() if key[0] in ('hV','hC') and u in key[1]]
                    # hC centered at u, and hV centered at u iff non-root.
                    assert len(gone)==comb(kt,2)+h[u]*kt
                    assert all(nb not in key[1] or key[0] not in ('hV','hC') for key in mb)
                assert sum(local_b.values())==sum(mb.values())-phi
                lphi+=bud*sum(local_b.values())
                for key,val in local_b.items(): levels[key[0]]+=bud*val
            for key,val in local.items(): levels[key[0]]+=val
        assert sum(levels.values())==lphi
        assert sum(L)+core==sum(Q(L[u],D[u])+L[u]+sum(Q(1,D[z]) for z in nbr[u]) for u in range(core))
        # Independent physical-vertex enumeration when small. The large cases
        # aggregate identical leaf attachments, retaining their exact multiplicity.
        if sum(L)+core<=100:
            full=par+[u for u in range(core) for _ in range(L[u])]
            adj=[[] for _ in full]
            for v,p0 in enumerate(full[1:],1): adj[v].append(p0); adj[p0].append(v)
            brute=Q(0)
            for v,p0 in enumerate(full):
                ll=L.copy(); pp=par
                if v<core: ll[v]+=1
                else: ll[p0]-=1; ll.append(1); pp=par+[p0]
                brute+=sum(Q(1,len(adj[j])) for j in adj[v])*(sum(snapshot(pp,ll)[0].values())-phi)
            assert brute==lphi
        res=lphi-Q(C['P'],C['Q'])*phi
        assert res*(1 if C['mode']==1 else -1)<=0,(name,res)
        print(f'{name}: n={core+sum(L)} core={core} drift residual/S sign={(res>0)-(res<0):+d}; all exact checks passed',flush=True)
        return moved,births

    def hub(k,L,h):
        return ([-1,0]+[1]*k,[1,L]+[1]*k) if h else ([-1]+[0]*k,[L]+[1]*k)

    cases=[('T3',[-1],[2]),('L4-interlocking',[-1,0,0,1,1,1,3,3,4,6,6],[1,2,1,1,1,1,1,1,1,1,1]),
           ('mixed-tail-edge',[-1,0,1,1,2,2,3,5],[d,3,d+3,d+7,1,2,1,1]),
           ('inner-to-tail',[-1,0,1,1],[d-2,d-4,1,1])]
    for k,L,h in ((kt,3,0),(kt,5,1),(kt+1,2,0),(kt+1,2,1),(kf,1,0),(kf+1,1,1),
                  (kf+2,0,0),(kf+2,lf+1,1),(3,lf,0),(3,lf+1,1),(0,lf+2,0)):
        cases.append((f'k{k}-L{L}-h{h}',*hub(k,L,h)))
    M=Counter(); N=Counter()
    for name,par,L in cases:
        m,b=checks(par,L,name); M.update(m); N.update(b)
    print('transfer categories:',sorted(M.items()))
    print('birth categories:',sorted(N.items()))
    assert {(x[0],x[1]) for x in M}=={('f','g'),('g','hV'),('g','hC'),('hV','qP4'),('hV','qYt'),('hV','qYb'),('hC','qYt'),('hC','qYb'),('hC','qS3')}
    assert {(x[0],x[1]) for x in N}=={('f','f'),('f','g'),('g','hV'),('g','hC'),('hV','qP4'),('hV','qYt'),('hV','qYb'),('hC','qYb'),('hC','qS3')}
    for dst in ('hV','hC'): assert M['g',dst,'tail-T2/Yref'] and N['g',dst,'tail']
    # Exact labelled growth law through size 8. This is a small independent
    # sanity check on the product's indexing, first step, and inequality signs.
    law={(-1,0,0):Q(1)}; factor=Q(1); lam=Q(C['P'],C['Q'])
    for size in range(3,9):
        mean=Q(0); nxt={}
        for par,pr in law.items():
            adj=[[] for _ in par]
            for v,p0 in enumerate(par[1:],1): adj[v].append(p0); adj[p0].append(v)
            mean+=pr*sum(len(row)>=2 for row in adj)
            if size<8:
                for v,row in enumerate(adj): nxt[par+(v,)]=pr*sum(Q(1,size*len(adj[j])) for j in row)
        assert sum(law.values())==1
        bound=(Q(146042017510,184124847315) if C['mode']==2 else Q(584179227080,115462837236))*factor
        assert (bound-mean)*(1 if C['mode']==1 else -1)>=0
        factor*=1+lam/size; law=nxt
    print(f'exact labelled finite-size product sanity check, n=3..8; E[C_8]={mean}')


if __name__=='__main__':
    run(load_cert(sys.argv[1]))
