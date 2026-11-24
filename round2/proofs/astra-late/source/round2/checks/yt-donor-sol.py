#!/usr/bin/env python3
"""Independent Yt table reconstruction, labelled controls, integer lifetime and ledger."""
import argparse
from fractions import Fraction as Q
from functools import cache
from hashlib import file_digest
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import ceil, floor, lcm
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
O=ROOT/'round2/occupation-neon'
F=O/'pressure-memory/fork-lifetime'
PINS={F/'old/README.md':'5012b327fa0ac25af9e9fb9617aa239410d907be4b456b5e92c8dc09976a44a5',
      F/'old/replay.py':'622a5ac57871d44c4852bab3eb2f616a7354ac977fb16ad594da48e9141b35a3',
      F/'replay.py':'db0fba60406930918139d59731e6cac545bc877955508aa8f819967def980999',
      O/'all-contrast/scan.py':'6f2a074e9a695b4969cb58d39376e2a1f8c17ab41b04a3b2ec775a18f05ca02d',
      O/'higher/extract.py':'9f40a50a3623756b277243d192f672e76c45b4b9784715ed4d0d8fb27d36131c',
      O.parent/'attacks/neon1/audit.py':'d548a9aafa899ecab2b00a7f6ad0a4e6f828600088780f0a15f08927f17f0bfe',
      O/'pressure-memory/hv-lifetime/release/inputs.sha256':'e1eb82077a9100e1f35b6d02732948198664a413346be5185681940f7d7f3bb8',
      O/'pressure-memory/hv-lifetime/release/release.sha256':'8db4750e9957b2a48c2525dadd6f2a033063220d7d864c2e10c8369d5d3bd91d',
      O/'pressure-memory/path-fast/low-release/README.md':'9d12d1158f2210eb28cafdd16dbbbbb9c9512f207dd1449e30aab9609a748e97',
      O/'pressure-memory/path-fast/low-release/release.sha256':'b004e6df46b8d63f1a6dbe937c2ed01af050130a38584fc3520aec09134eb474'}
BLOB=('9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
      '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41')
LO,HI=Q('.565673'),Q('.565775'); GRID=10**16; U=lcm(10**6,*range(2,17))
V=(2,0); PAR=(-1,0,1,1); CH=((1,),(2,3),(),()); Z=(0,0,0,0)
CAP=(16,5,5,5); K=(5,3,2,2)
PARTS=('mother','apex','trunk','tips')
EXPECTED={
 'mother':(54088,'92700142673180172981/2065673000000000000000000000','153409020386489/400000000000000000000',
                  '13401473874619/5000000000000000000000','73222219688163/200000000000000000000'),
 'apex':(377798,'106499397693923/2500000000000000000000','13560645764043/50000000000000000000',
                 '182096203280233/2000000000000000000000','8660368912643/20000000000000000000'),
 'trunk':(207878,'365035451093613/5000000000000000000000','15366877880081/80000000000000000000',
                  '222918995717350291137/5164182500000000000000000000','69215475548891/400000000000000000000'),
 'tips':(501548,'1434768825321273847891/10328365000000000000000000000','198547226598651/400000000000000000000',
                 '318421100779199/1000000000000000000000','164249049594217/200000000000000000000')}


def sha(p):
    with p.open('rb') as f:return file_digest(f,'sha256').hexdigest()


def load(n,p):
    s=spec_from_file_location(n,p);m=module_from_spec(s);s.loader.exec_module(m);return m


def setat(x,i,v):return x[:i]+(v,)+x[i+1:]


def box(c):return min(Q(0),c/2)


def raw(c,path):
    """Derive Yt from its own potential and the three connected deletions."""
    x=tuple((d-k-1,k,1) for d,k in path)
    if x[1] not in c.yt:return None
    def yt(xx):
        a,t,l,r=xx
        return int(c.yttab[c.par[a],c.yt[t],c.tp(l,r)]) if t in c.yt else 0
    def hv(a,t,b):return int(c.hvtab[c.par[a],c.hv[t],c.child[b]]) if t in c.hv else 0
    def hc(t,b,d):return int(c.hctab[c.hc[t],c.tp(b,d)]) if t in c.hc else 0
    def plus(a):return (a[0]+1,a[1],a[2])
    a,t,l,r=x; old=yt(x); base=-c.lam*old; dif=[]
    for i,y in enumerate(x):
        L,k,h=y;gain=yt(setat(x,i,plus(y)))-old
        nb=[j for j in range(4) if PAR[j]==i or PAR[i]==j]
        base+=(L+sum(Q(1,path[j][0]) for j in nb))*gain
        if L:base+=Q(L,path[i][0])*(yt(setat(x,i,(L-1,k+1,h)))-old)
        dif.append(gain)
    bch,bpar=(Q(int(rates[c.state[t]]),c.Sr) for rates in (c.rch,c.rpar))
    for present,removed in ((l,3),(r,2)):
        base+=(Q(1,path[removed][0])-bch)*(hv(a,plus(t),present)-hv(a,t,present))
    base+=(Q(1,path[0][0])-bpar)*(hc(plus(t),l,r)-hc(t,l,r))
    eps=1 if c.mode==2 else -1
    counts=(1,*(k-len(kids) for (_,k),kids in zip(path,CH)))
    slopes=tuple(Q(eps*g,c.S) if n else Q(0) for g,n in zip((dif[0],*dif),counts))
    signed=Q(eps*base,c.S)
    return signed,slopes


def geometry(x,atom):
    omitted=sum(k-len(CH[i]) for i,(_,k) in enumerate(x))
    size=4+sum(d-k-1 for d,k in x)+2*omitted
    cores=4+omitted
    if atom:
        i,d,k=atom
        assert x[i][1]>len(CH[i])
        size+=abs(d)-1+k-1;cores+=k
    return size,cores


def labelled():
    # Build a permanent apex with a mother and nominated biological children.
    # Check every legal chronology including extras and both exterior tip roles.
    ncases=0
    for carrier in (0,1,2,3):
        for ev in permutations('TLRXA'):
            if ev.index('T')>min(ev.index('L'),ev.index('R')):continue
            if carrier and ev.index('X')<ev.index('TLR'[carrier-1]):continue
            # Enumerate independently using physical parent IDs and direct combinations.
            parent={'a':'m','T':'a','L':'T','R':'T','A':'a','X':('a','T','L','R')[carrier]}
            tags={(('a',),None)};born={'a'}
            for w in ev:
                g=parent[w];born.add(w)
                assert g in born
                added=[]
                for (slots,atom) in tags:
                    n=len(slots)
                    if n<4 and slots[PAR[n]]==g:added.append((slots+(w,),atom))
                    roles=(2,3) if carrier>=2 else (carrier,)
                    if atom is None and any(i<n and slots[i]==g for i in roles):added.append((slots,w))
                assert len(added)==len(set(added)) and not tags.intersection(added)
                tags.update(added)
                direct=set()
                for t in [v for v in ev if v in born and parent[v]=='a']:
                    for b,c in combinations([v for v in ev if v in born and parent[v]==t],2):
                        ids=('a',t,b,c)
                        for e in [None,*[v for v in born if v not in ids and parent[v] in
                                          (ids[2:] if carrier>=2 else (ids[carrier],))]]:
                            direct.add((ids,e))
                assert {z for z in tags if len(z[0])==4}==direct
            ncases+=1
    # Actual minimal trees assembled from known census, including an aged donor
    # whose child/leaf identities make it pass a gate the unmarked crown fails.
    x=((4,2),(3,2),(2,0),(2,0)); atom=(0,5,1)
    parent=(-1,0,1,2,2,1,1,3,4,5,5,5,5,12)
    # Root 0 is the mother. Tagged cores are 1..4; outside donor 5 has
    # three leaf children and one promoted child 12, which has its own leaf.
    deg=[0]*len(parent)
    for v,u in enumerate(parent[1:],1):deg[u]+=1;deg[v]+=1
    states=lambda i:(deg[i],sum(deg[j]>=2 for j,u in enumerate(parent) if u==i))
    assert tuple(map(states,(1,2,3,4)))==x and states(5)==(5,1)
    assert geometry(x,())==(9,5) and geometry(x,atom)==(13,6)==geometry(x,(0,-5,1))
    assert len(parent)-1==13 and sum(d>=2 for d in deg[1:])==6
    assert Q(1,deg[5])==Q(1,5) and Q(1,deg[1])==Q(1,4)
    print('labelled chronology cases',ncases,'; physical donor-only geometry witness (9,5) -> (13,6)',flush=True)


def table_checks(S,C):
    rows=[S.Row(c) for c in C]; scanned=0
    slots=[[(d,k) for d in range(2,D+1) for k in range(len(CH[i]),min(K[i],d-1)+1)]
           for i,D in enumerate(CAP)]
    table={}
    for x in product(*slots):
        dsk=tuple((d-k-1,k,1) for d,k in x)
        both=[]
        for w in rows:
            ref=raw(w.c,x)
            assert (ref is not None)==w.present('Yt',dsk)
            if ref is None:continue
            a,p=w.coefficients('Yt',dsk)
            assert ref==(a,(p[0],*p[1::2]))
            assert not any(p[2::2])
            n=(1,*(k-len(CH[i]) for i,(_,k) in enumerate(x)))
            assert a+sum(j*box(c) for j,c in zip(n,ref[1]))>=0
            swap=raw(w.c,x[:2]+x[2:][::-1]);assert swap==(a,ref[1][:-2]+ref[1][-2:][::-1])
            both.append(ref[1]);scanned+=1
        assert len(both) in (0,2)
        if both:table[x]=tuple(both)
    print('exact signed affine row/port/tip-exchange comparisons',scanned,flush=True)
    return table


def extrema(lane,reward,clocks):
    # Maximize the signed fraction by enumerating every breakpoint corner.
    sign=1 if lane%2 else -1
    d=int(U*LO);n=sign*reward*U; extra=[(0,int(U*(HI-LO)))];
    for a,b,y in clocks:
        z=sign*y[lane]; n+=a*z;d+=a
        if a<b:extra.append((z,b-a))
    origin_n,origin_d=n,d
    best_n,best_d=n,d
    for z,h in sorted(extra,reverse=True):
        n+=z*h;d+=h
        if n*best_d>best_n*d:best_n,best_d=n,d
    tick=-((-best_n)//best_d)
    assert origin_n-tick*origin_d+sum(h*max(0,z-tick) for z,h in extra)<=0
    return sign*tick


def solve(table,part):
    carriers={'mother':(),'apex':(0,),'trunk':(1,),'tips':(2,3)}[part]
    state_count=0
    rc=[0]+[U//d for d in range(1,17)]

    @cache
    def toll(x,port):
        both=table.get(x)
        if both is None:return Z,Z
        a,b=(s[port] for s in both)
        ticks=lambda e,f:(floor(e*GRID),ceil(e*GRID),floor(f*GRID),ceil(f*GRID))
        return ticks(-box(a),-box(b)),ticks(a,b)

    def value(x,m,e):
        n=len(x)
        outstanding=tuple(sum(j>=n for j in CH[i]) for i in range(n))
        if any(d>CAP[i] or k+outstanding[i]>K[i] for i,(d,k) in enumerate(x)):return Z
        if carriers and not e:
            if not any((x[i][1]+outstanding[i]<K[i] if i<n else K[i]>len(CH[i])) for i in carriers):return Z
        if n==4 and (x[2],e and e[0]==2)<(x[3],e and e[0]==3):
            # Pick the reverse representative from the producer at ties and strict order.
            x=x[:2]+(x[3],x[2]);e=(5-e[0],*e[1:]) if e and e[0]>=2 else e
        return backward(x,m,e)

    @cache
    def backward(x,m,e):
        nonlocal state_count
        state_count+=1;n=len(x)
        r,s=Z,Z
        if n==4 and (not carriers or e):
            size,cores=geometry(x,e)
            if size>=21 or size>=13 and cores>=5:r,s=toll(x,e[0]+1 if e else 0)
        clocks=[]
        def add(lo,hi,z,w=Z):
            assert 0<=lo<=hi
            if hi:clocks.append((lo,hi,tuple(a+b for a,b in zip(z,w)) if w!=Z else z))
        for i,(d,k) in enumerate(x):
            children=[j for j in CH[i] if j<n];L=d-k-1
            inc=value(setat(x,i,(d+1,k)),m,e)
            internal=L*U+sum(rc[x[j][0]] for j in children)+ (rc[x[PAR[i]][0]] if i else 0)
            add(internal,internal,inc)
            if i==0:add(rc[m] if m>0 else 0,rc[abs(m)],inc,s if not carriers else Z)
            marked=e and i==e[0]
            if marked:
                a=e[1];add(rc[a] if a>0 else 0,rc[abs(a)],inc,s)
            add(0,(k-len(children)-bool(marked))*rc[2],inc)
            if L:
                bud=setat(x,i,(d,k+1)); vals=[value(bud,m,e)]
                if n<4 and PAR[n]==i:vals.append(value(bud+(V,),m,e))
                if not e and i in carriers:vals.append(value(bud,m,(i,2,0)))
                add(L*rc[d],L*rc[d],tuple(map(sum,zip(*vals))))
        if e and e[1]>0:
            i,d,k=e; L=d-k-1;low=L*U+rc[x[i][0]]
            add(low,low+k*rc[2],value(x,m,(i,3 if d==2 else -4,k)))
            if L:add(L*rc[d],L*rc[d],value(x,m,(i,d,k+1)))
        if m==2:
            low=rc[x[0][0]];add(low,low+rc[2],value(x,-3,e))
        return tuple(extrema(j,y,clocks) for j,y in enumerate(r))

    virgin,older=(value((V,),m,()) for m in (2,-3))
    weight=Q(1,2)/(LO+Q(3,2))
    result=tuple(op((older[j], older[j]+weight*(virgin[j]-older[j])))*scale/GRID
                 for j,(op,scale) in enumerate(((min,LO),(max,HI),(min,LO),(max,HI))))
    count=state_count;backward.cache_clear();toll.cache_clear()
    assert result==tuple(map(Q,EXPECTED[part][1:])) and count==EXPECTED[part][0],(part,count,result)
    print(part,count,'states',*(str(v) for v in result),flush=True)
    return result


def compose(results):
    sums=tuple(sum(a[i] for a in results.values()) for i in (0,2))
    assert sums==tuple(map(Q,('3092299274719269569703/10328365000000000000000000000',
                             '9405336866046559940037/20656730000000000000000000000')))
    den=tuple(map(Q,('1357493815758092773/1000000000000000000','273436004214812983/200000000000000000')))
    trials=tuple(map(Q,('.565677930574','.565770634792')))
    net=tuple(map(Q,('6817529487964063107239353589443/8925000000000000000000000000000000000',
                     '280520213267354027660117849251/425000000000000000000000000000000000')))
    hv=(Q(10361983749,10**16),Q(14352363253,10**16))
    old=tuple(t+s*(n+h)/d for t,s,n,h,d in zip(trials,(1,-1),net,hv,den))
    assert old==tuple(map(Q,('6853561875850472647665032629043939443/12115632305640977999025000000000000000',
                             '328740990573504279104102687021150749/581051508956477588875000000000000000')))
    old_p4=(Q(359297001126132137,10**24),Q(6963149521191559,25000000000000000000000))
    new_p4=(Q(23575773396079784114802186917,25500000000000000000000000000000000),
            Q(43746697941322995419,56100000000000000000000000))
    show=lambda v,b:f'{(ceil(v*10**15) if b else floor(v*10**15))//10**15}.{(ceil(v*10**15) if b else floor(v*10**15))%10**15:015d}'
    wide=tuple(map(Q,('6853566920645426224860150084809360393/12115632305640977999025000000000000000',
                      '10848445658578037268151132921697974717/19174699795563760432875000000000000000')))
    for name,base,delta in [('f7cc2e1',old,(0,0)),('f4f7ea3',wide,tuple(a-b for a,b in zip(new_p4,old_p4)))]:
        assert base==tuple(a+s*b/d for a,s,b,d in zip(old,(1,-1),delta,den))
        endpoints=tuple(a+s*w/d for a,s,w,d in zip(base,(1,-1),sums,den))
        assert endpoints==tuple(t+s*(n+h+b+w)/d for t,s,n,h,b,w,d in zip(trials,(1,-1),net,hv,delta,sums,den))
        print(name,'+ Yt exact',*(str(e) for e in endpoints),'outward',show(endpoints[0],0),show(endpoints[1],1),flush=True)
    print('exact Yt sums',*sums,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('payload',nargs=2,type=Path)
    p.add_argument('--parts',nargs='+',choices=PARTS,default=PARTS)
    p.add_argument('--ledger-only',action='store_true',help='recompose already verified receipt fractions')
    a=p.parse_args();assert __debug__
    for f,h in PINS.items():assert sha(f)==h,f
    if a.ledger_only:
        compose({k:tuple(map(Q,v[1:])) for k,v in EXPECTED.items()});return
    for f,h in zip(a.payload,BLOB):assert sha(f)==h,f
    b=load('fork_paths',F/'replay.py') # imports scanner only, not its solver
    import scan
    from extract import Cert
    C=tuple(Cert(p) for p in a.payload)
    assert [(c.mode,c.D,c.S,c.lam) for c in C]==[(2,18,10**11,LO),(1,18,10**11,HI)]
    assert b.LOW==LO and b.HIGH==HI and b.GRID==GRID
    labelled();table=table_checks(scan,C)
    results={part:solve(table,part) for part in a.parts}
    if set(results)==set(PARTS):compose(results)


if __name__=='__main__':main()
