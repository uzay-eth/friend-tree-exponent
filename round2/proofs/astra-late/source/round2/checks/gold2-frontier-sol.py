#!/usr/bin/env python3
"""Independent GOLD2 frontier checks: physical birth flux and signed repair."""
import argparse
import json
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import ceil, floor, lcm
from mmap import ACCESS_READ, mmap
from pathlib import Path
from runpy import run_path

R = Path(__file__).resolve().parents[1]
CHECKS = R/'checks'
G = R/'occupation-quintet-gold2'
FRONTIER = G/'host/frontier'
S = 2**40
GRID = 10**12
SIDES = ('lower','upper')
COMPACT = 'diagnostic/lower-sibling/'
COMPACT_PINS = {
    COMPACT+'price.py':'589e5d8d14f91621d0f5ad6bfecfdc00e4018a05bb3100416d0de8b22a4587f1',
    COMPACT+'compact.log':'c3d7d979fe82c4ac626cbdbef5f461411e966abd4e53ab35529a9cd5bb75c02a',
}
COMPACT_TICKS = dict(sibling=Q(16363762149897743918379249111,65029541606195200000000),
                     upstream=Q(1742793847212361294561574577,65029541606195200000000),
                     box=Q(26558583743732931982533,65029541606195200000000))
PINS = {
    CHECKS/'gold2-host-h3-sol.py':'ef62e72358be9f868e3a39ab996292f861b93ee264b1c7756c8364019a3697b7',
    G/'lifetime/pressure/frontier/replay.py':'1add7809bf6f509432b3931f2470e089eb8faf55ef0ca97a8828f18314ffa7df',
    G/'lifetime/pressure/replay.py':'a2d3a053ba785abad1a338bc929eae8c70cbdde367d4edce08dd93fc632a962c',
    G/'lifetime/check.py':'ddae6be5f70649fcd778f872215d5d0b9826aecea061eefc198463369f756eca',
}


def raw(path):
    data=path.read_bytes()
    assert sha256(data).hexdigest()==PINS[path],path
    return data


def load(path):
    spec=spec_from_file_location(path.stem,path)
    m=module_from_spec(spec)
    exec(compile(raw(path),str(path),'exec'),m.__dict__)
    return m


def tags(lp,limit,prior):
    """Rebuild the actual mother-degree birth-flux inequalities."""
    lo,hi=prior
    kernel={}
    for incoming,target,upper in product(range(2,limit+2),range(2,limit+1),(False,True)):
        disc=lo if upper else hi
        W={}
        for D in range(target,1,-1):
            for c in range(D,0,-1):
                b=Q(D-c,D)
                up=W.get((D+1,c),0)
                birth=b*((D==target)+W.get((D,c+1),0))
                rates=(D-c,D-c+Q(c-1,2)+Q(1,incoming))
                vals=[(birth+r*up)/(disc+b+r) for r in rates]
                w=Q((ceil(max(vals)*GRID) if upper else floor(min(vals)*GRID)),GRID)
                assert (w>=max(vals) if upper else w<=min(vals)) and w>=0
                W[D,c]=w
        kernel[incoming,target,upper]=W[2,1]
    for i,j in product(range(2,limit+2),range(2,limit+1)):
        assert 0<=kernel[i,j,False]<=kernel[i,j,True]<=1/lo
    keys=list(product(range(2,limit+1),(False,True)))
    ids={k:len(lp.rows)+i for i,k in enumerate(keys)}
    lp.rows.update({('tag',k):i for k,i in ids.items()})
    for i,(kind,terms,m) in enumerate(lp.cols):
        if kind!='occ' or not m[2]: continue
        D,k,b,*_=m
        s=min(D,limit+1)
        more=tuple((ids[j,up],z) for j,up in keys if (z:=b*((s==j)-kernel[s,j,up])))
        lp.cols[i]=(kind,terms+more,m)
    for key in keys:
        lp.col(('birth-tag-slack',key),{('tag',key):1 if key[1] else -1})
    return kernel


def signed_price(lp,prior,objective,scale=10**6):
    """Exact per-core lower bound, including signed zero/negative objectives."""
    from scipy.optimize import linprog
    lo,hi=prior
    cost=[Q(objective(m))/scale if kind=='occ' else Q(0) for kind,_,m in lp.cols]
    A,rhs=lp.matrix()
    assert lp.rhs=={'birth':1}
    y=linprog(list(map(float,cost)),A_eq=A,b_eq=rhs,method='highs-ipm',
              options=dict(dual_feasibility_tolerance=1e-9,primal_feasibility_tolerance=1e-9))
    if not y.success: y=linprog(list(map(float,cost)),A_eq=A,b_eq=rhs,method='highs-ds')
    assert y.success,y.message
    yn=[round(v*GRID) for v in y.eqlin.marginals]
    defects={}
    for (kind,terms,_),c in zip(lp.cols,cost,strict=True):
        defect=sum(a*yn[i] for i,a in terms)-GRID*c
        defects[kind]=max(defects.get(kind,0),defect)
    assert all(sum(a*yn[i] for i,a in terms)<=GRID*c+defects[kind]
               for (kind,terms,_),c in zip(lp.cols,cost,strict=True))
    eta=scale*(yn[lp.rows['birth']]-sum(defects.values())/lo)/GRID
    return min(lo*eta,hi*eta)


def controls(audit):
    """Compare every tag entry to the separate earlier marked-birth implementation."""
    model=audit.load(R/'occupation-neon/pressure-memory/birth-law/replay.py')
    host=audit.load(G/'host/reader.py')
    docs=host.receipt()
    prior=tuple(Q(docs[s]['certified_exact']) for s in SIDES)
    model.LO,model.HI,model.SAT=*prior,7
    first=model.build(24,sided=True)
    second=model.build(24,sided=True)
    kernel=tags(first,12,prior)
    old=load(G/'lifetime/pressure/frontier/replay.py')
    assert old.tagged_births(second,12,*prior)==kernel
    assert (first.rows,first.cols,first.rhs)==(second.rows,second.cols,second.rhs)
    assert (len(first.rows),len(first.cols))==(623,5564)
    assert signed_price(first,prior,lambda m:-3*m[2])==-3*prior[1]
    assert signed_price(first,prior,lambda m:3*m[2])==3*prior[0]
    print('PASS',len(kernel),'physical tag roots; 623 rows / 5564 columns; signed birth-row controls',flush=True)


def singleton_controls(audit):
    life=load(G/'lifetime/check.py')
    host=audit.load(G/'host/reader.py')
    for side in SIDES:
        with host.open_host(side) as w, mmap(w.p.fd,0,access=ACCESS_READ) as mm:
            t=life.tail.Snapshot(mm,w.golden)
            count=0
            for D in (2,3,24,25,251,252,253,3999,4000,4001,8000,10**8):
                for c in sorted({1,2,min(D,12),min(D,22),min(D,23),min(D,30),min(D,71)}):
                    A,B=life.row(t,D,c,w.cap_f)
                    assert (A,B)==w.slack_ticks(D,c)
                    for rho in (Q(0),Q(c,3),Q(c,2)):
                        D0=D-c
                        recv=rho if w.center(D,c) is None else c*w.center(D,c)
                        change=w.f(D+1,c)-w.f(D,c)
                        promote=(w.f(D,c+1)-w.f(D,c)+w.f(2,1)+w.g((D,c+1),(2,1))) if D0 else 0
                        direct=w.cap_f+w.sign*((D0+recv)*change+Q(D0,D)*promote-w.lam*w.f(D,c))
                        assert A+B*rho==direct>=0
                        count+=1
            print(side,'sparse singleton charge controls',count,flush=True)


def live_roots(w,pressure):
    life=pressure.life
    with mmap(w.p.fd,0,access=ACCESS_READ) as mm:
        t=life.tail.Snapshot(mm,w.golden)
        W={}
        life.price(t,7999,30,w.prior[1],w.cap_f,False,(701,W))
        roots,states,corners=pressure.live(t,w.cap_f,W,700,12,5,3,2,w.prior[1],False)
    assert (states,corners)==(2800083,17262)
    assert set(roots)=={(0,0)}|{(D,c) for D in range(2,6) for c in range(2,D+1)}
    print(w.side,'singleton states/corners',states,corners,flush=True)
    return roots


UPPER_SOURCES = {
    'edge18':('edges/upper/ceiling8.log',r'^EDGE18 (\d+/\d+) ',S),
    'fy_small':('edges/upper/absorb15.log',r'^TOTAL FY floor (\d+/\d+) ',S),
    'fy_wide':('payer/upper/README.md',r'^\| `251/12` \| `35/6` \| floor \| `(\d+/\d+)`',1),
}


def double_tail(prior,cap):
    lo,hi=Q(565713486,10**9),Q(565736987,10**9)
    assert lo<=prior[0]<=prior[1]<=hi and cap==463
    def hit(upper):
        nxt=[Q(1)]*46
        disc=lo if upper else hi
        for D in range(251,1,-1):
            cur=[Q(int(upper))]*46
            for c in range(min(44,D),0,-1):
                b=Q(D-c,D)
                vals=[(rate*nxt[c]+b*cur[c+1])/(disc+rate+b)
                      for rate in (D-c,D-Q(c,2))]
                cur[c]=Q(ceil(max(vals)*GRID) if upper else floor(min(vals)*GRID),GRID)
            nxt=cur
        return nxt[1]
    a,b=hit(False),hit(True)
    q=1/(1+lo)
    value=cap*max(Q(0),a-(10*b+q**10/(1-q))/252)*a/hi
    assert Q('.3606267611')<value<Q('.3606267612')
    return value


def corrected_upper_h3(frozen,payload):
    import ast
    import re
    source=frozen('higher/cohort.py')
    wrapper=frozen('host/frontier/replay-upper-h3.py')
    assert sha256(source).hexdigest()=='b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c'
    assert sha256(wrapper).hexdigest()=='49a4deeb6f8c3cf215ea7d093e2a2ab4c34d2a5af241bda30e0a3e3cf23ad0fe'
    old=b'vg-(-GRID*affine[3][i]//1)'
    assert source.count(old)==1
    executed=source.replace(old,b'vg+(-GRID*affine[3][i]//1)')
    log=frozen('host/frontier/upper-h3.log').decode().splitlines()
    assert log[0]==f'corrected cohort source {sha256(source).hexdigest()} executed {sha256(executed).hexdigest()}'
    prefix=f'upper original Gold2 {payload} namespace'
    assert log[1].startswith(prefix)
    opts=ast.parse('namespace'+log[1][len(prefix):],mode='eval').body
    assert not opts.args and {k.arg:ast.literal_eval(k.value) for k in opts.keywords}==dict(
        side='upper',min_mother=80,mother=480,mother_stock=12,middle=2,middle_stock=2,
        young=12,young_stock=3,sibling=2,upstream=2,outside_threshold=True,exclude_strict=False)
    assert 'literal labelled H3 rows 42' in log and 'PASS corrected original Gold2 upper H3' in log
    [norm]=[m[1] for line in log if (m:=re.match(r'^order floor physical hmy (\d+/\d+) ',line))]
    [ticks]=[m[1] for line in log if (m:=re.fullmatch(r'corrected upper mother-high ticks/core (\d+/\d+)',line))]
    assert 0<Q(ticks)==S*Q(norm)
    return Q(ticks)


def corrected_upper_increments(frozen):
    import ast
    import re
    source=frozen('post-fy/upper-frontier/frontier.py')
    wrapper=frozen('host/frontier/replay-upper-increments.py')
    assert sha256(source).hexdigest()=='1e2071a42f1b39c52c647f6a3b0f0205b4f12309d107aa8874e4b22afd7a9592'
    assert sha256(wrapper).hexdigest()=='480e62f6d24e83faa61f5a9b8d0f25693e6f29d07c773bdbf75a16bc9c937d2a'
    assert sha256(frozen('post-fy/replay.py')).hexdigest()=='d8a9a8a65d6b652151c41de18fa7a0f4e8aec4a6ed8e89395412e4ddd5bc3c23'
    old=b'self.round(w/width)'
    assert source.count(old)==1
    executed=source.replace(old,b'self.round(d*w/width)')
    log=frozen('host/frontier/upper-increments.log').decode().splitlines()
    assert log[0]==f'corrected increment source {sha256(source).hexdigest()} executed {sha256(executed).hexdigest()}'
    assert 'PASS corrected original Gold2 upper H3 increments' in log
    if log[-1]!='PASS corrected original Gold2 upper H3 increments':
        assert log[-1].endswith(' exit=0')
    out={}
    for name,cut,stock,low,c_min,m_min,fixtures in (
        ('middle3',120,10,80,2,3,True),('stock',79,12,19,9,2,False),
    ):
        i=log.index('corrected upper increment '+name)+1
        if fixtures:
            assert log[i].startswith('literal physical checks ')
            i+=1
        label='original upper Gold2 '
        assert log[i].startswith(label)
        assert ast.literal_eval(log[i][len(label):])==dict(host=5,stock=4,receiver=cut,rstock=stock,
            donor=7,tstock=3,old=low,last=cut,mother_min_stock=c_min,middle_min=m_min,child_min=2,
            ceiling=False,fixtures=fixtures,reserve=False)
        norm,ticks=re.match(r'^normalized (\d+/\d+) \S+ ticks (\d+/\d+) ',log[i+1]).groups()
        assert 0<Q(ticks)==S*Q(norm)
        out[f'h3_mother_high_{name}']=Q(ticks)
    return out


def corrected_upper_threshold(frozen):
    import ast
    import re
    source=frozen('post-fy/replay.py')
    wrapper=frozen('diagnostic/upper-correction/replay-threshold.py')
    assert sha256(source).hexdigest()=='d8a9a8a65d6b652151c41de18fa7a0f4e8aec4a6ed8e89395412e4ddd5bc3c23'
    assert sha256(wrapper).hexdigest()=='125b7b106273d8ef4796bf0183a98b081ccb6206652d57cc11cda93ef81c4c56'
    shift=b'shift=self.round(w/width) if width else 0'
    old=(b"if a.kind=='threshold' and a.side=='lower':\n"
         b'                    assert total==Q(1084423030952180046412413909,1704710415481443450880000000000000)')
    assert source.count(shift)==source.count(old)==1
    executed=source.replace(shift,b'shift=self.round(d*w/width) if width else 0').replace(
        old,b'if True:\n                    print("corrected exact ticks/core", total*s.scale, flush=True)')
    log=frozen('diagnostic/upper-correction/threshold-upper.log').decode().splitlines()
    assert log[0]==f'threshold source {sha256(source).hexdigest()} executed {sha256(executed).hexdigest()}'
    assert log[1].startswith('original Gold2 ')
    assert ast.literal_eval(log[1][len('original Gold2 '):])==dict(side='upper',kind='threshold',
        ceiling=False,limit=12,wing=32,arm=4,host=5,stock=4,receiver=79,rstock=8,donor=18,
        tstock=3,old=19,young=18,reserve=False,fixtures=True)
    assert log[2].startswith('literal physical checks ') and log[-1].endswith(' exit=0')
    ticks=Q(re.fullmatch(r'corrected exact ticks/core (\d+/\d+)',log[3])[1])
    assert ticks==S*Q(re.match(r'^H3 normalized (\d+/\d+) ',log[4])[1])
    return ticks


def compact_sibling(audit):
    """Rebuild the live-prefix recursion using literal corner enumeration."""
    from collections import Counter
    for name,pin in COMPACT_PINS.items():
        assert sha256((G/name).read_bytes()).hexdigest()==pin
    assert sha256((CHECKS/'gold2-h3-sol.py').read_bytes()).hexdigest()=='08458a2bbb3d0915165ed63bb264328f1f3c7be787d71e632658a8f0a077179d'
    assert sha256((G/'post-fy/replay.py').read_bytes()).hexdigest()=='d8a9a8a65d6b652151c41de18fa7a0f4e8aec4a6ed8e89395412e4ddd5bc3c23'
    base=run_path(str(CHECKS/'gold2-h3-sol.py'))
    writer=base['writer']
    host=audit.load(G/'host/reader.py')
    dims={'h':(79,8),'m':(5,4),'y':(5,3),'s':(3,3),'u':(3,3)}
    parent={'m':'h','y':'m','s':'h','h':'u'}
    orders={'sibling':('hmys','hmsy','hsmy'),'upstream':('uhmy',),'box':('hmy',)}
    over=(0,0); newborn=(2,1); grid=2**20
    # Independently enumerate labelled paths and their chronological orders in
    # every seven-vertex increasing tree. Later births have one prefix parent.
    counts={k:0 for k in orders}
    for tail in product(*(range(i) for i in range(1,7))):
        par=(None,*tail)
        for y in range(1,7):
            m=par[y]; h=par[m] if m else None
            if not h: continue
            counts['box']+=1
            u=par[h]
            if u:
                assert u<h<m<y
                counts['upstream']+=1
            for s in range(h+1,7):
                if s==m or par[s]!=h: continue
                observed=''.join(k for _,k in sorted(zip((h,m,y,s),'hmys')))
                assert orders['sibling'].count(observed)==1
                counts['sibling']+=1
    assert all(counts.values())
    print('PASS compact registration partitions',counts,flush=True)

    def increment(B,z):
        return max(0,z-B)-max(0,-B)
    for B,x,y in product(range(-3,8),range(5),range(5)):
        assert 0<=increment(B,x)+increment(B,y)<=increment(B,x+y)
    with host.open_host('lower') as w, writer.Snapshot('lower') as snap:
        identity=lambda st:(st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns,st.st_ctime_ns)
        assert identity(snap.stat)==identity(w.p.stat)
        e=writer.Rows(snap,w.p,w.golden)
        assert e.sig==-1 and e.hcap==408 and snap.scale==S
        literal_counts=Counter()
        for dh,ds,du in product((19,38,79),(2,3,5),(2,3,5)):
            base['literal'](w.p,e,(None,0,1,2,3,2,2),(2,du,dh,2,2,ds,3),literal_counts)
        assert literal_counts['negative H3'] and literal_counts['rejected donor-rate mutations']
        print('PASS compact literal charged rows',dict(literal_counts),flush=True)
        lo,hi=w.prior
        den=lcm(hi.denominator,*range(2,81))
        disc=hi.numerator*(den//hi.denominator)
        cache={}
        def row(h,m,y):
            key=h,m,y
            if key not in cache:
                A,(bm,bh,by)=e.triple(m,h,y)
                B=A+max(0,bh)*Q(h[1]-1,2)+max(0,bm)*Q(m[1]-2,2)+max(0,by)*Q(y[1]-1,2)
                assert B<=408
                cache[key]=B,bh
            return cache[key]
        terminal=None
        result={}
        for kind,group in orders.items():
            prices=[]
            for order in group:
                next_layer={}
                for length in range(len(order),0,-1):
                    roles=order[:length]; complete=length==len(order)
                    if complete and kind=='sibling' and terminal is not None:
                        next_layer={tuple(dict(zip('hmys',p))[r] for r in roles):v for p,v in terminal.items()}
                        continue
                    ns=tuple(tuple(j for j,t in enumerate(roles) if parent.get(t)==r or parent.get(r)==t)
                             for r in roles)
                    born_at=None if complete else roles.index(parent[order[length]])
                    axes=[]
                    for j,r in enumerate(roles):
                        D,C=dims[r]; required=max(1,len(ns[j])+int(j==0))
                        axes.append(([over] if r in 'su' else [])+[(d,c) for d in range(D,1,-1)
                                     for c in range(min(d,C),required-1,-1)])
                    table={}; corners=0
                    for p in product(*axes):
                        reward=0
                        if complete:
                            marks=dict(zip(roles,p)); h,m,y=(marks[r] for r in 'hmy')
                            if 19<=h[0] and h[0]>y[0] and m in e.mids:
                                B,b=row(h,m,y)
                                if kind=='box': reward=max(0,-B)
                                else:
                                    r='s' if kind=='sibling' else 'u'; q=marks[r]
                                    x=Q(1,q[0]) if q!=over else Q(1,4) if b>=0 else 0
                                    z=max(0,b/2)-b*x
                                    assert z>=0
                                    reward=increment(B,z)
                        n=floor(grid*reward)*den; d=disc; ports=[]
                        for i,(D,C) in enumerate(p):
                            r=roles[i]
                            if not D:
                                assert r in 'us'
                                # Overflow mark is frozen. For the sole missing
                                # upstream registration, zero lower-bounds its
                                # actual nonnegative birth rate.
                                continue
                            inside=sum(den//p[j][0] for j in ns[i] if p[j]!=over)
                            width=(C-len(ns[i]))*(den//2)+sum(den//4 for j in ns[i] if p[j]==over)
                            gain=(D+1,C) if r not in 'us' or D<dims[r][0] else over
                            vg=table.get(p[:i]+(gain,)+p[i+1:],0)
                            rate=(D-C)*den+inside
                            n+=rate*vg; d+=rate
                            if width: ports.append((width,vg))
                            if D>C:
                                post=p[:i]+((D,C+1),)+p[i+1:]
                                vb=table.get(post,0)
                                if i==born_at: vb+=next_layer.get(post+(newborn,),0)
                                rate=(D-C)*(den//D)
                                n+=rate*vb; d+=rate
                        samples=[(n+sum(a*v for a,v in chosen),d+sum(a for a,_ in chosen))
                                 for bits in product((0,1),repeat=len(ports))
                                 for chosen in [tuple(t for t,bit in zip(ports,bits) if bit)]]
                        v=min(a//b for a,b in samples)
                        assert v>=0 and all(a-b*v>=0 for a,b in samples)
                        table[p]=v; corners+=len(samples)
                    if complete and kind=='sibling': terminal=table
                    next_layer=table
                    print('compact',kind,roles,len(table),'states',corners,'exact corners',flush=True)
                prices.append(lo*Q(next_layer[(newborn,)],grid))
            result[kind]=sum(prices)
            assert result[kind]==COMPACT_TICKS[kind],(kind,result[kind])
            print('PASS compact',kind,'exact ticks/core',result[kind],flush=True)
        snap.check_unchanged()
    return result


def replay(audit,receipt,side,run_roots=False):
    import re
    inputs=receipt['frozen_inputs']
    def frozen(name):
        data=(G/name).read_bytes()
        assert sha256(data).hexdigest()==inputs[name],name
        return data
    for name in inputs: frozen(name)
    host=audit.load(G/'host/reader.py')
    name=G/'host/conditioned24.json'
    assert receipt['host_sha256']==audit.PINS[name] and receipt['original_receipt_sha256']==host.RECEIPT
    assert receipt['original_source']==host.SOURCE and receipt['lp_git_blob']=='7dd40ac41e5d95ae41843da51de436fd45b781f1'
    profile=json.loads(audit.raw(name))['results'][side]
    expected=receipt['results'][side]
    docs=host.receipt(); original=docs[side]
    prior=tuple(Q(docs[t]['certified_exact']) for t in SIDES)
    old_f,old_edge,_=audit.credits(side)
    ledger={}
    if side=='lower':
        ledger.update(edge18=old_edge['edge18'],fy_small=old_edge['absorb12'],fy_wide=old_edge['old_payer'])
        report=frozen('../checks/gold2-lower-h3-alt-sol.md').decode()
        frozen('../checks/gold2-lower-h3-alt-sol.py')
        for key,band,value in (
            ('h3_threshold','19..79',Q(5481957648397036499283,13005908321239040000000)),
            ('h3_seam','80..100',Q(229762209993424478073686403,16257385401548800000000)),
        ):
            [amount]=re.findall(r'^\| Mother `'+re.escape(band)+r'[^|]+\| `(\d+/\d+)`',report,re.M)
            assert Q(amount)==value
            ledger[key]=value
        if 'base_receipt_sha256' in receipt:
            for name,pin in COMPACT_PINS.items():
                assert sha256(frozen(name)).hexdigest()==pin
            log=frozen(COMPACT+'compact.log').decode()
            assert log.rstrip().endswith(' exit=0')
            del ledger['h3_threshold']
            for kind in ('sibling','upstream'):
                [v]=re.findall(r'^exact '+kind+r' ticks/core (\d+/\d+) ',log,re.M)
                assert Q(v)==COMPACT_TICKS[kind]
                ledger[f'h3_live_{kind}']=Q(v)
    else:
        for key,(name,pat,scale) in UPPER_SOURCES.items():
            [norm]=re.findall(pat,frozen(name).decode(),re.M)
            ledger[key]=scale*Q(norm)
        ledger['h3_threshold']=corrected_upper_threshold(frozen)
        if 'h3_mother_high' in expected['ledger']:
            ledger['h3_mother_high']=corrected_upper_h3(frozen,docs[side+'_sha256'])
        if 'h3_mother_high_middle3' in expected['ledger'] or 'h3_mother_high_stock' in expected['ledger']:
            ledger.update(corrected_upper_increments(frozen))
        ledger['edge_double_tail']=double_tail(prior,original['cap_g'])
    assert ledger=={k:Q(v) for k,v in expected['ledger'].items()}
    roots={tuple(map(int,k.split(','))):Q(v) for k,v in expected['singleton_roots'].items()}
    root_doc=json.loads(frozen(f'lifetime/pressure/frontier/marginals/{side}-floor.json'))
    assert (root_doc['side'],root_doc['bound'])==(side,'floor')
    assert root_doc['gold2_receipt_sha256']==host.RECEIPT
    assert root_doc['gold2_payload_sha256']==docs[side+'_sha256']
    assert tuple(map(Q,root_doc['prior']))==prior and Q(root_doc['discount'])==prior[1]
    assert Q(root_doc['overflow_ticks'])==0 and (root_doc['states'],root_doc['corners'])==(2800083,17262)
    assert root_doc['options']==dict(mother=5,live_cut=700,live_arity=12,tip=3,donors=2,tag_limit=12)
    for path,pin in PINS.items():
        if path.is_relative_to(G):
            name=str(path.relative_to(G))
            assert root_doc['sources'][f'{G.name}/{name}']==inputs[name]==pin
    assert roots=={tuple(map(int,k.split(','))):Q(v) for k,v in root_doc['roots'].items()}
    assert set(roots)=={(0,0)}|{(D,c) for D in range(2,6) for c in range(2,D+1)}
    assert all(v>=0 and (v*GRID).denominator==1 for v in roots.values())
    with host.open_host(side) as w:
        if run_roots:
            print(side,'original inode authenticated; starting singleton replay',flush=True)
            assert roots==live_roots(w,load(G/'lifetime/pressure/replay.py'))
        cells={tuple(map(int,k.split(','))):tuple(map(int,v)) for k,v in profile['cells'].items()}
        caps=audit.cap_cells(w)
        assert len(cells)==307 and all(h==caps[k] and 0<=h<=w.higher_cap and u<=w.global_high
                                       for k,(_,u,h) in cells.items())
        prior=w.prior; d=w.sign
        assert tuple(map(Q,profile['prior']))==prior==tuple(map(Q,expected['original_prior']))
        assert Q(profile['nominal'])==w.lam==Q(expected['nominal'])
        assert profile['payload_sha256']==w.docs[side+'_sha256'] and profile['cap_error']==w.golden['err_per_core']
        assert expected['payload_sha256']==w.docs[side+'_sha256']
        assert expected['higher_cap']==w.higher_cap and expected['h3_family_cap']==original['verify']['caps']['h3']
        assert expected['original_upper_denominator']==w.global_high
        model=audit.load(R/'occupation-neon/pressure-memory/birth-law/replay.py')
        def lp_for(bracket):
            model.LO,model.HI,model.SAT=*bracket,7
            obj=model.build(24,sided=True)
            tags(obj,12,bracket)
            assert (len(obj.rows),len(obj.cols))==(623,5564)
            return obj
        base=lp_for(prior)
        def reward(m):
            D,k,b,*_=m
            assert b==0 or D>5 or (D,k+2) in roots
            return b*roots.get((D,k+2),roots[0,0])
        F=signed_price(base,prior,reward)
        assert F==Q(root_doc['results']['control']['singleton_ticks'])
        E=w.golden['err_per_core']; B=sum(ledger.values())-E
        assert E==w.cap_f+max(0,w.golden['cap_g'])+w.higher_cap
        M=Q(profile['denominator'][1]); saved=Q(profile['saved_plus_reduced'])
        assert F>old_f and B+F>0 and 0<M<=w.global_high and 0<=saved<=w.higher_cap
        def trial(z,restricted=False):
            lo,hi=prior
            assert lo<=z<=hi
            bracket=((lo,z) if d==1 else (z,hi)) if restricted else prior
            a=d*(z-w.lam)
            assert a>0
            def obj(m):
                _,u,h=cells[m[:2]]
                return reward(m)+w.higher_cap-h-a*u
            return B+signed_price(lp_for(bracket) if restricted else base,bracket,obj)
        scalar=w.lam+d*(B+F)/w.global_high
        separate=w.lam+d*(B+F+saved)/M
        margin=trial(separate)
        assert margin>0
        uniform=separate+d*margin/M
        fixed=d*Q(ceil(d*uniform*10**12)+1,10**12)
        fixed_margin=trial(fixed,True)
        endpoint=fixed if fixed_margin>0 else uniform
        target=Q(receipt['target']); assert target==Q(241,426)
        target_margin=trial(target,True)
        exact=dict(singleton=F,replaced_singleton=old_f,error=E,net=B+F,d24_upper_denominator=M,
                   standalone_host_savings=saved,scalar_endpoint=scalar,separate_endpoint=separate,
                   uniform_margin=margin,uniform_endpoint=uniform,fixed_point=fixed,fixed_margin=fixed_margin,
                   final_endpoint=endpoint,target_margin=target_margin,inward_target_gap=d*(target-endpoint))
        exact.update(target_uniform_margin=trial(target),scalar_target_margin=B+F-d*(target-w.lam)*w.global_high)
        if side=='lower' and 'lower_99_175_margin' in expected:
            exact['lower_99_175_margin']=trial(Q(99,175),True)
            if 'base_receipt_sha256' in receipt:
                assert exact['lower_99_175_margin']>0 and endpoint>Q(99,175)
                print('PASS revised restricted 99/175 margin',exact['lower_99_175_margin'],flush=True)
        assert all(Q(expected[k])==v for k,v in exact.items())
        assert fixed_margin>0 and target_margin<0
        assert expected['outward']==audit.outward(endpoint,side)
        assert d*(endpoint-prior[0 if side=='lower' else 1])>0
        print(side,'endpoint',endpoint,expected['outward'],'exact fixed margin',fixed_margin,
              'restricted target margin',target_margin,flush=True)
        return endpoint


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--singleton-controls',action='store_true')
    cli.add_argument('--compact-sibling',action='store_true')
    cli.add_argument('--receipt-sha256',help='expected final frontier receipt SHA-256')
    cli.add_argument('--live-sibling-receipt',action='store_true')
    cli.add_argument('--side',choices=SIDES,action='append')
    cli.add_argument('--replay-singleton-roots',action='store_true',help='expensive full Bellman replay')
    args=cli.parse_args()
    assert __debug__
    for p in PINS: raw(p)
    audit=load(CHECKS/'gold2-host-h3-sol.py')
    if args.singleton_controls or not (args.receipt_sha256 or args.compact_sibling): controls(audit)
    if args.singleton_controls: singleton_controls(audit)
    if args.compact_sibling: compact_sibling(audit)
    if args.receipt_sha256:
        path=FRONTIER/('receipt-live-sibling.json' if args.live_sibling_receipt else 'receipt.json')
        data=path.read_bytes()
        assert sha256(data).hexdigest()==args.receipt_sha256
        receipt=json.loads(data)
        if args.live_sibling_receipt:
            base=(FRONTIER/'receipt.json').read_bytes()
            assert sha256(base).hexdigest()==receipt['base_receipt_sha256']
            assert receipt['results']['upper']==json.loads(base)['results']['upper']
        results={side:replay(audit,receipt,side,args.replay_singleton_roots) for side in args.side or SIDES}
        if set(results)==set(SIDES):
            assert results['lower']<Q(241,426)<results['upper'] and receipt['target_inside']
            print('PASS exact both-sign frontier receipt',args.receipt_sha256)
