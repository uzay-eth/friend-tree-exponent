#!/usr/bin/env python3
"""Cold receipt, exact physical-birth LP dual, and replacement arithmetic."""
import argparse
from collections import defaultdict
import ctypes as C
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
import json
from math import ceil, floor, lcm
from pathlib import Path
import random
import re
import subprocess

R2 = Path(__file__).resolve().parents[1]
PM = R2/'occupation-neon/pressure-memory'
FAST = PM/'path-fast'
PARTS = ('minimum','mother','0','1','2','3')
LO,HI = Q('.565673'),Q('.565775')
G = 10**18
PINS = {
    'path-fast/solver.cc':'af0df8e91f6ea7575ef77a050d8bf55be9745277105c3f75b19227d19ae8abfb',
    'path-fast/charges.hh':'f9edec256b6ef4cd07d4a222110265facfa7c8580b49b3b4922dbaa041b281e0',
    'path-fast/replay.py':'ef31abbde283d6b2bcc9dafbf31473c8bf87ea7e07340fce3dec787deafbf8a2',
    'path-fast/pool.py':'553fa8b1209a7c0110eab5f360a0c97c13455480765a19df95d1d3ae328bbf66',
    'path-fast/medium-entries.json':'30c8466ce23478b0a135b4838ab098bc36dbb85b66a84ed382b64acd55a6f76a',
    'path-fast/wide-head-entries.json':'ca4d2be94bf6f83275ab38938857b6303e6ebe99f4d83b4e63570ea480e5523e',
    'path-fast/wide-tail-entries.json':'fa4caa147b8f8e45e1b85f5a368267437eaec522f07b39019e0ec57e5e191a6a',
    'path-fast/wide-mother12-entries.json':'4b992c7a272d89e6445fa9e41b17b89ed7fe781ce98fced3117c1ef2057cda48',
    'donor_lift.py':'b9423487db066349cdb8bbd7a64380041c65fe2fae14d05907a4932775bd195b',
    'birth-law/replay.py':'92a99ef12d8567798c4ed4dc3beb6535d3f18ef3ad5028aff0ec84e3684a5ba8',
    'p4-donor-medium.log':'5302e1414a632177f1254001854adc6c171aeb2a0df559cd4924d67b01730994',
    'hv-lifetime/release/README.md':'b420ec4f03f8b314ac25873ad94fc68cccc5f298592a38d829f5a9b209870270',
    'combined/README.md':'5faa5879b10043993ff081506aad898ce02e62e33081ba2d1cc9732a3966ee31',
    'path_band.py':'3bb442ee42b93fa431ae00e8b9cbf0128fbb013c6f2cdca143a6ab49d3c2f642',
    'path.py':'28ce3034eb096c6ebc53c339a2ee6f2375cd72e26cd85b333e76072b299eaaed',
    '../../occupation-neon/higher/extract.py':'9f40a50a3623756b277243d192f672e76c45b4b9784715ed4d0d8fb27d36131c',
    '../../attacks/neon1/audit.py':'d548a9aafa899ecab2b00a7f6ad0a4e6f828600088780f0a15f08927f17f0bfe',
    '../../checks/stationary-p4-sol.py':'7432857b1b9fc05b7850cb9d2a6f22fc8200d22dc77f134bbd75eb69da2f0653',
}


def module(name,file):
    s=spec_from_file_location(name,file)
    m=module_from_spec(s); s.loader.exec_module(m)
    return m


def entry(name, dims, ars, cap, parts):
    f=FAST/(name+'-entries.json'); j=json.loads(f.read_text())
    assert (j['degrees'],j['arities'],j['side_cap'],j['mother_cap'],j['include_old_band'])==(list(dims),list(ars),3,cap,False)
    assert j['grid']==G and tuple(map(Q,j['discount']))==(LO,HI)
    assert j['channels']==['lower_lo','lower_hi','upper_lo','upper_hi']
    rows=j['conditional']; assert tuple(rows)==parts
    for r in rows.values():
        assert set(r)=={'overflow',*map(str,range(2,cap+1))}
        assert all(len(z)==4 and all(type(x) is int and 0<=x<2**60 for x in z)
                   and z[0]<=z[1] and z[2]<=z[3] for z in r.values())
    return rows


def component(r):
    return tuple((min(x[j] for x in r.values())*LO/G,max(x[j+1] for x in r.values())*HI/G) for j in (0,2))


def totals(rows):
    return tuple(tuple(sum(component(r)[i][j] for r in rows.values()) for j in (0,1)) for i in (0,1))


def medium(rows):
    lines=(PM/'p4-donor-medium.log').read_text().splitlines()
    assert len(lines)==22
    meta={}
    for i,p in enumerate(PARTS):
        head=lines[3*i+1].split(); assert head[:3]==['component',p,'censuses']
        meta[p]=tuple(int(head[j]) for j in (3,6,8))
        assert meta[p][2]==4*meta[p][1]
        stored=tuple(tuple(map(Q,s.split()[1:3])) for s in lines[3*i+2:3*i+4])
        assert component(rows[p])==stored
    print('PASS archived Python / native medium: all 24 endpoints;',sum(x[1] for x in meta.values()),'states,',sum(x[2] for x in meta.values()),'drift checks')
    return meta


def fresh(root, meta):
    for j,k in (('medium-native','medium'),('medium-ubsan','medium'),('wide-head','wide-head'),('wide-tail','wide-tail'),('wide-mother12','wide-mother12')):
        data=(root/(j+'.json')).read_bytes(); a=json.loads(data); a.pop('library_sha256')
        b=json.loads((FAST/(k+'-entries.json')).read_text()); b.pop('library_sha256')
        assert a==b, j
        print('PASS fresh entry replay',j,sha256(data).hexdigest())
    for line in (root/'medium-ubsan.log').read_text().splitlines():
        x=re.match(r'component (\S+) censuses (\d+) live states (\d+).* residuals (\d+) slots (\d+)',line)
        if x:
            p,a,b,c,d=x.groups(); assert tuple(map(int,(a,b,c)))==meta.pop(p) and b==d
    assert not meta
    for name,expected in (('wide-head',76893297),('wide-tail',64145862),('wide-mother12',25476804)):
        log=(root/(name+'.log')).read_bytes()
        rows=re.findall(rb' live states (\d+).* residuals (\d+) slots (\d+)',log)
        assert sum(int(a) for a,_,_ in rows)==expected and all(4*int(a)==int(b) and a==c for a,b,c in rows)
        assert len(rows)=={'wide-head':4,'wide-tail':2,'wide-mother12':1}[name]
        print('PASS fresh native totals',name,expected,4*expected,'log SHA-256',sha256(log).hexdigest())


def kernel():
    result=[]
    for m in range(2,14):
        # Compute all destination classes at once, using independent use of
        # the lower gain corner on denominators and zero beyond degree 34.
        values={}
        zero=(0,)*12
        for d in range(34,1,-1):
            for k in range(d-1,-1,-1):
                l=d-k-1; b=Q(l,d)
                grow=values.get((d+1,k),zero)
                prom=values.get((d,k+1),zero)
                values[d,k]=tuple(floor(min((G0*b*(j==min(d,13))+b*p+r*g)/(HI+b+r)
                                            for r in (Q(l),Q(l)+Q(k,2)+Q(1,m))))
                                  for j,g,p in zip(range(2,14),grow,prom))
        result.append(tuple(Q(v,G0) for v in values[2,0]))
    assert all(0<=sum(row)<=1 for row in result)
    return tuple(result)


G0=10**16


def dual(lp, c):
    from scipy.optimize import linprog
    A,rhs=lp.matrix()
    options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9}
    p=linprog(list(map(float,c)),A_eq=A,b_eq=rhs,method='highs-ipm',options=options)
    if not p.success:
        p=linprog(list(map(float,c)),A_eq=A,b_eq=rhs,method='highs-ds')
    assert p.success,p.message
    y=tuple(Q(round(x*10**12),10**12) for x in p.eqlin.marginals)
    defects=defaultdict(Q)
    products=[]
    for z,(kind,vec,_) in zip(c,lp.cols):
        gap=sum(a*y[i] for i,a in vec)-z
        defects[kind]=max(defects[kind],gap); products.append((kind,gap))
    assert all(gap<=defects[k] for k,gap in products)
    assert set(defects)=={'occ','edge','stubs','tail-up','palm-slack'}
    return sum(a*b for a,b in zip(rhs,y))-sum(defects.values())/LO,sum(defects.values())/LO


def lift_controls(lp):
    """Contract physical census/edge corners at the degree and arity seams."""
    by_name={idx:name for name,idx in lp.rows.items()}
    edge={desc:dict((by_name[i],q) for i,q in vec) for kind,vec,desc in lp.cols if kind=='edge'}
    def weights(a,b,x):
        assert a<=x<=b
        return ((a,Q(1)),) if a==b else ((a,(b-x)/(b-a)),(b,(x-a)/(b-a)))
    count=0
    for d,k in ((2,0),(3,1),(4,1),(12,4),(13,5),(18,8),(19,6),(20,7),(1000,23)):
        p,b=Q(1,d),Q(d-k-1,d); s=Q('.56572'); u=Q(k,3); v=Q(2,7)
        cut=(d,k) if d<=18 else (19,min(k,7)); a,f=cut
        sat=a==19 and f==7
        if sat:
            bp=[(bb,pp,w*z) for pp,w in weights(0,Q(1,19),p)
                for bb,z in weights(0,1-8*pp,b/(1-8*p)*(1-8*pp))]
        elif a==19:
            bp=[(1-(k+1)*pp,pp,w) for pp,w in weights(0,Q(1,19),p)]
        else:
            bp=[(b,p,Q(1))]
        out=defaultdict(Q)
        for (bb,pp,w),(ss,w2),(uu,w3),(vv,w4) in product(bp,weights(LO,HI,s),
                 weights(0,0 if sat else Q(k,2),0 if sat else u),weights(0,Q(1,2),v)):
            # The discount is not in the column metadata; recover its exact
            # corner from the independent mass row b−s.
            for name,x in cols_for(lp,by_name,(a,f,bb,pp,uu,vv),ss):
                out[name]+=w*w2*w3*w4*x
        cell=(d if d<=12 else 13 if d<=18 else 19,min(k,3))
        if sat:
            out['arity']+=k-7; out['start-n',cell]+=k-7
            out['rho']+=u; out['up',cell]+=u
        assert (out['birth'],out['mass'],out['arity'],out['rho'])==(b,b-s,k-1,u+v+s-1)
        assert (out['end-n',cell],out['end-p',cell],out['start-n',cell],out['start-p',cell],out['up',cell],out['down',cell])==(1,p,k,k*p,u,v)
        expect=defaultdict(Q,{cut:s,(2,0):-s}) if cut!=(2,0) else defaultdict(Q)
        if d<=18:
            g=d-k-1+u+v
            expect[cut]+=g; expect[d+1,min(k,7) if d==18 else k]-=g
        if b and not sat: expect[cut]+=b; expect[a,f+1]-=b
        assert {x:q for (kind,x),q in ((name,q) for name,q in out.items() if isinstance(name,tuple)) if kind=='age' and q}=={x:q for x,q in expect.items() if q}
        count+=1
    a,b=(13,3),(19,1); p,q=Q(1,15),Q(1,97)
    lifted=defaultdict(Q)
    for (x,w),(y,z) in product(weights(Q(1,18),Q(1,13),p),weights(0,Q(1,19),q)):
        for key,val in edge[a,b,x,y].items(): lifted[key]+=w*z*val
    assert tuple(lifted[k] for k in [('start-n',a),('start-p',a),('end-n',b),('end-p',b),('up',a),('down',b)])==(-1,-p,-1,-q,-q,-p)
    print('PASS exact physical LP column lifts:',count,'censuses and one mixed/overflow edge')


def cols_for(lp,names,desc,s):
    for kind,vec,m in lp.cols:
        if kind!='occ' or m!=desc: continue
        row=tuple((names[i],q) for i,q in vec)
        if dict(row).get('mass',Q())==desc[2]-s:
            return row
    raise AssertionError((desc,s))


def pool(lp,rows,cap):
    boxes=[]
    for j in (0,2):
        result=[]
        for lane,sign in ((0,1),(1,-1)):
            c=[]
            for kind,_,desc in lp.cols:
                if kind!='occ': c.append(Q()); continue
                d,_,b,_,_,_=desc
                c.append(sign*b*Q(rows[d if d<=cap else 0][j+lane],10**12))
            bound,repair=dual(lp,c)
            result.append(sign*bound/10**6)
            print('exact dual repair, cap/sign/lane',cap,j,lane,repair)
        corners=tuple(s*v for s,v in product((LO,HI),result))
        a=max(LO*min(r[j] for r in rows.values())/G,min(corners))
        b=min(HI*max(r[j+1] for r in rows.values())/G,max(corners))
        assert 0<a<=b
        boxes.append((a,b))
    return tuple(boxes)


def join(rows,refine=None):
    cap=12 if refine else 3
    joined={}
    for d in (0,*range(2,cap+1)):
        key='overflow' if d==0 else str(d)
        mother=refine[key] if refine else rows['mother'][key]
        joined[d]=tuple(mother[j]+sum(rows[p][key if d in (2,3) else 'overflow'][j] for p in PARTS if p!='mother') for j in range(4))
    return joined


def compose(price,old):
    frozen=(Q(6853561875850472647665032629043939443,12115632305640977999025000000000000000),
            Q(328740990573504279104102687021150749,581051508956477588875000000000000000))
    assert all(str(q) in (PM/'hv-lifetime/release/README.md').read_text() for q in frozen)
    den=(Q(1357493815758092773,10**18),Q(273436004214812983,2*10**17))
    assert old==(Q(359297001126132137,10**24),Q(6963149521191559,25000000000000000000000))
    end=tuple(v+eps*(new-was)/d for v,eps,new,was,d in zip(frozen,(1,-1),price,old,den))
    assert frozen[0]<end[0]<end[1]<frozen[1]
    ticks=floor(end[0]*10**15),ceil(end[1]*10**15)
    print('replacement increments:',*(a-b for a,b in zip(price,old)))
    print('exact f7cc2e1 replacement:',*end,'outward:',*(f'0.{v:015d}' for v in ticks))


BRIDGE = r'''
#include "solver.cc"
extern "C" int point(const Tables* t,const int* d,const int* k,int old,I* out) {
    Reward r; std::array<Census,4> s;
    for (int i=0;i<4;++i) s[i]={d[i],k[i]};
    try {
        int status = Charges(t,old)(s,r);
        if (status==1) {
            std::copy_n(r.minimum.data(),4,out);
            for(int i=0;i<5;++i) {
                std::copy_n(r.source[i].data(),4,out+4+8*i);
                std::copy_n(r.slope[i].data(),4,out+8+8*i);
            }
        }
        return status;
    } catch (...) { return 0; }
}
extern "C" int step(const Tables* t,const int* d,const int* k,int n,int part,int dy) {
    static const int db[]={16,9,9,16}, kb[]={4,4,4,3};
    static Solver solver(db,kb,3,12,t,true,nullptr);
    State s; s.n=n; s.y={dy,0};
    for (int i=0;i<n;++i) s.path[i]={d[i],k[i]};
    try { solver.part=part; return solver.live(s); }
    catch (...) { return -1; }
}
extern "C" I rounded(I a,I b,int upper) { return ratio(W(a),b,upper!=0); }
'''


def wire_check(blobs,root):
    # The bridge includes the exact native definitions without editing producer.
    # A single Solver instance is reused in the exit probe (part is set per call).
    so=root/'path-fast-wire.so'
    subprocess.run(['g++','-std=c++20','-O2','-Wall','-Wextra','-Werror','-fsanitize=undefined',
                    '-fno-sanitize-recover=all','-shared','-fPIC','-I',str(FAST),'-x','c++','-', '-o',str(so)],input=BRIDGE,text=True,check=True)
    n=module('native',FAST/'replay.py'); tables=tuple(n.original(b,s) for b,s in zip(blobs,('lower','upper')))
    U=lcm(10**9,*range(2,18))
    for t in tables:
        peaks=[max(abs(int(a.min())),abs(int(a.max()))) for a in (t.c.p4tab,t.c.hvtab)]
        # 4 gains at ≤16U, 4 promotions at ≤U, 2 transfers at ≤U,
        # and 16 directional floor copies at U/2; ∣delta∣≤2*peak.
        peak=max(peaks); bound=(1+2*(64+4+2+8))*peak*U*10**7
        assert peak<2**50 and bound<2**127
        print('authenticated full P4/hV tick peaks and charge bound bits:',*peaks,bound.bit_length(),flush=True)
    tab,indices=n.native_tables(tables)
    lib=C.CDLL(str(so)); dim=C.c_int*4
    lib.point.argtypes=(C.POINTER(n.Tables),n.INDEX,n.INDEX,C.c_int,C.POINTER(n.I))
    lib.step.argtypes=(C.POINTER(n.Tables),n.INDEX,n.INDEX,C.c_int,C.c_int,C.c_int)
    lib.rounded.argtypes=(n.I,n.I,C.c_int); lib.rounded.restype=n.I
    rng=random.Random(1609916)
    counts=defaultdict(int)
    ds=(16,9,9,16); ks=(4,4,4,3)
    samples=[tuple(zip(d,k)) for d,k in product(product((2,12,13,16),(3,8,9),(3,6,7,9),(2,15,16)),
                                               ((1,1,1,0),(1,3,1,0),(2,2,2,1),(4,4,4,3)))
             if all(x>y for x,y in zip(d,k))]
    for _ in range(2500):
        d=[rng.randrange(2,m+1) for m in ds]
        samples.append(tuple((m,rng.randrange(i<3,min(m,ks[i]+1))) for i,m in enumerate(d)))
    for path in samples:
        d,k=(dim(*z) for z in zip(*path))
        x=tuple((a-b-1,b,1) for a,b in path)
        off=sum(k)-3; size=sum(d)-sum(k)+2*off
        old=d[0]<=12 and d[1]<=8 and d[2]<=6 and d[3]<=16 and k[0]==1 and k[1]<=3 and k[2]==1 and k[3]==0
        include=(size>20 or size>12 and off>0) and not old and all(t.supported('P4',x) for t in tables)
        out=(n.I*44)()
        assert lib.point(tab,d,k,1,out)==(1 if include else -1),path
        counts['accepted' if include else 'rejected']+=1
        if include:
            mins=[]; rest=[[],[],[],[],[]]
            for t in tables:
                a,free=t.row('P4',x); a/=t.c.S; slope=[Q(0)]*5
                for i,parent,m,w in free:
                    j=0 if parent else i+1
                    assert m==(1 if parent else k[i]-int(i<3))
                    slope[j]=Q(w,t.c.S); a+=m*min(Q(0),slope[j])/2
                mins.extend((floor(a*G),ceil(a*G)))
                for j,z in enumerate(slope):
                    h=max(Q(0),-z)*G/2
                    rest[j].append((floor(h),ceil(h),floor(z*G),ceil(z*G)))
            want=mins+[v for (l,r) in rest for v in (*l[:2],*r[:2],*l[2:],*r[2:])]
            assert list(out)==want,(path,[(i,a,b) for i,(a,b) in enumerate(zip(out,want)) if a!=b])
        for j in range(4):
            e=list(d); e[j]=ds[j]+1
            assert lib.step(tab,dim(*e),k,4,-2,-1)==0
            e=list(k); e[j]=ks[j]+1
            assert lib.step(tab,d,dim(*e),4,-2,-1)==0
            counts['irreversible exits']+=2
    a=dim(5,2,2,2)
    for role,k,q in product(range(4),range(5),(-1,0,2)):
        if role==0 and q>=0 and k==0:
            continue
        demand=1+(role==0 and q<0)
        assert lib.step(tab,a,dim(k,0,0,0),1,role,q)==(k+demand<=4)
        counts['short-prefix demand']+=1
    for a,b,up in product((-2**60+1,-27,-3,-1,0,1,3,27,2**60-1),(1,2,3,19),(0,1)):
        assert lib.rounded(a,b,up)==(ceil if up else floor)(Q(a,b))
    assert lcm(10**6,*range(2,17))==9009000000
    print('PASS independent wide C++/Fraction wire and UBSan exit/rounding fixtures:',dict(counts),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--fresh',type=Path)
    p.add_argument('--wire',type=Path,nargs=2,metavar=('LOWER','UPPER')); p.add_argument('--scratch',type=Path)
    a=p.parse_args()
    assert __debug__
    for file,pin in PINS.items(): assert sha256((PM/file).read_bytes()).hexdigest()==pin,file
    if a.wire:
        assert a.scratch and a.scratch.is_dir()
        wire_check(a.wire,a.scratch)
    med=entry('medium',(9,7,6,10),(3,3,3,2),3,PARTS)
    meta=medium(med)
    dims,ars=(16,9,9,16),(4,4,4,3)
    head=entry('wide-head',dims,ars,3,PARTS[:4]); tail=entry('wide-tail',dims,ars,3,PARTS[4:]); rows=head|tail
    fine=entry('wide-mother12',dims,ars,12,('mother',))['mother']
    if a.fresh: fresh(a.fresh,meta)
    # The finer control is nested at both exact low bins and every former
    # overflow sub-bin. This also catches copying cap-3 instead of cap-12 mother.
    for key,z in fine.items():
        c=rows['mother'][key if key in ('2','3') else 'overflow']
        assert all(c[j]<=z[j]<=z[j+1]<=c[j+1] for j in (0,2)),key
    shared,refined=join(rows),join(rows,fine)
    assert all(refined[d][j]-shared[d if d in (2,3) else 0][j]==fine[str(d) if d else 'overflow'][j]-rows['mother'][str(d) if d in (2,3) else 'overflow'][j]
               for d in refined for j in range(4))
    print('wide exact uniform component sums:',totals(rows))
    bl=module('bl',PM/'birth-law/replay.py')
    source=module('birth_ref',R2/'checks/stationary-p4-sol.py')
    k=kernel(); assert k==source.class_kernel()
    lp=bl.build(18,True,k); assert (len(lp.rows),len(lp.cols))==(462,3820)
    lift_controls(lp)
    old=pool(lp,shared,3); new=pool(lp,refined,12)
    print('exact shared cap 3:',old); print('exact shared mother-only cap 12:',new)
    assert all(b[0]>=a[0] for a,b in zip(old,new))
    compose(tuple(v[0] for v in new),tuple(v[0] for v in totals(med)))


if __name__=='__main__': main()
