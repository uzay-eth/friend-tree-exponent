#!/usr/bin/env python3
"""Exact independent hV tail, physical atom, and installment checks."""
import argparse
from ast import literal_eval
from fractions import Fraction as Q
from functools import cache
from hashlib import file_digest, sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations, product
from math import ceil, comb, floor, prod
from pathlib import Path
import random
import re

R2 = Path(__file__).resolve().parents[1]
BASE = 'occupation-neon/pressure-memory/hv-lifetime'
PINS = {
    BASE+'/old/replay.py': 'bc92b36d485a4f3e8249dd670b017f1385296242f22a474fa9ced2de92110b35',
    BASE+'/replay.py': '0c4d5824c88005eba91184a38df28920527f4058ee4a0aa732fd368fbe49666a',
    'occupation-neon/higher/older/replay.py': '2e600622e108be8791912577adcbb9291e17e3412763abcf5e7ca5c03130fb2f',
    'occupation-neon/higher/extract.py': '9f40a50a3623756b277243d192f672e76c45b4b9784715ed4d0d8fb27d36131c',
    'attacks/neon1/audit.py': 'd548a9aafa899ecab2b00a7f6ad0a4e6f828600088780f0a15f08927f17f0bfe',
}
HASHES = ('9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
          '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41')
ROLES = (0,0,1,2)


def load(path):
    spec = spec_from_file_location('hv_old', R2/path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def embed(c, x):
    d = sum(x)
    if d < 18:
        return ((x,Q(1)),)
    k,h = min(x[1],18),x[2]
    z = Q(18,d)
    keys = (('inf',), *((f'T{i}',k,h) for i in (1,2,3)))
    return tuple((j,comb(3,i)*z**i*(1-z)**(3-i)) for i,j in enumerate(keys))


def contract(c, x, entry):
    return sum((prod(q for _,q in r)*int(entry(*(v for v,_ in r)))
                for r in product(*(embed(c,v) for v in x)) if all(q for _,q in r)),Q())


def hv(c, x):
    if x[1][1] > 17:
        return Q(0)
    return contract(c,x,lambda a,b,d: c.hvtab[c.par[a],c.hv[b],c.child[d]] if b in c.hv else 0)


def gain(x,i):
    l,k,h = x[i]
    return x[:i]+((l+1,k,h),)+x[i+1:]


def flags(c,x):
    a,b,d = x
    return a in c.p4 and b in c.p4, a in c.ybu and b in c.ybw, b in c.yt, b in c.p4 and d in c.p4


def classifier(old, tables):
    rng = random.Random(18947)
    targets = []
    for da,db,dc in product((9,10,17,18,19,48,10**30+7),(9,10,17,18,10**25+3),(9,17,18,10**40+11)):
        targets.append(((da,2),(db,2),(dc,1)))
    for k in (16,17,18,19,50):
        targets += [((10**60+3,k),(10**35+3,2),(10**55+7,k)),
                    ((10**35+3,2),(10**45+1,k),(18,1))]
    targets.extend((((2,1),(2,1),(2,0)),((17,1),(18,17),(17,0))))
    for _ in range(100):
        ds = [rng.randrange(2,1000) for _ in range(3)]
        targets.append(tuple((d,rng.randrange(int(i<2),min(d,24))) for i,d in enumerate(ds)))
    counts = {}
    for tab in tables:
        c = tab.c
        for path in targets:
            x = tuple((d-k-1,k,1) for d,k in path)
            a = hv(c,x)
            diffs = [hv(c,gain(x,i))-a for i in range(3)]
            status, slopes = [],[]
            for i, n, routed in zip(ROLES,(1,path[0][1]-1,path[1][1]-1,path[2][1]),flags(c,x)):
                z = tab.sign*diffs[i]/c.S
                label = ('dead' if path[1][1]>17 else 'absent' if not n else 'routed' if routed else
                         'positive' if z>0 else 'negative' if z<0 else 'zero')
                status.append(label)
                slopes.append(z if label in ('negative','positive') else Q())
                counts[label] = counts.get(label,0)+1
            assert tuple(status) == old.signs(tab,path), path
            assert tuple(slopes) == old.H.row(tab,path)[1:], path
    assert set(counts) >= {'absent','routed','dead','positive','negative'}
    print(f'PASS rational Bernstein seam/tails: {len(targets)} censuses per witness; {counts}')


@cache
def adjacency(t):
    a = [[] for _ in t]
    for x,p in enumerate(t[1:],1):
        a[x].append(p)
        a[p].append(x)
    return tuple(tuple(xs) for xs in a)


def state(t,i):
    a = adjacency(t)
    return (sum(len(a[x])==1 for x in a[i]),sum(t[x]==i and len(a[x])>1 for x in a[i]),int(i!=0))


def edge_transfer(c, x, i):
    def g(states):
        return contract(c,states,lambda a,b:c.g[c.state[a],c.state[b]])
    l,k,h = x[i]
    d = l+k+h
    if d < 18:
        return g(gain(x,i))-g(x)
    if k>17:
        return Q(0)
    inf, t1 = c.state['inf',],c.state['T1',k,h]
    z = Q(18,d)
    out = Q()
    for v,w in embed(c,x[1-i]):
        ix = c.state[v]
        if i:
            gap = int(c.g[ix,inf])-int(c.g[ix,t1]); ref = int(c.yw[c.par[v],k])
        else:
            gap = int(c.g[inf,ix])-int(c.g[t1,ix]); ref = int(c.yp[(1-h)*19+k,c.child[v]])
        out += w*(z*z*(1-z)*gap/6 + z**3*ref)
    return out


def center(c,x,par):
    _,k,h = x
    if sum(x)<18:
        return Q(int((c.rpar if par else c.rch)[c.state[x]]),c.Sr)
    return Q(int((c.trpar if par else c.trch)[(1-h)*19+min(k,18)]),c.Sr)


def physical_raw(tab,t,v,groups):
    """Arrow drift of the existing hV, plus four-core births and routed transfers."""
    c = tab.c
    a = adjacency(t)
    x = tuple(state(t,i) for i in v)
    w = hv(c,x)
    ds = [hv(c,gain(x,i))-w for i in range(3)]
    out = -c.lam*w
    for src, neighbors in enumerate(a):
        p = Q(1,len(neighbors))
        for dst in neighbors:
            if dst in v:
                out += p*ds[v.index(dst)]
            if src in v and len(a[dst])==1:
                i = v.index(src)
                u = t+(dst,)
                after = tuple(state(u,j) for j in (*v,dst))
                ch = [[1],[2],[],[]]; ch[i].append(3)
                # Only an entirely new four-body term is owned by this old hV.
                out += p*(hv(c,after[:3])-w+c.potential(ch,after,parts=True)[(6,5,4)[i]])
    for i, group, parent, routed in zip(ROLES,groups,(True,False,False,False),flags(c,x)):
        if routed:
            out += ds[i]*sum(center(c,x[i],parent)-Q(1,len(a[j])) for j in group)
    # g→hV transfer along either of the deleted tips.
    if x[1][1] <= 17:
        out += (Q(1,sum(x[0]))-center(c,x[1],True))*edge_transfer(c,x[1:],0)
        out += (Q(1,sum(x[2]))-center(c,x[1],False))*edge_transfer(c,x[:2],1)
    return tab.sign*out/c.S


def planted(t,a):
    incl = {a}
    for j in range(a+1,len(t)):
        if t[j] in incl:
            incl.add(j)
    return len(incl),sum(len(adjacency(t)[j])>1 for j in incl)


def make_tree(degrees, fan=(2,2,2), young=False):
    t = [-1,0,1,2,3,4,0,2,2]
    for i,(d,n) in enumerate(zip(degrees,fan),3):
        for j in range(n):
            y = len(t); t.append(i)
            child = len(t); t.extend([y]*(1 if young else 2+j))
            if not young:
                t.extend([child]*(1+(j+i)%3))
        t.extend([i]*(d-len(adjacency(tuple(t))[i])))
    return tuple(t)


def physical(old, tables):
    sizes, checked, mutations = 0,0,0
    for d in (10,17,18,20,21,48):
        for ds in ((d,4,4),(d,10,10),(d,18,19)):
            t = make_tree(ds)
            a = adjacency(t); v = (3,4,5)
            x = tuple(state(t,i) for i in v)
            path = tuple((sum(s),s[1]) for s in x)
            assert tuple(z[0] for z in path) == ds
            groups = ([2],*([j for j in a[i] if t[j]==i and len(a[j])>1 and j not in v] for i in v))
            assert all(len(a[j])>2 for group in groups for j in group)
            assert all(len({len(a[j]) for j in group})==len(group) for group in groups)
            sz, cores = planted(t,3)
            counts = sum(k for _,k in path)-2
            smin = sum(d-k for d,k in path)+2*counts
            assert sz >= smin and cores >= 3+counts
            assert d<21 or sz>=d>20
            for y in (None,(0,0), *( (len(a[j]),state(t,j)[1]) for j in sum((list(g) for g in groups[1:]),[]) )):
                plus, more = (y[0]+y[1]-2,y[1]) if y and y[0] else (0,0)
                s,c = smin+plus,3+counts+more
                assert sz>=s and cores>=c
                assert old.proven_outside(path,y) == (s>20 or (s>12 and c>=5))
                assert not old.proven_outside(path,y) or sz>20 or (sz>12 and cores>=5)
                sizes+=1
            for tab in tables:
                row = old.H.row(tab,path)
                coeff = row[1:]
                floor0 = row[0]+sum(len(g)*min(Q(),b/2) for g,b in zip(groups,coeff))
                terms = [[b*Q(1,len(a[j]))-min(Q(),b/2) for j in g] for g,b in zip(groups,coeff)]
                actual = physical_raw(tab,t,v,groups)
                assert floor0>=0 and all(z>=0 for g in terms for z in g)
                assert actual == floor0+sum(sum(g) for g in terms), (ds,actual,row)
                mutations += sum(g[-1]!=0 for g in terms[1:])
                checked+=1
    for d in (20,21):
        t = make_tree((d-2,2,2),fan=(0,0,0)); x=(3,4,5)
        path=tuple((len(adjacency(t)[i]),state(t,i)[1]) for i in x)
        assert planted(t,3)==(d,3) and old.proven_outside(path,None)==(d>20)
    for d in (20,21,48):
        t = make_tree((d,2,2),fan=(0,0,0))
        assert planted(t,3)==(d+2,3) and planted(t,3)[0]>=d
    assert mutations
    print(f'PASS labelled nonvirgin adversary: {checked} signed raw rows, {sizes} mask lifts, {mutations} detectable omitted donor terms')


def clock_checks(old):
    rng = random.Random(671)
    H = old.H
    for _ in range(128):
        vals = tuple(tuple(Q(rng.randrange(-100,101),rng.randrange(1,19)) for _ in range(3)) for _ in range(5))
        edges = [(Q(1,i+2),Q(3,i+2),tuple(rng.randrange(-10**20,10**20) for _ in range(6)),i+1)
                 for i in range(4)]
        bounds = []
        for j in range(6):
            actual = [(vals[0][j//2]*H.GRID+sum(r*(v[j]+vals[p][j//2]*H.GRID)
                        for r,(_,_,v,p) in zip(rates,edges)))/(s+sum(rates))
                      for s in (H.LO,H.HI) for rates in product(*((a,b) for a,b,_,_ in edges))]
            bounds.append(ceil(max(actual)) if j%2 else floor(min(actual)))
        assert old.control(vals,edges)==tuple(bounds)
    seen = 0
    for role in range(3):
        valid = []
        for order in permutations('ABCY'):
            if not order.index('A')<order.index('B')<order.index('C') or order.index('Y')<order.index('ABC'[role]):
                continue
            valid.append(order)
            reg = {(('A',),None):1}
            for lab in order[1:]:
                parent = {'B':'A','C':'B','Y':'ABC'[role]}[lab]
                new = dict(reg)
                for (spine,y),n in reg.items():
                    if len(spine)<3 and spine[-1]==parent:
                        z=(spine+(lab,),y); new[z]=new.get(z,0)+n
                    if y is None and len(spine)>role and spine[role]==parent:
                        z=(spine,lab); new[z]=new.get(z,0)+n
                reg=new
            assert reg[tuple('ABC'),'Y']==1 and all(y not in x and n==1 for (x,y),n in reg.items())
        assert len(valid)==3-role
        seen+=len(valid)
    # A physical nonfounder D=2 stock with a leaf has precisely two gain
    # arrows and one promotion; with a core child its promotion stock is zero.
    for t in ((-1,0,1,2,0),(-1,0,1,2,3,0)):
        a=adjacency(t); x=2
        assert len(a[x])==2
        lp=sum(Q(1,len(a[j])) for j in a[x])
        birth=sum(Q(1,len(a[x])) for j in a[x] if len(a[j])==1)
        if state(t,x)==(1,0,1):
            assert lp+birth>=Q(3,2) and birth==Q(1,2)
        else:
            assert state(t,x)==(0,1,1) and birth==0
    print(f'PASS {seen} legal atom birth orders, 128 exact Bellman corner boxes, and virgin/bare mother controls')


def archive_ticks(root,w):
    names = ('apex20','apex48','parent20','parent48','sides')
    hashes = ('f0731f373f334af28fce128ce1e22d6f24073b24c0ad08de421833624d9d4ac2',
              '9a4e008d839c0c6415763b96bd8139a0becac231fd6a4ec40f143e6d7d662795',
              'e1e726a748aa2213ea6da64cfc44412d4a91d6b959cde9da2c3a5357fb2f456b',
              '0a5611854b3ba0a44f29b826156532ace0997137f21d3ced7add2abc5444980c',
              '2f8c4f13c712f869c0a16dd4edc33dc289da3e42be3e0b6573aeeca1f1243be3')
    headers = ((20,2,'apex',4,170010),(48,21,'apex',4,463058),
               (20,2,'parent',2,18470),(48,21,'parent',2,49438),
               (11,2,'middle',2,76382),(11,2,'tip',2,85103))
    ticks = []
    for name,pin in zip(names,hashes):
        data = (root/(name+'.log')).read_bytes()
        assert sha256(data).hexdigest()==pin
        *rows, exit_line = data.decode().splitlines()
        assert exit_line.endswith('exit=0')
        if name.startswith('parent'):
            unconditional = tuple(map(int,rows.pop(0).split()[-2:]))
        for row in rows:
            d,host,role,cap,states = headers[len(ticks)]
            dims = '(11, 10, 10) [3, 2, 2]' if d==11 else f'({d}, 5, 5) [6, 2, 2]'
            actual, values = row.split(' outward 1e-16 ticks ')
            assert actual == (f'{dims} old min 8 atom {role} degree cap {cap} disjoint False '
                              f'ead947b complement True min apex {host} states {states}')
            low,up,_ = literal_eval(values)
            pair=(low[0],up[0])
            if role=='parent':
                pair=tuple(floor(w*x) for x in pair)
                assert pair==unconditional
            ticks.append(pair)
    print('PASS five hashed archive receipts, six role/cap/state rows, and both parent Palm reweightings')
    return ticks


def composition(old,archive):
    note = (R2/BASE/'old/EAD947B.md').read_text()
    ticks = [(int(a),int(b)) for a,b in re.findall(r'^\| (?:Apex child[^|]*|Parent[^|]*|Middle child|Tip child) \| [^|]* \| (\d+) \| (\d+) \|',note,re.M)]
    assert len(ticks)==6, ticks
    floors = tuple(Q(sum(col),10**16) for col in zip(*ticks))
    assert floors==(Q(10361983749,10**16),Q(14352363253,10**16))
    need = Q(1,10**6)
    assert min(floors)>need
    lo=old.H.LO
    w=1-Q(1,2)/(lo+Q(3,2))
    assert 0<w<1
    if archive:
        assert archive_ticks(archive,w)==ticks
    print('PASS hV floors and strict threshold margins:',*floors,*(a-need for a in floors),'Palm weight floor',w)
    denom=(Q(1357493815758092773,10**18),Q(273436004214812983,2*10**17))
    ead=(Q(2586244278310822397274177587149,4571936719109803018500000000000),
         Q(6574835699364880591816419305857,11621030179129551777500000000000))
    e911=(Q(6853552627779976665165032629043939443,12115632305640977999025000000000000000),
          Q(328741600548942531604102687021150749,581051508956477588875000000000000000))
    for pin,start,expected in (('ead947b',ead,(565678819949264,565769419622021)),
                               ('e911a1c',e911,(565679256596412,565769102233118))):
        a,b=(v+s*r/d for v,s,r,d in zip(start,(1,-1),floors,denom))
        assert start[0]<a<b<start[1] and (floor(a*10**15),ceil(b*10**15))==expected
        if pin=='ead947b':
            assert str(a) in note and str(b) in note
        print(f'PASS {pin} exact composition:',a,b,'outward',*(f'0.{k:015d}' for k in expected))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('lower',type=Path); parser.add_argument('upper',type=Path)
    parser.add_argument('--archive',type=Path,help='also verify the five pinned external producer receipts')
    args=parser.parse_args()
    assert __debug__
    for path,pin in PINS.items():
        assert sha256((R2/path).read_bytes()).hexdigest()==pin,path
    old=load(BASE+'/old/replay.py')
    tabs=[]
    for path,pin,mode in zip((args.lower,args.upper),HASHES,(2,1)):
        with path.open('rb') as f:
            assert file_digest(f,'sha256').hexdigest()==pin
        c=old.H.Cert(path)
        assert (c.mode,c.D,c.K,c.S)==(mode,18,17,10**11)
        tabs.append(old.H.Local(c))
    composition(old,args.archive)
    classifier(old,tabs)
    physical(old,tabs)
    clock_checks(old)


if __name__=='__main__':
    main()
