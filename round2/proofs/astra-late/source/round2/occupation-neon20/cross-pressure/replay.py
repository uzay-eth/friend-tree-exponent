"""Same-generator D20/D18 stationary charged pressure on common inner masks."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
import json
from pathlib import Path
import sys
from time import process_time, monotonic


HERE = Path(__file__).resolve().parent
OLD = HERE.parents[1]/'occupation-neon'
ROWS = HERE.parent/'rows'
sys.path.insert(0,str(OLD/'fringe'))
from charges import Charges, Cert
from extract import PAYLOAD_SHA

sys.path.insert(0,str(OLD/'embedding-flow'))
from solve import IncidenceModel, lookup_family, roles, vertex

LO, HI = Q('0.565681716887067'), Q('0.565766354357632')
ZERO = Q(0)
FAMILIES = ('f','g','hV')
PINS = {
    ROWS/'extract.py':'f9352f511c7eb38f94aee4cf6bbce4a3286ae0a4dcb13f577578773ab4dbf1fe',
    ROWS/'lower.json':'b59a2c20b31e7b09f88df436e94e75921d22ac424a759a58c3caf118f7692d4d',
    ROWS/'upper.json':'eeb67cf94792f0a6372f72b02214fa07656a166a65935997d06703a8d9ed85ad',
    OLD/'fringe/charges.py':'44ed53fc7a0fc9ad25d941cd4160fe43aa9d00da334adb4e7de91a9a1a7c97e8',
    OLD/'higher/extract.py':'9f40a50a3623756b277243d192f672e76c45b4b9784715ed4d0d8fb27d36131c',
}
DEN18 = {'lower':(Q(113465364604,4*10**11),Q('1.356833021311')),
         'upper':(Q(115462837236,4*10**11),Q('1.367169995210'))}
OUT = 10**6


def load(path,name):
    spec = spec_from_file_location(name,path)
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def affine(table,f,x,check_box=True):
    terms = table.terms if isinstance(table,Charges) else table.affine
    a,ports = terms(f,x)
    if check_box:
        assert a+sum(n*min(ZERO,b/2) for _,_,n,b in ports)>=0,(f,x,a,ports)
    return a,ports


def fixtures(new, tables):
    """Whole charged sum equals the same labelled-arrow generator for all four."""
    checked = 0
    for p,l in [((-1,0,1,2,3),(1,3,2,3,4)),((-1,0,1,1,2,2),(1,1,2,1,2,1)),
                ((-1,0,1,2),(1,15,17,18)),((-1,0,1,1,1),(1,23,1,2,3))]:
        par = new.audit.build(p,l)
        graph = new.audit.graph(par)
        children,x,core = new.audit.physical_view(par)
        assert len(core)==len(p)
        arrows = new.audit.arrows(par)
        for tab in tables:
            c = tab.c
            phi = c.potential(children,x)
            sign = 1 if c.mode==2 else -1
            whole = sign*(sum(r*(c.potential(*new.audit.physical_view(par+(i,))[:2])-phi)
                              for i,r in arrows.items())-c.lam*phi)/c.S
            rows = ZERO
            for f,ids in new.embeddings(p):
                a,free = affine(tab,f,tuple(x[i] for i in ids),False)
                for j,mother,k,b in free:
                    v = ids[j]
                    donors = (p[v],) if mother else tuple(i for i in children[v] if i not in ids)
                    assert len(donors)==k
                    a += b*sum((Q(1,len(graph[core[w]])) for w in donors),ZERO)
                assert a>=0
                rows += a
            assert whole==rows,(whole,rows,p,l)
            checked+=1
    print('PASS',checked,'whole-potential labelled-arrow/charged-sum equalities on a common generator',flush=True)


def add(*terms):
    row = defaultdict(Q)
    for x,w in terms:
        for i,v in x.items():
            row[i] += w*v
    return {i:v for i,v in row.items() if v}


def mask(model,tables,limits):
    costs = {f:[defaultdict(Q) for _ in tables] for f in FAMILIES}
    n = defaultdict(int)
    for level in range(1,4):
        for tree in model.levels[level]:
            f,ps = lookup_family(tree)
            if f not in FAMILIES:
                continue
            slots = tuple(vertex(tree,p) for p in ps)
            if max(x.degree for x in slots)>limits[level-1]:
                continue
            x = tuple((t.degree-t.arity-1,t.arity,1) for t in slots)
            n[f]+=1
            for tab,row in zip(tables,costs[f],strict=True):
                a,ports = affine(tab,f,x)
                model.mass(row,tree,a)
                for i,mother,count,b in ports:
                    col = model.port(tree,ps[i],mother)
                    assert col is not None
                    if b:
                        row[col] += b
    print('common owned mask rows',dict(n),flush=True)
    return costs


def named(model,cost):
    return {model.variables[i][0]:v for i,v in cost.items() if v}


class Bounds:
    def __init__(self,module,model,seconds):
        self.p = module
        self.model = model
        self.lp = module.Pricing(model)
        self.seconds = seconds

    def one(self,c,up=False):
        obj = {i:-v for i,v in c.items()} if up else c
        bound = self.lp.run(obj,'ipm',self.seconds,1,1)
        return -bound.value if up else bound.value

    def both(self,c):
        return self.one(c),self.one(c,True)

    def ratio_test(self,a,b,t,up=False):
        bound = self.one(add((a,1),(b,-t)),up)
        return bound<=0 if up else bound>=0


def render(label, q):
    print(label,*q,'decimals',*[float(u) for u in q],flush=True)


def rate_convert(band):
    values = [s*q for s,q in product((LO,HI),band)]
    return min(values),max(values)


def ratio(b,a):
    assert a[0]>0
    return min(p/q for p,q in product(b,a)),max(p/q for p,q in product(b,a))


def whole(side,lams,den,mask18,mask20):
    """Unbooked omitted pressure: use nonnegativity and common-mu Dynkin."""
    d18,d20 = den
    sign = 1 if side=='lower' else -1
    # All-rate inequalities: P_i=sign*(mu-lambda_i)*πΦ_i, Π_i≥masked_i≥lower_i.
    m18,m20 = max(ZERO,mask18[0]),max(ZERO,mask20[0])
    lo,hi = LO,HI
    for l,(_,d),m in zip(lams,den,(m18,m20),strict=True):
        if sign==1: lo=max(lo,l+m/d)
        else: hi=min(hi,l-m/d)
    assert lo<=hi
    def evaluate(d1,d2,up):
        # Difference of a linear upper bound and max of two lower bounds,
        # or the reversed expression. Include every interior kink exactly.
        l1,l2 = lams
        if up:
            candidates = [lo,hi]+[s for s in (l1+sign*m18/d1,) if lo<s<hi]
            return max(sign*(s-l2)*d2-max(m18,sign*(s-l1)*d1) for s in candidates)
        candidates = [lo,hi]+[s for s in (l2+sign*m20/d2,) if lo<s<hi]
        return min(max(m20,sign*(s-l2)*d2)-sign*(s-l1)*d1 for s in candidates)
    whole_pressure = evaluate(d18[1],d20[0],False),evaluate(d18[0],d20[1],True)
    omitted = tuple(max(ZERO,max(sign*(s-l)*d[1] for s in (lo,hi))-m)
                    for l,d,m in zip(lams,den,(m18,m20),strict=True))
    render(side+' whole pressure difference (global ceilings and common mu)',whole_pressure)
    render(side+' omitted-side pressure ceilings (D18,D20)',omitted)
    render(side+' omitted contrast uncertainty',(-omitted[0],omitted[1]))
    numer = tuple(sign*(s-lams[1])/(sign*(s-lams[0])) for s in (lo,hi))
    render(side+' whole D20/D18 charged pressure ratio (common mu)',
           (min(numer)*d20[0]/d18[1],max(numer)*d20[1]/d18[0]))


def ratio_search(engine,a,b,initial,rounds=4):
    lo,hi = initial
    # Use short rational probes; each tested sign is certified separately.
    for up in (False,True):
        good,bad = (hi,lo) if up else (lo,hi)
        for _ in range(rounds):
            t = Q(round((good+bad)*OUT/2),OUT)
            if not min(good,bad)<t<max(good,bad):break
            if engine.ratio_test(a,b,t,up):good=t
            else:bad=t
        if up:hi=good
        else:lo=good
    return lo,hi


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for d in (18,20):
        for s in ('lo','hi'):
            p.add_argument(f'--d{d}-{s}',type=Path,default=Path(f'b3f{d}_{s}.bin'))
    p.add_argument('--caps',nargs='+',default=['4:4:4:4'],help='nonincreasing degree caps by model order')
    p.add_argument('--share',nargs='+',type=int,choices=range(4),default=[0,3])
    p.add_argument('--seconds',type=float,default=15)
    p.add_argument('--families',nargs='+',choices=FAMILIES,default=FAMILIES)
    p.add_argument('--ratios',action='store_true')
    p.add_argument('--ratio-rounds',type=int,default=2)
    p.add_argument('--joint-only',action='store_true')
    p.add_argument('--four-way',action='store_true',help='also allocate a signed four-potential contrast jointly')
    p.add_argument('--pay',help='hold paid f:g:hV cut fixed across model cuts')
    a = p.parse_args()
    assert __debug__ and 0<a.seconds<=120 and a.ratio_rounds>=0
    start = process_time(),monotonic()
    for path,want in PINS.items():
        assert sha256(path.read_bytes()).hexdigest()==want,path
    new = load(ROWS/'extract.py','cross_d20_rows')
    pricing = load(OLD/'embedding-flow/pricing/replay.py','cross_pricing')
    from hashlib import file_digest
    tables = {}
    den = {}
    stamps = {}
    for name,tag in (('lower','lo'),('upper','hi')):
        oldpath = getattr(a,'d18_'+tag)
        stamps[oldpath] = new.stamp(oldpath)
        with oldpath.open('rb') as fp:
            assert file_digest(fp,'sha256').hexdigest()==PAYLOAD_SHA[name]
        c = Cert(oldpath)
        c2 = new.Cert(getattr(a,'d20_'+tag),name)
        eps = 1 if name=='lower' else -1
        assert c.mode==c2.mode==(2 if eps==1 else 1) and c2.lam-c.lam==eps*Q(7,10**6)
        data = json.loads((ROWS/(name+'.json')).read_bytes())
        assert data['payload_sha256']==new.HASH[name] and tuple(data['root_ticks'])==new.ROOT[name]
        assert data['source_commit']==new.COMMIT and data['source_blobs']==new.SOURCE
        tab = new.Rows(c2)
        for l,k,h,b,r,*_ in data['f_rows']:
            assert tab.first_body((l,k,h))==(Q(b),Q(r))
        den[name] = DEN18[name],(c2.c0,Q(data['recharge']['ceiling']))
        tables[name] = Charges(c),tab
        print(name,'authenticated D18/D20 and matched',len(data['f_rows']),'D20 exported f rows',flush=True)
    fixtures(new,[t for ts in tables.values() for t in ts])
    for cuts in a.caps:
        limits = tuple(map(int,cuts.split(':')))
        assert len(limits)==4 and 2<=min(limits)<=max(limits)<=16
        paid = tuple(map(int,a.pay.split(':'))) if a.pay else limits[:3]
        assert len(paid)==3 and all(2<=d<=c for d,c in zip(paid,limits[:3],strict=True))
        comparison = {}
        for sharing in a.share:
            model = IncidenceModel(limits[0],3,sharing,low=LO,high=HI,caps=limits[:max(3,sharing+1)])
            engine = Bounds(pricing,model,a.seconds)
            print('cuts',limits,'paid',paid,'shared through',sharing,'columns/equalities/inequalities',len(model.variables),
                  len(model.equalities),len(model.upper_rows),flush=True)
            contrasts={}
            components={}
            for side,ts in tables.items():
                fs = mask(model,ts,paid)
                parts = [([f],*fs[f]) for f in a.families]
                # Cumulative ownership: each connected embedding has one charge.
                if len(parts)>1:
                    parts.append((a.families,*(add(*((fs[f][i],1) for f in a.families)) for i in range(2))))
                for fam,c1,c2 in parts[-1:] if a.joint_only else parts:
                    label = side+' '+'+'.join(fam)
                    key = side,tuple(fam)
                    objnames = named(model,c1),named(model,c2)
                    assert comparison.setdefault(key,objnames)==objnames
                    diff = add((c2,1),(c1,-1))
                    direct = rate_convert(engine.both(diff))
                    render(label+' shared-variable D20−D18 per core',direct)
                    if len(fam)!=len(a.families):
                        continue
                    contrasts[side] = diff,direct
                    b1,b2 = (rate_convert(engine.both(c)) for c in (c1,c2))
                    components[side] = c1,c2,b1,b2
                    render(label+' masked D18',b1)
                    render(label+' masked D20',b2)
                    blind = b2[0]-b1[1],b2[1]-b1[0]
                    render(label+' independent-pressure subtraction',blind)
                    # Same physical source on both sides: the positive common
                    # normalization cancels in the stationary pressure ratio.
                    if a.ratios and b1[0]>0:
                        interval = ratio(b2,b1)
                        # Intersect with the separate enclosure on (D20-D18)/D18.
                        via_diff = ratio(direct,b1)
                        seed = max(interval[0],1+via_diff[0]),min(interval[1],1+via_diff[1])
                        bisect = ratio_search(engine,c2,c1,seed,a.ratio_rounds)
                        render(label+' masked D20/D18 ratio (independent)',interval)
                        render(label+' masked D20/D18 ratio (certified shared)',bisect)
                    whole(side,tuple(t.c.lam for t in ts),den[side],b1,b2)
            if a.four_way:
                (l,lb),(u,ub) = (contrasts[s] for s in ('lower','upper'))
                all_four = rate_convert(engine.both(add((l,1),(u,1))))
                render('four-potential shared contrast (sum of lower and upper charged pressure changes)',all_four)
                render('four-potential independent-side subtraction interval', (lb[0]+ub[0],lb[1]+ub[1]))
                if a.ratios:
                    l,u = (components[s] for s in ('lower','upper'))
                    total18 = tuple(x+y for x,y in zip(l[2],u[2],strict=True))
                    total20 = tuple(x+y for x,y in zip(l[3],u[3],strict=True))
                    if total18[0]>0:
                        blind = ratio(total20,total18)
                        shift = ratio(all_four,total18)
                        seed = max(blind[0],1+shift[0]),min(blind[1],1+shift[1])
                        render('all-four masked ratio (total D20 / total D18)',ratio_search(
                            engine,add((l[1],1),(u[1],1)),add((l[0],1),(u[0],1)),seed,a.ratio_rounds))
    for ts in tables.values():
        ts[1].c.check_unchanged()
    assert all(new.stamp(path)==stamp for path,stamp in stamps.items())
    print('PASS CPU/wall',round(process_time()-start[0],3),round(monotonic()-start[1],3),flush=True)


if __name__=='__main__':
    main()
