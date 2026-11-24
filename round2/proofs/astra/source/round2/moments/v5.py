#!/usr/bin/env python3
"""Physical occupation dual for the frozen nominal four-vertex witness."""

import argparse
import gzip
import hashlib
import json
import re
import subprocess
import sys
from bisect import bisect_left
from collections import defaultdict
from fractions import Fraction as F
from math import ceil, floor
from pathlib import Path

from dual import Solver, verify
from marginal_lift import HighDegrees
from occupation import HIGH, LOW, ONE, TAIL, ZERO, Occupation, add, axpy, coarse, combine
from proofcheck import IntegerProof
from resources import resources


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'round2/outputs'
sys.path.insert(0, str(ROOT/'round2/motif-cert'))

DEGREE = (8, 16, 32, 48, 71, 95, 127, 191, 255, 383, 511, 767, 1023, 1399)
CHILD = (0, 1, 2, 3, 7, 15, 31, 63)


def group(d, k, K):
    return d*(d-1)//2-1+k if d <= K else 10000+100*min(5, bisect_left(DEGREE, d))+bisect_left(CHILD, k)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def vector_hash(row):
    h=hashlib.sha256(b'moments-sparse-v1\n')
    for j,a in sorted(row.items()):
        if a: h.update(f'{j}={a}\n'.encode('ascii'))
    return h.hexdigest()


def matrix_hash(model):
    h=hashlib.sha256(b'moments-matrix-v1\n')
    def line(parts):
        h.update((json.dumps(parts,separators=(',',':'),default=str)+'\n').encode('ascii'))
    for name,cap in zip(model.names,model.cap):
        line(['column',name,str(cap)])
    for sense,rows in [('eq',model.equal),('ge',model.atleast)]:
        for name,row,rhs in rows:
            line([sense,name,str(rhs),vector_hash(row)])
    return h.hexdigest()


def prior_interval():
    """Verify the two wider V5 sheets before narrowing a stationary balance."""
    ends,guards=[],{}
    for side in ('lower','upper'):
        path=OUT/f'moments-v5-anchor-{side}-dual.json.gz'
        with gzip.open(path,'rt') as stream:
            w=json.load(stream)
        if w.get('use_prior') or w['side']!=side or 'shift' not in w['rows']:
            raise ValueError('ineligible stationary predecessor')
        subprocess.run([sys.executable,'-B',__file__,'check',path],check=True,stdout=subprocess.DEVNULL)
        with gzip.open(OUT/f'moments-v5-{side}-k{w["K"]}.tsv.gz','rt') as stream:
            head=next(stream).split()
        nominal=F(int(head[2]),int(head[3]))
        q=F(w['duals']['shift']['bound'])
        ends.append(nominal+(q if side=='lower' else -q))
        guards[path.name]=digest(path)
    if not LOW<=ends[0]<=ends[1]<=HIGH:
        raise ValueError('predecessor endpoints contradict the initial interval')
    return ends,guards


def analytic_caps(path, lam):
    import tail
    h, newborn, slope = tail.read_tail(path)
    q, x = tail.Q, tail.X
    a, b, c, e = ((h['f_p0'], h['f_p1'], 0, 0) if h['mode'] == 1 else
                  (0, h['f_p0'], -h['f_p1'], -h['f_p2']))
    w = a+b*q+c*q*q+e*q*x
    den = (1+x)*(1+x)
    gain = b*(1-q)*(1+x)+c*(1-q)*(2*q+(1+q)*x)+e*((1-2*q)*x-q*x*x)
    bud = -b*x+c*(-2*q*x+x*x)-e*x*x
    sign = 1 if h['mode']==1 else -1
    polys = [(rho*gain+(q*(bud+newborn+slope*x)-w*lam)*den)*sign for rho in (q,(1+q)*F(1,2))]
    xmax = F(1, h['DF'])
    result = []
    for i in range(16):
        hi = max(max(tail.bernstein(p, F(i,16), F(i+1,16), xmax)) for p in polys)
        result.append(ceil(hi if hi >= 0 else hi/(1+xmax)**2))
    return result


def lift_inputs(side, K, profile):
    def filename(name):
        return OUT / name
    side_num = 1 if side=='upper' else 2
    low = filename('mc-low-'+side+'-nominal.jsonl')
    signed_rows = [json.loads(z) for z in low.read_text().splitlines()]
    lam = F(signed_rows[0]['P'], signed_rows[0]['Q'])
    snapshot = filename('mc-'+side+'.bin')
    profile_tail = filename('mc-mbar-'+side+'-tail.tsv')
    higher = filename('mh-'+side+'-nominal.log')
    two = filename('occupation-two-'+side+'.txt')
    five = filename('occupation-five-6.txt' if side=='upper' else 'occupation-five-lower.txt')
    data = {}
    with gzip.open(profile, 'rt') as stream:
        header = list(map(int, next(stream).split()[1:]))
    assert header[0]==side_num and F(header[1], header[2])==lam
    match = re.search(r'higher_residual_per_core=(-?\d+)/(\d+)', higher.read_text())
    assert match and int(match[2])==header[3]
    data['higher'] = int(match[1])
    credits = {}
    # Our role names identify the donor. The historical two-site file names
    # the receiver: "parent vertex slack" receives from the younger child.
    for kind, receiver, txt in [('parent','child',two),('child','parent',two),('quartet',None,five)]:
        lines = [s for s in txt.read_text().splitlines() if s.startswith('mode '+str(side_num)+' ')
                 and (receiver is None or (' '+receiver+' vertex slack ') in s)]
        assert len(lines)==1
        p = re.search(r'newborn reserve (\d+)/(\d+) at discount 1133/2000', lines[0])
        assert p and int(p[2])==header[3]
        credits[kind] = int(p[1])
    data['reserves'] = credits
    data['higher_hashes'] = {p.name:digest(p) for p in (low,higher,two,five,profile_tail,profile)}
    data['snapshot_sha256'] = digest(snapshot)
    data['mbar_sha256'] = digest(filename('mc-mbar-'+side+'.i64'))
    data['analytic_f_caps'] = analytic_caps(snapshot, lam)
    data['analytic_mbar'] = [list(map(int,z.split())) for z in profile_tail.read_text().splitlines()[1:]]
    data['base_f'] = max([r['cap'] for r in signed_rows if r['family'].startswith('f_')]+data['analytic_f_caps'])
    data['base_g'] = next(r['cap'] for r in signed_rows if r['family']=='g')
    return data


def score_inputs(side,profile):
    """The unsigned denominator bind uses only the potential and comparison tables."""
    snapshot=OUT/f'mc-{side}.bin'
    maxima=OUT/f'mc-mbar-{side}.i64'
    tail_file=OUT/f'mc-mbar-{side}-tail.tsv'
    with gzip.open(profile,'rt') as st:
        header=next(st).split()
    if int(header[1])!=(1 if side=='upper' else 2):
        raise ValueError('wrong unsigned potential')
    return {'score_inputs':{p.name:digest(p) for p in (snapshot,maxima,tail_file,profile)},
            'analytic_mbar':[list(map(int,s.split())) for s in tail_file.read_text().splitlines()[1:]]}


class Application:
    def __init__(self, model, profile, ext, degree_tests=False, fy_path=None, score_path=None,floor_path=None,
                 tip_path=None,star_path=None,star_children=0,gap_limits=None,edge_score_path=None,score_only=False,
                 tail_score=False):
        self.m, self.ext = model, ext
        self.degree_tests = degree_tests
        vertices, cells, edge_bounds = {}, {}, {}
        with gzip.open(profile, 'rt') as stream:
            header = list(map(int, next(stream).split()[1:]))
            for line in stream:
                tag, *rest = line.split()
                z = list(map(int,rest))
                if tag=='v':
                    vertices[z[0],z[1]+1] = tuple(z[2:])
                elif tag=='b':
                    cells[tuple(z[:3])] = tuple(z[3:])
                elif tag=='e':
                    edge_bounds[z[0],z[1]]=tuple(z[2:])
                else:
                    raise ValueError(tag)
        side, pn, pd, self.unit, df, p0, p1, p2, kt = header
        self.lam, self.side = F(pn,pd), side
        for kb in range(len(CHILD)+1):
            kmin = 0 if kb==0 else CHILD[kb-1]+1
            kmax = CHILD[kb] if kb<len(CHILD) else None
            for qi in range(16):
                qlo,qhi = F(qi,16), F(qi+1,16)
                if kmax is not None:
                    qlo = max(qlo, 1-F(kmax+1,df))
                if qlo>qhi:
                    continue
                if side==1:
                    fmax = p0+p1*qhi
                else:
                    top = min(qhi, max(qlo, F(p0,2*p1)))
                    fmax = p0*top-p1*top*top
                bottom = max(intercept + max(0,F(slope,df)) for k,intercept,slope in ext['analytic_mbar']
                             if min(kmin,kt+1) <= k <= min(kmax if kmax is not None else kt+1,kt+1))
                cap=0 if score_only else ext['analytic_f_caps'][qi]
                cells[len(DEGREE),kb,qi] = (df,0,kmin,kmax,cap,ceil(fmax+bottom)+4,0,0)
        self.frow, self.grow, self.mbar, self.score_defect = {}, {}, {}, {}
        self.credit = {'parent': {}, 'child': {}}
        self.source = {}
        self.residual = model.lift.residual
        self.blocks = defaultdict(lambda: {w:{} for w in ('mass','slots','parent','child')})
        self.ranges = {}
        self.cell_state = []
        self.gap=None
        if gap_limits:
            from v5_surplus import SingletonGap
            self.gap=SingletonGap(gap_limits)
        self.tail_lines={}
        if tail_score:
            from score_tail import build_lines
            if not degree_tests: raise ValueError('affine score requires physical bin birth intensities')
            name='upper' if side==1 else 'lower'
            self.tail_lines=build_lines(OUT/f'mc-mbar-{name}.i64',header,model.lift.end,DEGREE,CHILD,
                                        cells,ext['analytic_mbar'])
        self._vertices(vertices)
        self._tail(cells)
        self._edges(edge_bounds)
        if self.gap: self.gap.build(self)
        if score_path:
            self._score_overlaps(score_path)
        self.fy_credit,self.fy_columns={},{}
        self.tip=self.star=None
        if fy_path:
            from v5_overlap import add_donors
            self.fy_credit,self.fy_columns=add_donors(self,fy_path)
        if tip_path:
            from v5_hub import TipEdges
            self.tip=TipEdges(self,tip_path)
        if star_path:
            from v5_star import HubStars
            self.star=HubStars(self,star_path,star_children)
        if edge_score_path: self._conditional_edges(edge_score_path)
        if degree_tests:
            self._degree_balance()
        if score_only:
            if fy_path or tip_path or star_path or floor_path or gap_limits:
                raise ValueError('the unsigned denominator model has no source variables')
            self.den=combine((1,self.mbar),(-1,self.score_defect))
            unit={j:1 for j in model.v.values()}
            self.rows={'M-per-core':(self.den,unit,-1),'M-per-birth':(self.den,model.nu,-1)}
        else: self._observables()
        if floor_path:
            self._lower_score(floor_path,header,cells)

    def record(self, gid, rlo, rhi, vectors):
        for key, expr in vectors.items():
            axpy(self.blocks[gid][key], expr)
        old = self.ranges.get(gid, (rlo, rhi))
        self.ranges[gid] = min(rlo,old[0]), max(rhi,old[1])

    def donor(self, lo, hi, slots, incoming, part):
        if lo > 0:
            axpy(self.credit[part], slots, F(lo,2*self.unit))
            axpy(self.credit[part], incoming, -F(lo,self.unit))
        elif hi < 0:
            axpy(self.credit[part], incoming, -F(hi,self.unit))

    def _vertices(self, vertices):
        m = self.m
        assert set(vertices)==set(m.finite)|set(m.lift.states)
        for s,(fc,sc,delta) in vertices.items():
            d,c=s
            if d<=m.cutoff:
                j=m.v[s]
                p={m.e[a,s]:m.inv[a] for a in m.parents}
                kids={m.e[s,w]:m.inv[w] for w in m.labels} if c>1 else {}
            else:
                j=m.lift.v[s];p={m.lift.p[s]:1}
                kids={m.lift.ch[s]:1} if c>1 else {}
            mass,slots={j:1},{j:c-1}
            gid=group(d,c-1,m.cutoff)
            fields=dict(mass=mass,slots=slots,parent=p,child=kids)
            self.record(gid,F(1,d),F(1,d),fields)
            if self.gap: self.gap.record(('v',s),gid,d,c,delta,delta,fields)
            self.frow[j]=F(fc,self.unit);self.mbar[j]=F(sc,self.unit)
            self.donor(delta,delta,mass,p,'parent')
            self.donor(delta,delta,slots,kids,'child')

    def _tail(self, cells):
        m, rem = self.m, self.residual
        totals = defaultdict(dict)
        for key, (dlo,dhi,klo,khi,fc,sc,dl,du) in sorted(cells.items()):
            _,_,qi=key
            rlo,rhi=F(1,dhi) if dhi else 0,F(1,dlo)
            mass,slots,p,kids=(m.variable(('v5-tail',key,name)) for name in ('mass','slots','parent','child'))
            m.ge(('k-min',key), {slots:1,mass:-klo})
            if khi is not None:
                m.ge(('k-max',key), {slots:-1,mass:khi})
            m.ge(('parent-cap',key),{mass:1,p:-2})
            m.ge(('child-cap',key),{slots:1,kids:-2})
            if self.degree_tests:
                b=m.variable(('v5-tail',key,'birth'))
                m.ge(('bin-birth-lo',key),{b:1,mass:-F(qi,16)})
                m.ge(('bin-birth-hi',key),{b:-1,mass:F(qi+1,16)})
                m.ge(('bin-degree-lo',key),{mass:1-rlo,b:-1,slots:-rlo})
                m.ge(('bin-degree-hi',key),{mass:rhi-1,b:1,slots:rhi})
                add(totals['birth'],b,1)
                self.cell_state.append((dlo,dhi,klo,khi,mass,b,p,kids))
                if self.tail_lines:
                    lines=self.tail_lines[key]
                    cap=max([F(0)]+[F(sc-a,self.unit) for _,a in lines])
                    gap=m.variable(('v5-unsigned-tail-excess',key),cap)
                    for slope,a in lines:
                        m.ge(('unsigned-tail-line',key,slope),{gap:1,mass:F(a-sc,self.unit),b:F(slope,self.unit)})
                    self.score_defect[gap]=1
            gid=group(dlo,klo,m.cutoff)
            assert dhi==0 or gid==group(dhi,klo,m.cutoff)
            fields=dict(mass={mass:1},slots={slots:1},parent={p:1},child={kids:1})
            self.record(gid,rlo,rhi,fields)
            if self.gap: self.gap.record(('b',key),gid,dlo,klo+1,dl,du,fields)
            self.frow[mass]=F(fc,self.unit);self.mbar[mass]=F(sc,self.unit)
            self.donor(dl,du,{mass:1},{p:1},'parent')
            self.donor(dl,du,{slots:1},{kids:1},'child')
            for name, j, invname in [('mass',mass,'inv'),('slots',slots,'slotinv'),
                                     ('parent',p,'parentinv'),('child',kids,'childinv')]:
                add(totals[name],j,1)
                add(totals[invname,'lo'],j,rlo)
                add(totals[invname,'hi'],j,rhi)
            add(totals['core-frac','lo'],mass,F(15-qi,16))
            add(totals['core-frac','hi'],mass,F(16-qi,16))
        for name in ('mass','slots','parent','child'):
            m.eq(('tail-all',name),combine((1,totals[name]),(-1,rem[name])))
        if self.degree_tests:
            m.eq('bin-total-birth',combine((1,totals['birth']),(-1,rem['mass']),(1,rem['inv']),(1,rem['slotinv'])))
        bands={name:rem[name] for name in ('inv','slotinv','parentinv','childinv')}
        bands['core-frac']=combine((1,rem['inv']),(1,rem['slotinv']))
        for name,known in bands.items():
            m.ge(('tail-rlo',name),combine((1,known),(-1,totals[name,'lo'])))
            m.ge(('tail-rhi',name),combine((-1,known),(1,totals[name,'hi'])))

    def _edges(self,bounds):
        m=self.m
        collapsed=defaultdict(dict)
        inc,out,rin,rout=(defaultdict(dict) for _ in range(4))
        bounds={key:z for key,z in bounds.items() if key[0] in self.blocks and key[1] in self.blocks}
        self.pair_columns={}
        for (a,b),(cap,defect) in sorted(bounds.items()):
            j=m.variable(('v5-edge',a,b))
            self.pair_columns[a,b]=j
            self.grow[j]=F(cap,self.unit);self.score_defect[j]=F(defect,self.unit)
            add(inc[b],j,1); add(out[a],j,1)
            for i,extreme in enumerate(('lo','hi')):
                add(rin[b,extreme],j,self.ranges[a][i]);add(rout[a,extreme],j,self.ranges[b][i])
            add(collapsed[a if a<10000 else -1,b if b<10000 else -1],j,1)
        for gid,fields in self.blocks.items():
            m.eq(('v5-in',gid),combine((1,inc[gid]),(-1,fields['mass'])))
            m.eq(('v5-out',gid),combine((1,out[gid]),(-1,fields['slots'])))
            for tag,mom,actual in [('in',rin,fields['parent']),('out',rout,fields['child'])]:
                m.ge(('v5-reciprocal-lo',tag,gid),combine((1,actual),(-1,mom[gid,'lo'])))
                m.ge(('v5-reciprocal-hi',tag,gid),combine((-1,actual),(1,mom[gid,'hi'])))
        for (a,b),j in m.e.items():
            ca=group(a[0],a[1]-1,m.cutoff) if a[0] else -1
            cb=group(b[0],b[1]-1,m.cutoff) if b[0] else -1
            add(collapsed[ca,cb],j,-1)
        for key,row in collapsed.items():
            m.eq(('v5-edge-overlap',key),row)

    def _score_overlaps(self,path):
        prices={}
        with gzip.open(path,'rt') as stream:
            for line in stream:
                tag,a,m,w,p=line.split()
                prices[tag,int(a),int(m),int(w)]=F(int(p),self.unit)
        def label(s):
            d,c=s
            return d*(d-1)//2-1+c-1 if d else -1
        for (a,m,b),j in self.m.t.items():
            add(self.score_defect,j,prices.get(('t',label(a),label(m),label(b)),0))
        for (a,m,_,b),j in self.m.q.items():
            add(self.score_defect,j,prices.get(('q',label(a),label(m),label(b)),0))
        for (m,a,b),j in self.m.ch.items():
            x,y=sorted((label(a),label(b)))
            add(self.score_defect,j,prices.get(('c',label(m),x,y),0))

    def _degree_balance(self):
        """Test (1-T/d)_+ 1{k>=h}; gains and threshold-crossing buds are nonnegative."""
        m=self.m
        exact=[]
        for s in [*m.finite,*m.lift.states]:
            d,c=s
            if d<=m.cutoff:
                v=m.v[s]
                inp={m.e[a,s]:m.inv[a] for a in m.parents}
                if c>1:
                    axpy(inp,{m.e[s,w]:m.inv[w] for w in m.labels})
            else:
                v=m.lift.v[s];inp={m.lift.p[s]:1}
                if s in m.lift.ch: inp[m.lift.ch[s]]=1
            exact.append((d,c-1,v,inp))
        for T in sorted({0,m.lift.end,72,128,256,512,1024,1400}):
            for h in (0,1,2,3,4,8,16):
                if T==0 and h==0:
                    continue
                value_lo,value_hi,flux_lo,flux_hi={},{},{},{}
                if h==0 and T<2:
                    axpy(flux_lo,m.nu,1-F(T,2))
                for d,k,v,inp in exact:
                    f=max(0,1-F(T,d))
                    b=F(d-k-1,d)
                    if k>=h:
                        delta=F(T,d*(d+1)) if d>=T else 0
                        add(value_lo,v,f)
                        add(flux_lo,v,delta*(d-k-1))
                        axpy(flux_lo,inp,delta)
                    if k+1==h:
                        add(flux_lo,v,b*f)
                value_hi.update(value_lo)
                flux_hi.update(flux_lo)
                for dl,dh,kl,kh,v,b,p,kids in self.cell_state:
                    f0=max(0,1-F(T,dl))
                    f1=max(0,1-F(T,dh)) if dh else 1
                    certain=kl>=h
                    possible=kh is None or kh>=h
                    if possible:
                        add(value_hi,v,f1)
                        low_d=max(T,dl)
                        dL=F(T,low_d+1) if not dh or dh>=T else 0
                        add(flux_hi,b,dL)
                        add(flux_hi,p,dL/low_d)
                        add(flux_hi,kids,dL/low_d)
                    if certain:
                        add(value_lo,v,f0)
                        dL=F(T,dh+1) if dh and dl>=T else 0
                        add(flux_lo,b,dL)
                        if dh:
                            add(flux_lo,p,dL/dh)
                            add(flux_lo,kids,dL/dh)
                    if kl<=h-1 and (kh is None or kh>=h-1):
                        add(flux_hi,b,f1)
                    if kl==kh==h-1:
                        add(flux_lo,b,f0)
                m.ge(('degree-test-low',T,h),combine((1,flux_hi),(-m.low,value_lo)))
                m.ge(('degree-test-high',T,h),combine((m.high,value_hi),(-1,flux_lo)))

    def _conditional_edges(self,path):
        by_triple={}
        with gzip.open(path,'rt') as stream:
            for line in stream:
                tag,a,b,value=line.split()
                a,b=int(a),int(b)
                price=F(int(value),self.unit)
                if tag=='e': add(self.score_defect,self.pair_columns[a,b],price)
                elif tag=='t': by_triple[a,b]=price
                else: raise ValueError(tag)
        for (a,s,b),j in self.m.t.items():
            if a[0]: continue
            mid=group(s[0],s[1]-1,self.m.cutoff)
            end=group(b[0],b[1]-1,self.m.cutoff) if b[0] else -1
            add(self.score_defect,j,by_triple.get((mid,end),0))

    def _observables(self):
        m,ext,S=self.m,self.ext,self.unit
        unit={j:1 for j in m.v.values()}
        self.den=combine((1,self.mbar),(-1,self.score_defect))
        self.net=combine((F(ext['reserves']['quartet'],S),m.nu),(-1,self.frow),(-1,self.grow),
                         (-F(ext['higher'],S),unit),(1,self.fy_credit))
        for role in self.credit:
            # A per-vertex gain jump is bounded by 32*S; the retained finite
            # high receiver has at most DF-1 slots. These generous caps affect
            # only the exact rounding correction to a dual.
            j=m.variable(('usable-singleton-credit',role), 32*1400)
            self.source[role]=j
            extra=self.gap.extra[role] if self.gap else {}
            m.ge(('singleton-lifetime',role),combine((1,{j:1}),(-F(ext['reserves'][role],S),m.nu),(-1,extra)))
            m.ge(('singleton-occupation',role),combine((1,{j:1}),(-1,self.credit[role])))
            add(self.net,j,1)
        base=sum(ext['reserves'].values())
        self.base=combine((F(base,S),m.nu),(-F(ext['base_f']+ext['base_g']+ext['higher'],S),unit))
        signed_singleton=combine((1,self.base),(F(ext['base_f'],S),unit),(-1,self.frow))
        signed_pairs=combine((1,signed_singleton),(F(ext['base_g'],S),unit),(-1,self.grow))
        self.rows={'shift':(self.net,self.den,1),'M-per-core':(self.den,unit,-1),
                   'M-per-birth':(self.den,m.nu,-1),'denominator-only-shift':(self.base,self.den,1),
                   'signed-singleton-shift':(signed_singleton,self.den,1),
                   'signed-row-shift':(signed_pairs,self.den,1),
                   'f-cap':(self.frow,unit,-1),'paired-edge-cap':(self.grow,unit,-1)}
        if self.fy_columns:
            self.rows['edge-donor-min']=self.fy_credit,unit,1
        if self.gap:
            self.rows['extra-singleton-min']=combine((1,self.gap.extra['parent']),(1,self.gap.extra['child'])),unit,1

    def _lower_score(self,path,header,cells):
        m=self.m
        lookup={key:j for j,key in enumerate(m.names)}
        low,tail_rows={},{}
        with gzip.open(path,'rt') as stream:
            for line in stream:
                tag,*fields=line.split(); z=list(map(int,fields))
                if tag=='v':
                    d,k,score=z;s=(d,k+1)
                    j=m.v[s] if d<=m.cutoff else m.lift.v[s]
                    low[j]=F(score,self.unit)
                elif tag=='b':
                    low[lookup['v5-tail',tuple(z[:3]),'mass']]=F(z[3],self.unit)
                elif tag=='a':tail_rows[z[0]]=z[1:]
                else:raise ValueError(tag)
        side,_,_,_,df,p0,p1,p2,kt=header
        for key,(_,_,kmin,kmax,*_) in cells.items():
            if key[0]!=len(DEGREE):continue
            qlo,qhi=F(key[2],16),F(key[2]+1,16)
            if kmax is not None:qlo=max(qlo,1-F(kmax+1,df))
            f=(p0+p1*qlo if side==1 else
               min(p0*q-p1*q*q for q in (qlo,qhi))-F(p2,df)*qhi)
            b=min(z+min(0,F(slope,df)) for k,(z,slope) in tail_rows.items()
                  if min(kmin,kt+1)<=k<=min(kt+1,kmax if kmax is not None else kt+1))
            low[lookup['v5-tail',key,'mass']]=F(floor(f+b)-4,self.unit)
        self.score_floor=low
        self.rows['M-per-birth-min']=low,m.nu,1
        self.rows['fixed-reserve-shift-ceiling']=self.base,low,-1


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=('solve','check'))
    ap.add_argument('witness',type=Path)
    ap.add_argument('--side',choices=('upper','lower'))
    ap.add_argument('--cutoff',type=int,default=4)
    ap.add_argument('--vertex-cutoff',type=int,default=32)
    ap.add_argument('--degree-tests',action='store_true')
    ap.add_argument('--edge-donor',action='store_true')
    ap.add_argument('--score-overlap',action='store_true')
    ap.add_argument('--lower-score',action='store_true')
    ap.add_argument('--tip-edges',action='store_true')
    ap.add_argument('--star-children',type=int,choices=(0,3,7),default=0)
    ap.add_argument('--use-prior',action='store_true')
    ap.add_argument('--original-score',action='store_true')
    ap.add_argument('--singleton-gap',action='store_true')
    ap.add_argument('--edge-score',action='store_true')
    ap.add_argument('--score-only',action='store_true')
    ap.add_argument('--coefficient-pin',type=Path,help='write exact coefficient and price hashes during a signed replay')
    ap.add_argument('--tail-score',action='store_true')
    ap.add_argument('--rows',nargs='+',default=['shift','M-per-core','denominator-only-shift','f-cap','paired-edge-cap'])
    args=ap.parse_args()
    resources()
    if args.action=='solve':
        if args.side is None:
            ap.error('choose a side for the proposal')
        sheet=dict(version=1,side=args.side,K=args.cutoff,M=args.vertex_cutoff,rows=args.rows,
                   degree_tests=args.degree_tests,edge_donor=args.edge_donor,
                   score_overlap=args.score_overlap,lower_score=args.lower_score,tip_edges=args.tip_edges or bool(args.star_children),
                   star_children=args.star_children,use_prior=args.use_prior,original_score=args.original_score,
                   singleton_gap=args.singleton_gap,edge_score=args.edge_score,score_only=args.score_only,
                   tail_score=args.tail_score,duals={})
    else:
        with gzip.open(args.witness,'rt') as stream:
            sheet=json.load(stream)
    side,K,M=sheet['side'],sheet['K'],sheet['M']
    profile=OUT/f'moments-v5-{side}-k{K}.tsv.gz'
    score_only=sheet.get('score_only',False)
    if score_only and (sheet.get('use_prior') or sheet.get('original_score')):
        raise ValueError('standalone nominal denominator needs the original growth interval')
    inputs=score_inputs(side,profile) if score_only else lift_inputs(side,K,profile)
    fy_path=OUT/f'moments-v5-{side}-FY.tsv.gz' if sheet.get('edge_donor') else None
    if fy_path:inputs['FY_sha256']=digest(fy_path)
    score_path=OUT/f'moments-v5-{side}-score-k{K}.tsv.gz' if sheet.get('score_overlap') else None
    if score_path:inputs['score_sha256']=digest(score_path)
    edge_score_path=OUT/f'moments-v5-{side}-edge-score-k{K}.tsv.gz' if sheet.get('edge_score') else None
    if edge_score_path:inputs['edge_score_sha256']=digest(edge_score_path)
    floor_path=OUT/f'moments-v5-{side}-floor-k{K}.tsv.gz' if sheet.get('lower_score') else None
    if floor_path:inputs['floor_sha256']=digest(floor_path)
    tip_path=OUT/f'moments-v5-{side}-tip.tsv.gz' if sheet.get('tip_edges') else None
    if tip_path:inputs['tip_sha256']=digest(tip_path)
    star_limit=sheet.get('star_children',0)
    star_path=OUT/f'moments-v5-{side}-star.tsv.gz' if star_limit else None
    if star_path:inputs['star_sha256']=digest(star_path)
    interval=(LOW,HIGH)
    if sheet.get('use_prior'):
        interval,inputs['stationary_predecessors']=prior_interval()
    if sheet.get('original_score'):
        inputs['original_score_sha256']=digest(ROOT/f'verify/motif-bracket/reference/{side}.bin.xz')
    gap_limits=None
    if sheet.get('singleton_gap'):
        from v5_surplus import source_limits
        gap_limits=source_limits((OUT/f'occupation-two-{side}.txt').read_text())
        inputs['singleton_source_geometry']=gap_limits
    if args.action=='solve':
        sheet['inputs']=inputs
    if sheet['version']!=1 or inputs!=sheet['inputs']:
        raise ValueError('input digest or nominal constants changed')
    m=Occupation(K,low=interval[0],high=interval[1]);HighDegrees(m,M)
    app=Application(m,profile,inputs,sheet.get('degree_tests',False),fy_path,score_path,floor_path,
                    tip_path,star_path,star_limit,gap_limits,edge_score_path,score_only,sheet.get('tail_score',False))
    if sheet.get('original_score'):
        from original_score import envelope
        old,old_hash=envelope(app,side)
        if old_hash!=inputs['original_score_sha256']: raise ValueError('original potential changed')
        app.rows['original-M-per-birth']=old,m.nu,-1
    solver=Solver(m,warm=True) if args.action=='solve' else None
    checker=IntegerProof(m)
    pin=None
    if args.coefficient_pin:
        if args.action!='check' or score_only:
            raise ValueError('coefficient pin requires a signed replay')
        pin=dict(price_sheet=args.witness.name,price_sha256=digest(args.witness),sigma=1 if side=='upper' else -1,
                 nominal=str(app.lam),unit=app.unit,growth_interval=[str(v) for v in interval],
                 columns=len(m.names),equalities=len(m.equal),inequalities=len(m.atleast),matrix_sha256=matrix_hash(m),
                 base_caps=[inputs[k] for k in ('base_f','base_g','higher')],reserve_ticks=inputs['reserves'],
                 objectives={})
    print(f'{side} joint D={K} vertex D={M}: {len(m.names)} columns, {len(m.equal)} equalities, '
          f'{len(m.atleast)} inequalities; nominal {app.lam}',flush=True)
    if score_only: print('Unsigned potential only. Physical occupancy rows on the original growth interval.',flush=True)
    for name in sheet['rows']:
        num,den,sign=app.rows[name]
        if solver:
            initial=F(8,10**6) if name.endswith('shift') else 2 if name=='original-M-per-birth' else 0
            sheet['duals'][name]=solver.ratio(num,den,sign,grid=10**12,initial=initial)
        w=sheet['duals'][name];q=F(w['bound'])
        obj=combine((sign,num),(-sign*q,den))
        a=verify(m,obj,w)
        if a<0 or str(a)!=w['margin'] or checker.evaluate(obj,w)!=a or w['direction']!=sign:
            raise ArithmeticError(('dual failed',name,a))
        if pin is not None:
            pin['objectives'][name]=dict(direction=sign,bound=str(q),verified_margin=str(a),
                                        numerator_sha256=vector_hash(num),denominator_sha256=vector_hash(den),
                                        signed_test_sha256=vector_hash(obj))
        print(f'{name} {float(q):.12f} = {q}; exact margin {float(a):.4g}',flush=True)
        if name=='shift':
            endpoint=app.lam+(q if side=='lower' else -q)
            print(f'physical endpoint {float(endpoint):.12f} = {endpoint}',flush=True)
        if solver:
            partial=dict(sheet,rows=list(sheet['duals']))
            tmp=args.witness.with_suffix('.tmp')
            with tmp.open('wb') as raw, gzip.GzipFile(fileobj=raw,mode='wb',filename='',mtime=0) as out:
                out.write((json.dumps(partial,separators=(',',':'))+'\n').encode())
            tmp.replace(args.witness)
    if pin is not None:
        args.coefficient_pin.write_text(json.dumps(pin,indent=2)+'\n')


if __name__=='__main__':
    main()
