"""Physical stationary occupation relaxation; build exactly, solve approximately."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path


def acc(table, key, j, x):
    if x:
        r = table[key]
        r[j] = r.get(j, 0) + x


def add(*terms):
    out = {}
    for scale, expr in terms:
        for j, x in expr.items():
            out[j] = out.get(j, 0) + scale*x
    return {j: x for j, x in out.items() if x}


def v(j, x=1):
    return {j: x}


class Model:
    def __init__(self, N, D3=0, plant=None):
        self.N, self.D3, self.plant = N, D3, plant
        self.names, self.cap, self.rows = [], [], []
        self.states = [(d,c) for d in range(2,N+1) for c in range(1,d+1)]
        self.states += [(0,c) for c in range(1,N+2)]  # D>N; last c is overflow c>=N+1
        self.p = {s: self.var(('p',s)) for s in self.states}
        self.edge_drift, self.flow, self.em = (defaultdict(dict) for _ in range(3))
        self.out, self.inc, self.ox, self.ix, self.cx, self.px = (defaultdict(dict) for _ in range(6))
        self.ext = {}
        self.b = {}
        self.build()

    def var(self, name, cap=1):
        self.names.append(name)
        self.cap.append(cap)
        return len(self.names)-1

    def bval(self, s):
        d,c = s
        return F(d-c, d) if not self.plant == 'promotion' else F(d-c, d+1)

    def xs(self, s):
        d,c = s
        return (F(1,d),) if d else (F(0),F(1,max(self.N+1,c)))

    def grow(self, s):
        d,c = s
        return (d+1,c) if d < self.N else (0,c)

    def promo(self, s):
        d,c = s
        return d, min(c+1,self.N+1)

    def row(self, name, expr, rhs=0, a=0, b=0):
        self.rows.append((name, expr, rhs, a, b))

    def eq(self, name, expr, rhs=0):
        self.row(('+',name),expr,rhs)
        self.row(('-',name),add((-1,expr)),-rhs)

    def station(self, name, pop, low, high=None):
        # low <= Af <= high, so low <= b f and a f <= high.
        self.row(('station-high',name),low,b=pop)
        self.row(('station-low',name),add((-1,high if high is not None else low)),a=pop)

    def build(self):
        N = self.N
        cond = defaultdict(dict)
        extall = defaultdict(dict)
        for p,c in product(self.states, repeat=2):
            if p[1] == 1:
                continue
            st = p,c
            for xp,xc in product(self.xs(p),self.xs(c)):
                j = self.var(('e',p,c,xp,xc))
                acc(self.em, st, j, 1)
                acc(self.out,p,j,1); acc(self.inc,c,j,1)
                acc(self.ox,p,j,xc); acc(self.ix,c,j,xp)
                acc(self.cx,p,j,xp); acc(self.px,c,j,xc)
                for i, (s, xi, xj) in enumerate(((p,xp,xc),(c,xc,xp))):
                    if s[0]:
                        acc(self.flow,(i,st),j,s[0]-s[1]+(0 if self.plant=='donor' and i==1 else xj))
                    bp = 1-s[1]*xi if not s[0] else self.bval(s)
                    if s != (0,N+1) and bp:
                        for sign, target in ((-1,st),(1,st[:i]+(self.promo(s),)+st[i+1:])):
                            acc(self.edge_drift,target,j,sign*bp)
                    # The overflow promotion changes neither endpoint bin.
                for label, key, x in (('out',p,xp),('in',c,xc)):
                    acc(cond,(label,key),j,x)
            for mode,s,mult in (('ancestor',p,1),('sibling',p,p[1]-2),('child',c,c[1]-1)):
                if not s[0] or mult == 0:
                    continue
                j = self.var(('external',mode,st),F(mult,2))
                self.ext[mode,st] = j
                self.row(('donor-cap',mode,st),add((1,v(j)),(-F(mult,2),self.em[st])))
                acc(extall,(mode,s),j,1)

        allp = {j:1 for j in self.p.values()}
        allk = {j:1 for e in self.out.values() for j in e}
        self.eq('total',allp,1)
        self.eq('children',allk,1)
        for s,j in self.p.items():
            d,c = s
            self.eq(('incoming',s),add((1,self.inc[s]),(-1,v(j))))
            pc = add((1,self.cx[s]),(1,self.px[s]))
            self.b[s] = add((1,v(j)),(-1,pc)) if not d else v(j,self.bval(s))
            krow = add((1,self.out[s]),(-(c-1),v(j)))
            xrow = add((1,cond['out',s]),(-(c-1),cond['in',s]))
            if s != (0,N+1):
                self.eq(('outgoing',s),krow)
                if not d:
                    self.eq(('reciprocal-outgoing',s),xrow)
            else:
                self.row('overflow-children',add((-1,krow)))
                self.row('overflow-child-reciprocal',add((-1,xrow)))
                self.row('overflow-child-complement',add((-self.xs(s)[-1],krow),(1,xrow)))
            if not d:
                self.row(('positive-promotion',s),add((-1,self.b[s])))
                self.row(('promotion-cap',s),add((1,self.b[s]),(-1,v(j))))
            if d:
                for mode,mult,source in (('ancestor',c-1,self.ix),('sibling',c-2,self.ox),('child',1,self.ox)):
                    if mode=='sibling' and mult<0:
                        continue
                    self.eq(('donor-sum',mode,s),add((1,extall[mode,s]),(-mult,source[s])))

        totalb = add(*((1,b) for b in self.b.values()))
        self.row('promotion-high',totalb,b=allp)
        self.row('promotion-low',add((-1,totalb)),a=allp)
        birth = (2,1)
        root = defaultdict(dict)
        def transfer(dest, prev, rate):
            root[dest] = add((1,root[dest]),(1,rate))
            root[prev] = add((1,root[prev]),(-1,rate))
        gains = defaultdict(dict)
        for s,j in self.p.items():
            d,c = s
            if d:
                g = gains[d] = add((1,gains[d]),(1,gr := add((1,v(j,d-c)),
                         (self.plant!='donor',self.ix[s]),(1,self.ox[s]))))
                transfer(self.grow(s),s,gr)
            if self.b[s] and c < N+1 and (not d or c < d):
                transfer(self.promo(s),s,self.b[s])
                k = self.promo(s),birth
                self.edge_drift[k] = add((1,self.edge_drift[k]),(1,self.b[s]))
        overflow = 0,N+1
        k = overflow,birth
        self.edge_drift[k] = add((1,self.edge_drift[k]),(1,self.b[overflow]))
        root[birth] = add((1,root[birth]),(1,totalb))
        for s,j in self.p.items():
            self.station(('p',s),v(j),root[s])
        for d,g in gains.items():
            self.station(('crossing',d),{j:1 for (di,_),j in self.p.items() if di==0 or di>d},g)
        for st,pop in self.em.items():
            for i,s in enumerate(st):
                if not s[0]:
                    continue
                modes = ('ancestor','sibling') if i==0 else ('child',)
                gr = add((1,self.flow[i,st]),*((1,v(self.ext[m,st])) for m in modes if (m,st) in self.ext))
                dest = st[:i]+(self.grow(s),)+st[i+1:]
                self.edge_drift[dest] = add((1,self.edge_drift[dest]),(1,gr))
                self.edge_drift[st] = add((1,self.edge_drift[st]),(-1,gr))
        for st,pop in self.em.items():
            self.station(('edge',st),pop,self.edge_drift[st])
        self.reciprocal(cond, gains)

    def reciprocal(self, cond, gains):
        # H=1[D>N]/D. Gain decrement is (x-c*x^2+rho*x^2)/(1+x).
        # Its numerator is nonnegative; dropping rho gives a lower bound,
        # and rho<=c/2 yields an upper bound. Linearize with x<=U.
        lo, hi, mass = [], [], []
        U = F(1,self.N+1)
        for s in self.states:
            if s[0]:
                continue
            x = cond['in',s]
            cx = add((1,self.cx[s]),(1,self.px[s]))
            mass.append((1,x))
            # x-c*x^2 >= x-U*c*x. Also <=x; adding rho*x^2 never exceeds x.
            lo.extend(((1/(1+U),x),(-U/(1+U),cx)))
            hi.append((1,x))
        boundary = F(1,self.N+1),gains[self.N]
        self.station('tail-reciprocal',add(*mass),add(boundary,(-1,add(*hi))),add(boundary,(-1,add(*lo))))

    def exact_rows(self,a,b):
        for name,e,r,low,high in self.rows:
            yield name,add((1,e),(a,low or {}),(-b,high or {})),r


def solve(model,a,b,method='highs-ipm', force_elastic=False):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix
    vals, cols, ip, rhs = [], [], [0], []
    for _, e, r, low, high in model.rows:
        for j,x in add((1,e),(a,low or {}),(-b,high or {})).items():
            vals.append(float(x)); cols.append(j)
        ip.append(len(vals)); rhs.append(float(r))
    A = csr_matrix((vals,cols,ip), shape=(len(rhs),len(model.names)))
    rhs = np.array(rhs)
    obj = np.zeros(A.shape[1])
    options = {'primal_feasibility_tolerance':1e-8,'dual_feasibility_tolerance':1e-8,
               'ipm_optimality_tolerance':1e-10}
    if not force_elastic:
        from types import SimpleNamespace
        from scipy.optimize._highspy import _core as core
        equal, skip = [], set()
        for i, (name, *_) in enumerate(model.rows[:-1]):
            if isinstance(name,tuple) and (name[0]=='+' or (name[0]=='station-high' and a==b)):
                delta = A[i]+A[i+1]
                if not delta.nnz and rhs[i]+rhs[i+1]==0:
                    equal.append(i)
                    skip.update((i,i+1))
        keep = sorted(set(range(len(rhs)))-skip.difference(equal))
        mat = A[keep].tocsc()
        lp = core.HighsLp()
        lp.num_col_, lp.num_row_ = mat.shape[1],mat.shape[0]
        lp.a_matrix_.num_col_, lp.a_matrix_.num_row_ = lp.num_col_,lp.num_row_
        lp.a_matrix_.format_ = core.MatrixFormat.kColwise
        lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = mat.indptr,mat.indices,mat.data
        lp.col_cost_, lp.col_lower_, lp.col_upper_ = obj, np.zeros(len(obj)), np.full(len(obj),np.inf)
        lp.row_upper_ = rhs[keep]
        eq = set(equal)
        lp.row_lower_ = [rhs[i] if i in eq else -np.inf for i in keep]
        h = core._Highs()
        opts = core.HighsOptions()
        opts.output_flag, opts.threads = False, 1
        opts.solver = {'highs-ipm':'ipm','highs-ds':'simplex','pdlp':'pdlp'}[method]
        opts.run_crossover = 'off'
        for k,value in options.items(): setattr(opts,k,value)
        h.passOptions(opts); h.passModel(lp); h.run()
        x = np.asarray(h.getSolution().col_value)
        if len(x)==len(obj) and np.isfinite(x).all():
            residual = rhs-A@x
            violation = max(0.,-min(residual),-min(x))
            if violation <= 1e-8:
                return SimpleNamespace(x=x,fun=0.,ineqlin=SimpleNamespace(residual=residual),original_violation=violation)
        if method=='pdlp':
            raise RuntimeError(f'PDLP did not produce a checked primal: {h.getModelStatus()}')
    # Uniform row elasticization: its dual is a Farkas candidate iff optimum>0.
    from scipy.sparse import hstack
    aug = hstack([A, -np.ones((len(rhs),1))], format='csc')
    result = linprog(np.r_[obj,1],A_ub=aug,b_ub=rhs,
                     bounds=(0,None),method=method,options=options)
    if not result.success:
        raise RuntimeError(result.message)
    return result


def save(path, model, a,b,result):
    from verify_farkas import fingerprint
    y = -result.ineqlin.marginals
    cert = {'N':model.N,'D3':model.D3,'slice':[str(a),str(b)],
            'fingerprint':fingerprint(model.exact_rows(a,b),model.cap),
            'dual':[[i,str(F(float(x)).limit_denominator(10**15))] for i,x in enumerate(y) if x>1e-14]}
    Path(path).write_text(json.dumps(cert,separators=(',',':'))+'\n')


def cli(construct, default_motif=0):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cutoff',type=int,default=8)
    ap.add_argument('--motif',type=int,default=default_motif)
    ap.add_argument('--prior',nargs=2,default=['.5650','.5665'])
    ap.add_argument('--slices',type=int,default=5)
    ap.add_argument('--certificate-prefix')
    ap.add_argument('--method',choices=['highs-ipm','highs-ds','pdlp'],default='highs-ipm')
    ap.add_argument('--probe',nargs=2,help='check one fixed slice, without slicing the prior')
    ap.add_argument('--elastic',action='store_true',help='always solve the elastic LP')
    args = ap.parse_args()
    if args.slices<1 or args.cutoff<2:
        ap.error('positive slice count and degree cutoff >= 2 required')
    from verify_farkas import check
    from time import perf_counter
    start = perf_counter()
    m = construct(args.cutoff,args.motif)
    print(f'N={m.N} D3={m.D3} variables={len(m.names)} rows={len(m.rows)} build={perf_counter()-start:.3f}s',flush=True)
    lo,hi = map(F,args.prior)
    def run(a,b,tag):
        if not 0<a<=b:
            ap.error('slice must have 0 < lower <= upper')
        start = perf_counter()
        r = solve(m,a,b,args.method,args.elastic)
        margin = None
        if r.fun>1e-10:
            y = ((i,F(float(-x)).limit_denominator(10**15)) for i,x in enumerate(r.ineqlin.marginals) if x< -1e-14)
            margin = check(m.exact_rows(a,b),m.cap,dict(y))[0]
            if margin<0 and args.certificate_prefix:
                save(f'{args.certificate_prefix}-{tag}.json',m,a,b,r)
        max_violation = getattr(r,'original_violation',max(0.,r.fun-min(r.ineqlin.residual)))
        print(f'[{float(a):.9f},{float(b):.9f}] elastic={r.fun:.5g} row_violation<={max_violation:.3g} exact='
              f'{float(margin) if margin is not None else "--"} wall={perf_counter()-start:.3f}s',flush=True)
        return r, margin is not None and margin<0
    if args.probe:
        run(*map(F,args.probe),'probe')
        return
    if lo==hi:
        run(lo,hi,'point')
        return
    l,l_bad = run(lo,lo,'lo0'); h,h_bad = run(hi,hi,'hi0')
    step = (hi-lo)/args.slices
    trial = 0
    while l_bad and lo+step<=hi:
        trial += 1
        _,bad = run(lo,lo+step,f'lo{trial}')
        if not bad:
            step /= 2
            if step<F(1,10**6): break
        else: lo+=step
    step = (F(args.prior[1])-F(args.prior[0]))/args.slices
    trial = 0
    while h_bad and hi-step>=lo:
        trial += 1
        _,bad = run(hi-step,hi,f'hi{trial}')
        if not bad:
            step /= 2
            if step<F(1,10**6): break
        else: hi-=step
    print(f'prior after exact exclusions: [{lo}, {hi}]',flush=True)
    if l.fun==0 and h.fun==0:
        pvars = [(abs(l.x[j]-h.x[j]), s, l.x[j], h.x[j]) for s,j in m.p.items()]
        ext = sum(abs(l.x[j]-h.x[j]) for j in m.ext.values())
        print('endpoint p L1=',sum(x[0] for x in pvars),'ext donor L1=',ext,
              'tail p L1=',sum(x[0] for x in pvars if not x[1][0]),
              'largest p gaps=',sorted(pvars, reverse=True)[:4],flush=True)
        groups = defaultdict(lambda:[0,0.,0.])
        for j, name in enumerate(m.names):
            group = groups[name[0]]
            group[0] += 1
            group[1] += abs(l.x[j]-h.x[j])
            group[2] += max(l.x[j],h.x[j])
        for k,(n,difference,scale) in groups.items():
            print(f'endpoint group {k:16} n={n:5d} L1={difference:.6g} mass_scale={scale:.6g}',flush=True)


if __name__=='__main__':
    cli(Model)
