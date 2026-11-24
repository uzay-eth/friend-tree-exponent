#!/usr/bin/env python3
"""Independent inverse search: integer slabs, high precision confirmation, local controls."""
import os
for var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[var] = '1'

import argparse
from dataclasses import dataclass
from functools import cache
from itertools import combinations
from math import ceil, floor, gcd, log, exp
import time

import mpmath as mp
import numpy as np


IR_LO = '0.482266918323368718268897622913549185185240982589329346611091'
IR_HI = '0.482266918323368718268897622913549185185240982589329346611092'
NCONTROL = 512
SEED = 68133021


@dataclass(frozen=True)
class Fn:
    name: str
    f: object
    # Upper bounds on |f|, |f'| and |f''| for lo <= x < 1.
    bounds: object

    @cache
    def at(self, x):
        return self.f(x)


def fixed(a, b, c):
    return lambda lo: (a, b, c)


def power(k):
    return Fn(f'x^{k}' if k > 1 else 'x', lambda x: x**k, fixed(1, k, k*(k-1)))


X = power(1)
PWR = [None, X] + [power(k) for k in range(2, 9)]


def times_power(fn, k):
    def bounds(lo):
        a, b, c = fn.bounds(lo)
        return a, b+k*a, c+2*k*b+k*(k-1)*a
    return Fn(f'{PWR[k].name}*{fn.name}', lambda x: x**k * fn.at(x), bounds)


def fz(z):
    bound = max(1, exp(-z))
    return Fn(f'F[{z:g}]', lambda x: mp.hyp1f1(1, x+2, -z)/(x+1),
              fixed(bound, bound, 2*bound))


F = fz(1)
FVAR = [fz(z) for z in (-1, .5, 2, 3)]
GP = Fn('Gamma(1+x)', lambda x: mp.gamma(1+x), fixed(2, 3, 8))
G = Fn('Gamma(x)', mp.gamma, lambda lo: (1/lo+1, 1/lo**2+1, 2/lo**3+2))
GVAR = [GP]
for z in (.5, 1, 2):
    GVAR.extend((Fn(f'lowerGamma(1+x,{z:g})',
                    lambda x, z=z: mp.gammainc(1+x, 0, z), GP.bounds),
                 Fn(f'upperGamma(1+x,{z:g})',
                    lambda x, z=z: mp.gammainc(1+x, z, mp.inf), GP.bounds)))

HVAR = []
for b, z in ((2, -1), (2, 1), (3, 2)):
    HVAR.append(Fn(f'1F1(1+x;{b};{z})',
                   lambda x, b=b, z=z: mp.hyp1f1(1+x, b, z),
                   fixed(exp(abs(z)), abs(z)*exp(abs(z)), z*z*exp(abs(z)))))
for b, z in ((1, 1), (2, .5), (2, 2)):
    HVAR.append(Fn(f'U(1,{b}+x,{z:g})',
                   lambda x, b=b, z=z: mp.hyperu(1, b+x, z),
                   fixed(1/z+1/z**2, 1/z**2+2/z**3, 2/z**3+6/z**4)))
for z in (-.5, .5, .75):
    HVAR.append(Fn(f'2F1(1,1+x;2;{z:g})',
                   lambda x, z=z: mp.hyp2f1(1, 1+x, 2, z),
                   fixed(1/(1-abs(z)), abs(z)/(1-abs(z))**2, 2*z*z/(1-abs(z))**3)))
for z in (-.5, .5):
    HVAR.append(Fn(f'Mellin[1/(1-({z:g})t)]',
                   lambda x, z=z: mp.hyp2f1(1, 1+x, 2+x, z)/(x+1), fixed(2, 2, 4)))
for z in (.5, .75):
    HVAR.append(Fn(f'Mellin[1_{{t>{z:g}}}]',
                   lambda x, z=z: (1-mp.power(z, x+1))/(x+1), fixed(1, 1, 2)))

EX = Fn('exp(x)', mp.exp, fixed(exp(1), exp(1), exp(1)))
EM = Fn('exp(-x)', lambda x: mp.exp(-x), fixed(1, 1, 1))
LX = Fn('log(x)', mp.log, lambda lo: (-log(lo), 1/lo, 1/lo**2))
PSI = Fn('psi(x)', lambda x: mp.digamma(x), lambda lo: (1/lo+1, 1/lo**2+2, 2/lo**3+3))
SPARSE = [X, EX, EM, LX]
for name, base, ub in (('2', 2, 2), ('3', 3, 3), ('pi', None, 3.15)):
    SPARSE.append(Fn(f'{name}^x', lambda x, base=base: mp.power(base or mp.pi, x),
                     fixed(ub, ub*log(ub), ub*log(ub)**2)))
SPARSE += [G, PSI, Fn('cos(pi*x)', lambda x: mp.cospi(x), fixed(1, 3.15, 10)),
           Fn('sin(pi*x)', lambda x: mp.sinpi(x), fixed(1, 3.15, 10))]

LVAR = [Fn('x*log(x)', lambda x: x*mp.log(x), lambda lo: (1, 1-log(lo), 1/lo))]
for k in (-2, -1, -.5, .5, 1, 2):
    LVAR.extend((Fn(f'x^({k:g}*x)', lambda x, k=k: x**(k*x),
                    lambda lo, k=k: (exp(abs(k)), exp(abs(k))*abs(k)*(1-log(lo)),
                                     exp(abs(k))*(k*k*(1-log(lo))**2+abs(k)/lo))),
                 Fn(f'x*exp({k:g}*x)', lambda x, k=k: x*mp.exp(k*x),
                    fixed(exp(max(k, 0)), exp(max(k, 0))*(1+abs(k)),
                          exp(max(k, 0))*(2*abs(k)+k*k)))))
for c, branch in ((.5, 0), (1, 0), (2, 0), (-.25, 0), (-.25, -1)):
    LVAR.append(Fn(f'W[{branch}]({c:g}*x)', lambda x, c=c, k=branch: mp.lambertw(c*x,k),
                   lambda lo: (8+2*abs(log(lo)), 3/lo, 6/lo**2)))


def constant(name, function):
    def bounds(lo):
        return float(abs(function()))+1e-12, 0, 0
    return Fn(name, lambda x: function(), bounds)


CONSTANTS = [constant('pi', lambda: mp.pi), constant('e', lambda: mp.e),
             constant('log(2)', lambda: mp.log(2)), constant('log(3)', lambda: mp.log(3)),
             constant('EulerGamma', lambda: mp.euler), constant('zeta(3)', lambda: mp.zeta(3)),
             constant('Catalan', lambda: mp.catalan), constant('1/pi', lambda: 1/mp.pi),
             constant('W(1)', lambda: mp.lambertw(1)), constant('phi', lambda: mp.phi)]
CONSTANTS += [constant(f'sqrt({n})', lambda n=n: mp.sqrt(n))
              for n in (2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19)]


@cache
def row_count(hh, mode):
    """Exact primitive equation count, including bounded constant term."""
    ans = 0
    for d in range(1, min(hh)+1):
        # Mobius inversion. All modes require at least one nonzero coefficient.
        n, mu, p = d, 1, 2
        while p*p <= n:
            if n % p == 0:
                n //= p
                if n % p == 0:
                    mu = 0
                    break
                mu = -mu
            p += 1
        if n > 1:
            mu = -mu
        h0, *hs = (h//d for h in hh)
        if mode == 'poly' or mode == 'field':
            v = hs[0]
            for h in hs[1:]:
                v *= 2*h+1
        elif mode == 'kernel':
            h = hs[0]
            v = ((2*h+1)**3-1)//2 * (2*h+1)**2
        elif mode == 'sparse':
            v = hs[0]
            for h in hs[1:]:
                v *= 2*h
        else:
            raise ValueError(mode)
        ans += mu*v*(2*h0+1)
    return ans


@cache
def coefficient_rows(hs, mode):
    axes = [np.arange(-h,h+1,dtype=np.int16) for h in hs]
    if mode != 'kernel':
        axes[0] = np.arange(1,hs[0]+1,dtype=np.int16)
    if mode == 'sparse':
        axes[1:] = [a[a != 0] for a in axes[1:]]
    rows = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1,len(hs))
    if mode == 'kernel':
        a,b,c = rows[:,:3].T
        rows = rows[(a>0) | ((a==0)&(b>0)) | ((a==0)&(b==0)&(c>0))]
    return rows


@dataclass
class Template:
    basis: tuple
    hh: tuple
    mode: str
    label: str = ''

    def __post_init__(self):
        self.count = row_count(self.hh, self.mode)


def template(basis, h, mode, label=''):
    return Template(tuple(basis), (h,)*(1+len(basis)), mode,
                    label or ','.join(f.name for f in basis))


def kernel_templates(fs, h):
    return [template([f,times_power(f,1),times_power(f,2),X,PWR[2]], h,
                     'kernel', f.name) for f in fs]


def families():
    for d,h in ((2,75),(3,24),(4,11),(5,7),(6,5)):
        yield f'Algebraic d={d}', [template(PWR[1:d+1][::-1],h,'poly')]
    fields = []
    for a,b in ((2,13),(2,3),(2,5),(3,5),(3,7),(5,13)):
        basis = [X]+[constant(f'sqrt({n})',lambda n=n:mp.sqrt(n)) for n in (a,b,a*b)]
        fields.append(Template(tuple(basis),(9,16,9,9,9),'field',f'Q(sqrt({a}),sqrt({b}))'))
    yield 'Multiquadratic', fields
    yield 'HKRZ F[1]', kernel_templates([F], 7)
    yield 'HKRZ other z', kernel_templates(FVAR, 5)
    yield 'Gamma kernels', kernel_templates(GVAR, 5)
    yield 'Kummer U / hypergeom / Mellin', kernel_templates(HVAR, 4)
    yield 'Sparse elementary + special', [template(b,4 if n<=4 else 2 if n==5 else 1,'sparse')
                                           for n in range(1,len(SPARSE)+1) for b in combinations(SPARSE,n)]
    yield 'Affine two constants', [Template((X,a,b),(10,14,10,10),'sparse',
                                           f'{a.name}+{b.name}')
                                  for a,b in combinations(CONSTANTS,2)]
    yield 'Self powers and Lambert', [template([X,f],32,'sparse',f.name) for f in LVAR]
    fs = list(dict.fromkeys([F,*FVAR,*GVAR,*HVAR,*SPARSE[1:]]))
    yield 'Affine fixed points', [template([X,f],28,'sparse',f.name) for f in fs]


@dataclass
class Window:
    x: object
    radius: object


@dataclass
class Hit:
    root: object
    formula: str
    sub: str = ''


def distinct(hits):
    unique = []
    for hit in sorted(hits, key=lambda h: h.root):
        if not unique or abs(hit.root-unique[-1].root)>mp.mpf('1e-54'):
            unique.append(hit)
    return unique


def polynomial_exact_degree(aa):
    # Only invoked for the actual target/control/nearest-window roots. Gauss's
    # lemma delegates exact Z[x] irreducibility, not a numerical guess.
    from sympy import Poly, Symbol
    return Poly.from_list(aa, Symbol('x')).is_irreducible


def equation(t, aa):
    return ' '.join([str(aa[-1])] + [f'{a:+d}*{b.name}' for a,b in zip(aa,t.basis) if a]) + ' = 0'


def confirm(t, aa, w):
    if t.mode == 'poly' and not polynomial_exact_degree(aa):
        return [], 0
    terms = [(a,b) for a,b in zip(aa,t.basis) if a]

    def f(x):
        return aa[-1]+mp.fsum(a*b.f(x) for a,b in terms)

    # Constant-coefficient linear templates have an exact one-step formula.
    if t.basis[0] == X and all(b.bounds(.4)[1] == 0 for b in t.basis[1:]):
        r = w.x-f(w.x)/aa[0]
        return ([r] if abs(r-w.x)<=w.radius else []), 0

    def solve(lo, hi, flo, fhi):
        if not flo:
            return lo
        if not fhi:
            return hi
        try:
            root = mp.findroot(f, (lo,hi), tol=mp.mpf('1e-142'), maxsteps=22)
            if lo <= root <= hi and abs(f(root))<mp.mpf('1e-69'):
                return root
        except (ValueError, ZeroDivisionError):
            pass
        for i in range(265):
            m = (lo+hi)/2
            fm = f(m)
            if not fm or hi-lo < mp.mpf('1e-77'):
                return m
            if (fm > 0) == (flo > 0):
                lo,flo = m,fm
            else:
                hi,fhi = m,fm
        return (lo+hi)/2

    # Derivative bounds rule out hidden pairs inside tiny windows. Report
    # intervals where neither monotonicity nor the absence of a root is clear.
    def walk(lo,hi,flo,fhi,depth):
        mid=(lo+hi)/2
        rad=(hi-lo)/2
        fm=f(mid)
        lower=np.nextafter(float(lo),0)
        bounds=[(abs(a),b.bounds(lower)) for a,b in terms]
        m1,m2=(sum(a*bd[k] for a,bd in bounds)*(1+1e-14) for k in (1,2))
        if abs(fm)>rad*m1+mp.mpf('1e-85'):
            return [],0
        der = mp.diff(f,mid)
        if abs(der)>rad*m2:
            if not flo or not fhi or (flo>0)!=(fhi>0):
                return [solve(lo,hi,flo,fhi)],0
            return [],0
        if abs(fm)>abs(der)*rad+m2*rad**2/2+mp.mpf('1e-85'):
            return [],0
        if depth==12:
            if not flo or not fhi or (flo>0)!=(fhi>0):
                return [solve(lo,hi,flo,fhi)],1
            return [],1
        a,n=walk(lo,mid,flo,fm,depth+1)
        b,m=walk(mid,hi,fm,fhi,depth+1)
        return a+b,n+m

    lo,hi = w.x-w.radius,w.x+w.radius
    return walk(lo,hi,f(lo),f(hi),0)


def scan_template(t, windows, hitlists):
    lo = np.nextafter(float(min(w.x-w.radius for w in windows)),0)
    radii = np.array([float(w.radius) for w in windows])
    values = np.array([[float(f.at(w.x)) for w in windows] for f in t.basis])
    amp,lip,_ = np.array([f.bounds(lo) for f in t.basis]).T*(1+1e-14)
    rows = coefficient_rows(t.hh[1:],t.mode)
    uncertain=0
    for start in range(0,len(rows),2048):
        block=rows[start:start+2048]
        z=-(block.astype(float)@values)
        ab=np.abs(block)
        eps=2e-13*(1+ab@amp)
        tol=(ab@lip)[:,None]*radii+eps[:,None]
        n=np.rint(z)
        mask=(np.abs(n)<=t.hh[0]) & (np.abs(z-n)<=tol)
        rr,cc=np.nonzero(mask | (tol>=.5))
        for i,j in zip(rr,cc):
            aa=[int(y) for y in block[i]]
            if tol[i,j]<.5:
                acs=(int(n[i,j]),)
            else:
                acs=range(max(-t.hh[0],ceil(z[i,j]-tol[i,j])),
                          min(t.hh[0],floor(z[i,j]+tol[i,j]))+1)
            for a0 in acs:
                if gcd(*aa,a0)!=1:
                    continue
                roots,u=confirm(t, aa+[a0], windows[j])
                uncertain+=u
                for root in roots:
                    hitlists[j].append(Hit(root,equation(t,aa+[a0]),t.label))
    return uncertain


def fwindows(center,sigma):
    r=3*sigma
    near=min(mp.mpf('0.000004'), center/4,(1-center)/4)
    lower=max(r*mp.mpf('1.01'),center-mp.mpf('.045'))
    upper=min(1-r*mp.mpf('1.01'),center+mp.mpf('.045'))
    a=max(mp.mpf(0),center-2*r-lower)
    b=max(mp.mpf(0),upper-center-2*r)
    # When windows are so wide that nearby disjoint controls do not exist, use
    # the available full interval and make the overlap explicit in the output.
    rng=np.random.default_rng(SEED)
    uu=[mp.mpf(float(u)) for u in rng.random(NCONTROL)]
    if a+b:
        ctrl=[lower+u*(a+b) if u*(a+b)<a else upper-(1-u)*(a+b) for u in uu]
    else:
        lo=max(r,(center-mp.mpf('.045')))
        hi=min(1-r,(center+mp.mpf('.045')))
        ctrl=[lo+u*(hi-lo) for u in uu]
        print('WARNING: target/control windows can overlap at this sigma')
    ir=(mp.mpf(IR_LO)+mp.mpf(IR_HI))/2
    return [Window(center,r),Window(ir,mp.mpf('1e-50')),Window(center,near)]+[Window(x,r) for x in ctrl]


def scan_rational(windows,H=10000):
    hitlists=[[] for _ in windows]
    # In (0,1): sum phi(q), q=2..H.
    phi=list(range(H+1))
    for p in range(2,H+1):
        if phi[p]==p:
            for q in range(p,H+1,p):
                phi[q]-=phi[q]//p
    qv=np.arange(2,H+1)
    for j,w in enumerate(windows):
        lows=np.ceil(float(w.x-w.radius)*qv-3e-12).astype(int)
        highs=np.floor(float(w.x+w.radius)*qv+3e-12).astype(int)
        for i in np.flatnonzero(lows<=highs):
            q=int(qv[i])
            for p in range(max(1,int(lows[i])),min(q-1,int(highs[i]))+1):
                if gcd(p,q)==1:
                    root=mp.mpf(p)/q
                    if abs(root-w.x)<=w.radius:
                        hitlists[j].append(Hit(root,f'{p}/{q}'))
    return sum(phi[2:]),hitlists


@cache
def expression_catalogue(H):
    # Every expression x^y, x/y and log_x(y) is tried as written and multiplied
    # by p/q, 1<=p,q<=H, gcd(p,q)=1. Also x^(p/q) and x^(-p/q).
    qs=np.array([(p,q) for p in range(1,H+1) for q in range(1,H+1) if gcd(p,q)==1])
    cc=[(f.name,f.at(mp.mpf(0))) for f in CONSTANTS]
    matrices=[[],[],[]]
    meta=[]
    for i,(an,a) in enumerate(cc):
        for j,(bn,b) in enumerate(cc):
            if i==j:
                continue
            matrices[0].append(float(a)**float(b))
            matrices[1].append(float(a)/float(b))
            matrices[2].append(log(float(b))/log(float(a)))
            meta.append((an,a,bn,b))
    block=np.array(matrices)[:,:,None]*(qs[:,0]/qs[:,1])[None,None,:]
    pp=np.array([float(c) for _,c in cc])[None,:,None]**(np.array([1,-1])[:,None,None]*(qs[:,0]/qs[:,1])[None,None,:])
    catalogues=[]
    for vals in (block,pp):
        mid=vals.ravel()
        ii=np.flatnonzero((mid>=-4e-13)&(mid<=1+4e-13))
        ix=ii[np.argsort(mid[ii])]
        catalogues.append((vals.shape,mid[ix],ix))
    return cc,meta,qs,catalogues,block.size+pp.size


def explicit_expressions(windows,H=84):
    cc,meta,qs,catalogues,n=expression_catalogue(H)
    hitlists=[[] for _ in windows]
    for w,hits in zip(windows,hitlists):
        lo=float(w.x-w.radius)-4e-13
        hi=float(w.x+w.radius)+4e-13
        for which,(shape,vs,ix) in enumerate(catalogues):
            l,r=np.searchsorted(vs,(lo,hi),side='left')
            for k,i,j in zip(*np.unravel_index(ix[l:r],shape)):
                p,q=map(int,qs[j]); fac=mp.mpf(p)/q
                if which==0:
                    an,a,bn,b=meta[i]
                    value=(a**b if k==0 else a/b if k==1 else mp.log(b)/mp.log(a))*fac
                    form=(f'{an}^{bn}',f'{an}/{bn}',f'log_{an}({bn})')[k]
                    text=f'({p}/{q})*({form})'
                else:
                    an,a=cc[i]; p=p if k==0 else -p
                    value=a**(mp.mpf(p)/q)
                    text=f'{an}^({p}/{q})'
                if abs(value-w.x)<=w.radius:
                    hits.append(Hit(value,text))
    return n,hitlists


def extra_nearest(center, windows, hitlists, run):
    radius=windows[2].radius
    cap=min(mp.mpf('.003'),center/4,(1-center)/4)
    while len(distinct(hitlists[2]))<3 and radius<cap:
        radius=min(8*radius,cap)
        bigger=run(Window(center,radius))
        hitlists[2]=bigger
    return radius


def fmt(h,w,sigma):
    return f'{mp.nstr(h.root,20)} ({mp.nstr((h.root-w.x)/sigma,8)} sigma): {h.formula}'


def summarize(name, n, windows, hitlists, sigma, uncertain, radius, temps=None):
    h,ir,near,*ctrl=map(distinct,hitlists)
    counts=np.array([len(c) for c in ctrl]); empirical=(1+sum(counts>=len(h)))/(1+len(ctrl))
    hist=','.join(f'{k}:{v}' for k,v in zip(*np.unique(counts,return_counts=True)))
    print(f'ROW | {name} | {n} | {len(h)} | {counts.mean():.6f} | {counts.max()} | {empirical:.6f} | {len(ir)} | {hist} | {uncertain}',flush=True)
    for hit in sorted(near,key=lambda q:abs(q.root-windows[0].x))[:3]:
        print('  near',fmt(hit,windows[0],sigma))
    if len(near)<3:
        print('  nearest search radius',mp.nstr(radius,8),'returned only',len(near))
    eq_counts={t.label:t.count for t in temps} if temps else {}
    for hit in h:
        print(f'  TARGET (family n={n}; template n={eq_counts.get(hit.sub,n)}):',fmt(hit,windows[0],sigma))
    for hit in ir:
        print(f'  IR-0 (family n={n}): {mp.nstr(hit.root,75)}:',hit.formula)
    return name,n,len(h),str(empirical),int(counts.max()),len(ir)


def pslq_ir(x):
    ff=F.at(x)
    ee=mp.exp(x)
    g=mp.gamma(x)
    lx=mp.log(x)
    pp=[x**k for k in range(9)]
    gp=lambda t: mp.gamma(t)*mp.digamma(t)
    cases=[('powers 0..2',pp[:3],10**18,lambda t:[t**k for k in range(3)]),
           ('powers 0..4',pp[:5],10**11,lambda t:[t**k for k in range(5)]),
           ('powers 0..6',pp[:7],10**8,lambda t:[t**k for k in range(7)]),
           ('powers 0..8',pp,10**6,lambda t:[t**k for k in range(9)]),
           ('1,x,x²,F,xF',[1,x,x*x,ff,x*ff],10**11,lambda t:[1,t,t*t,F.f(t),t*F.f(t)]),
           ('1,x,x²,F,xF,x²F',[1,x,x*x,ff,x*ff,x*x*ff],10**9,
            lambda t:[1,t,t*t,F.f(t),t*F.f(t),t*t*F.f(t)]),
           ('1,x,e^x,xe^x,e^-x',[1,x,ee,x*ee,1/ee],10**11,
            lambda t:[1,t,mp.exp(t),t*mp.exp(t),mp.exp(-t)]),
           ('1,x,log(x),xlog(x)',[1,x,lx,x*lx],10**14,
            lambda t:[1,t,mp.log(t),t*mp.log(t)]),
           ('1,x,Gamma(x),Gamma\'(x)',[1,x,g,gp(x)],10**14,
            lambda t:[1,t,mp.gamma(t),gp(t)])]
    # Coefficients are returned by PSLQ only when |a_i| < maxcoeff. Tolerance
    # allows the certified IR bracket's 5e-61 radius at the stated heights.
    print('IR PSLQ: working dps=95; strict maxcoeff; maxsteps=20000')
    for label,vs,height,basis in cases:
        tol=mp.mpf('1e-58')*height
        rel=mp.pslq(mp.matrix(vs),tol=tol,maxcoeff=height,maxsteps=20000)
        if rel:
            def f(t):
                return mp.fsum(a*b for a,b in zip(rel,basis(t)))
            residual=f(x)
            lo,hi=map(mp.mpf,(IR_LO,IR_HI))
            ends=[f(t) for t in (lo,hi)]
            print(f'  {label}; bound={height}; tol={mp.nstr(tol)}; approximate relation={rel}; residual={mp.nstr(residual,15)}')
            print('    endpoint residuals:', ', '.join(mp.nstr(v,15) for v in ends))
            if label.startswith('powers'):
                m=sum(k*(k-1)*abs(a) for k,a in enumerate(rel))
                m1=abs(mp.diff(f,x))+m*(hi-lo)/2
                if abs(residual)>m1*(hi-lo)/2+mp.mpf('1e-80'):
                    print('    rejected: |residual| exceeds the maximum change across the certified bracket',
                          mp.nstr(m1*(hi-lo)/2,15))
            der=mp.diff(f,x)
            if der:
                r=mp.findroot(f,(x,x-residual/der),tol=mp.mpf('1e-151'),verify=False)
                print('    numerical root minus IR midpoint:',mp.nstr(r-x,20))
        else:
            print(f'  {label}; bound={height}; tol={mp.nstr(tol)}; no relation returned')


def main():
    start=time.monotonic()
    mp.mp.dps=95
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('VALUE'); parser.add_argument('SIGMA')
    args=parser.parse_args()
    x,s=map(mp.mpf,(args.VALUE,args.SIGMA))
    if not (0<3*s<min(x,1-x)):
        parser.error('require SIGMA > 0 and [VALUE-3*SIGMA,VALUE+3*SIGMA] strictly inside (0,1)')
    ws=fwindows(x,s)
    print(f'VALUE={args.VALUE}; SIGMA={args.SIGMA}; radius={mp.nstr(3*s,15)}; controls={NCONTROL}; seed={SEED}')
    print(f'IR-0 midpoint={mp.nstr(ws[1].x,62)}; radius=1e-50')
    print('ROW | family | candidate descriptions tried | target distinct roots | ctrl mean | ctrl max | empirical upper-tail p | IR distinct roots | ctrl histogram(count:windows) | unresolved subintervals',flush=True)
    rows=[]
    for name,calc in [('Rational',scan_rational),('Expressions in constants',explicit_expressions)]:
        n,hh=calc(ws)
        rad=extra_nearest(x,ws,hh,lambda w:calc([w])[1][0])
        rows.append(summarize(name,n,ws,hh,s,0,rad))
    for name,temps in families():
        hits=[[] for _ in ws]
        n=sum(t.count for t in temps)
        unresolved=0
        if name=='Affine two constants':
            for t in temps:
                t.label='constant pair'
        for t in temps:
            unresolved+=scan_template(t,ws,hits)
        def rerun(w):
            nonlocal unresolved
            h=[[]]
            for t in temps:
                unresolved+=scan_template(t,[w],h)
            return h[0]
        radius=extra_nearest(x,ws,hits,rerun)
        rows.append(summarize(name,n,ws,hits,s,unresolved,radius,temps))
    pslq_ir(ws[1].x)
    print(f'Tests={len(rows)}; Bonferroni alpha per row={.05/len(rows):.8f}; wall={time.monotonic()-start:.3f}s')


if __name__=='__main__':
    main()
