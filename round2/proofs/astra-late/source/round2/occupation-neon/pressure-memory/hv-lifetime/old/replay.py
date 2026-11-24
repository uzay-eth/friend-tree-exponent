#!/usr/bin/env python3
"""Older hV: exact external signs and compressed physical lifetime prices."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import cache
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import ceil, floor, lcm
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('hv_lifetime', HERE.parent/'replay.py')
H = module_from_spec(spec)
spec.loader.exec_module(H)
ZERO = H.ZERO
Z = H.NO_CREDIT
GRID = H.GRID
SCALE = lcm(H.LO.denominator,H.HI.denominator,*range(2,129))
SLOPES = ('parent', 'apex-child', 'middle-child', 'tip-child')


def route(e, x, port):
    A,B,C=x; c=e.c
    return (A in c.p4 and B in c.p4, A in c.ybu and B in c.ybw, B in c.yt,
            B in c.p4 and C in c.p4)[port]


def knots(c,x):
    if sum(x)<c.D:
        return (x,), (1,)
    d=sum(x); a=c.D; b=d-a
    return c.tail(*x),(b**3,3*a*b*b,3*a*a*b,a**3)


def derivative_sign(e,x,i):
    """Exact homogeneous numerator of the derivative (up to positive factors)."""
    c=e.c
    y=x[i]; z=(y[0]+1,*y[1:]); d=sum(y)
    keys, before=knots(c,y); ks, after=knots(c,z)
    B=d**3 if d>=c.D else 1; A=(d+1)**3 if d+1>=c.D else 1
    diff=Counter(dict(zip(ks,(v*B for v in after))))
    diff.subtract(dict(zip(keys,(v*A for v in before))))
    entries=[tuple(zip(*knots(c,s))) if j!=i else tuple(diff.items()) for j,s in enumerate(x)]
    val=0
    for (u,a),(v,b),(w,c) in product(*entries):
        if a and b and c: val+=a*b*c*e.weight('hV',(u,v,w))
    return (val>0)-(val<0)


def signs(e,path):
    x=tuple((d-k-1,k,1) for d,k in path)
    if not e.live('hV',x):return ('dead',)*4
    out=[]
    for port,(i,n) in enumerate(zip((0,0,1,2),(1,path[0][1]-1,path[1][1]-1,path[2][1]))):
        if not n:out.append('absent')
        elif route(e,x,port):out.append('routed')
        else:out.append({-1:'negative',0:'zero',1:'positive'}[e.sign*derivative_sign(e,x,i)])
    return tuple(out)


def verify_atom_split(E):
    paths=(((6,2),(4,2),(3,1)),((10,3),(9,2),(9,2)),((40,5),(2,1),(13,2)))
    for x in paths:
        counts=(1,x[0][1]-1,x[1][1]-1,x[2][1]); rates=tuple(tuple(Q(1,j+2) for j in range(n)) for n in counts)
        for e in E:
            a,*coeff=H.row(e,x)
            baseline=a+sum(n*min(0,c/2) for n,c in zip(counts,coeff))
            assert baseline>=0
            score=baseline+sum(c*p-min(0,c/2) for c,rs in zip(coeff,rates) for p in rs)
            assert score==a+sum(c*sum(rs) for c,rs in zip(coeff,rates))
    # The two legal birth orders (side, path) and (path, side) reach one
    # registration apiece; no promotion can label the same child twice.
    for births in ('sp','ps'):
        prefixes={(None,None)}
        for j in range(2):
            prefixes|={(j,m) for s,m in prefixes if s is None}|{(s,j) for s,m in prefixes if m is None}
        assert len([1 for s,m in prefixes if None not in (s,m) and births[s]=='s' and births[m]=='p'])==1


@cache
def score(E,path,objective):
    a,b=(H.row(e,path) for e in E)
    n=(1,path[0][1]-1,path[1][1]-1,path[2][1])
    if objective!='reward':
        aa,bb=(sum(min(0,s*m/2) for s,m in zip(c[1:],n)) for c in (a,b))
        if objective=='excess':
            a,b=(-aa,*a[1:]),(-bb,*b[1:])
        elif objective in ('apex-demand','tip-demand','middle-demand','parent-demand'):
            j=('parent-demand','apex-demand','middle-demand','tip-demand').index(objective)
            a,b=((abs(c[j+1])*n[j]/2,Q(0),Q(0),Q(0),Q(0)) for c in (a,b))
        else:
            a,b=((sum(abs(s)*m/2 for s,m in zip(c[1:],n)),Q(0),Q(0),Q(0),Q(0)) for c in (a,b))
    return tuple((p,q,p-q) for p,q in zip(a,b))


@cache
def rate(q):
    assert SCALE%q.denominator==0
    return q.numerator*(SCALE//q.denominator)


def control(values,edges):
    den=lcm(*(q.denominator for row in values for q in row))
    base=rate(H.LO); swing=rate(H.HI)-base
    nums=[v.numerator*(den//v.denominator)*SCALE*GRID for v in values[0] for _ in (0,1)]
    choices=[[(0,swing)] for _ in nums]
    cs=[ZERO]+[tuple(v.numerator*(den//v.denominator)*GRID for v in row for _ in (0,1)) for row in values[1:]]
    for a,b,successor,port in edges:
        A=rate(a); B=rate(b)-A
        if not (A or B): continue
        base+=A
        for j,(s,v) in enumerate(zip(cs[port],successor)):
            f=s+v*den
            nums[j]+=A*f
            if B:choices[j].append((f,B))
    out=[]
    for j,(num,opts) in enumerate(zip(nums,choices)):
        total=base
        for f,b in sorted(opts,reverse=bool(j%2)):
            if (f*total<=num if j%2 else f*total>=num):break
            num+=f*b;total+=b
        d=total*den
        out.append(-((-num)//d) if j%2 else num//d)
    return tuple(out)


@cache
def atom_score(E,x,role):
    a,b=(H.row(e,x)[role+2] for e in E)
    # The atom's sole existing send is paid as an affine edge credit.
    # Subtract the corresponding one-donor directional floor only.
    z,w=max(Q(0),-a)/2,max(Q(0),-b)/2
    return ((z,w,z-w), *((a,b,a-b) if i==role else Z for i in range(-1,3)))


def outside_frozen(size, cores):
    """Complement of the ead947b hV SIZE/CORE reserve, before interior checks."""
    return size > 20 or (size > 12 and cores >= 5)


def proven_outside(path, atom):
    """Exact predicate on guaranteed lower counts; no inferred full fringe."""
    n = sum(k for _, k in path)-2
    size = sum(d-k for d,k in path)+2*n
    if atom is not None and atom[0]:
        # Increase the selected side branch's two-vertex minimum using its
        # distinct own children; each was born with two vertices.
        d,k=atom
        size+=d+k-2
        n+=k
    return outside_frozen(size,3+n)


def enclosure(E,caps,ks,minimum,track,objective,disjoint,atom_role=None,atom_cap=2,frozen=False,host_min=2,parent_max=Q(1,2)):
    zs=(Z,)*5

    @cache
    def F(x,v=0,b=0,atom=None):
        length=len(x)
        if any(d>caps[i] or k>ks[i] for i,(d,k) in enumerate(x)) or (length<3 and x[-1][1]==ks[length-1]):
            return ZERO
        values=zs
        if (length==3 and x[0][0]>=host_min and max(d for d,k in x)>=minimum
                and (not disjoint or H.paid(x,True)) and (not frozen or proven_outside(x,atom))):
            if atom_role is None:
                values=score(E,x,objective)
            elif atom is not None or atom_role==-1:
                values=atom_score(E,x,atom_role)
        edges=[]
        for i,(d,k) in enumerate(x):
            l=d-k-1
            g=F(x[:i]+((d+1,k),)+x[i+1:],v,b,atom)
            fixed=Q(l)+sum((Q(1,x[j][0]) for j in (i-1,i+1) if 0<=j<length),Q(0))
            edges.append((fixed,fixed,g,0))
            if i==0:edges.append((Q(0),parent_max,g,1))
            n=k-(i+1<length)
            if i==track:
                exact=Q(v+b,2); upper=exact+Q(n-v-b,3)
                assert n>=v+b
            else:exact,upper=Q(0),Q(n,2)
            if i==atom_role and atom is not None:
                assert n>=1
                dd,_=atom
                send=Q(1,dd) if dd else Q(0)
                edges.append((send,send if dd else Q(1,atom_cap+1),g,i+2))
                upper-=Q(1,2)
            edges.append((exact,upper,g,i+2 if atom_role is None else 0))
            if l:
                post=x[:i]+((d,k+1),)+x[i+1:]
                continuation=F(post,v+int(i==track),b,atom)
                if i==atom_role and atom is None:
                    # Separate new atom branch: this newborn is a side child.
                    new=F(post,v,b,H.VIRGIN)
                    continuation=tuple(a+c for a,c in zip(continuation,new))
                if i==length-1 and length<3:
                    # In the extending prefix this promotion's newborn is
                    # the selected child, never a side donor.
                    extension=F(post+(H.VIRGIN,),v,b,atom)
                    continuation=tuple(a+c for a,c in zip(continuation,extension))
                r=Q(l,d);edges.append((r,r,continuation,0))
        if v:
            g=v*(1+Q(1,x[track][0])); r=Q(v,2)
            edges.extend(((g,g,F(x,v-1,b),0),(r,r,F(x,v-1,b+1),0)))
        if b:
            g=Q(b,x[track][0]);edges.append((g,g+Q(b,2),F(x,v,b-1),0))
        if atom is not None and atom[0]:
            dd,kk=atom; l=dd-kk-1
            g=l+Q(1,x[atom_role][0])
            after=(dd+1,kk) if dd<atom_cap else (0,0)
            edges.append((g,g+Q(kk,2),F(x,v,b,after),0))
            if l:
                r=Q(l,dd)
                edges.append((r,r,F(x,v,b,(dd,kk+1)),0))
        return control(values,edges)

    start=F((H.VIRGIN,)); n=F.cache_info().currsize;F.cache_clear()
    return [(floor(min(start[j]*H.LO,start[j]*H.HI)*10**16/GRID),
             ceil(max(start[j+1]*H.LO,start[j+1]*H.HI)*10**16/GRID)) for j in (0,2,4)],n


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lower',type=Path);p.add_argument('upper',type=Path)
    p.add_argument('--caps',nargs=3,type=int,default=(18,7,7));p.add_argument('--arity',nargs=3,type=int,default=(3,2,2));
    p.add_argument('--old',type=int,default=8);p.add_argument('--track',type=int,nargs='+',choices=(0,1,2,3),default=[0,1],
                   help='0 free, 1 tip, 2 apex, 3 middle')
    p.add_argument('--objective',nargs='+',choices=('reward','excess','width','apex-demand','tip-demand','middle-demand','parent-demand'),default=['reward'])
    p.add_argument('--disjoint',action='store_true');p.add_argument('--sign-only',nargs=6,type=int,metavar='N')
    p.add_argument('--ladder',type=int,nargs='+',help='vary common cap (keep each coordinate no larger than --caps)')
    p.add_argument('--atoms',nargs='+',choices=('parent','apex','middle','tip'));p.add_argument('--atom-cap',type=int,choices=(2,3,4),default=2)
    p.add_argument('--frozen',action='store_true',help='exclude ead947b size/core reserve using proved focal-fringe counts')
    p.add_argument('--host-min',type=int,default=2,help='minimum paid apex degree; use disjoint bands')
    a=p.parse_args();assert all(2<=d<=128 for d in a.caps) and all(1<=k<=17 for k in a.arity)
    for name,f in zip(('lower','upper'),(a.lower,a.upper)):
        with f.open('rb') as s:digest=hashlib.file_digest(s,'sha256').hexdigest()
        if H.PAYLOAD_SHA[name]!=digest:p.error(f'{name} hash mismatch')
    E=tuple(H.Local(H.Cert(f)) for f in (a.lower,a.upper))
    verify_atom_split(E)
    if a.sign_only:
        path=tuple(zip(a.sign_only[::2],a.sign_only[1::2]));assert all(0<=k<d for d,k in path) and all(k for d,k in path[:2])
        for e in E:
            result=signs(e,path);print(dict(zip(SLOPES,result)),flush=True)
            if e.live('hV',tuple((d-k-1,k,1) for d,k in path)):
                ss=H.row(e,path)[1:]
                assert all(z==0 if r in ('absent','routed','zero') else (z>0)==(r=='positive') and z!=0 for r,z in zip(result,ss))
                print('exact slopes',*(str(s) for s in ss),flush=True)
        return
    for cap in a.ladder or [max(a.caps)]:
        d=tuple(min(cap,v) for v in a.caps)
        if a.atoms:
            for name in a.atoms:
                role=('parent','apex','middle','tip').index(name)-1
                result,n=enclosure(E,d,a.arity,a.old,-1,'excess',a.disjoint,role,a.atom_cap,a.frozen,a.host_min,
                                   Q(1,3) if role==-1 else Q(1,2))
                if role==-1:
                    # Only retain the Palm births with mother D≥3. At D=2
                    # a birth requires virgin stock; its age hazard is at
                    # least 3/2 and a virgin never returns. Hence the mass
                    # of D=2 mothers is at most 1/[2(LO+3/2)].
                    w=1-Q(1,2)/(H.LO+Q(3,2)); assert w>0
                    assert result[0][0]>=0 and result[1][0]>=0
                    print('parent D_birth>=3 unconditional floors, outward 1e-16 ticks',*(floor(w*l) for l,_ in result[:2]),flush=True)
                print(d,a.arity,'old min',a.old,'atom',name,'degree cap',a.atom_cap,'disjoint',a.disjoint,
                      'ead947b complement',a.frozen,'min apex',a.host_min,'states',n,'outward 1e-16 ticks',result,flush=True)
            continue
        for t in a.track:
            for o in a.objective:
                result,n=enclosure(E,d,a.arity,a.old,(-1,2,0,1)[t],o,a.disjoint,frozen=a.frozen,host_min=a.host_min)
                print(d,a.arity,'old min',a.old,'memory',t,o,'disjoint',a.disjoint,'states',n,'outward 1e-16 ticks',result,flush=True)


if __name__=='__main__':main()
