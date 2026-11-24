"""Direct physical reward and Bellman probes for the Gold2 lower singleton."""

import ast
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from json import loads
from mmap import mmap, ACCESS_READ
from pathlib import Path
from struct import unpack_from
from types import SimpleNamespace

R=Path(__file__).resolve().parents[1]
P=R/'occupation-quintet-gold2/lifetime/check.py'
sp=spec_from_file_location('lifetime',P)
a=module_from_spec(sp); sp.loader.exec_module(a)


def reward(t,d,c,receive):
    def f(d,c):
        if d<4000: return t.f(d,c)
        q=Q(d-c,d)
        return t.B*q+t.C*q*q+t.E*q/d
    L=d-c
    x=f(d,c)
    gain=f(d+1,c)-x
    bud=0
    if L:
        if d<252:
            edge=t.g((d,c+1) if c<22 else (d,'over'),(2,1))
        else:
            edge=(1-Q(252,d))*t.g(('inf',0),(2,1))+Q(252,d)*t.g(('tail',min(c+1,23)),(2,1))
        bud=f(d,c+1)-x+t.f(2,1)+edge
    if d<252 and c<=22:
        # Stored common reciprocal per core port, NOT the physical receive.
        center=Q(unpack_from('<q',t.raw,512+8*t.ids[d,c])[0],2**30)
        receive=c*center
    assert t.s==-1
    return 55+(L+receive)*gain+Q(L,d)*bud-t.lam*x


def main():
    gold=loads(a.tail.RECEIPT.read_text()); rec=gold['lower']
    lo,hi=(Q(gold[k]['certified_exact']) for k in ('lower','upper'))
    with (a.tail.BASE/'quintet_lower.bin').open('rb') as fp,mmap(fp.fileno(),0,access=ACCESS_READ) as raw:
        t=a.tail.Snapshot(raw,rec)
        pos=neg=zero=0
        for d,c in product((2,19,251,252,253,429,430,431,3999,4000,7999,8000),(1,2,3,17,18,19,22,23)):
            if c>d: continue
            A,B=a.row(t,d,c,55)
            for r in (Q(0),Q(c,5),Q(c,2)):
                assert reward(t,d,c,r)==A+B*r
            pos+=B>0; neg+=B<0; zero+=B==0
        print('PASS direct charge vs writer intercept/slope at seams:',pos,neg,zero,flush=True)

        W={}
        a.price(t,7999,30,hi,55,False,(7999,W))
        W[8000]=a.suffix(t,8000,30,hi,55,False)
        broad=0
        for d,c in product((251,252,3999,4000,7999),(1,18,22,23,30)):
            for r in (Q(c,7),Q(c,3),Q(c,2)):
                L=d-c
                value=(a.GRID*reward(t,d,c,r)+(L+r)*W[d+1][c]+Q(L,d)*W[d][c+1])/(hi+L+r+Q(L,d))
                assert W[d][c]<=value
                broad+=1
        print(f'PASS {broad} full broad quotients at layout, analytic and terminal seams',flush=True)
        m=a.neighbour_algorithm(lo,hi,gold)
        seen=[]
        def probe(d,c,mom,tip,read,v,lookup):
            if (d,c,mom,tip) not in ((430,18,2,(3,1)),(430,18,-3,(4,0)),(430,1,2,(0,0)),(252,18,2,(3,2)),(251,18,2,(3,1))): return
            j,i=tip; L=d-c
            # One concrete physical donor completion: all old untracked core
            # donors degree two, unknown mother degree three, unknown child four.
            # The old mother has her own degree-two parent and otherwise leaves.
            # The exact latest child has i core slots, one of which is focal.
            rr=Q(c-1-int(j>0),2)+Q(1,abs(mom))+(Q(1,j) if j else 0)
            nxt=read(d+1,c,mom,tip); bud=read(d,c+1,mom,(2,1))
            events=[(L+rr,nxt),(Q(L,d),bud)]
            if mom>0: events.append((Q(1,d)+Q(1,2),read(d,c,-3,tip)))
            if i:
                events.append((j-i+Q(1,d)+Q(i-1,2),read(d,c,mom,(4,0))))
                if j>i: events.append((Q(j-i,j),read(d,c,mom,(j,i+1))))
            physical=(a.GRID*reward(t,d,c,rr)+sum(rate*w for rate,w in events))/(hi+sum(rate for rate,w in events))
            assert v<=physical and lookup==max(v,W[d][c])
            # The fallback is a separate global subsolution; it need not pass
            # this row with the discarded, smaller conditional proposal.
            seen.append((d,c,mom,tip,a.row(t,d,c,55)[1],physical-v))

        # Observe the solver's actual candidate before applying its fallback.
        ns=m.conditional.__globals__.copy(); ns['probe']=probe
        tree=ast.parse((R/'occupation-quintet/lifetime/conditioned/replay.py').read_bytes())
        f=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='conditional')
        loop=next(n for n in f.body if isinstance(n,ast.For))
        for _ in range(3): loop=next(n for n in loop.body if isinstance(n,ast.For))
        loop.body.extend(ast.parse('probe(d,c,m,t,read,v,lookup[d,c,m,t])').body)
        exec(compile(ast.fix_missing_locations(ast.Module([f],[])),'<instrumented>', 'exec'),ns)
        def row(d,c):
            x,y=a.row(t,d,c,55);return x-55,y
        roots,n,corners=ns['conditional'](SimpleNamespace(side='lower',row_ticks=row),W,430,18,2,3,False)
        assert len(seen)==5 and any(s[4]!=0 for s in seen)
        z=m.birth_mix(roots,False)
        assert z==Q(106709227048431615194065682065030463,290238887678415021500000000000)
        assert a.fmt(lo*z,False,9)=='207990.215716345'
        print('PASS five direct physical Bellman quotients; d,c,mother,child,reward-slope:')
        for d,c,mom,tip,slope,_ in seen: print(d,c,mom,tip,slope)
        print('PASS exact conditioned mixed value and outward 207990.215716345 ticks/core',n,corners)


if __name__=='__main__': main()
