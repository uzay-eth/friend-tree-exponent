#!/usr/bin/env python3
"""Exact Gold2 cross-check using literal ordinary-tree event enumeration."""
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from os import close
from pathlib import Path
from subprocess import check_output

HERE=Path(__file__).resolve().parent
R2=HERE.parents[1]

def load(name,path):
    spec=spec_from_file_location(name,path)
    m=module_from_spec(spec);spec.loader.exec_module(m)
    return m

gold=load('gold2_wire',R2/'attacks/neon2-quintet/gold2/physical/check.py')
literal=load('literal',R2/'occupation-quintet/diagnostic/exact.py')


def main():
    examples=(
        (4000,5,2,3,2,2,3,5,65,2,253,251,252),
        (3999,6,3,3,2,11,2,253,38,2,3,2,4001),
        (4001,5,3,5,252,2,2,48,37,2,255,6,252),
        (251,5,2,4,2,2,5,12,8,2,3,2,3),
    )
    for side in ('lower','upper'):
        path=gold.ROOT/f'quintet_{side}.bin';p=gold.Gold(path)
        assert p.h['P']=={'lower':56571368,'upper':56573686}[side]
        gold.transfer(p)  # All 20 signed knots, interpolation and ordered FY cancellation.
        used=set();fy_slots=0
        for i,ds in enumerate(examples):
            f=HERE/f'fixture-gold2-{i}.tsv'
            if i==3:literal.dump(f,(2,*ds),(-1,*(j+1 for j in literal.PAR)))
            else:literal.dump(f,ds)
            counts,w,raw,ch,moves,fy,checks,n=literal.evaluate(p,f)
            out=check_output([HERE/'census',path,f],text=True).splitlines()
            seen=set()
            for line in out:
                tag,*a=line.split('\t')
                if tag=='tree':assert max(abs(Q(z)-sum(ch.values())/p.h['S']) for z in a[5:])<Q('3e-11')
                if tag=='family':
                    k=a[0];assert int(a[1])==counts[k];seen.add(k)
                    for z,t in zip(a[2:5],(w[k],raw[k],ch[k])):assert abs(Q(z)-t/p.h['S'])<Q('3e-11')
                if tag=='move':
                    k,l=(gold.old.NAMES[int(j)] for j in a[:2]);v=moves[k,l]/p.h['S']
                    assert abs(Q(a[2])-v)<Q('3e-11')
                    if v:used.add((k,l))
            assert seen==set(gold.old.NAMES)
            cpp=sum(Q(line.split('\t')[-1]) for line in out if line.startswith('bin\t1\t'))
            assert abs(cpp-fy/p.h['S'])<Q('3e-11')
            fy_slots+=checks['slots']
            print(side,i,'N',n,'counts',dict(counts),'FY',dict(checks),'total exact',sum(ch.values())/p.h['S'])
        assert len(used)==8 and fy_slots>0
        print(side,'all signed transfer channels',sorted(used),'FY actual slots',fy_slots)
        close(p.fd)


if __name__=='__main__':main()
