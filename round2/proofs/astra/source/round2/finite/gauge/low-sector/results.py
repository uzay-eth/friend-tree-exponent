"""Recompute enclosure comparisons from the raw rounded rational log lines."""
from ast import literal_eval
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import ceil, comb, floor
from pathlib import Path
import re

OUT=Path(__file__).resolve().parents[3]/'outputs'
I=r'\[([\d.]+),([\d.]+)\]'


@cache
def log(name):
    result={}
    for line in (OUT/f'finite-low-{name}.log').read_text().splitlines():
        if m:=re.search(r'labelled multicore rows (\d+)',line):
            result['rows']=int(m[1])
        if line.startswith('label-erased low shapes by core count '):
            result['count']=literal_eval(line[line.index('{'):])
        if m:=re.match(r'wall_seconds (\S+) user_seconds (\S+) sys_seconds (\S+) max_rss_kib (\d+) exit 0',line):
            result['wall'],result['cpu'],result['sys']=map(F,m.groups()[:3])
            result['rss']=int(m[4])
        if m:=re.match(r'column (\S+) seed (\S+)',line):
            row=result[m[2],m[1]]={}
        if m:=re.match('lower endpoint '+I+' upper endpoint '+I,line):
            row['lower'],row['upper']=tuple(map(F,m.groups()[:2])),tuple(map(F,m.groups()[2:]))
        if m:=re.match('total width '+I,line):
            row['total']=tuple(map(F,m.groups()))
        if m:=re.match(r'width (\S+) '+I,line):
            row[m[1]]=tuple(map(F,m.groups()[1:]))
    assert 'rss' in result
    return result


def groups(row,which):
    return sum(v[0] for k,v in row.items() if which(k))


def count_shapes(kmax):
    nums=[]
    for k in range(1,kmax+1):
        @cache
        def multisets(size,children,i):
            if not children or not size:
                return int(not children and not size)
            if i>size or i>len(nums):
                return 0
            return sum(comb(nums[i-1]+j-1,j)*multisets(size-i*j,children-j,i+1)
                       for j in range(min(children,size//i)+1))
        nums.append(sum((7-max(2,1+m))*multisets(k-1,m,1) for m in range(min(5,k-1)+1)))
    return nums


def radius_interval(matrix):
    a,b=matrix[0]
    c,d=matrix[1]
    lo,hi=max(a,d),max(a+b,c+d,max(a,d))
    for _ in range(100):
        mid=(lo+hi)/2
        if (mid-a)*(mid-d)<b*c:
            lo=mid
        else:
            hi=mid
    return lo,hi


def asdec(q):
    return f'{float(q):.9f}'


if __name__=='__main__':
    counts=count_shapes(6)
    assert counts[:5]==[5,25,185,1530,13915]
    assert all(log('five')['count'][k]==counts[k-1] for k in range(1,6))
    print('unmarked complete low counts k=1..6',counts)
    cases=(('two-level',None),('three-timing','two-level'),('four-timing','three-timing'),
           ('root-high-three','four-timing'),('below-high-three','four-timing'),
           ('high-three','four-timing'),('five','four-timing'),('five-high-four','four-timing'))
    for name,base in cases:
        data=log(name)
        row=data['B','B']
        assert groups(row,lambda k:k not in ('lower','upper','total'))<=row['total'][1]
        fields=(name,data['rows'],float(data['wall']),float(data['cpu']),data['rss'],asdec(row['total'][0]))
        if base:
            prev=log(base)
            saved=prev['B','B']['total'][0]-row['total'][1]
            fields+=(asdec(saved),asdec(1000*saved/(data['rows']-prev['rows'])),
                     asdec(data['cpu']/(100*saved)))
        print('cost',*fields)
    old=log('three-timing')['B','B']
    four=log('four-timing')['B','B']
    assert old['low_4'][0]>F('.7388153629')
    for label,pred in (('five low',lambda k:k=='low_5'),
                       ('one high four',lambda k:k in ('one_root_4','one_below_4')),
                       ('all remaining',lambda k:k not in ('low_4','low_5','one_root_4','one_below_4','total','upper','lower'))):
        print('four-core width transfer',label,asdec(groups(four,pred)-groups(old,pred)))
    a=groups(four,lambda k:k in ('one_root_3','one_below_3'))
    b=four['lower'][0]+groups(four,lambda k:k.startswith(('low_','several_')))
    c=four['lower'][0]+groups(four,lambda k:k not in ('lower','upper','total','one_root_3','one_below_3','low_5'))
    assert a>F('1.39960671') and b>F('1.0169') and c>F('1.2033')
    print('oracle lower bounds: three-core one-high alone, all one-high and star exact, both dominant classes exact',
          *map(asdec,(a,b,c)))
    for s,left,right in (('.565615','four-left-a','four-left-b'),('.565809','four-ledger','four-ledger')):
        L=[[log(src)[seed,col]['lower'][0] for src,col in ((left,'A'),(right,'B'))] for seed in ('A','B')]
        H=[[log(src)[seed,col]['upper'][1] for src,col in ((left,'A'),(right,'B'))] for seed in ('A','B')]
        lo,hi=radius_interval(L)[0],radius_interval(H)[1]
        assert lo<1<hi
        print('no-prior full radius',s,floor(lo*10**8),ceil(hi*10**8),'units 1e-8')
        test=log('four-test-left' if s=='.565615' else 'four-test-right')
        test_hi=max(5*test['A','.2']['upper'][1],test['B','.2']['upper'][1])
        assert test_hi<hi
        print('including direct g=(1/5,1), no-prior radius',s,
              floor(lo*10**8),ceil(test_hi*10**8),'units 1e-8')
    new=log('four-ledger')['B','B']
    orig=(OUT/'finite-low-four-first.log').read_text()
    m=re.search(r'true row 1/16 \[([^,]+),([^\]]+)\]',orig)
    assert max(abs(F(m[1])-new['lower'][0]),abs(F(m[2])-new['upper'][1]))<F(2,10**12)
    print('independent backward/forward four-core endpoints agree to 2e-12')
    chamber=defaultdict(lambda:[0,0])
    total_states=0
    for name in OUT.glob('finite-low-bellman-chamber-*.log'):
        content=name.read_text()
        total_states+=int(re.search(r'part 0 of 1 states (\d+)',content)[1])
        assert 'exit 0' in content
        for kind,pos,mode,low,up in re.findall(r'(gain|extra) (root|below) (\w+) lower (\d+) upper (\d+)',content):
            for key in ((kind,pos,mode),(kind,'all',mode)):
                chamber[key][0]+=int(low)
                chamber[key][1]+=int(up)
    assert total_states==518
    for (kind,pos,mode),value in sorted(chamber.items()):
        print('chamber',kind,pos,mode,'lower and upper savings units 1e-12',*value)
    new_lo=new['lower'][0]+F(chamber['gain','all','actual'][0],10**12)
    new_hi=new['upper'][1]-F(chamber['gain','all','actual'][1],10**12)
    assert new_lo>=F('.314726046818') and new_hi<=F('3.140363448468')
    print('four-sector with actual joint chamber',new_lo,new_hi)
    audit=(OUT/'finite-low-lower-model.log').read_text()
    numbers=re.findall(r'root (\S+) upper (\d+) / 1e12',audit)
    assert len(numbers)==2
    for (tag,v),limit in zip(numbers,(F(2,3),F(1))):
        assert F(int(v),10**12)<F('.882016')*limit
    print('homogeneous lower-comparison contraction at .565615 < .882016')
