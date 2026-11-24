#!/usr/bin/env python3
"""Independent expected row and envelope counters from only a Claim 6 header."""
from itertools import product
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import struct
import sys

with Path(sys.argv[1]).open('rb') as f:
    h=struct.unpack('<Q63q',f.read(512))
assert h[0]==0x4654325F55355631
(df,d2,k2,d3m,k3m,d3e,k3e,d4i,k4i,sp4,ds3,ks3,d4e,k4e,
 d5i,k5i,sp5,d5o,k5o,d5c,k5c,sb5,ds4,ks4,d5e,k5e)=h[10:36]


def inner(dcut,kcut,minimum):
    return [(d,c) for d in range(2,dcut) for c in range(minimum,min(d,kcut)+1)]


def endpoint(d,k):
    return len(inner(d,k,1))+max(0,d-k-1)+k+2


g2,g3,g4,g5=[endpoint(d,k) for d,k in ((d2,k2),(d3e,k3e),(d4e,k4e),(d5e,k5e))]
m3=len(inner(d3m,k3m,2)); s3=len(inner(ds3,ks3,3)); s4=len(inner(ds4,ks4,4))
p4=[(a,b) for a,b in product(inner(d4i,k4i,2),repeat=2) if a[0]+b[0]<=sp4]
p5=[(a,b,c) for a,b,c in product(inner(d5i,k5i,2),repeat=3) if a[0]+b[0]+c[0]<=sp5]
b5=[(a,b) for a,b in product(inner(d5o,k5o,3),inner(d5c,k5c,2)) if a[0]+b[0]<=sb5]
fin=sum(d+1 for d in range(2,df)); cent=len(inner(d2,k2,1))
rows={'f_in':cent,'f_tab':fin-cent,'g':comb(g2+1,2)-1,'h3':m3*comb(g3+1,2),
      'p4':sum(g4*g4 if a<b else comb(g4+1,2) for a,b in p4 if a<=b),
      's3':s3*comb(g4+2,3),
      'p5':sum(g5*g5 if a<c else comb(g5+1,2) for a,b,c in p5 if a<=c),
      'b5':len(b5)*comb(g5+1,2)*g5,'s4':s4*comb(g5+3,4)}
env={'1d':m3*g2+len(p4)*g3+len(p5)*g4,
     '2d':s3*comb(g3+1,2)+len(b5)*(g4*g4+comb(g4+1,2)),
     '3d':s4*comb(g4+2,3)}
print('signed mode, nominal',h[1],h[2],h[3])
print('endpoint and block counts',(g2,g3,g4,g5),(m3,len(p4),s3,len(p5),len(b5),s4))
print(json.dumps({'rows':rows,'envelope_checks':env},sort_keys=True))
if len(sys.argv)==3:
    result=json.loads(Path(sys.argv[2]).read_bytes())
    result=result.get('lower' if h[1]==2 else 'upper',result)
    record=result.get('verify',result)
    assert record['mode']==h[1] and (record['P'],record['Q'])==h[2:4]
    assert record['rows']==rows and record['envelope_checks']==env
    print('PASS actual coverage counters')
    if 'tail' in result:
        t=result['tail']; w=record['comp_fin']; caps=record['caps']; a=record['alpha_ticks']
        low=min(w['min_W_c0'],w['min_W'],t['w_c0_min_ticks'],
                t['w_c1_min_ticks']+w['g_far_min']-a,t['w_c2_min_ticks']+w['g_far_min'],
                t['w_min_ticks']+w['g_far_min']+a)
        high=max(w['max_W_c0'],w['max_W'],t['w_max_ticks']+w['g_far_max'])
        err=max(0,max(caps['f_in'],caps['f_tab'],t['cap_ticks'])+max(0,caps['g'])+record['higher_per_core'])
        assert 0<low<=high and (low,high,err)==(result['m_min'],result['m_max'],result['err_per_core'])
        edge=Q(h[2],h[3])+(1 if h[1]==1 else -1)*Q(err,low)
        assert edge==Q(result['certified_exact'])
        print('PASS independent signed composition',low,high,err,edge)
