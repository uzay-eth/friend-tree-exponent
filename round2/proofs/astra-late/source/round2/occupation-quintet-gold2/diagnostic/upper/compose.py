#!/usr/bin/env python3
"""Exact upper-only extension of the frozen D24 uniform-prior ledger."""
import ast
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
src=HERE.parents[1]/'host/h3-joint/receipt.json'
raw=src.read_bytes()
assert sha256(raw).hexdigest()=='ae6f5cde2ef710b4196c3c445a3786580e5ccc9dfd9206b2962f60628866c6ed'
p=json.loads(raw)['results']['upper']
log=(HERE/'mother80-480.log').read_text()
assert 'exit=0' in log and 'literal labelled H3 rows 42' in log
assert 'min_mother=80, mother=480, mother_stock=12, middle=2, middle_stock=2, young=12, young_stock=3' in log
assert p['payload_sha256'] in log
floor,ticks=map(Q,re.search(r'EXACT upper physical H3 complement (\S+) ticks/core (\S+)',log).groups())
assert ticks==floor*2**40 and ticks>0
M=Q(p['d24_denominator'][1]);t=p['uniform_trial']
assert t['prior']==p['original_prior'] and not t['restricted'] and t['proved']
old=Q(t['point'])-Q(t['margin'])/M
assert old==Q(p['uniform_endpoint'])
new=Q(t['point'])-(Q(t['margin'])+ticks)/M
q=Q(241,426); scale=10**15
assert q<new<Q(p['final_endpoint'])

def outward(endpoint):
    out=(endpoint*scale).__ceil__()
    return f'{out//scale}.{out%scale:015d}'

def frontier_ticks(name, *, old, last, stock, min_stock, min_middle):
    path=HERE.parents[1]/'post-fy/upper-frontier'/name
    starts, data, status=path.read_text().splitlines()
    assert starts.startswith('original upper Gold2 ') and status.endswith(' exit=0')
    assert ast.literal_eval(starts.removeprefix('original upper Gold2 '))==dict(
        host=5, stock=4, receiver=last, rstock=stock, donor=7, tstock=3,
        old=old, last=last, mother_min_stock=min_stock, middle_min=min_middle,
        child_min=2, ceiling=False, fixtures=False, reserve=False)
    units, value=map(Q,re.fullmatch(r'normalized (\S+) \S+ ticks (\S+) \S+ states \d+ paid \d+ CPU seconds \S+',data).groups())
    assert units*2**40==value and value>0
    return value

# The low band is beyond the frozen arity cap of 8. The high band is
# beyond the frozen degree cap of 79 and excludes our middle D=2.
stock=frontier_ticks('stock-floor.log',old=19,last=79,stock=12,min_stock=9,min_middle=2)
middle=frontier_ticks('middle3-floor.log',old=80,last=120,stock=10,min_stock=2,min_middle=3)
assert stock+middle==Q('13266167036667732636019911/364805903750000000000')
joint=Q(t['point'])-(Q(t['margin'])+ticks+stock+middle)/M
assert joint==new-(stock+middle)/M and q<joint<new
print('floor ticks/core',ticks,float(ticks))
print('exact upper endpoint',new)
print('outward upper endpoint',outward(new))
print('certified inward step relative to previous fixed endpoint',float(Q(p['final_endpoint'])-new))
print('remaining gap to 241/426',float(new-q))
print('additional sufficient fixed-denominator numerator (normalized/core)',float((new-q)*M/2**40))
print('corrected frontier ticks/core',stock+middle,float(stock+middle))
print('joint floor ticks/core',ticks+stock+middle,float(ticks+stock+middle))
print('joint exact upper endpoint',joint)
print('joint outward upper endpoint',outward(joint))
print('joint remaining gap to 241/426',float(joint-q))
print('joint additional sufficient fixed-denominator numerator (normalized/core)',float((joint-q)*M/2**40))
