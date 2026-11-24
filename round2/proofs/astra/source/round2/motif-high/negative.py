#!/usr/bin/env python3
"""Check both frozen endpoints using private copies. No witness is opened for writing."""
import argparse
from fractions import Fraction as Q
import hashlib
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile

from oracle import Certificate, HEADER

NAME = ('V','C','P','F','Y','T')
TARGETS = {'upper': (56581165,100000000), 'lower': (5656453,10000000)}
EXPECTED = {'upper': dict(zip(NAME,(167,125,35,31,64,87))),
            'lower': dict(zip(NAME,(5647,904,705,985,1723,50)))}


def run(args, expect, contains):
    out = subprocess.run(list(map(str,args)),text=True,capture_output=True)
    msg = out.stdout + out.stderr
    if out.returncode != expect or contains not in msg:
        raise RuntimeError(f"unexpected exit {out.returncode}, expected {expect}: {msg}")
    return msg.strip()


def put(path, offset, n):
    with open(path,'r+b') as f:
        f.seek(HEADER.size+8*offset)
        f.write(struct.pack('<q',n))


def mutation_indices(c):
    g3, g4 = c.levels[3],c.levels[4]
    u, v, b = (0,1),(0,2),(0,3)
    # Test orientations, near-knot successors, INF, T1, overflow, and repeat
    # states on DISTINCT children. Physical internal slots must pass support.
    for f,states,row in (
        ('V',(None,u,None),(g3.inf+2,g3.id[u],g3.id[(g3.d0-2,0)])),
        ('C',(v,None,None),(g3.id[v],g3.inf,g3.inf+1)),
        ('P',(None,u,u,None),(g4.inf+g4.kt+2,g4.id[u],g4.id[u],g4.id[(g4.d0-2,0)])),
        ('F',(v,None,u,None),(g4.id[v],g4.inf+1,g4.id[u],g4.inf)),
        ('Y',(None,v,None,None),(g4.inf+2,g4.id[v],g4.inf+1,g4.inf+1)),
        ('T',(b,None,None,None),(g4.id[b],g4.inf,g4.inf,g4.inf+g4.kt+2))):
        a = c.active(f,states)
        assert a >= 0
        off = c.offset(f,a,row)
        yield f,c.sections[f][0]+off,c.tables[f][off]


def tail_control(path, fixture, checker):
    shutil.copyfile(fixture,path)
    c = Certificate(fixture)
    g = c.levels[3]
    u,up,w = (0,1),(1,1),(1,0)
    a = c.active('V',(None,up,None))
    off = c.offset('V',a,(g.inf,g.id[up],g.id[w]))
    put(path,c.sections['V'][0]+off,c.sign*c.S)
    rows = [f'V {i} {g.id[u]} {g.id[w]}' for i in (g.inf,g.inf+2)]
    out = subprocess.run([str(checker),str(path),'--query'],input='\n'.join(rows)+'\n',
                         text=True,capture_output=True,check=True)
    inf,knot = map(Q,out.stdout.splitlines())
    # D_parent=2D0, k_parent=1, u=(0,1), w=(1,0). The only nonzero
    # signed jump is u gain, of amplitude S(1-t), with t=1/2.
    d0 = g.d0
    physical = (Q(1,2*d0)+Q(1,2))*Q(c.S,2)
    good = (inf+knot)/2
    bad = inf/2  # At T1, naive same-corner credit uses the zero T1 jump.
    assert inf == Q(c.S,2) and knot == Q(c.S,d0)
    assert physical-bad == good-physical == Q(c.S,4*d0)
    print(f"quadratic-tail negative control: naive shortfall={physical-bad} units; checked surplus={good-physical}",flush=True)


def digest(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(2**20),b''): h.update(block)
    return h.hexdigest()


def run_side(side,path,checker,scratch):
    c = Certificate(path)
    p,q = TARGETS[side]
    before = digest(path)
    print(f"{side} sha256={before} target={p}/{q}",flush=True)
    shutil.copyfile(path,scratch)
    # Frozen-file maxima are independently re-read from the unmodified source.
    for f,ofs,old in mutation_indices(c):
        base = [checker,path,'--pq',p,q,'--family',f]
        baseline = run(base,0,'max_ceil=')
        peak = int(baseline.split('max_ceil=')[1].split()[0])
        if peak != EXPECTED[side][f]: raise RuntimeError(f"frozen {side} {f} maximum changed to {peak}")
        put(scratch,ofs,old-c.sign*c.S)
        try:
            mutated = run([checker,scratch,'--pq',p,q,'--family',f,'--budget',peak],2,'exceeds budget')
            new_max = int(mutated.split('max_ceil=')[1].split()[0])
            print(f"family={f} baseline={peak} planted={new_max} exit=2",flush=True)
        finally: put(scratch,ofs,old)
    # One independent source-envelope and one arithmetic preflight mutation.
    ofs = c.sections['eyv1'][0]; old = c.data[ofs]
    put(scratch,ofs,-1)
    run([checker,scratch,'--family','Y'],1,'negative source envelope radius')
    put(scratch,ofs,old)
    print("negative cross-level radius rejected exit=1",flush=True)
    ofs = c.sections['T'][0]; old = c.data[ofs]
    put(scratch,ofs,2**62)
    run([checker,scratch,'--family','T'],1,'outside arithmetic range')
    put(scratch,ofs,old)
    print("integer outside proven arithmetic range rejected exit=1",flush=True)
    for plant,text in [('transfer','centered transfer mismatch'),('birth','newborn booking mismatch'),
                       ('credit','physical row exceeds exact corner majorant')]:
        cmd = ['python3','-B',Path(__file__).with_name('physical_check.py'),path,
               '--check',checker,'--pq',p,q,'--nmax',4,'--plant',plant]
        msg = run(cmd,2,text)
        print(f"{plant}: {msg.splitlines()[-1]} exit=2",flush=True)
    assert digest(scratch) == before == digest(path)
    print("restored private copy and original SHA-256 agree",flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=Path('round2/outputs'))
    args = ap.parse_args()
    checker = args.out/'mh-check'
    with tempfile.TemporaryDirectory(prefix='private-mutation-',dir=Path(__file__).parent) as d:
        scratch = Path(d)/'private.bin'
        tail_control(scratch,args.out/'mc-fixture.bin',checker)
        for side in TARGETS: run_side(side,args.out/f'mc-{side}.bin',checker,scratch)
    print('PASS',flush=True)


if __name__ == '__main__': main()
