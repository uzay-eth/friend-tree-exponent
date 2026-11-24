#!/usr/bin/env python3
"""Rebind private raw mutations, then require numerical rejection by stock verifier."""
import argparse
import hashlib
import os
from pathlib import Path
import shlex
import struct
import subprocess
import time

MAX_W = 16 << 40
MAX_E = 32 << 40


def digest(path):
    with path.open('rb') as fp:
        return hashlib.file_digest(fp, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('snapshot', type=Path)
    p.add_argument('verifier', type=Path)
    p.add_argument('offsets', type=Path)
    p.add_argument('scratch', type=Path)
    a = p.parse_args()
    size, *offsets = map(int, subprocess.check_output([a.offsets, a.snapshot], text=True).split())
    assert size == a.snapshot.stat().st_size
    with a.snapshot.open('rb') as f:
        vals = []
        for offset in offsets:
            f.seek(offset)
            vals.append(struct.unpack('<q', f.read(8))[0])
    old, up, mid, rad = vals
    assert abs(up-old-mid)<=rad
    print('offsets bytes',offsets,'raw old/up/mid/rad',vals,'base lhs/rhs',abs(up-old-mid),rad,flush=True)
    for label, idx, cap in [('g',0,MAX_W),('env23-mid',2,MAX_E)]:
        def lhs(v):
            return abs(up-v-mid) if idx==0 else abs(up-old-v)
        new=max((-cap,cap),key=lhs)
        assert lhs(new)>rad and new!=vals[idx]
        dst=a.scratch/f'{label}.bin'
        cmd=['cp','--reflink=auto','--',str(a.snapshot),str(dst)]
        print(shlex.join(cmd),flush=True)
        try:
            subprocess.run(cmd,check=True)
            assert os.stat(dst).st_ino != os.stat(a.snapshot).st_ino
            with dst.open('r+b') as f:
                f.seek(offsets[idx]); f.write(struct.pack('<q',new))
            # This is the NEW expected digest. The numerical verifier has no
            # input digest option; its caller normally checks this binding.
            bound=digest(dst)
            print(label,'offset',offsets[idx],'old',vals[idx],'new',new,'new sha256',bound,
                  'lhs',lhs(new),'rhs',rad,flush=True)
            cmd=[str(a.verifier),str(dst),'--threads','1']
            print(shlex.join(cmd),flush=True)
            now=time.monotonic()
            run=subprocess.run(cmd,capture_output=True,text=True,timeout=70)
            print('exit',run.returncode,'seconds',round(time.monotonic()-now,2),'stdout',repr(run.stdout),
                  'stderr',repr(run.stderr),flush=True)
            assert run.returncode==2 and run.stderr.strip()=='ERROR: env23 verification failed' and not run.stdout
        finally:
            dst.unlink(missing_ok=True)


if __name__=='__main__':main()
