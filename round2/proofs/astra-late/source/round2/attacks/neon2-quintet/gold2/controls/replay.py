#!/usr/bin/env python3
"""Clone the authenticated lower snapshot, swap two env23 joins, run stock verifier once."""
from pathlib import Path
import hashlib
import os
import struct
import subprocess
import sys
import tempfile


def require(ok, what):
    if not ok:
        raise RuntimeError(what)


def memory_kib():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            return int(line.split()[1])
    raise RuntimeError('MemAvailable unavailable')


def main(base, planner, temp):
    src = base/'quintet_lower.bin'
    size = 6344902168
    require((base/'hashes.log').read_text().splitlines()==[
        'quintet_lower.bin: OK','quintet_upper.bin: OK'], 'both upstream hashes required')
    require(src.stat().st_size==size, 'source size')
    print('available_kib',memory_kib(),flush=True)
    require(memory_kib()>=12*1024**2, 'need at least 12 GiB available')
    require(len(os.sched_getaffinity(0))==1 and os.nice(0)>=15, 'require single CPU and nice >=15')
    proc = subprocess.run([planner,src],check=True,stdout=subprocess.PIPE,text=True,timeout=10)
    patches = [tuple(map(int,line.split())) for line in proc.stdout.splitlines()]
    require(patches==[(6141599992,-93711837683,-175951565211),
                      (6141600000,-175951565211,-93711837683)], 'unexpected patch plan')
    with tempfile.TemporaryDirectory(prefix='gold2-env23-',dir=temp) as d:
        dst = Path(d)/'lower-env23-swap.bin'
        subprocess.run(['cp','--reflink=always','--',src,dst],check=True,timeout=15)
        require(not os.path.samefile(src,dst) and dst.stat().st_nlink==1, 'private inode required')
        dst.chmod(0o600)
        with dst.open('r+b') as f:
            require(not os.path.samefile(src,dst) and os.fstat(f.fileno()).st_nlink==1, 'private inode required')
            for off,old,new in patches:
                f.seek(off)
                require(struct.unpack('<q',f.read(8))[0]==old,'preimage')
                f.seek(off)
                f.write(struct.pack('<q',new))
                print('offset',off,'old',old,'new',new,flush=True)
            f.flush()
            require(os.fstat(f.fileno()).st_size==size,'mutant size')
        with dst.open('rb') as f:
            digest = hashlib.file_digest(f,'sha256').hexdigest()
        print('mutant_sha256',digest,flush=True)
        require(memory_kib()>=12*1024**2, 'need at least 12 GiB before verifier')
        run = subprocess.run([base/'verify_quintet',dst,'--threads','1'],capture_output=True,text=True,timeout=150)
        print('verifier_exit',run.returncode,'stdout',repr(run.stdout),'stderr',repr(run.stderr),flush=True)
        require(run.returncode==2 and run.stdout=='' and run.stderr=='ERROR: env23 verification failed\n',
                'expected numerical env23 rejection')


if __name__=='__main__':
    main(*map(Path,sys.argv[1:]))
