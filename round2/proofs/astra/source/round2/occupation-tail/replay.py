#!/usr/bin/env python3
"""Rebuild the physical witnesses and bound the additional discounted occupation."""
import argparse
from contextlib import redirect_stdout
from fractions import Fraction
from io import StringIO
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(HERE.parent/"occupation"))
from replay import after_base


def compile_cpp(source, destination):
    subprocess.run(["c++","-std=c++20","-O3","-Wall","-Wextra","-fno-fast-math",
                    "-fsanitize=undefined","-fno-sanitize-recover=all",source,
                    "-lgmpxx","-lgmp","-o",destination], check=True)


def run(args, work):
    base = work/"base"
    snapshots = args.snapshots.resolve() if args.snapshots else base
    cmd = [sys.executable,"-B",HERE.parent/"motif-cert/replay.py","--work",base]
    if args.snapshots:
        cmd.extend(["--snapshots",snapshots])
    subprocess.run(cmd,cwd=ROOT,check=True)
    baseline = StringIO()
    with redirect_stdout(baseline):
        after_base(base,snapshots,work)
    print(baseline.getvalue().rstrip(),flush=True)
    summary = re.findall(r"newborn lifetime ceiling (\d+)/1099511627776\s+"
                         r"(lower|upper) fraction ([0-9]+/[0-9]+)",baseline.getvalue())
    if [side for _,side,_ in summary]!=["lower","upper"]:
        raise ArithmeticError("missing physical baseline")

    solver = work/"affine"
    compile_cpp(HERE/"solver.cpp",solver)
    for y,side,start in summary:
        label = "1" if side=="upper" else "2"
        file = snapshots/f"mc-{side}.bin"
        command = [solver,file,label,"63","4","1400,72,56,56,24,16,24,10,7,6",
                   "40,16,10,14,9,7,8,5,4,3"]
        out = subprocess.check_output(command,cwd=ROOT,text=True)
        print(side,out.rstrip(),flush=True)
        record = re.fullmatch(r"mask 63 at 1133/2000 rows 55259 510900 21172384 10109560 "
                              r"positive sources(?: \d+){4} newborn (\d+)/1099511627776\s*",out)
        if record is None:
            raise ArithmeticError("missing additional-occupation rows")
        v = int(record.group(1))
        right = side=="upper"
        if v!=(17458365 if right else 15289028):
            raise ArithmeticError("additional occupation differs")
        end = Fraction(start)+(-1 if right else 1)*Fraction(v,int(y))
        scaled = end*10**10
        tick = -(-scaled//1) if right else scaled//1
        if tick!=(5657979941 if right else 5656581227):
            raise ArithmeticError("unexpected physical endpoint")
        print(side,"fraction",end,"outward",f"0.{tick:010d}",flush=True)
        subprocess.run([sys.executable,"-B",HERE/"physical.py",solver,file,label],cwd=ROOT,check=True)
    print("Outward physical interval [0.5656581227, 0.5657979941].",flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--snapshots",type=Path,help="read and verify an existing frozen pair")
    args = ap.parse_args()
    if not __debug__:
        raise RuntimeError("enable Python assertions")
    if hasattr(os,"sched_getaffinity"):
        os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    outputs = HERE.parent/"outputs"
    outputs.mkdir(exist_ok=True)
    with TemporaryDirectory(dir=outputs,prefix="affine-") as tmp:
        run(args,Path(tmp))


if __name__=="__main__":
    main()
