#!/usr/bin/env python3
"""Rebuild both signed witnesses, check all rows, then compute the physical interval."""
import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S = 2**40
LOW = Fraction(113, 200)
DEGREE = "192,48,16,9,6"
CORES = "30,14,9,6,5"
PAIR_DEGREE = "1399,1024,16,5,4"
PAIR_CORES = "48,24,8,5,4"


def require(condition, reason):
    if not condition:
        raise ArithmeticError(reason)


def capture(cmd):
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)
    return proc.stdout


def match_quantity(expr, text):
    tail = re.search(expr, text)
    require(tail is not None, "missing fixed-point quantity")
    value, scale = map(int, tail.groups())
    require(value >= 0 and scale == S, "invalid value or units")
    return value


def after_base(base, snapshots, work):
    exe = work/"reserve"
    subprocess.run(["c++","-std=c++17","-O3","-Wall","-Wextra","-fno-fast-math",
                    "-fsanitize=undefined","-fno-sanitize-recover=all",HERE/"five_reserve.cpp","-o",exe],
                   cwd=ROOT, check=True)
    ends = {}
    for side in ("lower","upper"):
        frozen = json.loads((base/f"{side}-verified.json").read_text())
        physical = snapshots/f"mc-{side}.bin"
        p, q = frozen["trial"]
        lam = Fraction(p,q)
        e = frozen["error_per_core"]
        require(0 < frozen["min_weight"] < frozen["max_weight"] and e >= 0, "raw comparison")
        v = []
        for suffix in ([DEGREE,CORES], [PAIR_DEGREE,PAIR_CORES,"1133","2000","parent"],
                       [PAIR_DEGREE,PAIR_CORES,"1133","2000","child"]):
            out = capture([exe,physical,*suffix])
            print(out.rstrip(), flush=True)
            v.append(match_quantity(r"newborn reserve (\d+)/(\d+) at discount 1133/2000\s*$",out))
        out = capture([sys.executable,"-B",HERE/"lifetime.py",base/f"mc-mbar-{side}.i64",
                       base/f"{side}-comparison.json"])
        print(out.rstrip(), flush=True)
        y = match_quantity(r"newborn lifetime ceiling (\d+)/(\d+)\s*$",out)
        require(y > 0 and sum(v) > e/LOW, "the signed improvement does not follow")
        improvement = (sum(v)-e/LOW)/y
        ends[side] = lam + (improvement if side == "lower" else -improvement)
        print(side, "fraction", ends[side], flush=True)
        print(capture([sys.executable,"-B",HERE/"check_five.py",exe,physical]).rstrip(), flush=True)
    decimal_scale = 10**10
    lo, hi = ends["lower"], ends["upper"]
    require(Fraction(565615,10**6) < lo <= hi < Fraction(565809,10**6), "location of interval")
    lower_tick = lo.numerator*decimal_scale//lo.denominator
    upper_tick = -(-hi.numerator*decimal_scale//hi.denominator)
    require((lower_tick, upper_tick) == (5656523639,5658045676), "changed reported endpoints")
    print(f"Outward physical interval [0.{lower_tick:010d}, 0.{upper_tick:010d}].")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshots",type=Path,help="verify the two frozen binaries from this directory")
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError("run with Python assertions enabled")
    if hasattr(os,"sched_getaffinity"):
        os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    output = ROOT/"round2/outputs"
    output.mkdir(exist_ok=True)
    with TemporaryDirectory(prefix="occupation-",dir=output) as tmp:
        work = Path(tmp)
        base = work/"base"
        arguments = [sys.executable,"-B",HERE.parent/"motif-cert/replay.py","--work",base]
        snapshots = args.snapshots.resolve() if args.snapshots else base
        if args.snapshots:
            arguments.extend(["--snapshots",snapshots])
        subprocess.run(arguments,check=True,cwd=ROOT)
        after_base(base,snapshots,work)


if __name__ == "__main__":
    main()
