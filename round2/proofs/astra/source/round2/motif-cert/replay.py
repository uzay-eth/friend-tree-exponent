#!/usr/bin/env python3
"""Build and replay the integer witness. Only private work files are created."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import mkdtemp

from tail import verify as tail_bound

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
HIGH=HERE.parent/"motif-high"
M=json.loads((HERE/"certificate.json").read_text())
UNIT=M["S"]
SIDES=("upper","lower")


def require(test,message):
    if not test: raise ValueError(message)


def digest(path):
    with path.open("rb") as f: return hashlib.file_digest(f,"sha256").hexdigest()


def call(cmd,log=None,live=False):
    if live:
        with log.open("w") as out:
            subprocess.run(cmd,check=True,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT)
        return log.read_text()
    proc=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True)
    if log: log.write_text(proc.stdout+proc.stderr)
    if proc.returncode:
        raise ValueError(f"exit {proc.returncode}: {cmd}; {proc.stderr[-800:]} {proc.stdout[-800:]}")
    return proc.stdout


def compile_tools(work):
    for name,expected in M["source_sha256"].items():
        require(digest(HERE.parent/"motif"/name)==expected,"proposal/layout source digest mismatch: "+name)
    cxx=["c++","-std=c++20","-O3","-Wall","-Wextra"]
    units=(("export",["cc","-std=gnu11","-O3","-ffp-contract=off","-fno-fast-math",HERE/"export.c","-lm"]),
           ("envelopes",cxx+[HERE/"envelopes.cpp"]),("low",cxx+[HERE/"low.cpp"]),
           ("comparison",cxx+[HERE/"comparison.cpp"]),("high",cxx+[HIGH/"check.cpp","-lgmpxx","-lgmp"]))
    for name,cmd in units: call(cmd+["-o",work/name],work/(name+"-compile.log"))


def build(work,side,variant="full"):
    tag=side if variant=="full" else variant
    path=work/f"mc-{tag}.bin"
    require(not path.exists(),"existing output; choose a fresh rebuild directory")
    part=work/f"mc-{tag}.part"
    call([work/"export",side,variant,part],work/f"build-{tag}.log",True)
    call([work/"envelopes","build",part],work/f"build-{tag}-envelopes.log")
    call([work/"envelopes","check",part],work/f"preflight-{tag}-envelopes.log")
    if variant=="full":require(digest(part)==M[side]["sha256"],"regenerated witness digest mismatch")
    os.link(part,path);part.unlink()
    return path


def measure(path,work,side,pq,label,profile=False):
    p,q=pq;fr=F(p,q)
    env=call([work/"envelopes","check",path],work/f"{label}-envelopes.log").strip()
    n1,n2=M["envelopes"]
    require(env==f"envelopes 1D={n1} 2D={n2} pass","incomplete envelope verification")
    low=call([work/"low","scan",path,str(p),str(q)],work/f"{label}-low.jsonl")
    caps={};seen=set()
    for line in low.splitlines():
        row=json.loads(line);fam=row["family"]
        require(fam not in seen and fam in ("f_inner","f_table","g") and row["rows"]==M["rows"][fam] and
                (row["P"],row["Q"])==(p,q),"incomplete/mismatched low scan")
        caps[fam]=row["cap"];seen.add(fam)
    require(seen=={"f_inner","f_table","g"},"missing low row family")
    analytic=tail_bound(path,fr)
    caps["f_analytic"]=analytic["cap"]
    (work/f"{label}-analytic.json").write_text(json.dumps(analytic)+"\n")
    txt=call([work/"high",path,"--pq",str(p),str(q)],work/f"{label}-high.log")
    require(f"mode={side.upper()} P={p} Q={q} S={UNIT} Sr={M['S_r']}" in txt.splitlines(),"high scan target")
    high_seen=set()
    for fam,n,cap in re.findall(r"^family=([VCPFYT]) rows=(\d+) .*?max_ceil=(-?\d+) ",txt,re.M):
        require(fam not in high_seen and int(n)==M["rows"][fam],"incomplete high family")
        caps[fam]=int(cap);high_seen.add(fam)
    require(high_seen==set("VCPFYT"),"missing high families")
    nums=re.findall(r"^higher_residual_per_core=(\d+)/(\d+)$",txt,re.M)
    require(len(nums)==1 and int(nums[0][1])==UNIT,"missing high aggregation")
    higher=int(nums[0][0])
    args=[work/"comparison",path]
    if profile:args+=["--profile",work/("mc-mbar-"+side)]
    comp=json.loads(call(args,work/f"{label}-comparison.json"))
    require(comp["comparison"] and comp["unit"]==UNIT and comp["finite_rows"]==M["rows"]["f_inner"]+M["rows"]["f_table"] and
            0<comp["lower"]<=comp["upper"],"comparison failure")
    error=max(caps["f_inner"],caps["f_table"],caps["f_analytic"])+caps["g"]+higher
    # Take the nonnegative part for the base bound; a signed negative cap is still valid.
    base=fr+(1 if side=="upper" else -1)*F(max(0,error),comp["lower"])
    return {"trial":[p,q],"caps":caps,"higher_per_core":higher,"error_per_core":error,
            "min_weight":comp["lower"],"max_weight":comp["upper"],"base":str(base)}


def test(side,path,work,fixture):
    expected=M[side];p,q=expected["trial"]
    args=[sys.executable,"-B"]
    oracle=call(args+[HERE/"low_oracle.py",path,"--target",f"{p}/{q}","--low",work/"low"],work/f"{side}-oracle.log")
    tree=call(args+[HERE/"tree_check.py",path,"--target",f"{p}/{q}","--cap",str(expected["error_per_core"]),
                    "--lo",str(expected["min_weight"]),"--hi",str(expected["max_weight"]),
                    "--low",work/"low","--high",work/"high"],work/f"{side}-tree.log")
    hcheck=call(args+[HIGH/"physical_check.py",path,"--check",work/"high","--pq",str(p),str(q),"--nmax","7"],work/f"{side}-high-physical.log")
    mut=call(args+[HERE/"mutations.py",side,path,fixture,work],work/f"{side}-mutations.log",True)
    print(oracle.strip(),tree.strip(),hcheck.splitlines()[-1],mut.splitlines()[-1],sep="\n",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--snapshots",type=Path,help="verify existing mc-upper.bin and mc-lower.bin; default regenerates both")
    ap.add_argument("--work",type=Path,help="a new or empty private work directory, default round2/outputs/mc-replay-*")
    ap.add_argument("--side",choices=SIDES)
    ap.add_argument("--checks-only",action="store_true",help="run all universal arithmetic, skipping independent tree and mutation controls")
    a=ap.parse_args()
    out=HERE.parent/"outputs";out.mkdir(exist_ok=True)
    if a.work:
        work=a.work.resolve();work.mkdir(parents=True,exist_ok=True)
        require(not any(work.iterdir()),"existing work directory is not empty")
    else:work=Path(mkdtemp(prefix="mc-replay-",dir=out))
    # Restrict independent proof, build, logs, profiles and mutation copies to user-selected work.
    os.environ.update({"TMPDIR":str(work),"OMP_NUM_THREADS":"1","PYTHONDONTWRITEBYTECODE":"1"})
    if hasattr(os,"sched_getaffinity"):
        cpu=min(os.sched_getaffinity(0));os.sched_setaffinity(0,{cpu})
    print("work",work,flush=True)
    compile_tools(work)
    fixture=build(work,"upper","fixture") if not a.checks_only else None
    for side in ([a.side] if a.side else SIDES):
        path=(a.snapshots.resolve()/f"mc-{side}.bin") if a.snapshots else build(work,side)
        expected=M[side]
        require(path.stat().st_size==M["bytes"] and digest(path)==expected["sha256"],"frozen witness size or digest mismatch")
        result=measure(path,work,side,expected["trial"],side,True)
        for k in ("caps","higher_per_core","error_per_core","min_weight","max_weight"):
            require(result[k]==expected[k],f"{side} arithmetic result differs for {k}")
        (work/f"{side}-verified.json").write_text(json.dumps(result,indent=2)+"\n")
        print(f"{side} S={UNIT} cap={result['error_per_core']} min={result['min_weight']} max={result['max_weight']} "
              f"base={result['base']} universal arithmetic PASS",flush=True)
        if not a.checks_only:test(side,path,work,fixture)
    print("PASS",flush=True)


if __name__=="__main__":
    try:main()
    except (ValueError,ArithmeticError,OSError,subprocess.CalledProcessError) as e:
        print("REJECTED",e,file=sys.stderr);sys.exit(2)
