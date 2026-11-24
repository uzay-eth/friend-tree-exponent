"""Force independent acceptance predicates to reject incorrect low and high witnesses."""
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil
import struct
import subprocess
import sys
from low_oracle import Local, probes
from replay import HERE, M, digest
from tail import verify as tail_bound
import negative as high_negative
import corner_audit as high_reference


def must_reject(fn,description):
    try: fn()
    except ValueError as e:
        print("rejected",description,str(e),flush=True)
    else: raise ArithmeticError("accepted planted defect: "+description)


def limit(value,cap):
    if value>cap: raise ValueError(f"row {value} > budget {cap}")


def run_expected(cmd,code,reason):
    pr=subprocess.run(cmd,text=True,capture_output=True)
    if pr.returncode!=code or reason not in pr.stdout+pr.stderr:
        raise ArithmeticError(f"unexpected result {pr.returncode} from {cmd}: {pr.stdout[-150:]} {pr.stderr[-150:]}")
    print("rejected",reason,"exit",code,flush=True)


def patch(path,offset,fmt,value,action):
    with path.open("r+b") as f:
        f.seek(offset); saved=f.read(struct.calcsize(fmt))
        try:
            f.seek(offset);f.write(struct.pack(fmt,value));f.flush()
            action()
        finally:
            f.seek(offset);f.write(saved);f.flush()


def run(side,original,fixture,work):
    m=M[side];target=F(*m["trial"])
    K=m["caps"];S=M["S"];low=work/"low";envelope=work/"envelopes";high=work/"high"
    c=Local(original); g=c.l2; sgn=c.sign
    with TemporaryDirectory(dir=work,prefix="mutate-") as tmp:
        scratch=Path(tmp)/"copy.bin";shutil.copyfile(original,scratch)
        def put(section,i,value,fn):
            patch(scratch,248+8*(c.sections[section][0]+i),"<q",value,fn)
        def row(tag,a,b,budget,pq=target):
            k=(tag,a,b);limit(probes(c,low,scratch,pq,[k])[k],budget)
        v=int(c.f((1,0)))
        put("f",0,v-sgn*S,lambda:must_reject(lambda:row("v",1,0,K["f_inner"]),"singleton"))
        a,b=g.id[0,1],g.id[1,0]
        off=a*g.n+b
        put("E",off,c.tables["E"][off]-sgn*S,
            lambda:must_reject(lambda:row("g",a,b,K["g"]),"edge"))
        L=c.par["DH"]-2;s=(L,1);u=g.id[s];nk=len(c.grid)
        old=c.data[c.sections["fy"][0]+u*nk]
        put("fy",u*nk,old+sgn*S,
            lambda:must_reject(lambda:row("g",g.inf,u,K["g"]),"FY incoming knot"))
        x=c.Sr//2 if sgn*(c.f((2,0))-c.f((1,0)))>0 else 0
        put("rpar",g.id[1,0],x,lambda:must_reject(lambda:row("v",1,0,K["f_inner"]),"vertex center"))
        trial=target-F(sgn,100)
        must_reject(lambda:row("v",1,0,K["f_inner"],trial),"stricter exponent")
        # Analytic parameters are not certified by a finite sample or metadata hash.
        p0=c.coef[0]//2 if side=="upper" else c.coef[0]+S//3
        patch(scratch,40,"<q",p0,lambda:must_reject(
            lambda:limit(tail_bound(scratch,target)["cap"],K["f_analytic"]),"analytic f tail"))
        for section,msg in [("e230","1D envelope failure"),("eyv0","2D envelope failure")]:
            v=c.data[c.sections[section][0]]
            put(section,0,v+1,lambda:run_expected([envelope,"check",scratch],2,msg))
        put("E",g.inf*g.n+g.inf,1,lambda:run_expected([low,"scan",scratch],2,"nonzero g(INF,INF)"))
        patch(scratch,0,"<Q",0,lambda:run_expected([envelope,"check",scratch],2,"invalid frozen geometry"))
        # D0_2 field offset: 8+4*4+6*8=72, then DF occupies four bytes.
        patch(scratch,76,"<i",2**30,lambda:run_expected([envelope,"check",scratch],2,"unsupported degree geometry"))
        if digest(scratch)!=m["sha256"]:raise ArithmeticError("mutation copy restoration")
        base=[sys.executable,"-B",HERE/"tree_check.py",original,"--target",str(target),"--cap",str(m["error_per_core"]),
              "--lo",str(m["min_weight"]),"--hi",str(m["max_weight"]),"--low",low,"--high",high,"--nmax","3"]
        run_expected(base+["--plant","fy"],1,"complete physical drift transfer identity failed")
        run_expected(base+["--plant","double-gap"],1,"credited physical row exceeds corners")
        high_negative.tail_control(scratch,fixture,high)
        high_reference.make_table(fixture,scratch,1 if side=="upper" else 2)
        high_reference.run(scratch,high)
        high_negative.run_side(side,original,high,scratch)
    if digest(original)!=m["sha256"]: raise ArithmeticError("original witness digest changed")
    print(f"{side} independent numerical, geometry, and physical mutation controls PASS",flush=True)


if __name__=="__main__":
    side,original,fixture,work=sys.argv[1:]
    run(side,*map(Path,(original,fixture,work)))
