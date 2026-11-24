#!/usr/bin/env python3
"""Cold audit of the eligibility of the additional source APIs.

Verifies, against the two frozen snapshots only:
  * exceptional FY knots and the discrete-grid support condition,
  * frozen positive allocations bound the maximized residual (floors nonnegative),
  * paid row and donor credits are disjoint summands (negative controls),
  * degree-72 true seam (patched-header control),
  * high-role permutation invariance at symmetric frozen shapes,
  * physical gains from the frozen tables and their integer floors.
"""

import argparse
from fractions import Fraction as Q
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

REPO = Path(__file__).resolve().parents[2]
SNAP = {"upper": REPO / "round2/outputs/mc-upper.bin",
        "lower": REPO / "round2/outputs/mc-lower.bin"}
FCAP = {"upper": 154457, "lower": 95343}

CPP = r'''
#include "round2/motif-cert/reserve/reserve.hpp"
#include "round2/motif-high/occupation/credit.hpp"
#include <cstdlib>
#include <iostream>
#include <array>
using namespace mh;
int main(int,char**argv){
  const char* file=argv[1]; int side=std::atoi(argv[2]);
  mc::Snapshot hdr(file);
  std::cout<<"GRID "<<hdr.H.D0_2<<" "<<hdr.H.DF<<" "<<hdr.H.nknots<<" "<<hdr.H.Q_K<<"\nKNOOTS";
  for(int j=0;j<hdr.H.nknots;j++) std::cout<<" "<<hdr.H.knot_num[j];
  std::cout<<"\n";
  motif_reserve::Scores sc(file);
  std::cout<<"SEAM 72 3 "<<sc.edge_residual_ceiling({72,2},{3,1})<<"\n";
  std::cout<<"SEAM 71 3 "<<sc.edge_residual_ceiling({71,2},{3,1})<<"\n";
  for(int D=2;D<=80;D++) std::cout<<"SING "<<D<<" 2 "<<sc.singleton_floor({D,2})<<"\n";
  int shown=0;
  for(int Dp: {6,8,12,24,48,72}) for(int Dc: {2,3,4,5}) {
    if(Dc>Dp) continue;
    for(auto slot: {motif_reserve::Slot::ParentOfParent,
                    motif_reserve::Slot::OtherChildOfParent,
                    motif_reserve::Slot::ChildOfChild})
      for(int dy: {2,3,4,7,16,72}) for(int pooled=0;pooled<2;pooled++) {
        auto don=pooled?motif_reserve::Donor::pool(dy):motif_reserve::Donor::exact(dy);
        long long df=sc.edge_donor_floor({Dp,2},{Dc,1},slot,don);
        if(!df) continue;
        auto ce=sc.edge_residual_ceiling({Dp,2},{Dc,1});
        auto ru=sc.edge_unused_floor({Dp,2},{Dc,1});
        std::cout<<"EDGE "<<Dp<<" "<<Dc<<" "<<(int)slot<<" "<<dy<<" "<<pooled
                 <<" "<<ce<<" "<<ru<<" "<<df<<"\n";
        if(++shown>=8) goto done_edges;
      }
  }
  done_edges:
  mh::Snapshot snap(file);
  mh_occ::Source src(file,side==1?56581165:5656453,side==1?100000000:10000000);
  auto level=[&](Family f)->const Cert4vLevel&{return f==mh::Edge?snap.c.L2:snap.level(f);};
  auto point=[&](const Cert4vLevel&g,int D,int c){int L=D-c,k=c-1;return D<g.D0?c4_lvl_inner(&g,L,k):c4_lvl_tail(&g,k);};
  auto val=[&](Family f,const mh_occ::Members&s,int n){
    const auto&g=level(f); std::array<int,4> x{};
    for(int j=0;j<n;j++) x[j]=point(g,s[j].degree,s[j].core_degree);
    int a=f==mh::Edge?0:snap.active(f,x);
    if(a<0) return 0LL;
    return f==mh::Edge?(long long)snap.sign*snap.weight[f][(size_t)x[0]*snap.c.L2.NG+x[1]]
                       :(long long)snap.value(f,a,x);
  };
  shown=0;
  for(Family f: {mh::V,mh::C}) {
    const auto&form=mh::shapes[f]; int child[4]={0,0,0,0};
    for(int i=1;i<form.n;i++) child[form.parent[i]]++;
    auto sh=mh_occ::Shape::canonical(f);
    for(int d0: {4,6,8,12,16}) for(int d1: {3,4,6,8}) for(int d2: {3,4,6,8}) {
      int deg[4]={d0,d1,d2,4};
      mh_occ::Members s{};
      for(int j=0;j<form.n;j++){int pos=sh.role[j]; s[pos]={deg[j],1+child[j]};}
      if(src.row(sh,s.data())==0) continue;
      for(int i=0;i<form.n;i++) {
        auto s2=s; s2[i].degree++;
        long long g=val(f,s2,form.n)-val(f,s,form.n);
        auto dir=(i==0)?mh_occ::Direction::Parent:mh_occ::Direction::Child;
        long long ex=src.triple_donor(sh,s.data(),i,dir,7,false);
        long long po=src.triple_donor(sh,s.data(),i,dir,7,true);
        if(!ex && !po) continue;
        std::cout<<"GAIN "<<(int)f<<" "<<d0<<" "<<d1<<" "<<d2<<" "<<i<<" "<<g<<" "<<ex<<" "<<po<<"\n";
        if(++shown>=8) goto done_gain;
      }
    }
  }
  done_gain:
  for(Family f: {mh::C,mh::Y,mh::T}) {
    const auto&form=mh::shapes[f]; int child[4]={0,0,0,0};
    for(int i=1;i<form.n;i++) child[form.parent[i]]++;
    auto sh=mh_occ::Shape::canonical(f);
    mh_occ::Shape sw; sw.family=f;
    sw.role = (f==mh::Y)? std::array<int,4>{0,1,3,2} : std::array<int,4>{0,2,1,3};
    int deg[4]={8,5,11,4};
    mh_occ::Members a{},b{};
    for(int j=0;j<form.n;j++){a[sh.role[j]]={deg[j],1+child[j]}; b[sw.role[j]]={deg[j],1+child[j]};}
    std::cout<<"PERM "<<(int)f<<" "<<src.row(sh,a.data())<<" "<<src.row(sw,b.data())<<"\n";
  }
  return 0;
}
'''


def qfloor(q):
    return q.numerator // q.denominator


def donor_floor(d, degree, pooled):
    p = Q(1, degree)
    if pooled:
        lo = Q(0)
        hi = Q(1, degree)
        val = min(Q(max(d, 0), 2) - d * x for x in (lo, hi))
        return qfloor(val)
    return qfloor(Q(max(d, 0), 2) - d * p)


def parse(out):
    data = {"GRID": None, "KNOOTS": None, "SING": [], "EDGE": [], "GAIN": [], "PERM": [], "SEAM": []}
    for line in out.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] in data and p[0] not in ("GRID", "KNOOTS"):
            data[p[0]].append(tuple(map(int, p[1:])))
        elif p[0] == "GRID":
            data["GRID"] = tuple(map(int, p[1:]))
        elif p[0] == "KNOOTS":
            data["KNOOTS"] = tuple(map(int, p[1:]))
    return data


def grid_audit(data):
    d0, df, nk, qk = data["GRID"]
    knots = data["KNOOTS"]
    assert d0 == 72 and df == 1400 and qk == 10080 and nk == len(knots)
    real = [k == 0 or qk % k == 0 for k in knots]
    bad = [(Q(k, qk), j) for j, (k, r) in enumerate(zip(knots, real)) if not r]
    assert {str(f) for f, _ in bad} == {"5/12", "11/24"}, bad
    assert qk % knots[1] == 0 and Q(knots[1], qk) == Q(1, d0)
    supported = 0
    for D in range(2, d0 + 1):
        j = 0
        while j + 2 < nk and D * knots[j + 1] < qk:
            j += 1
        on = (D * knots[j] == qk) or (D * knots[j + 1] == qk)
        if not on:
            assert real[j] and real[j + 1], (D, j)
            supported += 1
    print(f"[grid] knots {nk}, exceptional exactly {[str(f) for f,_ in bad]}, "
          f"{supported} interior degrees supported between real knots; seam knot 1/{d0} exact")


def audit(side, work, patched=False):
    exe = work / ("probe_patched" if patched else "probe")
    src = work / "probe.cpp"
    if not patched:
        src.write_text(CPP)
    if patched:
        text = (REPO / "round2/motif-cert/reserve/reserve.hpp").read_text()
        needle = "return {at(l2().INF_ID)-at(from),x.D+1}; // D=d2, so pre-gain t=1"
        assert needle in text
        shutil.copytree(REPO / "round2/motif-cert", work / "motif-cert", dirs_exist_ok=True)
        (work / "motif").mkdir(exist_ok=True)
        shutil.copy(REPO / "round2/motif/cert4v_format.h", work / "motif/cert4v_format.h")
        (work / "motif-cert/reserve/reserve.hpp").write_text(
            text.replace(needle, needle.replace("x.D+1", "x.D")))
        cpp = CPP.replace('"round2/motif-cert/reserve/reserve.hpp"',
                          f'"{work / "motif-cert/reserve/reserve.hpp"}"')
        src.write_text(cpp)
    subprocess.run(["c++", "-std=c++20", "-O2", "-Wall", "-Wextra", "-fno-fast-math",
                    "-I" + str(REPO), str(src), "-lgmpxx", "-lgmp", "-o", str(exe)],
                   check=True, cwd=REPO)
    out = subprocess.check_output([str(exe), str(SNAP[side]), "1" if side == "upper" else "2"],
                                  cwd=REPO, text=True)
    return parse(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("side", nargs="?", choices=("upper", "lower"), default="upper")
    args = ap.parse_args()
    if hasattr(os, "sched_getaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    with TemporaryDirectory(prefix="hy4-extra-", dir=REPO / "round2/outputs") as td:
        work = Path(td)
        data = audit(args.side, work)
        print(f"=== {args.side} ===")
        grid_audit(data)
        fcap = FCAP[args.side]
        assert all(0 <= f <= fcap + 5000 for _, _, f in data["SING"])
        print(f"[caps] {len(data['SING'])} singleton floors nonnegative, fcap={fcap}; "
              f"min={min(f for _,_,f in data['SING'])}, max={max(f for _,_,f in data['SING'])}")

        assert data["EDGE"] and data["GAIN"] and data["PERM"]
        overlap = 0
        for Dp, Dc, slot, dy, pooled, ce, ru, df in data["EDGE"]:
            assert ru == max(0, -ce)
            assert df >= 0
            wrong = max(0, -(ce - df))
            assert wrong > ru
            overlap += wrong - ru
        print(f"[reserve] {len(data['EDGE'])} donor queries: unused=max(0,-residual), all floors >=0; "
              f"recomputing the row with the actual donor gap would double-claim {overlap} extra units, "
              f"so paid U-psi and -Rmax are disjoint in the code")

        exact = pooled_n = 0
        for f, d0, d1, d2, i, g, ex, po in data["GAIN"]:
            if ex:
                assert ex == donor_floor(g, 7, False), (f, g, ex)
                exact += 1
            if po:
                assert po == donor_floor(g, 7, True), (f, g, po)
                pooled_n += 1
            if g % 7:
                assert ex != -(-g // 7) or po != -(-g // 7)
        print(f"[oracle] raw physical gains from frozen tables; {exact} exact and {pooled_n} pooled "
              f"donor floors match independently computed floor(J(d,1/7))")

        for f, a, b in data["PERM"]:
            assert a == b, (f, a, b)
        print("[roles] role swaps at symmetric frozen shapes C, Y and T give identical row values")

        patched = audit(args.side, work, patched=True)
        tr = {d: ce for d, _, ce in data["SEAM"]}
        al = {d: ce for d, _, ce in patched["SEAM"]}
        assert tr[71] == al[71] and tr[72] != al[72], (tr, al)
        print(f"[seam] patched header (D instead of D+1) leaves D=71 edge ceiling {tr[71]} "
              f"unchanged but changes D=72 from {tr[72]} to {al[72]}; true seam uses 73")
    print("EXTRA SOURCES COLD AUDIT: PASS")


if __name__ == "__main__":
    main()
