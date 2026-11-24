#!/usr/bin/env python3
"""Exact tensor Bernstein enclosure of the analytic vertex residual, including q=0, x=0."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse
import json
import struct


def header(fp):
    """Read the frozen layout independently of its C implementation."""
    data = fp.read(248)
    vals = struct.unpack("<Q4i6q44i", data)
    names = ("magic mode P Q reserved0 S S_r f_p0 f_p1 f_p2 alpha_lo_num "
             "DF D0_2 KT_2 D0_3 KT_3 DHU KHU DCU KCU D0_4 KT_4 "
             "DP4_1 DP4_2 DP4_S KP4_1 KP4_2 DF4_1 DF4_2 DF4_S KF4_1 KF4_2 "
             "DY4 KY4 DT4 KT4 Q_K nknots").split()
    h = dict(zip(names, vals))
    if h["magic"] != 0x4654325F34563031 or h["mode"] not in [1,2]:
        raise ValueError("invalid certificate header")
    return h


def read_tail(path):
    with Path(path).open("rb") as fp:
        h = header(fp)
        d, k = h["D0_2"], h["KT_2"]
        inner = d*(d-1)//2-1
        ng = inner+1+k+2
        nf = h["DF"]*(h["DF"]-1)//2-1
        start = (5+3*h["nknots"])*inner
        def item(ix):
            fp.seek(248+8*ix)
            return struct.unpack("<q",fp.read(8))[0]
        fnb = item(start)  # D=2,k=0
        gbase = start+nf
        gi = item(gbase+ng*inner)
        # After a promotion k' >= 1. Its newborn is the D=2,k=0 corner.
        newborn = [item(gbase+ng*(inner+1+j)) for j in range(1,k+2)]
        dg = d*((max if h["mode"]==1 else min)(newborn)-gi)
    return h,fnb+gi,dg


class Poly(dict):
    def __add__(self, other):
        other = poly(other)
        return Poly({k:self.get(k,0)+other.get(k,0) for k in self.keys() | other.keys()})
    __radd__ = __add__
    def __neg__(self):
        return Poly({k:-v for k,v in self.items()})
    def __sub__(self, other):
        return self+(-poly(other))
    def __rsub__(self, other):
        return poly(other)+(-self)
    def __mul__(self, other):
        ans = Poly()
        for (i,j),v in self.items():
            for (k,l),w in poly(other).items():
                key = (i+k,j+l)
                ans[key] = ans.get(key,0)+v*w
        return ans
    __rmul__ = __mul__
    def at(self,q,x):
        return sum(v*q**i*x**j for (i,j),v in self.items())


def poly(v):
    return v if isinstance(v,Poly) else Poly({(0,0):v})


Q = Poly({(1,0):1})
X = Poly({(0,1):1})


def bernstein(poly_, q0, q1, xmax):
    """Power to tensor Bernstein after q=q0+(q1-q0)*t, x=xmax*y."""
    p = {k:v for k,v in poly_.items() if v} or {(0,0):0}
    nq = max(i for i,j in p)
    nx = max(j for i,j in p)
    power = {}
    for (i,j),v in p.items():
        for a in range(i+1):
            key = a,j
            power[key] = power.get(key,0)+v*comb(i,a)*q0**(i-a)*(q1-q0)**a*xmax**j
    return [sum(c*F(comb(i,a),comb(nq,a))*F(comb(j,b),comb(nx,b))
                for (a,b),c in power.items() if a<=i and b<=j)
            for i in range(nq+1) for j in range(nx+1)]


def verify(path, target=None, panels=16):
    h,C,d = read_tail(path)
    lam = target if target is not None else F(h["P"],h["Q"])
    upper = h["mode"]==1
    a,b,c,e = ((h["f_p0"],h["f_p1"],0,0) if upper else
               (0,h["f_p0"],-h["f_p1"],-h["f_p2"]))
    # f = a+b*q+c*q^2+e*q*x, x=1/D. Work in integer weight units.
    w = a+b*Q+c*Q*Q+e*Q*X
    den = (1+X)*(1+X)
    gain = b*(1-Q)*(1+X)+c*(1-Q)*(2*Q+(1+Q)*X)+e*((1-2*Q)*X-Q*X*X)
    bud = -b*X+c*(-2*Q*X+X*X)-e*X*X
    cap = None
    xmax = F(1,h["DF"])
    for extreme,rho in enumerate((Q,(1+Q)*F(1,2))):
        p = (rho*gain+(Q*(bud+C+d*X)-w*lam)*den)*(1 if upper else -1)
        # Independent exact oracle at actual integer states, including the first seam.
        for D in (h["DF"], h["DF"]+1, 3*h["DF"]+17):
            x = F(1,D)
            for L in (0,1,D//3,D-1):
                q = L*x
                r = L if extreme==0 else F(D+L,2)
                direct = r*(w.at(F(L+1,D+1),F(1,D+1))-w.at(q,x))
                direct += q*(w.at(q-x,x)-w.at(q,x)+C+d*x)-lam*w.at(q,x)
                if p.at(q,x) != (1 if upper else -1)*direct*(1+x)**2:
                    raise ArithmeticError("gain/bud polynomial oracle")
        for j in range(panels):
            top = max(bernstein(p,F(j,panels),F(j+1,panels),xmax))
            bound = top if top>=0 else top/(1+xmax)**2
            cap = bound if cap is None else max(cap,bound)
    return {"family":"f_analytic","rows":"D>=DF","DF":h["DF"],"panels":panels,
            "cap":-(-cap.numerator//cap.denominator),"cap_rational":str(cap),"P":lam.numerator,"Q":lam.denominator}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("snapshot")
    p.add_argument("--target",type=F)
    p.add_argument("--panels",type=int,default=16)
    args=p.parse_args()
    if args.panels<1: p.error("at least one panel is required")
    print(json.dumps(verify(args.snapshot,args.target,args.panels),sort_keys=True))
