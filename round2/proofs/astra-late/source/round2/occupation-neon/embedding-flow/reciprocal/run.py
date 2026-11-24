"""Solve the physical inverse-degree moment relaxation with rational dual repair."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
from model import Model
from solve import original, optimize


def inner_charge(model, sources, names, theta):
    """Original first-body and edge contrast on the complete explicitly named degrees."""
    if model.threshold>17:
        raise ValueError("inner D18 extraction supports H<=17")
    cost=defaultdict(Q)
    for family, rows in (("f",model.roots),("g",model.edges)):
        if family not in names:
            continue
        for cells in rows:
            if any(c.degree==0 for c in cells):
                continue
            labels=tuple((c.degree-c.arity-1,c.arity,1) for c in cells)
            e=(0,)*len(cells)
            for table,sign in zip(sources,(1,-theta)):
                a,ports=table.row(family,labels)
                scale=sign/Q(table.c.S)
                model.add_observation(cost,cells,e,scale*a)
                for i,mother,_,slope in ports:
                    model.add_observation(cost,cells,e,scale*slope,
                                          kind=(i,"mother" if mother else "children"))
    return cost


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threshold",type=int,default=4)
    parser.add_argument("--arity",type=int,default=3)
    parser.add_argument("--order",type=int,default=1)
    parser.add_argument("--slabs",type=int,default=1)
    parser.add_argument("--no-edge-balances",action="store_true")
    parser.add_argument("--no-sharing",action="store_true")
    parser.add_argument("--objective",nargs="+",default=["birth-virgin","mixed-edge"],
                        choices=["birth-virgin","mixed-edge","old-mother-birth","f","g","fg"])
    parser.add_argument("--payloads",nargs=2,type=Path)
    parser.add_argument("--theta",type=Q,default=Q(1))
    parser.add_argument("--method",choices=["highs-ipm","highs-ds","highs"],default="highs-ipm")
    args=parser.parse_args()
    if args.no_edge_balances and any(n in ("g","fg") for n in args.objective):
        parser.error("g uses live exterior fluxes on the old edge")
    if args.no_edge_balances and args.no_sharing and "mixed-edge" in args.objective:
        parser.error("mixed-edge requires retained edge columns")
    model=Model(args.threshold,args.arity,args.order,args.slabs,
                edge_balances=not args.no_edge_balances,share=not args.no_sharing)
    print("threshold, arity, test degree, Bernstein degree, slabs",
          model.threshold,model.arity,model.m,model.n,model.slabs,flush=True)
    print("columns",len(model.variables),"balances and other equalities",len(model.equalities),
          "upper rows",len(model.upper_rows),flush=True)
    sources=None
    for name in args.objective:
        paired=name in ("f","g","fg")
        if paired:
            if not args.payloads:
                parser.error("a D18 objective requires both original payloads")
            if sources is None:
                sources=tuple(original(path,side) for path,side in zip(args.payloads,("lower","upper")))
                print("authenticated both physical D18 payloads",flush=True)
            c=inner_charge(model,sources,name,args.theta)
        else:
            c=model.objective(name)
        cuts,errors=[],[]
        for maximize in (False,True):
            val,err=optimize(model,c,maximize,method=args.method,scale=10**7 if paired else 10)
            cuts.append(val)
            errors.append(err)
        if paired or name=="mixed-edge":
            vals=[v*s for v in cuts for s in (model.low,model.high)]
            cuts=min(vals),max(vals)
            unit="per core"
        else:
            unit="per birth"
        print(name,unit,*cuts,"decimal",*[float(x) for x in cuts],
              "repair before normalization",*[float(x) for x in errors],flush=True)


if __name__=="__main__":
    main()
