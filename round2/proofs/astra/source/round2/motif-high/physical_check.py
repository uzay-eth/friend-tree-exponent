#!/usr/bin/env python3
"""Replay the h=1 dynamics on concrete core trees and compare with exact corner rows."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import product
import subprocess

from oracle import ARITY, CUTS, FAMILIES, INF, NB, Certificate, Tree, bud, degree, gain, inv, newborn_extensions, replace


def trees(cert, nmax):
    # Increasing parent lists exhaust every rooted core tree shape at these sizes.
    # Leaf multiplicities decorate the physical core tree; no leaf-rate sampling.
    d3, d4 = cert.par['D3'],cert.par['D4']
    for n in range(1,nmax+1):
        for seq in product(*(range(1,i) for i in range(2,n+1))):
            p = (-1,0)+seq
            k = [p.count(i) for i in range(n+1)]
            for tag in range(3):
                leaves = [0]
                for i in range(1,n+1):
                    if tag == 0: l = (i+i//2+n)%3
                    elif tag == 1: l = 1
                    elif k[i] == 0 or i == 1:
                        d = (d4-1,d4,d4+1,d3-1,d3,d3+1,2*d3+3)[(n+3*i)%7]
                        l = d-1-k[i]
                    else: l = (k[i]+i)%3
                    leaves.append(max(l,1-k[i]))
                yield Tree(p,leaves)
    # Each of the maximal shapes has simultaneous tail donors, endpoint
    # promotions and mixed exact/internal cutoffs in separate decorated trees.
    for parents, ends in [
            ((-1,0,1,2,3), (1,4)),           # vertical four
            ((-1,0,1,1,3), (2,4)),           # fork
            ((-1,0,1,2,2), (1,3,4)),         # parent and two children
            ((-1,0,1,1,1), (2,3,4)),         # trident
            ((-1,0,1,1,2,4,4), (3,5,6))]:
        for mode in range(8):
            leaves = [0]+[1]*(len(parents)-1)
            for j,x in enumerate(ends):
                d = (d4-1,d4,d4+1,d3-1,d3,d3+1,2*d3+1,1009)[(j+mode)%8]
                leaves[x] = max(1,d-1-parents.count(x))
            yield Tree(parents,leaves)
        # Tail all the ends at the same time, well away from both knots.
        leaves = [0]+[1]*(len(parents)-1)
        for j,x in enumerate(ends): leaves[x] = 2*d3+1+j*97
        yield Tree(parents,leaves)
    # High actual k on an END lies beyond all clamped categories; all of its
    # real outside children remain present in both jump rates and motif counts.
    p = [-1,0,1,2,3]+[4]*(cert.par['K3']+3)
    for big_l in (0,d3,d3*3):
        l = [0]+[1]*(len(p)-1)
        l[4] = big_l
        yield Tree(p,l)
    # Promotion can take a clamped category to overflow at either knot.
    for k in sorted({cert.par['K3'],cert.par['K4']}):
        p = [-1,0,1,2]+[3]*k
        yield Tree(p,[0,1,1,d3+3]+[1]*k)


def physical_case(t, c, records, counts, plant):
    old = t.active(c)
    phi = sum((w for _,_,_,w in old.values()),Q(0))
    direct = Q(0)
    newborns = {}
    fresh = len(t.p)
    for x in range(1,fresh):
        direct += t.r[x]*(t.after(x,False).potential(c)-phi)
        if t.b[x]:
            new = t.after(x,True).active(c)
            direct += t.b[x]*(sum((w for _,_,_,w in new.values()),Q(0))-phi)
            for vs,(_,_,_,w) in new.items():
                if fresh in vs: newborns[x,vs] = w

    debit, credit, booked = {}, {}, {}
    total = Q(0)
    for x in range(1,fresh):
        if t.b[x]:
            w = c.weight('E',(bud(t.s[x]),NB))
            total += t.b[x]*w
            booked[x,frozenset((x,fresh))] = w
    # The high-degree singleton row uses the total uncentered input. This
    # identity is independent of every four-core outside incidence below.
    fstart,_ = c.sections['f']
    for x,s in enumerate(t.s):
        if s is INF or not c.par['D2'] <= degree(s) < c.par['DF']-1: continue
        d, k = degree(s),s[1]
        id0, id1 = d*(d-1)//2-1+k, d*(d+1)//2-1+k
        delta = c.sign*(c.data[fstart+id1]-c.data[fstart+id0])
        nb = [t.p[x]]+t.ch[x]
        individual = [Q(max(delta,0),2)-inv(t.s[y])*delta for y in nb]
        whole = Q(max(delta,0)*(k+1),2)-sum(inv(t.s[y]) for y in nb)*delta
        assert sum(individual) == whole and min(individual) >= 0
        counts['f_credit_incidence'] += len(individual)
        counts['f_credit_positive'] += sum(q > 0 for q in individual)
    rows = []
    for support,(f,verts,st,w) in old.items():
        row = Q(0)
        outside_loss = Q(0)
        exact_four = ARITY[f] == 4 and all(z is not INF and degree(z) < c.par['D4'] for z in st)
        for j,x in enumerate(verts):
            if x == 0: continue
            d = c.weight(f,replace(st,j,gain(st[j])))-w
            row += t.r[x]*d
            if t.b[x]:
                post = replace(st,j,bud(st[j]))
                row += t.b[x]*(c.weight(f,post)-w)
                if ARITY[f] < 4:
                    larger,new = newborn_extensions(f,post,j)
                    if c.active(larger,new) >= 0:
                        nb = c.weight(larger,new)
                        key = (x,support|{fresh})
                        if plant != 'birth' or larger != 'F':
                            assert key not in booked, "duplicate physical birth"
                            booked[key] = nb
                            row += t.b[x]*nb
            if ARITY[f] == 4:
                if exact_four:
                    for y in [t.p[x]]+t.ch[x]:
                        if y == 0 or y in support: continue
                        q = Q(max(c.sign*d,0),2)-inv(t.s[y])*c.sign*d
                        assert q >= 0
                        outside_loss += q
                        counts['four_credit_incidence'] += 1
                        counts['four_credit_positive'] += q > 0
                continue
            for y in [t.p[x]]+t.ch[x]:
                if y in support: continue
                extended = support|{y}
                if extended not in old: continue
                larger = old[extended][0]
                if f == 'E':
                    cent = 4 if larger == 'C' else (2 if y == t.p[x] else 3)
                else:
                    cent = 0 if y == t.p[x] else 1
                key = (support,extended)
                assert key not in debit
                debit[key] = (inv(t.s[y])-c.center(cent,st[j]))*d
                row -= debit[key]
        for donor,contact,small,keep,cent in CUTS.get(f,[]):
            if plant == 'transfer' and f == 'F' and donor == 1: continue
            kept = tuple(st[k] for k in keep)
            key = (frozenset(verts[k] for k in keep),support)
            assert key[0] in old and key not in credit
            d = c.weight(small,replace(kept,keep.index(contact),gain(st[contact])))-c.weight(small,kept)
            credit[key] = (inv(st[donor])-c.center(cent,st[contact]))*d
            row += credit[key]
        total += row
        if f != 'E':
            rows.append((f,st,c.sign*(row-c.lam*w)+outside_loss*(2 if plant == 'credit' else 1)))
    if booked != newborns: raise ValueError("physical newborn booking mismatch")
    if debit != credit: raise ValueError("physical centered transfer mismatch")
    if total != direct: raise ValueError("decomposed drift differs from direct physical updates")
    counts['trees'] += 1
    counts['events'] += sum(1+bool(b) for b in t.b[1:])
    counts['births'] += len(booked)
    counts['transfers'] += len(debit)
    counts['nonzero_transfers'] += sum(v != 0 for v in debit.values())
    for f,st,r in rows:
        key = f,st
        if key not in records or r > records[key]: records[key] = r
        counts[f] += 1
        ntail = sum(s is not INF and degree(s) > c.par['D'+str(ARITY[f])] for s in st)
        counts['tail'+str(ntail)] += 1


def audit_rows(c, records, binary, path, counts):
    keys = sorted({(f,x) for f,st in records for x,_ in c.corners(f,st)})
    # The query path always returns a reduced exact rational, bypassing integer ceiling.
    inp = ''.join(f+' '+' '.join(map(str,x))+'\n' for f,x in keys)
    out = subprocess.run([binary,path,'--query','--pq',str(c.lam.numerator),str(c.lam.denominator)],
                         input=inp,text=True,capture_output=True,check=True)
    values = [Q(q) for q in out.stdout.splitlines()]
    assert len(values) == len(keys)
    lookup = dict(zip(keys,values))
    gaps = {}
    for (f,st),r in records.items():
        slack = sum((a*lookup[f,x] for x,a in c.corners(f,st)),Q(0))-r
        if slack < 0: raise ValueError(f"physical row exceeds exact corner majorant: {f} {st} excess={-slack}")
        if f not in gaps or slack < gaps[f]: gaps[f] = slack
        counts['strict'] += slack > 0
    assert set(gaps) == set(FAMILIES[1:])
    print(f"unique_physical_rows={len(records)} exact_corner_queries={len(keys)} strict_majorants={counts['strict']}")
    print("smallest_majorant_minus_physical_units", ' '.join(f'{f}={gaps[f]}' for f in FAMILIES[1:]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('file')
    ap.add_argument('--check',default='round2/outputs/mh-check')
    ap.add_argument('--nmax',type=int,default=6)
    ap.add_argument('--pq',type=int,nargs=2,metavar=('P','Q'))
    ap.add_argument('--plant',choices=('transfer','birth','credit'))
    args = ap.parse_args()
    c = Certificate(args.file)
    if args.pq: c.lam = Q(*args.pq)
    print(f"mode={'UPPER' if c.sign == 1 else 'LOWER'} target={c.lam}")
    counts, records = Counter(), {}
    for i,t in enumerate(trees(c,args.nmax)):
        try: physical_case(t,c,records,counts,args.plant)
        except ValueError as ex:
            print(f"case={i} parents={t.p} leaf_counts={t.l} {ex}")
            return 2
    print("physical identities", ' '.join(f"{k}={counts[k]}" for k in
              ('trees','events','births','transfers','nonzero_transfers')))
    print("physical row occurrences", ' '.join(f"{k}={counts[k]}" for k in FAMILIES[1:]))
    print("physical tail counts", ' '.join(f"{n}={counts['tail'+str(n)]}" for n in range(4)))
    print("retained credits", ' '.join(f"{k}={counts[k]}" for k in
                  ('four_credit_incidence','four_credit_positive','f_credit_incidence','f_credit_positive')))
    try: audit_rows(c,records,args.check,args.file,counts)
    except ValueError as ex:
        print(str(ex))
        return 2
    print("PASS")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
