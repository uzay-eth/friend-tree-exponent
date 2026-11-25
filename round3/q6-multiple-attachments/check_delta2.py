"""Portable export of the Q6 six-row exact finite checker (local-host pin stub removed;
original hashes recorded in ../PROVENANCE.md). All mathematical assertions unchanged."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product


def graph(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def raw_draws(adj, common, distinct=False):
    # Enumerate explicit uniform mediators and incident-edge indices, with copies.
    n = len(adj)
    law = Counter()
    for a in range(n):
        for b in ([a] if common else range(n)):
            size = len(adj[a])
            pairs = permutations(range(size), 2) if distinct else product(range(size), range(len(adj[b])))
            denom = (n if common else n*n)*size*((size-1) if distinct else len(adj[b]))
            for i,j in pairs:
                law[adj[a][i], adj[b][j]] += F(1,denom)
    assert sum(law.values()) == 1
    return law


def compare(adj, model):
    n = len(adj)
    p = [sum(F(row.count(v), n*len(row)) for row in adj) for v in range(n)]
    law = raw_draws(adj, model.startswith("A"), model=="A-simple")
    if model in ("B-rej-pair","B-rej-seq"):
        distinct_sum = sum(mass for (a,b),mass in law.items() if a!=b)
        law = {(a,b):mass/(distinct_sum if model=="B-rej-pair" else 1-p[a])
               for (a,b),mass in law.items() if a!=b}
        assert sum(law.values())==1
    for v,t in enumerate(p):
        native = simple = F(0)
        for (a,b),prob in law.items():
            hit = int(a==v)+int(b==v)
            native += (hit if model.endswith("-multi") else min(1,hit))*prob
            simple += min(1,hit)*prob
        if model in ("A-simple","B-multi","A-multi"):
            assert native==2*t
        if model=="A-simple":
            assert simple==2*t
        if model in ("B-simple-col","B-multi"):
            assert simple==2*t-t*t
        if model=="A-multi":
            rows = [F(row.count(v),len(row)) for row in adj]
            assert simple-t==sum(x*(1-x) for x in rows)/n
            pure = sum(len(row)==2 and row[0]==row[1] for row in adj)
            assert 0 <= simple-t <= F(n-pure,4*n)
        if model=="B-rej-pair":
            assert simple==2*t*(1-t)/(1-sum(x*x for x in p))
        if model=="B-rej-seq":
            assert simple==t*(1+sum(x/(1-x) for i,x in enumerate(p) if i!=v))
        if not model.endswith("-multi"):
            assert simple==native
    return len(p)


def main():
    diamond = graph(4,[(0,1),(0,2),(0,3),(1,2),(1,3)])
    mixed = graph(5,[(0,1),(0,1),(0,2),(0,3),(2,3),(2,4),(3,4)])
    edge, parallel = graph(2,[(0,1)]),graph(2,[(0,1)]*2)
    path = graph(3,[(0,1),(1,2)])
    cases = {"A-simple":[diamond], "A-multi":[mixed,parallel], "B-multi":[mixed,edge,parallel],
             "B-simple-col":[diamond,path,edge], "B-rej-pair":[diamond,path,edge], "B-rej-seq":[diamond,path,edge]}
    for model,graphs in cases.items():
        print(model, "finite old-label checks", sum(compare(adj,model) for adj in graphs))
    floor, cfin, lam, delta = F(10,49), F(40866367,200000000), F(49,50), F(1,20000)
    g, j = 1/(1-delta),1-5*delta
    tight = g*lam+2*(g-j)
    printed = g*lam+(g-j)*2*(floor+F(46,25))/cfin
    assert floor<cfin and tight<printed<F(99,100)
    assert printed==F(39488432850800,40047037208017)
    assert sum(range(3,339))==57288 and 338-3+1==336
    # Arithmetic comparison to the already certified seam bound in Astra's frozen replay.
    assert F(802252931359625003,30625000000000000000000)>F("0.00002619601")
    print("Global floor", floor, "tight lambda", tight, "printed lambda", printed)
    print("PASS")


if __name__=="__main__":
    main()
