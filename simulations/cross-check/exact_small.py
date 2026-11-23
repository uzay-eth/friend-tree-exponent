#!/usr/bin/env python3
"""Independent exact E[C_n] for the random friend tree, by full enumeration of unlabelled trees.

Written cold for the simulation cross-check (Claude Opus 5.5). It shares no code with
verify/motif-bracket/. State = canonical string of an unrooted tree (AHU encoding rooted at the
centre, or the smaller of the two encodings for a bicentral tree). From T_2 (one edge), each step
picks a uniform vertex V and a uniform neighbour W of V and attaches a new leaf to W.

Usage: exact_small.py NMAX [EXACT_FILE]
Prints n, number of reachable trees, E[C_n] (exact), and whether it equals column 3 of EXACT_FILE.
"""
import sys
from fractions import Fraction

def centres(adj):
    n = len(adj)
    deg = [len(a) for a in adj]
    layer = [v for v in range(n) if deg[v] <= 1]
    removed = len(layer)
    while removed < n:
        nxt = []
        for v in layer:
            for w in adj[v]:
                deg[w] -= 1
                if deg[w] == 1:
                    nxt.append(w)
        removed += len(nxt)
        layer = nxt
    return layer

def encode(adj, root, parent):
    # iterative AHU encoding
    order, stack = [], [(root, parent)]
    par = {root: parent}
    while stack:
        v, p = stack.pop()
        order.append(v)
        for w in adj[v]:
            if w != p:
                par[w] = v
                stack.append((w, v))
    code = {}
    for v in reversed(order):
        kids = sorted(code[w] for w in adj[v] if w != par[v])
        code[v] = "(" + "".join(kids) + ")"
    return code[root]

def canon(adj):
    cs = centres(adj)
    if len(cs) == 1:
        return encode(adj, cs[0], -1)
    a, b = cs
    # bicentral: encode as an edge
    ea, eb = encode(adj, a, b), encode(adj, b, a)
    return "E" + min(ea + eb, eb + ea)

def decode(code):
    # rebuild an adjacency list from a canonical code
    adj = []
    def build(s, i, parent):
        # s[i] == '('
        v = len(adj); adj.append([])
        if parent >= 0:
            adj[v].append(parent); adj[parent].append(v)
        i += 1
        while s[i] == '(':
            i = build(s, i, v)
        return i + 1
    if code[0] == "E":
        body = code[1:]
        # split the two balanced halves
        depth = 0
        for j, ch in enumerate(body):
            depth += 1 if ch == '(' else -1
            if depth == 0:
                break
        left, right = body[:j + 1], body[j + 1:]
        build(left, 0, -1)
        r = len(adj)
        build(right, 0, -1)
        adj[0].append(r); adj[r].append(0)
    else:
        build(code, 0, -1)
    return adj

def main():
    nmax = int(sys.argv[1])
    ref = {}
    if len(sys.argv) > 2:
        for line in open(sys.argv[2]):
            f = line.split()
            if len(f) >= 3 and f[0].isdigit():
                ref[int(f[0])] = (int(f[1]), Fraction(f[2]))
    dist = {canon([[1], [0]]): Fraction(1)}
    n = 2
    allok = True
    while True:
        ec = Fraction(0)
        for code, pr in dist.items():
            adj = decode(code)
            ec += pr * sum(1 for a in adj if len(a) >= 2)
        line = f"{n} {len(dist)} {ec}"
        if n in ref:
            ok = ref[n] == (len(dist), ec)
            allok &= ok
            line += "  MATCH" if ok else f"  MISMATCH ref={ref[n][0]} {ref[n][1]}"
        print(line, flush=True)
        if n == nmax:
            break
        new = {}
        for code, pr in dist.items():
            adj = decode(code)
            nv = len(adj)
            # W = w has probability sum over neighbours v of w of 1/(n deg v)
            for w in range(nv):
                pw = sum(Fraction(1, len(adj[v])) for v in adj[w]) / nv
                adj2 = [list(a) for a in adj]
                adj2.append([w]); adj2[w].append(nv)
                k = canon(adj2)
                new[k] = new.get(k, Fraction(0)) + pr * pw
        dist = new
        n += 1
    print("ALL_MATCH" if allok and ref else ("NO_REFERENCE" if not ref else "SOME_MISMATCH"))

if __name__ == "__main__":
    main()
