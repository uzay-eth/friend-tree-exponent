#!/usr/bin/env python3
"""Cold Gold2 audit: ordinary arrows, labelled registrations, and exact ledger."""
from collections import Counter
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
import json
from os import close
from pathlib import Path


ROUND = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = spec_from_file_location(name, ROUND / path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


wire = load('audit_wire', 'attacks/neon2-quintet/gold2/physical/check.py')
writer = load('audit_writer', 'occupation-quintet-gold2/post-fy/replay.py')
ordinary = wire.old


def continuation(path, marks):
    caps = ((79, 8), (5, 4), (18, 3))
    if any(d > D or c > C for (d, c), (D, C) in zip(map(marks.get, path), caps)):
        return 0
    # Deliberately asymmetric: post-birth marks and terminal order matter.
    return len(path) + sum((i + 1) * (d + 10*c) for i, (d, c) in enumerate(map(marks.get, path)))


def registration(parents, marks):
    out = set()
    for v in marks:
        chain = (v,)
        while len(chain) <= 3 and chain[0] != 0:
            out.add(chain)
            chain = (parents[chain[0]], *chain)
    return out


def literal(p, e, parents, degrees, counts):
    tree = ordinary.physical_tree([(u, v) for v, u in enumerate(parents) if v], degrees)
    core, marks = ordinary.physical_core(tree)
    pre = ordinary.subsets(p, core, marks)
    paths = registration(parents, marks)
    total = sum(continuation(path, marks) for path in paths)
    nxt = len(tree)
    lam = e.snap.lam
    rows = {path: lam * pre['h3', tuple(sorted(path))][0]
            for path in paths if len(path) == 3 and continuation(path, marks)
            and marks[path[0]][0] >= 19}
    rates = Counter()
    changes = {}
    reg_drift = Q(0)
    for u, nbrs in enumerate(tree):
        for v in nbrs:
            # One literal source->recipient arrow, evaluated in the ordinary graph.
            rate = Q(1, len(nbrs))
            counts['arrows'] += 1
            if v not in changes:
                after = list(tree) + [{v}]
                after[v] = tree[v] | {nxt}
                ac, am = ordinary.physical_core(after)
                succ = ordinary.subsets(p, ac, am)
                birth = v not in marks
                ap = dict(enumerate(parents))
                if birth:
                    ap[v] = u
                full = registration(ap, am)
                seeded = {(v,)} | {(*path, v) for path in paths if len(path) < 3 and path[-1] == u} if birth else set()
                assert full - paths == seeded and paths <= full
                if birth:
                    counts['seeded prefixes'] += len(seeded)
                    counts['cutoff exits'] += sum(bool(continuation(path, marks)) and not continuation(path, am) for path in paths)
                delta = sum(continuation(path, am) for path in paths | seeded) - total
                assert delta == sum(continuation(path, am) - continuation(path, marks) for path in paths) + sum(continuation(path, am) for path in seeded)
                changes[v] = am, succ, delta
            am, succ, delta = changes[v]
            reg_drift += rate * delta
            rates[v] += rate
            for path in rows:
                mother, middle, tip = path
                motif = frozenset(path)
                key = ('h3', tuple(sorted(path)))
                if v in motif:
                    jump = succ.get(key, (0,))[0] - pre[key][0]
                    rows[path] -= rate * jump
                    if u in core[v] - motif and any((kind, tuple(sorted(motif | {u}))) in pre for kind in ('p4', 's3')):
                        rows[path] += (rate - center(p, marks[v])) * jump
                        counts['outgoing'] += 1
                    if v == middle and u in motif:
                        keep, = motif - {middle, u}
                        edge = ('g', tuple(sorted((middle, keep))))
                        rows[path] -= (rate - center(p, marks[v])) * (succ[edge][0] - pre[edge][0])
                        counts['incoming'] += 1
                elif v not in marks and u in motif:
                    added = sum(succ.get((kind, tuple(sorted(motif | {v}))), (0,))[0] for kind in ('p4', 's3'))
                    rows[path] -= rate * (succ.get(key, (0,))[0] - pre[key][0] + added)
                    counts['nonzero newborn'] += bool(added)

    aggregate = Q(0)
    for u, (d, c) in marks.items():
        gains = d-c + sum(Q(1, marks[v][0]) for v in core[u])
        assert rates[u] == gains
        aggregate += gains * changes[u][2]
        leaves = tree[u] - core[u]
        assert sum(rates[v] for v in leaves) == Q(d-c, d)
        aggregate += sum(rates[v] * changes[v][2] for v in leaves)
    assert aggregate == reg_drift

    for path, h in rows.items():
        r, m, t = path
        aff, slopes = e.triple(marks[m], marks[r], marks[t])
        observed = aff + sum(slope * sum(Q(1, degrees[v]) for v in core[u] - set(path))
                             for slope, u in zip(slopes, (m, r, t)))
        assert observed == h <= e.hcap
        assert h + max(0, -h) == max(0, h)
        counts['negative H3'] += h < 0
        counts['H3 rows'] += 1
        # Plant a donor-denominator error on the first middle gain transfer.
        g = ('g', tuple(sorted((m, t))))
        mutation = (Q(1, degrees[r]) - Q(1, degrees[m])) * (changes[m][1][g][0] - pre[g][0])
        counts['rejected donor-rate mutations'] += observed != h + mutation


def center(p, mark):
    return Q(p.word(p.ends[2].ids[mark]), p.h['S_r'])


def terminal_order():
    ext = writer.ext
    engine = ext.Pass.__new__(ext.Pass)
    rs = ext.states(7, 4)
    ts = ext.states(6, 3)
    ri, rm = ext.axes(rs)
    ti, _ = ext.axes(ts)
    engine.nt = len(ts)
    engine.t0 = ti[(2, 1)]
    calls = {}

    def step(marks, adj, gs, bs):
        calls[marks] = bs
        return 100*marks[0][0] + marks[0][1]

    engine.step = step
    pb = list(range(1, len(rs)*len(ts)+2))
    edge_gain = [0]*(len(rs)+1)
    edge_bud = [400+i for i in range(len(rs)+1)]
    engine.edge((5, 3), rs, rm, edge_gain, edge_bud, pb, True)
    for r in rs:
        i = ri[r]
        assert calls[(5, 3), r][0] == edge_bud[i] + pb[i*len(ts)+ti[(2, 1)]]
    source = [700+i for i in range(len(rs))] + [0]
    engine.singleton(rs, rm, [source[ri.get((d, c+1), -1)] for d, c in rs])
    for d, c in rs:
        j = ri.get((d, c+1))
        assert calls[((d, c),)][0] == (100*d+c+1+source[j] if j is not None else 0)
    print('ordered terminal controls:', 2*len(rs))


def knots(p, e):
    count = 0
    assert e.knots == p.kn and len(p.kn) == 20
    for host, pay, donor in product(((34, 11), (80, 3), (150, 3), (251, 17)),
                                    ((2, 1), (5, 3), (6, 2), (34, 7)), (2, 3, 5, 12, 251, 252, 1000)):
        d, c = host
        assert host not in p.ins['3m']
        small = pay[0] <= 5
        large = 6 <= pay[0] <= 34 and d <= 150 and c <= 10 and donor <= 5
        assert not (small and large)
        start = p.n2 + 20*p.ends[2].ids[host]
        z = p.g((d+1, c), pay) - p.g(host, pay)
        assert z == e.g((d+1, c), pay) - e.g(host, pay)
        vals = tuple(p.word(start+j) for j in range(20))
        scores = tuple(a-x*z for x, a in zip(p.kn, vals))
        x = Q(1, donor)
        i = next(j for j in range(19) if p.kn[j] <= x <= p.kn[j+1])
        w = (x-p.kn[i])/(p.kn[i+1]-p.kn[i])
        s = (1-w)*scores[i] + w*scores[i+1]
        assert s == e.fy(host, x)-x*z <= max(scores)
        if donor >= 12:
            cut = Q(1, 12)
            assert s <= max(e.fy(host, y)-y*z for y in (*[k for k in p.kn if k < cut], cut))
        for k, a in zip(p.kn, vals):
            assert e.fy(host, k) == a
            count += 1
    return count


def ledger(doc):
    lo, hi = (Q(doc[s]['certified_exact']) for s in ('lower', 'upper'))
    assert Q('.565713486333') <= lo <= hi <= Q('.565736986679')
    cert = doc['lower']
    scale, target = 2**40, Q(99, 175)
    h = Q(1084423030952180046412413909, 1704710415481443450880000000000000)
    small = Q(39526146756920663734311, 109951162777600000000000000000)
    large = Q(1046834122639903380009, 21990232555520000000000000000)
    e, m, lam = cert['err_per_core'], cert['m_max'], Q(cert['nominal'])
    assert (e, m, lam) == (19217, 1808603193693, Q(7071421, 12500000))
    required = e + (target-lam)*m
    credit = scale*(h+small+large)
    lower = lam + (credit-e)/m
    assert credit > required == Q(97537456765729, 87500000) and lower > target
    print('required ticks:', required, 'surplus ticks:', credit-required, float(credit-required))
    print('new lower:', lower, float(lower), 'strict gap:', float(lower-target))


def main():
    docs = json.loads((wire.ROOT/'source'/wire.CERT/'certificate_gold2.json').read_bytes())
    cert = docs['lower']
    assert cert['verify']['sha256'] == wire.DIGEST['lower']
    p = wire.Gold(wire.ROOT/'quintet_lower.bin')
    counts = Counter()
    try:
        with writer.Snapshot('lower') as s:
            e = writer.Rows(s, p, cert)
            assert (e.sig, e.hcap) == (-1, 408)
            parents = (None, 0, 1, 2, 2, 2, 3, 3, 1, 1, 1, 1, 1, 1)
            for ds in ((7, 19, 5, 18, 18, 2, 2, 2, 5, 6, 22, 2, 4, 2),
                       (3, 79, 4, 3, 3, 2, 2, 2, 5, 6, 22, 2, 4, 2),
                       (3, 38, 5, 4, 4, 2, 2, 2, 5, 6, 22, 2, 4, 2),
                       (3, 18, 5, 4, 4, 2, 2, 2, 5, 6, 22, 2, 4, 2)):
                literal(p, e, parents, ds, counts)
            literal(p, e, (None, 0, 1, 2, 2, 3), (4, 65, 4, 5, 5, 2), counts)
            for key in ('outgoing', 'incoming', 'nonzero newborn', 'cutoff exits', 'negative H3', 'rejected donor-rate mutations'):
                assert counts[key], key
            print('literal:', dict(counts), 'FY knot slots:', knots(p, e), 'sparse words:', len(p.reads))
            s.check_unchanged()
    finally:
        close(p.fd)
    terminal_order()
    ledger(docs)


if __name__ == '__main__':
    main()
