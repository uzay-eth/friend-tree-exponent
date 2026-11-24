"""Independent physical-arrow and budget fixtures for the two frozen Q3 tails.

Run from the Astra root: taskset -c 12 nice -n 10 python3 -IB <this file>.
No writer code is used to generate the trees, physical rates or signed numerators.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha1
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / 'open/age-order/global-transport'
PINS = {
    'beta-clocks/README.md': 'de9eb915b26439eae7435bc3b3d8e10da5434b06',
    'beta-clocks/check.py': 'dbf17d8c48924078682545fab8f3614a9b1e091a',
    'beta-clocks/ledger.py': 'c4bfb0cb24279a0d615c52dd4089f9b0822d38c1',
    'sixteenth/README.md': '4693cf726cd415902da3223509f345de7dd2dd5b',
    'sixteenth/check.py': '44df8d01d2a27455326e0dd351a05580740ffd41',
    'sixteenth/replay.log': '4873f5ce59c9a902bfbdad5341cd5cb9551ecf0e',
    'quarter/check.py': '5dab19ad0224b96b63bcaf0a8915c623609609c7',
    'birth-window/sharp/POISSON-SPLICE.md': '2cfb2eadecfc04cf50cf69f0b5ec30c02593ce90',
    'birth-window/sharp/poisson-splice.py': 'c7afbae24076d75bfb534c6c6a5c4f7b3f27708f',
}


def tree(parents):
    adj = [set() for _ in range(len(parents)+1)]
    for v, u in enumerate(parents, 1):
        adj[u].add(v)
        adj[v].add(u)
    return adj


def arrows(adj):
    for u, ns in enumerate(adj):
        for v in ns:
            yield u, v, Q(1, len(ns))


def birth(adj, v):
    n = len(adj)
    return [ns | ({n} if u == v else set()) for u, ns in enumerate(adj)] + [{v}]


def physical():
    rows = Counter()
    for n in range(3, 8):
        for parents in product(*(range(i) for i in range(1, n))):
            adj = tree(parents)
            for u, ns in enumerate(adj):
                d = len(ns)
                leaves = {v for v in ns if len(adj[v]) == 1}
                for j in range(1, min(len(leaves), d-1)+1):
                    for chosen in combinations(sorted(leaves), j):
                        gain = loss = tracked_gain = Q(0)
                        for sender, target, w in arrows(adj):
                            nxt = birth(adj, target)
                            after = {v for v in chosen if len(nxt[v]) == 1}
                            if target == u:
                                gain += w
                                tracked_gain += w*(sender in chosen)
                                after.add(n)
                                assert len(nxt[u]) == d+1 and len(after) == j+1
                            elif len(after) != j:
                                assert target in chosen and sender == u
                                assert len(nxt[u]) == d and len(after) == j-1
                                loss += w
                                if j == 1:
                                    assert len(nxt[target]) == 2 and target in nxt[u]
                                    rows['last loss'] += 1
                        assert j <= gain <= d and tracked_gain == j and loss == Q(j,d)
                        assert loss <= Q(j,j+1) and j+Q(j,j+1) <= n
                        rows['positive stock'] += 1
                for v in ns:
                    k = len(adj[v])
                    if k < 2:
                        continue
                    gain = ret = grow = Q(0)
                    for sender, target, w in arrows(adj):
                        nxt = birth(adj, target)
                        up = len(nxt[u]) - d
                        side = len(nxt[v]) - k
                        assert up + side <= 1
                        gain += w*up
                        grow += w*side
                        if sender == v and up:
                            ret += w
                        if up:
                            assert nxt[n] == {u}
                    assert d+k <= n and ret == Q(1,k) and ret <= gain <= d
                    assert Q(1,d) <= grow <= k-1+Q(1,d) and gain+grow <= n
                    rows['specified return'] += 1
                    rows['fresh degree-one return'] += d == 1
                if d != 1:
                    continue
                v, = ns
                k = len(adj[v])
                j = sum(len(adj[w]) == 1 for w in adj[v])-1
                classes = Counter()
                for sender, target, w in arrows(adj):
                    nxt = birth(adj, target)
                    if target == u:
                        assert sender == v and len(nxt[v]) == k and len(nxt[u]) == 2
                        classes['promote'] += w
                    else:
                        h = len(nxt[v])
                        i = sum(len(nxt[w]) == 1 for w in nxt[v])-1
                        classes[h-k,i-j] += w
                assert classes['promote'] == Q(1,k)
                assert classes[0,-1] == Q(j,k)
                assert j+1 <= classes[1,1] <= Q(k+j+1,2)
                assert sum(classes.values()) == n
                rows['upper pending stock'] += 1
    print('PASS literal physical arrows:', dict(rows), flush=True)


def serial_and_parallel(adj, M, U, L):
    """Full two-jump arrow scan with a signed synthetic conditional payoff.

    Assign 0 to noninteracting terms (which cancel), ±g to serial, and
    choose the unfavorable legal sign/price on each adjacent parallel tree.
    """
    n = len(adj)
    d = list(map(len, adj))
    r = [sum((Q(1,d[v]) for v in ns), Q(0)) for ns in adj]
    s = [r[u]/(d[u]*(d[u]+1)) for u in range(n)]
    g = M-U[2,0]
    edges = {}
    for u, ns in enumerate(adj):
        for v in ns:
            if s[u] <= s[v]:
                continue
            par = birth(birth(adj,u),v)
            ku, kv = len(par[u]), len(par[v])
            ju = sum(len(par[w]) == 1 for w in par[u])-1
            edges[u,v] = max(0,U[ku,ju]-L[kv])
    total = Counter()
    for _, u, a in arrows(adj):
        for _, v, z in arrows(birth(adj,u)):
            value = (g if v == n else edges.get((u,v), 0)-edges.get((v,u), 0))
            if value:
                total[u,v] += a*z*value
    cost = sum((s[u]-s[v])*v0 for (u,v),v0 in edges.items())
    assert sum(total.values()) == g*sum(d[u]*s[u] for u in range(n))-cost
    return s, edges, sum(total.values())


def budgets():
    six = run_path(str(T/'sixteenth/check.py'))
    b = 20
    x = Q(1,16)
    M,L = six['lower'](b,x,x,80,9,12)
    U = six['tips'](b,x,x,32,18)
    g = M-U[2,0]
    fixture = [ [0]*(b-2), [i for i in range(b-2)],
                [0,1,1]+[0]*8+[2]*7,
                [0,1,2,3,4,4]+[0]*6+[5]*6,
                [0,1,2,3,4,5]+[0]*6+[6]*6 ]
    for parents in fixture:
        adj = tree(parents)
        s, edges, score = serial_and_parallel(adj,M,U,L)
        assert len(adj) == b-1 and score > 0
        # Leaving out the serial arrows must not pass even on a star.
        assert score - g*sum(len(ns)*s[u] for u,ns in enumerate(adj)) <= 0
        assert all(len(adj[u]) >= 1 and len(adj[v]) >= 1 for u,v in edges)
    print(f'PASS {len(fixture)} full two-arrow star/long-chain/high-degree fixtures; missing serial control rejected', flush=True)

    beta = run_path(str(T/'beta-clocks/check.py'))
    ledger = run_path(str(T/'beta-clocks/ledger.py'))
    for lo,hi,tar in beta['BOXES']:
        lo,hi = Q(lo),Q(hi)
        a,B,f,ft = map(Q,tar)
        cap = beta['price_caps'](lo)
        M = beta['product_floor'](lo,hi)
        P = ledger['prepare'](M,cap,lo)
        assert ledger['certify'](P,(a,B,f,ft))
        sites,bills,sparse,rich,high = P
        assert any(bill > 0 for bill in bills) and not ledger['certify'](P,(a,Q(1,10**12),f,ft))
        assert (a-f)*16+sparse > 2*a+B and rich > 2*a+B and high/33 <= B
        # Stress the analytic inequalities beyond the cutoff at literal core/leaf
        # ratios (both receive ends). The proof covers every integer d, not just this set.
        g,u2 = M-cap(2,lo),cap(2,lo)
        for d in (32,33,63,64,127,128,1024,10**6):
            for l in {0,1,d//2-1,d//2,d//2+1,d-1,d}:
                c = d-l
                for R in (0,Q(c,2)):
                    s = (l+R)/Q(d*(d+1))
                    t = Q(1,2*d)
                    left = (g-high)*d*s+l*(g*t-u2*max(0,t-s))
                    assert s*high <= B
                    if l == d:
                        assert left > 0
                    elif l <= d/2:
                        assert left >= sparse and left+c*(a-f) > 2*a+B
                    else:
                        assert left >= rich and left+c*(a-f) > 2*a+B
        print(f'PASS beta exact tariff [{lo},{hi}], {len(sites)} retained rows, base-reduction control and large-degree seams',flush=True)


if __name__ == '__main__':
    for rel, pin in PINS.items():
        raw = (T/rel).read_bytes()
        assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest() == pin, rel
    print(f'PASS {len(PINS)} audited path pins', flush=True)
    physical()
    budgets()
