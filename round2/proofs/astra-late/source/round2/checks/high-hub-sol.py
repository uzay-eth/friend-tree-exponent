#!/usr/bin/env python3
"""Cold high101 checks: live coefficients, seam multiplicity, and exact ledger."""

import argparse
from fractions import Fraction as Q
from hashlib import file_digest
from math import comb
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
O = ROOT / 'round2/occupation-neon'
H = O / 'high-hub-limit'
BASE = O / 'pressure-memory/fork-lifetime/yt-release'
COMMIT = 'aa511af202d0d64cd1a08f900e0516965d9db952'
PINS = {
    H/'bulk.py': '4104ee4c72a8c9a0579fb5b7d24a84bc7764ca1e10ae14e4799cb7e49f943c66',
    H/'check.py': 'e4c6b2aa2f53370c0b44790ff1e6c1e030c25368d03b41c9a395628ae3586798',
    H/'physical_check.py': 'a1a1e7f0b45f587b4a98e762ea3f6edb3cc424275382d607ff760c560c0875a4',
    H/'all-arms/exact_receipt.py': 'fed0796b2e14153b56d5d53975d5940ed5c183744573fe9e95f00111ed22a39c',
    O/'fringe/charges.py': '44ed53fc7a0fc9ad25d941cd4160fe43aa9d00da334adb4e7de91a9a1a7c97e8',
    O/'higher/extract.py': '9f40a50a3623756b277243d192f672e76c45b4b9784715ed4d0d8fb27d36131c',
    O/'higher/older/replay.py': '2e600622e108be8791912577adcbb9291e17e3412763abcf5e7ca5c03130fb2f',
    ROOT/'round2/sol-kernel-verifier/verify.py': '3b2e573568c9453e408fb401a5ac1e103da71f7caf071a62b28f73575e662c38',
    ROOT/'round2/attacks/neon1/audit.py': 'd548a9aafa899ecab2b00a7f6ad0a4e6f828600088780f0a15f08927f17f0bfe',
    ROOT/'round2/finite/residual_entries.py': '59f4ef6f915d6d7d932076dfe9170d7e6e5fa5d819cf3d0c11baf7932baa1277',
    BASE/'README.md': '22d4529edad049bae698d08c1cbdedbefeff902e456e2b56e96d4e5bdc808d4e',
    BASE/'verify.py': 'fac19a0ddbc53b1dfe77bb799f596e5b1a70c123f8bcd89fe6d5b1652693f6d8',
    O/'interpolate/edge/higher/transport/replay.py': 'a92a5e830ac6ed092eeb2be1c73709244f66aef9df1a227cf969785467a7e870',
}
M, K, CAP = 101, 12, 7
LO, HI = Q('.565673570322'), Q('.565775')
J = (2, 1, ())
Z = (Q(0),)*3
GRID = 10**23
DEN = Q('1.357493815758092773'), Q('1.367180021074064915')
D18_TRIAL = Q('.565673'), Q('.565775')
EDGE_TRIAL = Q('.565677930574'), Q('.565770634792')
BASE_NUM = (
    Q('21110394937523758224695119010393/8925000000000000000000000000000000000')
    + Q('3092299274719269569703/10328365000000000000000000000'),
    Q('36416704104097467168533889025283/14025000000000000000000000000000000000')
    + Q('9405336866046559940037/20656730000000000000000000000'),
)


def digest(path):
    with path.open('rb') as f:
        return file_digest(f, 'sha256').hexdigest()


def physical_registrations():
    # Check two equivalent representations, keyed by identity, with equal
    # shapes, several births at each seam degree, and permanent arm/arity exits.
    def history(count):
        yield ('gain', 100)
        for arm in range(count):
            yield ('birth', arm)
        if count:
            yield ('size', 0)
        yield ('gain', 101)
        for arm in range(count, K+2):
            yield ('birth', arm)
            yield ('size', arm)
        yield ('gain', 10**70)

    steps = 0
    for n in range(K+2):
        d, size, registered, entrance = 2, {}, set(), set()
        for event, i in history(n):
            if event == 'gain':
                if d < M <= i and len(size) <= K:
                    entrance = {a for a in size if size[a] <= CAP}
                d = i
            elif event == 'birth':
                size[i] = 2
                registered.add(i)
                if d >= M:
                    assert i not in entrance
                    entrance.add(i)
            else:
                size[i] = CAP+1
            direct = {a for a, s in size.items() if s <= CAP and d >= M and len(size) <= K}
            for ledger in (registered, entrance):
                assert direct == {a for a in ledger if size[a] <= CAP and d >= M and len(size) <= K}
            steps += 1
    print('PASS labelled seam/birth/size/arity control events', steps, flush=True)


def receipt(path, sign):
    header, *lines = path.read_text().splitlines()
    assert header == f'witness {sign} M={M} T={M} K={K} cap={CAP} LO={LO} HI={HI}'
    result = {}
    for line in lines:
        side, role, *fields = line.split()
        *parts, total = map(Q, fields)
        assert len(parts) == 3 and 0 < sum(parts) == total and min(parts) > 0
        assert (side, role) not in result
        result[side, role] = tuple(parts)
    assert set(result) == {(s, r) for s in ('floor', 'ceiling') for r in ('all', 'old', 'young')}
    for a, b, c, d, e, f in zip(*(result[s, r] for s in ('floor', 'ceiling') for r in ('all', 'old', 'young'))):
        assert max(a, b+c) <= min(d, e+f)
    return result


def envelope(p, k, upper):
    if p is None:
        return Z
    inf, rows = p
    aa, bb = rows[k]
    ans = []
    for f, a, b in zip(inf, aa, bb):
        assert a[0] == f and b[:2] == (0, 0)
        ends = []
        for u in (Q(0), Q(k, 2)):
            num = [a[j]+u*b[j]-(comb(3,j)*f if j <= 3 else 0) for j in range(1, 8)]
            ends.extend(sum(num[t]*Q(comb(i,t), M**t*comb(6,t)) for t in range(i+1)) for i in range(7))
        # 1/(1+p)^3 is at most one; only the adverse sign is needed.
        z = f + (max(Q(0), *ends) if upper else min(Q(0), *ends))/M
        ans.append(z if upper else max(Q(0),z))
    return tuple(ans)


def independent(blob, sign, expected):
    # Use pinned polynomial extraction and arm topology; do not call bulk.run,
    # stage_bounds, or band_price. Verify extra true-degree points directly.
    sys.path.insert(0, str(H))
    from check import Rewards, arms, v

    states, arcs = arms(CAP)
    r = Rewards(blob, states)
    assert (r.c.sign, r.c.c.D, r.c.c.K, r.c.c.S, r.c.c.lam) == (sign, 18, 17, 10**11, D18_TRIAL[sign == -1])
    polys = {x:r.polynomials(x) for x in states}
    checks = 0
    for x, p in polys.items():
        for d, k in ((101,1),(101,12),(102,12),(151,6),(10**71+37,12)):
            row = r.row(x,d,k)
            want = (Z,Z) if p is None else tuple(
                tuple(sum(c*Q(1,d)**i for i,c in enumerate(f))/((1+Q(1,d))**3) for f in a)
                for a in p[1][k])
            assert row == want, (sign,d,k,x)
            checks += 1
    print('PASS true-degree polynomial controls', sign, checks, flush=True)
    data = {x:(Q(v.size(x))-Q(1,v.degree(x)), v.plus(x), tuple(arcs[x].items())) for x in states}
    for upper in (False, True):
        s, norm = (LO, HI) if upper else (HI, LO)
        extreme = max if upper else min
        def rounded(vals):
            v0 = extreme(vals)*GRID
            n, d = v0.as_integer_ratio()
            return Q(-(-n//d) if upper else n//d, GRID)

        def qsolve(immediate, fixed, channels):
            return tuple(rounded((q+sum(rate*vec[i] for rate,vec in edges))/(fixed+sum(rate for rate,_ in edges))
                                 for edges in channels) for i,q in enumerate(immediate))

        arm = {}
        host = {}
        for k in range(K,-1,-1):
            beta = (Q(M-k-1,M), Q(1))
            if k < K:
                future, newborn = host[k+1], arm[k+1,J]
                host[k] = tuple(rounded(b*(a+z)/(s+b) for b in beta) for a,z in zip(future,newborn))
            else:
                host[k] = Z
            if k:
                for x in states:
                    intrinsic, plus, edges = data[x]
                    acc = tuple(q+sum(rate*arm[k,y][i] for y,rate in edges)
                                for i,q in enumerate(envelope(polys[x], k, upper)))
                    obs = arm.get((k,plus), Z) if upper else Z
                    nxt = arm.get((k+1,x), Z)
                    arm[k,x] = qsolve(acc,s+intrinsic,
                                      [((b,nxt),(p,obs)) for b in beta
                                       for p in ((Q(0),Q(1,M)) if upper else (Q(1,M),))])
        # Drive a different traversal through the degree slabs.
        ahead, hosts = arm, {k:(v0,Z,v0) for k,v0 in host.items()}
        for d in range(M-1,1,-1):
            here, cur = {}, {}
            p = Q(1,d)
            for k in range(min(K,d-1),-1,-1):
                leaf = d-k-1
                b = leaf*p
                if k:
                    for x in states:
                        intrinsic, plus, edges = data[x]
                        plus_v, birth_v = here.get((k,plus),Z), here.get((k+1,x),Z)
                        a = tuple(sum(rate*here[k,y][i] for y,rate in edges)+p*plus_v[i]+b*birth_v[i]
                                  for i in range(3))
                        g0 = leaf+Q(1,v.degree(x))
                        here[k,x] = qsolve(a,s+intrinsic+p+b,
                                           [((g,ahead[k,x]),) for g in (g0,g0+Q(k,2))])
                new = here.get((k+1,J),Z)
                resolv = cur.get(k+1,(Z,Z,Z))
                cur[k] = tuple(qsolve(tuple(b*(x+y) for x,y in zip(cont, new if j < 2 else Z)),s+b,
                                      [((g,hosts[k][j]),) for g in (Q(leaf),leaf+Q(k+1,2))])
                               for j,cont in enumerate(resolv))
            ahead,hosts = here,cur
        for role, vec in zip(('all','old','young'), hosts[0]):
            want = expected['ceiling' if upper else 'floor',role]
            got = tuple(norm*q for q in vec)
            assert got == want, (sign,upper,role,got,want)
        print('PASS independent exact Bellman', sign, 'ceiling' if upper else 'floor', '9 endpoints', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('lo',type=Path)
    p.add_argument('hi',type=Path)
    p.add_argument('receipts',type=Path)
    p.add_argument('--solve',choices=('lo','hi'))
    args = p.parse_args()
    for path, want in PINS.items():
        assert digest(path) == want, path
    physical_registrations()
    base = subprocess.check_output([sys.executable,'-EBu',str(BASE/'verify.py'),'--payloads',
                                    str(args.lo.resolve()),str(args.hi.resolve())],cwd=ROOT,text=True)
    assert 'PASS original D18 payload SHA-256s' in base
    bounds = [receipt(args.receipts/f'{k}.log',sign) for k,sign in (('lo',1),('hi',-1))]
    for i,side in enumerate(('lower','upper')):
        sign = 1-2*i
        prev = EDGE_TRIAL[i]+sign*BASE_NUM[i]/DEN[i]
        assert f'exact {side} {prev}\n' in base
        floor, ceil = [sum(bounds[i][direction,'all']) for direction in ('floor','ceiling')]
        endpoint = prev+sign*floor/DEN[i]
        assert endpoint == EDGE_TRIAL[i]+sign*(BASE_NUM[i]+floor)/DEN[i]
        a = endpoint*10**15
        n,d = a.as_integer_ratio()
        tick = n//d if i == 0 else -(-n//d)
        print(side, 'floor', floor, 'masked ceiling',ceil,'exact',endpoint,
              'outward',f'{tick//10**15}.{tick%10**15:015d}',flush=True)
    if args.solve:
        i = int(args.solve=='hi')
        independent((args.lo,args.hi)[i],1-2*i,bounds[i])


if __name__ == '__main__':
    main()
