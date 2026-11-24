"""Independent exact arithmetic for the D0=20 review; read-only on inputs."""
from fractions import Fraction as F
from pathlib import Path
import mmap
import struct
import sys

import b3f18_sol_trees as trees
import tree_identity_check as decoder

# Reuse the two independent Python evaluators, but use a lazy little-endian
# backing store, so no 11 GB copy is allocated just to examine selected rows.
original_load = decoder.load_cert


def lazy_load(path):
    with open(path, 'rb') as fp:
        data = mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ)

    class Rows:
        def fromfile(self, fp, count):
            self.offset, self.count = fp.tell(), count
            fp.seek(8 * count, 1)
            assert fp.tell() <= len(data)

        def __getitem__(self, i):
            assert 0 <= i < self.count, (i, self.count)
            return struct.unpack_from('<q', data, self.offset + 8*i)[0]

    old = decoder.array
    from types import SimpleNamespace
    decoder.array = SimpleNamespace(array=lambda _: Rows())
    try:
        return original_load(path)
    finally:
        decoder.array = old


def product_audit():
    lo, hi = F(565680, 10**6), F(565768, 10**6)
    a, b = F(146041684821, 184219674916), F(4*146044513514, 115454594604)
    assert lo == F(7071, 12500) and hi == F(70721, 125000)
    assert hi-lo == F(88, 10**6) < F(1, 10**4)
    assert a > F('0.792758') and b < F('5.059808')
    print(f'lo={lo}; hi={hi}; width={hi-lo}; strict decimal margins: '
          f'{a-F("0.792758")}, {F("5.059808")-b}')
    law, pl, ph = {(-1,0,0): F(1)}, F(1), F(1)
    for n in range(3, 9):
        new = {}
        mean = F(0)
        for par, prob in law.items():
            adj = [[] for _ in par]
            for v, p in enumerate(par[1:], 1):
                adj[v].append(p); adj[p].append(v)
            mean += prob * sum(len(ns) >= 2 for ns in adj)
            if n < 8:
                for w, ns in enumerate(adj):
                    new[par+(w,)] = prob * sum(F(1, n*len(adj[v])) for v in ns)
        assert sum(law.values()) == 1
        assert a*pl <= mean <= b*ph
        print(f'n={n}: {a*pl} <= {mean} <= {b*ph}')
        law, pl, ph = new, pl*(1+lo/n), ph*(1+hi/n)


if __name__ == '__main__':
    product_audit()
    for path in sys.argv[1:]:
        C = lazy_load(path)
        side = 'lo' if C['mode'] == 2 else 'hi'
        assert Path(path).name == f'b3f20_{side}.bin'
        names = 'D0 KT LF KF DP4 KP4 DYT KYT DYB KYB DS3 KS3 S S_r P Q'.split()
        hdr = tuple(C[x] for x in names)
        assert hdr == (20,19,16000,200,12,9,11,8,12,9,10,7,10**11,10**9,
                       565680 if side == 'lo' else 565768,10**6)
        pot = [C[x].count for x in ('f','g','hV','hC','qP4','qYt','qYb','qS3')]
        total = (C['qS3'].offset-96)//8 + pot[-1]
        assert (sum(pot),total,96+8*total) == (1398838433,1398860120,11190881056)
        print(f'{side}: {dict(zip(names,hdr))}; potential lengths={pot}; '
              f'potential sum={sum(pot)}; total={total}', flush=True)
        # Substitute the D0=20 comparison constants into the prior referee's
        # parameterized tree, charge and transfer checker without editing it.
        source = Path(trees.__file__).read_text()
        for a, b in zip('113465364604 184124847315 146042017510 115462837236 195820675182 146044806770 584179227080'.split(),
                        '113472636536 184219674916 146041684821 115454594604 200854054323 146044513514 584178054056'.split()):
            assert a in source
            source = source.replace(a, b)
        gl = {'__file__': trees.__file__, '__name__': 'audit_d20'}
        exec(compile(source, trees.__file__, 'exec'), gl)
        gl['run'](C)
        decoder.load_cert = lazy_load
        original_argv = sys.argv
        sys.argv = [decoder.__file__, path, '30']
        try:
            decoder.main()
        except SystemExit as exc:
            assert exc.code == 0
        finally:
            sys.argv = original_argv
