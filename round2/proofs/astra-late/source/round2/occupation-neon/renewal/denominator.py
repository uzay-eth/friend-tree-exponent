"""Extract integer upper charges for young vertices from pinned D18 payloads."""
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'attacks' / 'neon1'))
from audit import Cert


def charges(c, w):
    iw = c.child[w]
    # Charge the same motifs to the same vertices as the GMP comparison.
    # k<=1 has no apex charges. Yb contributions at a common bottom vertex
    # come from at most k_apex-1 siblings; preserve this multiplier by row.
    assert w[1] <= 1
    return [int(c.f[c.vi(*w)]),
            max(int(c.g[c.state[p], c.state[w]]) for p in c.par),
            max(0, int(c.hvtab[..., iw].max())),
            max(0, int(c.p4tab[..., iw].max())),
            max(0, max((u[1]-1)*int(c.ybtab[i, ..., iw].max())
                       for u, i in c.ybu.items()))]


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('certs', nargs='+', type=Path)
    for path in p.parse_args().certs:
        c = Cert(path)
        assert (c.D, c.K, c.LF, c.KF, c.S) == (18, 17, 14000, 180, 10**11)
        for w in ((1, 0, 1), (2, 0, 1), (3, 0, 1), (4, 0, 1),
                  (0, 1, 1), (1, 1, 1), (2, 1, 1)):
            a = charges(c, w)
            print(c.mode, *w, *a, sum(a), flush=True)
