#!/usr/bin/env python3
"""Convert between the pickle written by xcharge.py and a plain-text certificate.

Text format (one item per line, values as C99 hex floats, i.e. exact dyadic rationals):
  # mode upper|lower
  # Dcap <int>   # c1 <int>   # cmax <int>
  # starts <block starts of L, space separated>     (f is constant in L on [start_i, start_{i+1}-1]; the last
                                                      block is [last start, infinity))
  # class <index> <l> <k>      (non-TOP classes)     # TOP <index>
  # value <float>              (the value at which the LP was solved; informational)
  f <L> <c> <S>  <hex>         detailed state: L = block start, c, S = class indices joined by ',' ('-' if c = 0)
  F <L> <c>  <hex>             coarse state
  G <hex>                      beyond state (upper mode)
  k <b> <hex>                  kappa(b)
  nu <a> <b> <hex>             nu(a, b)
  Kp <hex>                     (upper mode; informational, the verifier recomputes max(0, max_b kappa(b)))

usage: xdump.py totext CERT.pkl OUT.txt | topkl CERT.txt OUT.pkl
"""
import sys, pickle
from xcharge import Scheme


def totext(src, dst):
    d = pickle.load(open(src, 'rb'))
    S = Scheme(*d['params'])
    with open(dst, 'w') as o:
        o.write(f"# mode {d['mode']}\n# Dcap {S.Dcap}\n# c1 {S.c1}\n# cmax {S.cmax}\n")
        o.write("# starts " + ' '.join(map(str, S.starts)) + "\n")
        for i, (l, k) in enumerate(S.classes):
            o.write(f"# class {i} {l} {k}\n")
        o.write(f"# TOP {S.TOP}\n# value {d['value']!r}\n")
        for key, v in d['sol'].items():
            h = float(v).hex()
            if key == 'G':
                o.write(f"G {h}\n")
            elif key == 'Kp':
                o.write(f"Kp {h}\n")
            elif key[0] == 'f' and len(key) == 4:
                o.write(f"f {key[1]} {key[2]} {','.join(map(str, key[3])) if key[3] else '-'} {h}\n")
            elif key[0] == 'f':
                o.write(f"F {key[1]} {key[2]} {h}\n")
            elif key[0] == 'k':
                o.write(f"k {key[1]} {h}\n")
            elif key[0] == 'nu':
                o.write(f"nu {key[1]} {key[2]} {h}\n")
            else:
                raise ValueError(key)


def topkl(src, dst):
    hdr = {}; sol = {}; starts = None
    for line in open(src):
        p = line.split()
        if not p:
            continue
        if p[0] == '#':
            if p[1] == 'starts':
                starts = [int(x) for x in p[2:]]
            elif p[1] in ('mode', 'Dcap', 'c1', 'cmax', 'value'):
                hdr[p[1]] = p[2]
            continue
        v = float.fromhex(p[-1])
        if p[0] == 'G':
            sol['G'] = v
        elif p[0] == 'Kp':
            sol['Kp'] = v
        elif p[0] == 'f':
            S = () if p[3] == '-' else tuple(int(x) for x in p[3].split(','))
            sol[('f', int(p[1]), int(p[2]), S)] = v
        elif p[0] == 'F':
            sol[('f', int(p[1]), int(p[2]))] = v
        elif p[0] == 'k':
            sol[('k', int(p[1]))] = v
        elif p[0] == 'nu':
            sol[('nu', int(p[1]), int(p[2]))] = v
    Lspec = starts if starts != list(range(len(starts))) else starts[-1]
    pickle.dump(dict(mode=hdr['mode'], params=(int(hdr['Dcap']), int(hdr['c1']), Lspec, int(hdr['cmax'])),
                     value=float(hdr['value']), sol=sol), open(dst, 'wb'))


if __name__ == '__main__':
    {'totext': totext, 'topkl': topkl}[sys.argv[1]](sys.argv[2], sys.argv[3])
