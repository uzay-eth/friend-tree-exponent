#!/usr/bin/env python3
"""Widths are exact dyadic differences of the logged binary64 eigenvalues."""
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
import sys

TAGS = 'reference zeroA zeroB zeroC clipA clipB clipAB clipC zeroABC coupled'.split()

def row(tag):
    ends, time, rss, last = [], [], [], []
    for mode in (2, 1):
        p = Path(f'{tag}_m{mode}.log')
        if not p.exists():
            return None
        lines = p.read_text().splitlines()
        end = [r for r in lines if r.startswith('FLOAT ')]
        tm = [r for r in lines if r.startswith('TIME ') and r.endswith('exit=0')]
        if len(end) != 1 or len(tm) != 1:
            return None
        ends.append(end[0].split('lam=')[1])
        time.append(tm[0].split('wall=')[1].split()[0])
        rss.append(tm[0].split('peak_KiB=')[1].split()[0])
        last.append([float(r.split('lam=')[1].split()[0]) for r in lines if r.startswith('SWEEP ')])
    assert all(len(l)==18 for l in last)
    width = Fraction(float(ends[1])) - Fraction(float(ends[0]))
    with localcontext() as ctx:
        ctx.prec = 80
        value = Decimal(width.numerator) / Decimal(width.denominator)
    shift = max(abs(l[-1]-l[-2]) for l in last)
    return f'| {tag} | {ends[0]} | {ends[1]} | {value} | {"/".join(time)} | {"/".join(rss)} | {shift:.2g} |'


if __name__ == '__main__':
    if sys.argv[1:] != ['--report']:
        for tag in sys.argv[1:] or TAGS:
            print(row(tag) or f'{tag}: pending')
    else:
        path = Path('REPORT.md')
        doc = path.read_text()
        rows = {tag: val for tag in TAGS if (val := row(tag))}
        table = ['| Case | λ lo | λ hi | Exact binary64 width | Wall lo/hi (s) | Peak lo/hi (KiB) | Last shift |',
                 '|:--|--:|--:|--:|--:|--:|--:|', *rows.values()]
        if len(rows) < len(TAGS):
            table += ['', 'Pending pairs: ' + ', '.join(t for t in TAGS if t not in rows) + '.']
        def width(tag):
            return Decimal(rows[tag].split('|')[4].strip())
        for t in ('zeroA','zeroB','zeroC','clipA','clipB','clipAB','clipC','zeroABC'):
            if 'reference' in rows and t in rows:
                table += ['', f'`{t}` narrows the reference width by {(1-width(t)/width("reference"))*100:.3f}%.']
        if 'coupled' in rows:
            w = width('coupled')
            base = Decimal('0.00013360681999630142')
            pure = Decimal('0.00012800250560285331')
            gain = (1-w/base)*100
            relation = 'below' if w < Decimal('0.0001') else 'at or above'
            table += ['', f'**Coupled endpoint-only BIG: {gain:.3f}% width reduction versus the earlier no-L5 baseline.** '
                          f'The absolute reduction is {(base-w)/(base-pure):.2f}× that of the earlier pure-inner L5 run. '
                          f'The measured width is {relation} `1.00e-4` at D0=12. '
                          'This result does not establish a D0=20/22 width.']
        if 'coupled' in rows and 'clipAB' in rows:
            gap = width('coupled')-width('clipAB')
            table += ['', f'**Coupled minus optimistic A+B width: {gap:.9g}.** '
                          f'The optimistic A+B run narrows the no-L5 baseline by '
                          f'{(1-width("clipAB")/Decimal("0.00013360681999630142"))*100:.3f}%.']
        for marker, body in [('RESULTS', '\n'.join(table)),
                             ('DERIVATION', Path('report_notes.md').read_text().strip())]:
            before, tail = doc.split(f'<!-- BEGIN {marker} -->\n', 1)
            _, after = tail.split(f'<!-- END {marker} -->', 1)
            doc = before + f'<!-- BEGIN {marker} -->\n' + body + f'\n<!-- END {marker} -->' + after
        path.write_text(doc)
