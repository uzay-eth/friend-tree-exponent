#!/usr/bin/env python3
"""Row/memory counts and coverage with fixed and native pure-inner exclusions."""
from pathlib import Path
from runpy import run_path

base = run_path(str(Path(__file__).resolve().parent.parent / 'mixed_l5_study/count_variants.py'))
cnt, shapes, size = (base[k] for k in ('cnt', 'shapes', 'size'))
for d in (6, 7, 8):
    for b in (1, 2):
        rows = sum(cnt(s, 'apex', d, b) for s in shapes if size(s) == 5)
        row_bytes = 8*rows
        print(f'D={d} buckets={b} rows={rows} row_bytes={row_bytes} '
              f'payload_D0=12:{row_bytes+8958236} D0=20:{row_bytes+8961884} D0=24:{row_bytes+8964284}')

source = (Path(__file__).resolve().parent.parent / 'mixed_l5_study/count_variants.py').read_text()
source = source[source.index('def count_l4_coverage():'):source.index("if __name__ == '__main__':")]
source = source.replace('def count_l4_coverage():', 'def count_l4_coverage(cap, exclude):')
source = source.replace('d<6, d<12 and k+h>a+h_inside, 6<=d<12, deg>1',
                        'd<cap, d<12 and k+h>a+h_inside, cap<=d<12, deg>1, d<exclude')
source = source.replace('if all(f[0] for f in slots):', 'if all(f[4] for f in slots):')
source = source.replace('for small,_,_,it in slots', 'for small,_,_,it,_ in slots')
source = source.replace('for small, ext, midrange, it in slots', 'for small, ext, midrange, it, _ in slots')
source = source.replace("print(f'{kind} mixed row/receiver incidences: {dict(c)}')",
                        "print(f'cap={cap} exclude_all_below={exclude} {kind}: {dict(c)} ({c[\"supported\"]/c[\"total\"]:.6%})')")
exec(source, base)
for cap, exclude in ((6, 6), (7, 6), (7, 7)):
    base['count_l4_coverage'](cap, exclude)
