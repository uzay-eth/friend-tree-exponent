#!/usr/bin/env python3
"""Build a Markdown table from complete paired diagnostic logs only."""
from pathlib import Path
import re


def extract(path):
    if not path.exists():
        return None
    s = path.read_text()
    lam = re.search(r'^FLOAT .* lam=(\S+)$', s, re.M)
    time = re.search(r'^TIME wall=(\S+) user=(\S+) sys=(\S+) peak_KiB=(\d+) exit=0$', s, re.M)
    lams = [float(x) for x in re.findall(r'^SWEEP .* lam=(\S+)', s, re.M)]
    if not lam or not time or len(lams) < 3:
        return None
    wall, user, sys, rss = time.groups()
    return float(lam[1]), float(wall), float(user) + float(sys), int(rss), abs(lams[-1] - lams[-3])


print('| D0 | Case | λ lo | λ hi | Width | Reduction | Wall lo/hi (s) | Peak RSS lo/hi (KiB) | Max last-two-sweep shift |')
print('|---:|:---|---:|---:|---:|---:|:---|:---|---:|')
for d in (12, 14, 16, 18, 20, 24, 25, 30):
    base = None
    for case in ('base', 'center', 'd5-6', 'd5-7', 'd5-8', 'h-center', 'l4-large',
                 'l5-6', 'l5-6-no-inner-slack', 'l5-6-no-mixed-slack', 'l5-7', 'l5-8'):
        hi, lo = [extract(Path(f'd{d}_{case}_m{m}.log')) for m in (1, 2)]
        if not hi or not lo:
            continue
        width = hi[0] - lo[0]
        if case == 'base':
            base = width
        gain = f'{100 * (1 - width/base):.3f}%' if base else '—'
        print(f'| {d} | {case} | {lo[0]:.17g} | {hi[0]:.17g} | {width:.17g} | {gain} | {lo[1]:.2f}/{hi[1]:.2f} | {lo[3]}/{hi[3]} | {max(lo[4], hi[4]):.2g} |')

for name in ('sample_build', 'sample_verify'):
    p = Path(name + '.log')
    if not p.exists():
        continue
    s = p.read_text()
    timer = re.search(r'^TIME .*$', s, re.M)
    if not timer or (name == 'sample_build' and 'Wrote sample.bin' not in s):
        continue
    verdict = re.search(r'^=== VERDICT: .*$', s, re.M)
    print(f'\n`{name}`: `{timer[0]}`' + (f' — `{verdict[0]}`' if verdict else '') + '.')
    if name == 'sample_verify' and 'exit=0' in timer[0] and verdict and 'PASS' in verdict[0]:
        cpu = sum(float(v) for v in re.findall(r'(?:user|sys)=(\S+)', timer[0]))
        rows = sum(map(int, re.findall(r'checked (\d+) states:', s)))
        print(f'\nGMP benchmark: {rows:,} printed row slots; {cpu:.2f} CPU seconds; '
              f'{cpu*1e9/rows:.1f} ns per printed slot (includes startup and comparison).')
        print('\nApproximate CPU hours per endpoint by original table size:')
        for label, n in [('D0=20 production', 1398838433), ('D0=22 production caps', 2155055415), ('D0=24 base', 1748953263),
                         ('D0=24 large', 3346876551), ('D0=25 large', 4060908685)]:
            print(f'- {label}: {cpu*n/rows/3600:.2f} h')

def width(case):
    hi, lo = [extract(Path(f'd12_{case}_m{m}.log')) for m in (1, 2)]
    return hi[0] - lo[0] if hi and lo else None


base, real = width('base'), width('l5-6')
if base is not None and real is not None:
    gain = 1 - real/base
    print(f'\nReal clipped L5 (D5=6): width change {gain:.3%} narrower versus D0=12 baseline '
          f'({"above" if gain > .3 else "below"} the 30% decision threshold).')
    for case, label in [('l5-6-no-inner-slack', 'all-inner'), ('l5-6-no-mixed-slack', 'mixed')]:
        candidate = width(case)
        if candidate is not None:
            print(f'- Suppressing {label} residual: width falls by {real - candidate:.17g} '
                  f'({1 - candidate/real:.3%} of the real L5 width; non-rigorous sensitivity).')
