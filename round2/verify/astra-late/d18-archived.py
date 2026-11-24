#!/usr/bin/env python3
"""Check the D18 high101 receipt and its endpoint arithmetic, without raw rows."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parents[2] / 'proofs/astra-late/source'
audit = run_path(str(root / 'round2/checks/high-hub-sol.py'))
release = audit['H'] / 'release'
manifest = {path: sha for sha, path in (line.split() for line in
            (release / 'release.sha256').read_text().splitlines())}
audit['physical_registrations']()
for i, side in enumerate(('lower', 'upper')):
    sign = 1 - 2*i
    path = release / (('lo', 'hi')[i] + '.log')
    assert audit['digest'](path) == manifest[str(path.relative_to(root))]
    bounds = audit['receipt'](path, sign)
    previous = audit['EDGE_TRIAL'][i] + sign*audit['BASE_NUM'][i]/audit['DEN'][i]
    endpoint = previous + sign*sum(bounds['floor', 'all'])/audit['DEN'][i]
    assert endpoint == (
        Q('14157243308580160981686611486065459691229739/25026934531690315946179968825000000000000000'),
        Q('44818646818747647275382808904477397643909921/79217319301603159409316399750000000000000000'),
    )[i]
    print(side, 'archived exact high101 endpoint', endpoint)
