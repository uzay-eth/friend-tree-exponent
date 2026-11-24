#!/usr/bin/env python3
"""Run an unchanged checker against an isolated tree with original source data."""
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory


def main():
    bundle = Path(__file__).resolve().parents[2] / 'proofs/astra-late'
    archive = bundle / 'authenticated'
    if len(sys.argv) < 2:
        raise SystemExit('usage: with-originals.py CHECKER.py [ARGS...]')
    script = Path(sys.argv[1])
    if script.is_absolute() or '..' in script.parts or script.suffix != '.py':
        raise SystemExit('checker must be a relative Python path in the source tree')
    with TemporaryDirectory(prefix='astra-late-') as scratch:
        source = Path(scratch) / 'source'
        shutil.copytree(bundle / 'source', source)
        for line in (archive / 'SHA256SUMS').read_text().splitlines():
            expected, name = line.split('  ', 1)
            relative = Path(name)
            if relative.is_absolute() or '..' in relative.parts or relative.suffix != '.md':
                raise SystemExit(f'invalid source-data path: {name}')
            data = (archive / relative.with_suffix('.md.txt')).read_bytes()
            if hashlib.sha256(data).hexdigest() != expected:
                raise SystemExit(f'archived source-data digest mismatch: {name}')
            dest = source / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
        return subprocess.call([sys.executable, '-IB', str(source / script), *sys.argv[2:]], cwd=source)


if __name__ == '__main__':
    raise SystemExit(main())
