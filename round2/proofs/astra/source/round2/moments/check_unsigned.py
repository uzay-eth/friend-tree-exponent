"""Replay a nominal denominator from its declared source and data files."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
from tempfile import TemporaryDirectory


def main():
    root=Path(__file__).resolve().parents[2]
    out=root/'round2/outputs'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sheet',type=Path)
    parser.add_argument('--empty-root',action='store_true',help='replay a hard-linked copy of only the declared payload')
    args=parser.parse_args()
    path=args.sheet.resolve()
    with gzip.open(path,'rt') as stream:
        sheet=json.load(stream)
    side,K=sheet['side'],sheet['K']
    if (side not in ('lower','upper') or not sheet.get('score_only') or sheet.get('use_prior')
            or set(sheet['rows'])!={'M-per-core','M-per-birth'}):
        raise ValueError('expected two standalone unsigned objectives on the original interval')
    allowed={path,out/f'mc-{side}.bin',out/f'mc-mbar-{side}.i64',out/f'mc-mbar-{side}-tail.tsv'}
    for group in ('','score-','edge-score-'):
        allowed.add(out/f'moments-v5-{side}-{group}k{K}.tsv.gz')
    seen=set()
    base=Path(__file__).resolve().parent
    sources={base/f'{name}.py' for name in ('check_unsigned','v5','occupation','marginal_lift','dual','proofcheck','resources')}
    if sheet.get('tail_score'): sources.add(base/'score_tail.py')
    source_seen={Path(__file__).resolve()}
    if args.empty_root:
        cache=base/'build'
        cache.mkdir(exist_ok=True)
        with TemporaryDirectory(dir=cache,prefix='unsigned-') as tmp:
            clean=Path(tmp)
            for f in sources|allowed:
                dest=clean/('round2/outputs' if f in allowed else 'round2/moments')/f.name
                dest.parent.mkdir(parents=True,exist_ok=True)
                os.link(f,dest)
            command=[sys.executable,'-B','round2/moments/check_unsigned.py',f'round2/outputs/{path.name}']
            subprocess.run(command,cwd=clean,check=True)
        print(f'Empty-root replay used exactly {len(sources)+len(allowed)} files.',flush=True)
        return

    def check_file(event,args):
        if event!='open' or not isinstance(args[0],str): return
        f=Path(args[0]).resolve()
        if f==path or f.is_relative_to(out):
            if f not in allowed: raise PermissionError(f'outside unsigned replay inputs: {f.name}')
            seen.add(f)
        elif f.is_relative_to(root):
            if f.suffix=='.pyc': raise FileNotFoundError('read the pinned Python source')
            if f not in sources: raise PermissionError(f'outside unsigned replay sources: {f.name}')
            source_seen.add(f)

    sys.addaudithook(check_file)
    verifier=Path(__file__).with_name('v5.py')
    sys.argv=[str(verifier),'check',str(path)]
    runpy.run_path(str(verifier),run_name='__main__')
    if seen!=allowed: raise ValueError('unsigned replay did not consume its declared inputs')
    if source_seen!=sources: raise ValueError('source inventory differs from actual imports')
    print(f'Isolated replay passed with {len(sources)} repository Python sources, the price sheet '
          'and six nominal potential data files.',flush=True)
    print('Verified price sheet SHA-256 '+hashlib.sha256(path.read_bytes()).hexdigest(),flush=True)


if __name__=='__main__':
    main()
