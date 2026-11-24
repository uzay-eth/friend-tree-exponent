"""Pin the band release, with separate exact replays of the two splice corners."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha1, sha256
from pathlib import Path
from runpy import run_path

ROUND = Path(__file__).resolve().parents[6]
EXT = ROUND / 'open/age-order/global-transport/annealed/extension'
FABLE = ROUND / 'attacks/fable5'


def pin(path, expected):
    data = path.read_bytes()
    assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest() == expected, path


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=20)
    parser.add_argument('--stop', type=int, default=40)
    args = parser.parse_args()
    assert 20 <= args.start < args.stop <= 40
    for name, digest in (
        ('third-physical/band/README.md', 'cb7af0fb9f25c802d7c240a8a9abd8dc6ae853e6'),
        ('third-physical/band/check.py', 'e2aa75a98bf62981a91e33c3b3664e2f55d248c6'),
        ('third-physical/README.md', '27c6b6c7a4aa559e166835f3adc17fd189c989bb'),
        ('third-physical/check.py', 'e39f8bfdf66f19e34a59247fa3649ccaf884358b'),
        ('first-sender/joint/README.md', '01e03b7e09596745e0a65ac789f98e71d26b4941'),
        ('first-sender/joint/check.py', 'd3b4ad5bc87760a843712d9af43af76450a07a37'),
        ('first-sender/refined/README.md', 'd41bde01411ec05690c48774f0a9efd9a0316fd2'),
        ('first-sender/refined/check.py', '4959deec3b3e94eb323dcf372f55939599e728c9'),
        ('first-sender/check.py', '5680c8e937ebafdec09c2e193e6b36fe47ac84b1'),
    ):
        pin(EXT / name, digest)
    manifest = FABLE / 'release.sha256'
    pin(manifest, '8387da535626a8bf0de45a3a4a60e6c08ca62ebe')
    for line in manifest.read_text().splitlines():
        digest, name = line.split()
        assert sha256((ROUND.parent / name).read_bytes()).hexdigest() == digest, name
    print('PASS band, physical inputs and archived Fable audit pins', flush=True)

    band = run_path(str(EXT / 'third-physical/band/check.py'))
    if args.start == 20:
        assert band['serial'](20) > Q('.2124')
        print('PASS b=20 serial coefficient > .2124', flush=True)
    cuts = tuple(map(Q, ('1/3', '41/120', '7/20', '3/8', '2/5', '1/2')))
    for b in range(args.start, args.stop):
        points = cuts if b <= 21 else (cuts[0], cuts[-1])
        for left, right in zip(points, points[1:]):
            ok, mode = band['ledger'](b, Q('.2124' if b == 20 else '.211'), left, right)
            assert ok, (b, left, right, mode)
            if b == 20:
                print(f'PASS b=20 [{left},{right}]: {mode}', flush=True)
        if b != 20:
            print(f'PASS b={b} [1/3,1/2]', flush=True)

    joint = run_path(str(EXT / 'first-sender/joint/check.py'))
    c = joint['coefficient'](40, 256, 24)
    assert c == Q(2138962243701651396480805099, 10240000000000000000000000000)
    assert c > Q(261,1250)
    joint['corner'](40, Q(1,3), Q(261,1250), Q(61,100))
    print(f'PASS b>=40 corner: {c} > 261/1250; paired debit strictly smaller', flush=True)


if __name__ == '__main__':
    main()
