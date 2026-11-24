"""Gold2 denominator exclusion using local potential ceilings and four physical rewards."""

import argparse
import hashlib
import importlib.util
import json
import re
from fractions import Fraction
from math import ceil, floor
from pathlib import Path

Q = Fraction
GOLD2 = Path(__file__).resolve().parents[2]
ROUND = GOLD2.parent
SCALE = 2**40
TARGET = Q(99, 175)
RECEIPT_PIN = '5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'
CELLS_PIN = '3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790'


def read_pinned(path, digest):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f'changed input {path}')
    return data


def match_fraction(data, pattern):
    values = re.findall(pattern, data.decode(), flags=re.MULTILINE)
    if len(values) != 1:
        raise ValueError(f'expected exactly one match for {pattern}')
    return Q(values[0])


def rounded_payments():
    singleton = read_pinned(GOLD2 / 'lifetime/README.md',
                            '17ad57b53161e61517a8ea776fa536abec39a76f967a1c57fe044b14fb418a5c')
    small = read_pinned(GOLD2 / 'edges/absorb12.log',
                        'eded191d098bb2c9296c5d7d2861a051bc303bce122953c518ebf1981d22d521')
    large = read_pinned(GOLD2 / 'payer/README.md',
                        '0020b89ff3d30dcf654225d7f37626d62c8f21b60312b7b2b342036faac6904f')
    sibling = read_pinned(GOLD2 / 'diagnostic/lower-sibling/compact.log',
                          'c3d7d979fe82c4ac626cbdbef5f461411e966abd4e53ab35529a9cd5bb75c02a')
    values = [
        match_fraction(singleton, r'^\| Conditioned floor \| \*\*([\d.]+)\*\*'),
        SCALE * match_fraction(small, r'^TOTAL (\d+/\d+) '),
        SCALE * match_fraction(large, r'Exact normalized value: `(\d+/\d+)`'),
        match_fraction(sibling, r'^exact sibling ticks/core (\d+/\d+) '),
    ]
    floors = [floor(v) for v in values]
    assert floors == [207990, 395261, 52341, 251635]
    return sum(floors)


def occupation_saving(cut, ceiling, local_upper, discount):
    """Floor the stopped lifetime of ceiling minus the host's local upper bound."""
    values = {}
    for degree in range(cut, 1, -1):
        for cores in range(degree, 0, -1):
            leaves = degree - cores
            saving = ceiling - local_upper[degree, cores]
            assert saving >= 0
            promotion = Q(leaves, degree)
            after_gain = values.get((degree + 1, cores), 0)
            after_promotion = values.get((degree, cores + 1), 0)

            def bellman_terms(receive):
                gain = leaves + receive
                numerator = saving + gain * after_gain + promotion * after_promotion
                return numerator, discount + gain + promotion

            ends = [bellman_terms(Q(0)), bellman_terms(Q(cores, 2))]
            current = min(floor(num / den) for num, den in ends)
            values[degree, cores] = current
            assert current >= 0
            for num, den in ends + [bellman_terms(Q(cores, 3))]:
                assert num >= den * current
    assert len(values) == cut * (cut + 1) // 2 - 1
    return values[2, 1]


def original_host_bounds(local_upper):
    path = GOLD2 / 'host/reader.py'
    read_pinned(path, '6be2e8aa5835ad0f1087179af5d84664d1176ed3ab071ca64f58af8d628fc910')
    spec = importlib.util.spec_from_file_location('gold2_host_repair', path)
    host = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(host)
    with host.open_host('lower') as original:
        for (degree, cores), upper in local_upper.items():
            assert original.comparison_ticks(degree, cores)[1] == upper
    print(f'PASS rebuilt {len(local_upper)} local ceilings from the authenticated Gold2 lower payload')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild-host', action='store_true')
    args = parser.parse_args()
    document = json.loads(read_pinned(ROUND / 'attacks/neon2-quintet/gold2/production.json', RECEIPT_PIN))
    cells = json.loads(read_pinned(GOLD2 / 'host/conditioned24.json', CELLS_PIN))
    assert cells['original_receipt_sha256'] == RECEIPT_PIN
    data = cells['results']['lower']['cells']
    local_upper = {(d, c): int(data[f'{d},{c-1}'][1])
                   for d in range(2, 25) for c in range(1, d + 1)}
    if args.rebuild_host:
        original_host_bounds(local_upper)

    lower = document['lower']
    lo, hi = (Q(document[side]['certified_exact']) for side in ('lower', 'upper'))
    ceiling, error, nominal = lower['m_max'], lower['err_per_core'], Q(lower['nominal'])
    assert ceiling == 1808603193693 and error == 19217 and nominal == Q(7071421, 12500000)
    net = rounded_payments() - error
    assert net == 888010 and net > 0
    for cut in (6, 18, 24):
        lifetime = occupation_saving(cut, ceiling, local_upper, hi)
        mean_ceiling = ceil(ceiling - lo * lifetime)
        endpoint = nominal + Q(net, mean_ceiling)
        margin = net - (TARGET - nominal) * mean_ceiling
        assert margin > {6: 4700, 18: 25700, 24: 27680}[cut]
        print(f'D<={cut}: V(2,1)={lifetime}, M<={mean_ceiling}, lower={endpoint}, '
              f'margin at 99/175={margin} ticks/core')
    right = Q(142, 251)
    assert TARGET < endpoint < Q(241, 426) < hi < right
    assert right.numerator * TARGET.denominator - right.denominator * TARGET.numerator == 1
    inside = {Q(a, b) for b in range(1, 427)
              for a in range(ceil(endpoint * b), floor(hi * b) + 1)}
    assert inside == {Q(241, 426)}
    print(f'PASS four integer reward floors, stationary potential bounds and exact '
          f'Farey exclusion through denominator 425; upper={hi}')


if __name__ == '__main__':
    main()
