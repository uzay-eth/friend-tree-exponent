"""Original Gold2 lower witness, three physical occupation floors, and Farey arithmetic."""
from fractions import Fraction
from hashlib import sha256
import json
from math import ceil, floor
from pathlib import Path
import re

Q = Fraction
GOLD2 = Path(__file__).resolve().parents[1]
ROUND = GOLD2.parent
SCALE = 2**40


def pinned(path, digest):
    data = path.read_bytes()
    assert sha256(data).hexdigest() == digest, path
    return data.decode()


def normalized_floor(path, digest, pattern):
    matches = re.findall(pattern, pinned(path, digest), re.MULTILINE)
    assert len(matches) == 1, (path, matches)
    return SCALE * Q(matches[0])


def rationals_in_interval(lower, upper, max_denominator):
    return {Q(a, b) for b in range(1, max_denominator + 1)
            for a in range(ceil(b * lower), floor(b * upper) + 1)}


def main():
    receipt = json.loads(pinned(
        ROUND / 'attacks/neon2-quintet/gold2/production.json',
        '5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'))
    lower, upper = receipt['lower'], receipt['upper']
    original, ceiling = Q(lower['certified_exact']), Q(upper['certified_exact'])
    nominal = Q(lower['nominal'])
    error = lower['cap_f'] + lower['cap_g'] + lower['higher_per_core']
    assert error == lower['err_per_core'] == 19217
    assert original == nominal - Q(error, lower['m_min'])
    assert ceiling == Q(upper['nominal']) + Q(upper['err_per_core'], upper['m_min'])
    assert 0 < lower['m_min'] < lower['m_max'] == 1808603193693

    h3 = normalized_floor(
        GOLD2 / 'post-fy/threshold-lower.historical.log',
        '39d0d02575cb05df74776466d787cf8f8ef90399003ee0709e2e7c0e9d5e9191',
        r'^H3 normalized orientations .* total (\d+/\d+) ')
    small = normalized_floor(
        GOLD2 / 'edges/absorb12.log',
        'eded191d098bb2c9296c5d7d2861a051bc303bce122953c518ebf1981d22d521',
        r'^TOTAL (\d+/\d+) ')
    older = normalized_floor(
        GOLD2 / 'payer/README.md',
        '0020b89ff3d30dcf654225d7f37626d62c8f21b60312b7b2b342036faac6904f',
        r'Exact normalized value: `(\d+/\d+)`')
    assert older == Q(1046834122639903380009, 20000000000000000)
    assert h3 > 699000 and small > 395000 and older > 52000

    target, right = Q(99, 175), Q(142, 251)
    for total, expected in (
        (699000 + 395000 + 52000, Q(12789408689335247753, 22607539921162500000)),
        (h3 + small + older, Q(126905728962161948969119589499374909399,
                                224328301684105796116575000000000000000)),
    ):
        net = total - error
        assert net > 0
        endpoint = nominal + Q(net, lower['m_max'])
        margin = net - (target - nominal) * lower['m_max']
        assert endpoint == expected and margin >= Q(2737543234271, 87500000)
        assert target < endpoint <= ceiling < right
        print('lower endpoint', endpoint, 'margin at 99/175 in ticks/core', margin)

    assert right - ceiling == Q(3697635659207, 56786772132550000000)
    assert right.numerator * target.denominator - target.numerator * right.denominator == 1
    assert rationals_in_interval(original, ceiling, 425) == {target}
    assert not rationals_in_interval(endpoint, ceiling, 425)
    assert rationals_in_interval(endpoint, ceiling, 426) == {Q(241, 426)}
    print('upper endpoint', ceiling)
    print('PASS three SHA-256-pinned physical floors and original signed receipt; '
          'all reduced denominators through 425 excluded; 241/426 lies inside')


if __name__ == '__main__':
    main()
