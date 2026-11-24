"""Verify the whole positive-power splice, the founder constant and beta telescopes."""
from fractions import Fraction as Q
from hashlib import sha1
from pathlib import Path

HUB = Path(__file__).resolve().parents[2]
PINS = {
    'all-powers/README.md': 'd0d8654bffacdf989e4a177c7174ffa549cd5100',
    'all-powers/check.py': '66b5049cd803efc502152f7dd4e4df0a6616e324',
    'mean-ledger/fractional/README.md': 'e813512d9b75df5b2273e2a726bc845e57261d66',
    'mean-ledger/fractional/check.py': '76bbdfc337be47cf8a0a91138ee5d191264017ab',
    'mean-ledger/fractional/star_bound.py': 'd0ff19f9a7400e75b1b9522b0b671a7fafe0bb4b',
    'fractional-middle/README.md': '5921ba52b80a3cfc8ccba1b16b6a77f66d7cfcfa',
    'fractional-middle/check.py': 'd8111c15376ad27f2fe0a583cabf8192f925ec88',
    'small-powers/README.md': 'c2ee87a3a4c3ec15595063abde1ec810261f0b9a',
    'small-powers/check.py': '40659d732606294bb3a10608cbfd69bfc0201c5d',
}
C = Q(9, 5000)


def main():
    for name, wanted in PINS.items():
        data = (HUB / name).read_bytes()
        assert sha1(f'blob {len(data)}\0'.encode() + data).hexdigest() == wanted, name

    boundaries = [Q(j, 512) for j in range(513)]
    assert boundaries[0] == 0 and boundaries[473] == Q(473, 512) and boundaries[-1] == 1
    assert min(Q(33, 1000), Q(1, 40), Q(406542520567, 2500000000000),
               Q(2085176012479, 12000000000000), Q(265041348327659, 2012097656250000),
               Q(1892360505093, 10**15), Q(22890923, 1058158080)) > C

    # A concave log ratio on [0,1] is bounded below by the smaller endpoint.
    leaf_hub_cap, leaf_mean_cap = Q(197, 305), Q(387, 1000)
    e_upper = Q(11, 4)
    assert min(2 / (e_upper * leaf_hub_cap),
               3 / (e_upper**2 * leaf_mean_cap)) > Q(51, 50)
    assert min(leaf_hub_cap, leaf_mean_cap) / 50 > C

    powers = (Q(1, 10**18), Q(1, 512), Q(9, 128), Q(1, 3), Q(1, 2),
              Q(473, 512), Q(1), Q(3, 2), Q(3), Q(7), Q(401))
    checks = 0
    for p in powers:
        beta = 1 / (1 + p)                  # B(1,p+1)
        step = 1 / ((p + 1) * (p + 2))    # Gamma(p+1)Gamma(1)/Gamma(3+p)
        total, first = Q(0), beta
        for n in range(1, 81):
            next_beta = n * beta / (n + p + 1)
            assert (beta - next_beta) / (p + 1) == step
            total += step
            assert (first - next_beta) / (p + 1) == total
            beta, step = next_beta, n * step / (n + p + 2)
            checks += 1
    print(f'PASS {len(PINS)} source pins, 512 fractional boxes joined, '
          f'all seed/star and founder bounds, {checks} exact beta telescopes')


if __name__ == '__main__':
    main()
