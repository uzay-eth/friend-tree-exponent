"""Exact checks for the real-power splice and its two endpoint calculations."""
from fractions import Fraction as Q
from hashlib import sha1
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PINS = {
    "mellin/README.md": "9e75bd41e8a204b100a086d21e9fea76ad844fb0",
    "mellin/check.py": "9786bbdf7230ee86b9aa518d421d9fb471866130",
    "mellin/lower-1-3.tsv": "defa8793b129e4b3c3f5ef6f23fd4033fc4ff8fc",
    "mean-ledger/real/README.md": "3562e4601bf35825e317cb77bed6f694506c0e79",
    "mean-ledger/real/check.py": "27234fd45649df783a6fa5378fc43da3c2c1017d",
}
for name, expected in PINS.items():
    data = (ROOT / name).read_bytes()
    assert sha1(f"blob {len(data)}\0".encode() + data).hexdigest() == expected

c = Q(9, 5000)
assert min(Q(1, 40), Q(1892360505093, 10**15), Q(22890923, 1058158080)) > c
cases = 0
for p in (Q(1), Q(129,128), Q(9,8), Q(2), Q(3), Q(7,2), Q(7), Q(51,2), Q(1000)):
    # B(1,p+1)=1/(p+1); both beta branches then use rational recurrences.
    beta, increment = 1/(p+1), 1/((p+1)*(p+2))
    original = beta
    accumulated = Q()
    assert (1-Q(1,2)) * Q(1,2) * 2/((p+1)*(p+2)) == increment/2
    for n in range(1, 80):
        next_beta = n*beta/(n+p+1)
        assert (beta-next_beta)/(p+1) == increment
        accumulated += increment
        assert accumulated == (original-next_beta)/(p+1)
        beta = next_beta
        increment *= n/(n+p+2)
        cases += 1

print("PASS five writer pins, three uniform margins, nine founder checks and",
      cases, "exact beta telescopes")
