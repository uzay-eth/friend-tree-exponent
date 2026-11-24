import contextlib
import io
import pathlib
import re
import sys
from fractions import Fraction

m = int(sys.argv[1])
src = pathlib.Path('round2/attacks/q3/hub_order_exact.py').read_text()
# Instrument the original display only; the original enumeration and arithmetic are unchanged.
src = re.sub(r'\{float\((lo|hi|d|e)\):\+?\.6f\}', r'{\1}', src)
s = io.StringIO()
with contextlib.redirect_stdout(s):
    exec(compile(src, 'hub_order_exact.py', 'exec'), {})
dp = {}
for line in pathlib.Path(sys.argv[2]).read_text().splitlines():
    t = line.split()
    if t[0] == 'm':
        use = int(t[1]) == m and t[3] == '0'
    elif use:
        dp[t[0], int(t[1])] = tuple(map(Fraction, t[4:]))
compared = 0
for line in s.getvalue().splitlines():
    match = re.match(r'\s+p_(\d+)( - p_\d+)? in \[([^,]+), ([^\]]+)\]', line)
    if match:
        u = int(match[1]); d = bool(match[2]); q = tuple(map(Fraction, match.groups()[2:]))
        ref = (Fraction(0), Fraction(0)) if d and u == 1 else dp['d' if d else 'p', max(u, 2)]
        assert q == ref, (line, ref)
        compared += 1
assert compared == 15
print(f'm={m}: all 30 exact endpoints match (8 marginals and 7 differences, including symmetry).')
