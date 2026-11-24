import csv
import math
import statistics as st
from collections import Counter, defaultdict

mu = 0.56572366  # numerical input, not a certified value
groups = defaultdict(list)
for filename, count in [("sample7.csv", 100), ("sample6.csv", 400)]:
    with open(filename) as f:
        rows = list(csv.DictReader(f))
    sizes = ([100000, 1000000, 10000000] if filename == "sample7.csv"
             else [100000, 1000000])
    assert [(int(x['run']), int(x['n'])) for x in rows] == [
        (i, n) for i in range(count) for n in sizes
    ]
    for row in rows:
        groups[int(row["n"])].append({k: int(v) for k, v in row.items()})

for n, rows in sorted(groups.items()):
    vals = [r["L"] for r in rows]
    assert all(r["L"] <= r["M2"] <= r["L"] + 2 for r in rows)
    base = math.log(n) / math.log(math.log(n))
    mean = st.mean(vals)
    se = st.stdev(vals) / math.sqrt(len(vals))
    print(f"n={n} R={len(vals)} hist={dict(sorted(Counter(vals).items()))} "
          f"mean={mean:.3f} SE={se:.3f} median={st.median(vals):g} "
          f"leading={mu*base:.3f} ratio={mean/(mu*base):.3f} "
          f"loglog={math.log(math.log(n)):.3f}")
    print("  max-component means: " + ", ".join(
        f"N={q}: {st.mean(r[f'M{q}'] for r in rows):.3f}" for q in [2, 3, 4]
    ))
