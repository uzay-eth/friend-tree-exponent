"""Evaluate the checked potential after each possible attachment in every small history."""
import argparse
from fractions import Fraction
from pathlib import Path
import subprocess
import sys

BANDED = Path(__file__).resolve().parents[1] / "banded-upper"
sys.path.insert(0, str(BANDED / "independent"))
from test_banded_audit import BandedHistoryTree, BandedSpec


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("lookup", type=Path)
    cli.add_argument("recipe", type=Path)
    cli.add_argument("--max-size", type=int, default=9)
    args = cli.parse_args()
    if not 3 <= args.max_size <= 12:
        cli.error("choose a maximum size from 3 through 12")
    geometry = BANDED / "geometry.txt"
    bands = [list(map(int, row.split())) for row in geometry.read_text().splitlines() if row.strip()]
    spec = BandedSpec([(upper, caps) for upper, *caps in bands])
    start = BandedHistoryTree(spec, [[1, 2], [0], [0]])
    states = {}
    weights = []
    checked = 0

    def visit(history, collect):
        nonlocal checked
        def potential(h):
            total = 0
            for u in h.core_vertices():
                state = h.state(u)
                if collect:
                    states.setdefault(state, len(states))
                else:
                    total += weights[states[state]]
            return total

        phi = potential(history)
        children = [history.step(w) for w in range(len(history.adj))]
        after = [potential(h) for h in children]
        if not collect:
            # Sum over ringing vertices. Each of its neighbours is equiprobable.
            drift = sum(Fraction(sum(after[w] - phi for w in adj), len(adj))
                        for adj in history.adj)
            if drift * 9983 > 5657 * phi:
                raise ValueError(f"drift failure at {history.adj}: {drift} vs {phi}*5657/9983")
            checked += 1
        if len(history.adj) < args.max_size:
            for h in children:
                visit(h, collect)

    visit(start, True)
    requests = []
    for kind, degree, core_degree, labels in states:
        assert kind == "FIN"
        fields = [degree, core_degree]
        for _, d, k, e, j in labels:
            fields.extend((d, k, e, j))
        requests.append(" ".join(map(str, fields)))
    run = subprocess.run([args.lookup.resolve(), geometry, args.recipe],
                         input="\n".join(requests) + "\n", text=True,
                         stdout=subprocess.PIPE, check=True)
    for row in run.stdout.splitlines():
        if row.startswith("weight "):
            weights.append(int(row.split()[1]))
        else:
            print(row)
    if len(weights) != len(states):
        raise ValueError("incomplete lookup")
    visit(start, False)
    print(f"PASS banded: {checked} attachment histories through n={args.max_size}, "
          f"{len(states)} distinct weight queries, failures=0")


if __name__ == "__main__":
    main()
