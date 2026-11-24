"""Test the fixed-size drift at both exact endpoints on the full small-tree law."""

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import file_digest
from mmap import ACCESS_READ, mmap
from pathlib import Path
import sys


PINNED_SIDES = {
    "lower": (
        "81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986",
        Fraction(2825181729593748861, 4994114995300000000),
        99_882_299_906, 1_808_571_640_916, 19_186, 2,
    ),
    "upper": (
        "8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383",
        Fraction(439124781075521909, 776183857175000000),
        155_236_771_435, 1_808_494_729_937, 19_783, 1,
    ),
}


def unrooted_key(neighbours):
    def rooted(vertex, previous):
        return tuple(sorted(rooted(w, vertex) for w in neighbours[vertex] if w != previous))

    return min(rooted(v, -1) for v in range(len(neighbours)))


def enumerate_sizes(last_size):
    seed = ((1,), (0, 2), (1,))
    key = unrooted_key(seed)
    layers = [{key: seed}]
    transitions = []
    # Every attachment destination has positive probability. Grouping successor trees
    # by isomorphism preserves all directed-arrow multiplicities.
    for size in range(3, last_size):
        layer = {}
        rows = {}
        for key, tree in layers[-1].items():
            row = defaultdict(Fraction)
            uniform_recipient = defaultdict(Fraction)
            for recipient, adjacent in enumerate(tree):
                successor = list(tree)
                successor[recipient] = adjacent + (size,)
                successor.append((recipient,))
                successor = tuple(successor)
                dest = unrooted_key(successor)
                layer.setdefault(dest, successor)
                row[dest] += sum(Fraction(1, size * len(tree[donor])) for donor in adjacent)
                uniform_recipient[dest] += Fraction(1, size)
            assert sum(row.values()) == sum(uniform_recipient.values()) == 1
            rows[key] = row, uniform_recipient
        transitions.append(rows)
        layers.append(layer)
    assert [len(layer) for layer in layers] == [1, 2, 3, 6, 11, 23, 47, 106, 235, 551][:last_size - 2]
    return layers, transitions


def core_state(tree):
    core = [v for v, neighbours in enumerate(tree) if len(neighbours) >= 2]
    index = {v: i for i, v in enumerate(core)}
    adjacency = [[index[w] for w in tree[v] if w in index] for v in core]
    return adjacency, [len(tree[v]) for v in core]


def check_side(reader, side, layers, transitions):
    _, endpoint, min_ticks, max_ticks, error_ticks, mode = PINNED_SIDES[side]
    direction = 1 if side == "lower" else -1
    assert reader.mode == mode
    assert reader.lam - direction * Fraction(error_ticks, min_ticks) == endpoint
    values = []
    for layer in layers:
        current = {}
        for key, tree in layer.items():
            adjacency, degrees = core_state(tree)
            ticks = reader.eval_phi(adjacency, degrees)
            assert min_ticks * len(adjacency) <= ticks <= max_ticks * len(adjacency)
            current[key] = ticks
        values.append(current)

    key, = layers[0]
    seed_value = values[0][key]
    law = {key: Fraction(1)}
    product = Fraction(1)
    wrong_rows = []
    checked = 0
    for index, rows in enumerate(transitions):
        size = index + 3
        following = defaultdict(Fraction)
        factor = 1 + endpoint / size
        for key, (physical, uniform_recipient) in rows.items():
            required = factor * values[index][key]
            next_ticks = values[index + 1]
            mean = sum(prob * next_ticks[dest] for dest, prob in physical.items())
            assert direction * (mean - required) >= 0, (side, size, key)
            wrong = sum(prob * next_ticks[dest] for dest, prob in uniform_recipient.items())
            if direction * (wrong - required) < 0:
                wrong_rows.append((size, key))
            for dest, prob in physical.items():
                following[dest] += law[key] * prob
            checked += 1
        assert set(following) == set(layers[index + 1]) and sum(following.values()) == 1
        product *= factor
        mean = sum(prob * values[index + 1][dest] for dest, prob in following.items())
        assert direction * (mean - seed_value * product) >= 0
        law = following
    return checked, len(wrong_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pinned_source", type=Path, help="Neon-2 clean checkout at b04cc50")
    parser.add_argument("lower", type=Path)
    parser.add_argument("upper", type=Path)
    parser.add_argument("--last-size", type=int, choices=range(5, 13), default=12)
    args = parser.parse_args()
    sys.path.insert(0, str(args.pinned_source.resolve() / "round2/neon2/quintet-cert"))
    from cert_check import PySnapshot

    layers, transitions = enumerate_sizes(args.last_size)
    total_faults = 0
    for side, path in (("lower", args.lower), ("upper", args.upper)):
        with path.open("rb") as stream:
            assert file_digest(stream, "sha256").hexdigest() == PINNED_SIDES[side][0]
            with mmap(stream.fileno(), 0, access=ACCESS_READ) as snapshot:
                reader = PySnapshot(raw=snapshot)
                checked, faults = check_side(reader, side, layers, transitions)
                reader.buf.release()
            print(f"{side}: {checked} conditional fixed-size rows, {sum(map(len, layers))} "
                  f"potential comparisons, {len(transitions)} product bounds, "
                  f"{faults} violations with uniform-recipient mutation", flush=True)
            total_faults += faults
    assert total_faults > 0
    print("PASS exact physical law and endpoint products through size", args.last_size)


if __name__ == "__main__":
    main()
