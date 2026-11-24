"""Hostwise comparison of the authenticated unoriented five-vertex potential."""

from collections import defaultdict
from fractions import Fraction
from math import ceil, comb, floor
from pathlib import Path
import hashlib
import json
import os
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
R2 = HERE.parents[1]
sys.path.insert(0, str(R2 / "attacks/neon2-quintet/physical"))
import check as physical

RECEIPT_HASH = "d26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb"
SNAPSHOT_HASH = {
    "lower": "81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986",
    "upper": "8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383",
}


def authenticated_receipt(path):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != RECEIPT_HASH:
        raise ValueError("the production receipt does not match the reviewed bytes")
    return json.loads(data)


def extrema(values):
    return int(values.min()), int(values.max())


class Comparison:
    """Original allocation to hosts, with half each four-path at its interior ends."""

    def __init__(self, path, side, receipt):
        self.snapshot = physical.Payload(path)
        self.h = self.snapshot.h
        self.mode = self.h["mode"]
        self.side = side
        self.scale = self.h["S"]
        self.global_low = receipt["m_min"]
        self.global_high = receipt["m_max"]
        self.alpha = receipt["verify"]["alpha_ticks"]
        self.family_caps = receipt["verify"]["caps"]
        self.higher_cap = receipt["higher_per_core"]
        if self.mode != (2 if side == "lower" else 1):
            self.close()
            raise ValueError("snapshot side disagrees with its receipt")
        stat = os.fstat(self.snapshot.fd)
        with os.fdopen(os.dup(self.snapshot.fd), "rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
            self.raw = np.memmap(stream, dtype="<i8", mode="r", offset=512)
        if digest != SNAPSHOT_HASH[side] or digest != receipt["verify"]["sha256"]:
            self.close()
            raise ValueError("snapshot digest disagrees with the pinned production input")
        self.identity = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        self.path = path
        self.g_min, self.g_max = self.edge_extrema()
        self.higher = self.higher_extrema()

    def close(self):
        raw = getattr(self, "raw", None)
        if raw is not None:
            raw._mmap.close()
        os.close(self.snapshot.fd)

    def check_identity(self):
        stat = os.stat(self.path)
        actual = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        if actual != self.identity:
            raise ValueError("snapshot changed during the replay")

    def edge_extrema(self):
        layout = self.snapshot
        count = layout.ends[2].n
        start = layout.offset["g"]
        minimum = np.full(count, np.iinfo(np.int64).max, dtype=np.int64)
        maximum = np.full(count, np.iinfo(np.int64).min, dtype=np.int64)
        for i in range(count):
            row = self.raw[start:start + i + 1]
            np.minimum(minimum[:i+1], row, out=minimum[:i+1])
            np.maximum(maximum[:i+1], row, out=maximum[:i+1])
            lo, hi = extrema(row)
            minimum[i] = min(minimum[i], lo)
            maximum[i] = max(maximum[i], hi)
            start += i + 1
        assert start == layout.offset["h3"]
        return minimum, maximum

    def higher_extrema(self):
        layout = self.snapshot
        intervals = defaultdict(lambda: [0, 0])
        branch_factors = defaultdict(int)
        counts = {j: layout.ends[j].n for j in (3,4,5)}
        read = defaultdict(int)

        def add(family, state, coefficient, bounds):
            cell = intervals[family, state]
            cell[0] = min(cell[0], coefficient * bounds[0])
            cell[1] = max(cell[1], coefficient * bounds[1])
            branch_factors[family, state] = max(branch_factors[family, state], coefficient)

        for family, space, level, degree in (("h3","3m",3,2),("s3","S3",4,3),("s4","S4",5,4)):
            length = comb(counts[level] + degree - 1, degree)
            for state, index in layout.ins[space].items():
                start = layout.offset[family] + index * length
                add(family, state, 1, extrema(self.raw[start:start+length]))
                read[family] += length

        inside = {name: tuple(states) for name, states in layout.ins.items()}
        for family, blocks in layout.block.items():
            count = counts[4 if family == "p4" else 5]
            for key, start in blocks.items():
                if family == "b5":
                    length = comb(count + 1, 2) * count
                else:
                    length = comb(count + 1, 2) if key[0] == key[-1] else count * count
                bounds = extrema(self.raw[start:start+length])
                read[family] += length
                if family == "p4":
                    left, right = (inside["4i"][i] for i in key)
                    add(family, left, right[1] - 1, bounds)
                    add(family, right, left[1] - 1, bounds)
                elif family == "p5":
                    left, center, right = (inside["5i"][i] for i in key)
                    add(family, center, (left[1] - 1) * (right[1] - 1), bounds)
                else:
                    hub, arm = inside["5o"][key[0]], inside["5c"][key[1]]
                    add(family, hub, arm[1] - 1, bounds)

        hosts = defaultdict(lambda: [0, 0])
        caps = defaultdict(int)
        for (family, state), bounds in intervals.items():
            c = state[1]
            if family in ("h3", "p4", "p5"):
                mult = comb(c, 2)
            elif family == "b5":
                mult = c * comb(c-1, 2)
            else:
                mult = comb(c, 3 if family == "s3" else 4)
            hosts[state][0] += mult * bounds[0]
            hosts[state][1] += mult * bounds[1]
            caps[state] += mult * branch_factors[family, state] * max(0, self.family_caps[family])
        self.weight_counts = dict(read)
        self.host_cap = {state: min(value, self.higher_cap) for state, value in caps.items()}
        return hosts

    def edge_bounds(self, d, c):
        ends = self.snapshot.ends[2]
        cells = np.minimum(c, ends.k + 1)
        anchor = min(d, ends.d)
        indices = np.array([ends.ids[anchor, i] for i in range(1, min(d, ends.k + 1) + 1)])
        lo = self.g_min[indices[cells - 1]]
        hi = self.g_max[indices[cells - 1]]
        if d > ends.d:
            # The numerators stay below 2^56 on the finite table at d<2000.
            lo = ((d - ends.d) * self.g_min[ends.inf] + ends.d * lo) // d
            hi = -(-((d - ends.d) * self.g_max[ends.inf] + ends.d * hi) // d)
        return np.minimum(lo, 0), np.maximum(hi, 0)

    def finite(self, d):
        """Return recharged lower and unrecharged upper for c=1,...,d."""
        if not 2 <= d < self.h["DF"]:
            raise ValueError("degree is outside the finite singleton table")
        start = self.snapshot.offset["f"] + (d-2) * (d+3) // 2
        f = self.raw[start+1:start+d+1]
        c = np.arange(1, d+1, dtype=np.int64)
        g_lo, g_hi = self.edge_bounds(d, c)
        lo = f + g_lo + self.alpha * (c-2)
        hi = f + g_hi
        if d < self.h["d3m"]:
            for k in range(2, min(self.h["K3m"], d) + 1):
                bounds = self.higher.get((d, k), (0, 0))
                lo[k-1] += bounds[0]
                hi[k-1] += bounds[1]
        return lo, hi

    def far(self, c, *, saturated=False):
        """Joint envelope on D>=DF at fixed c, or every c'>=c."""
        h = self.h
        df = h["DF"]
        ends = self.snapshot.ends[2]
        q0 = Fraction(0) if saturated else max(Fraction(0), 1 - Fraction(c, df))
        if self.mode == 1:
            assert h["f1"] >= 0
            low_f = h["f0"] + q0 * h["f1"]
            high_f = h["f0"] + h["f1"]
        else:
            assert h["f1"] > 0 and h["f2"] >= 0
            def quadratic(q):
                return q * (h["f0"] - h["f1"] * q)
            vertex = min(Fraction(1), max(q0, Fraction(h["f0"], 2*h["f1"])))
            low_f = min(quadratic(q0), quadratic(1)) - Fraction(h["f2"], df)
            high_f = quadratic(vertex)
        kap = min(c, ends.k + 1)
        classes = range(kap, ends.k+2) if saturated else (kap,)
        low_inf = int(self.g_min[ends.inf])
        high_inf = int(self.g_max[ends.inf])
        weight = Fraction(ends.d, df)
        low_t = min(int(self.g_min[ends.ids[ends.d, k]]) for k in classes)
        high_t = max(int(self.g_max[ends.ids[ends.d, k]]) for k in classes)
        low_g = min(0, low_inf, low_inf + weight*(low_t - low_inf))
        high_g = max(0, high_inf, high_inf + weight*(high_t - high_inf))
        lower = max(self.global_low, floor(low_f + low_g + self.alpha*(c-2)))
        upper = min(self.global_high, ceil(high_f + high_g))
        return lower, upper

    def verify_finite_extrema(self, receipt):
        lows, highs = [], []
        for d in range(2, self.h["DF"]):
            lo, hi = self.finite(d)
            lows.append(int(lo.min()))
            highs.append(int(hi.max()))
        expected = receipt["verify"]["comp_fin"]
        assert min(lows) == expected["min_W"] and max(highs) == expected["max_W"]
        assert min(lows) >= self.global_low and max(highs) <= self.global_high
        return min(lows), max(highs)

