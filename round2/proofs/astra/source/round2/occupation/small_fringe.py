#!/usr/bin/env python3
"""Lower bounds on the occupation of motif surplus during a newborn's life."""

import argparse
import lzma
import shutil
import sys
import tempfile
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "verify" / "tree-drift"))
from motif_identity import Cert

UNIT = 10**15
NEW = (1, ())
NB = (1, 0, 1)
U = F(1133, 2000)
L = F(113, 200)


def load_reference(mode):
    name = "upper" if mode == 1 else "lower"
    src = REPO / "verify" / "motif-bracket" / "reference" / (name + ".bin.xz")
    with tempfile.NamedTemporaryFile() as file, lzma.open(src) as stream:
        shutil.copyfileobj(stream, file)
        file.flush()
        return Cert(file.name)


def state(tree):
    return tree[0], len(tree[1]), 1


def up(s):
    return s[0] + 1, s[1], s[2]


def bud(s):
    return s[0] - 1, s[1] + 1, s[2]


def gain(tree):
    return tree[0] + 1, tree[1]


@cache
def nverts(tree):
    return 1 + tree[0] + sum(nverts(ch) for ch in tree[1])


@cache
def inv(tree):
    return F(1, sum(state(tree)))


@cache
def intrinsic(tree):
    """A jump acts at one target. Parent input at the root is separate."""
    leaves, children = tree
    p = inv(tree)
    out = Counter()
    input_rate = leaves + sum(map(inv, children), F())
    if input_rate:
        out[gain(tree)] += input_rate
    if leaves:
        out[(leaves - 1, tuple(sorted(children + (NEW,))))] += leaves * p
    for i, child in enumerate(children):
        subrow = Counter(intrinsic(child))
        subrow[gain(child)] += p
        others = children[:i] + children[i + 1:]
        for new, rate in subrow.items():
            out[(leaves, tuple(sorted(others + (new,))))] += rate
    assert sum(out.values()) == nverts(tree) - p
    return out


class Apex:
    """Sum each certified motif surplus once at its shallowest vertex."""

    def __init__(self, cert):
        self.cert = cert
        self.delta = cache(self.delta)
        self.lookup = cache(self.lookup)
        self.vertex = cache(self.vertex)
        for key in ("F", "G", "HV", "HC"):
            setattr(cert, key, cache(getattr(cert, key)))

    def lookup(self, name, s):
        c = self.cert
        return F(getattr(c, name)[c.inner[s]], c.Sr)

    def delta(self, family, states, i):
        c = self.cert
        fn = getattr(c, family)
        base = fn(*states)
        s = states[i]
        def change(move):
            return fn(*(states[:i] + (move(s),) + states[i + 1:])) - base
        return change(up), change(bud) if s[0] else F()

    def vertex(self, s):
        c = self.cert
        assert c.isinner(s) and s[2] == 1
        kept = self.lookup("rpar", s) + s[1] * self.lookup("rch", s)
        expr = (s[0] + kept) * (c.F(up(s)) - c.F(s))
        if s[0]:
            expr += F(s[0], sum(s)) * (c.F(bud(s)) - c.F(s) + c.F(NB) + c.G(bud(s), NB))
        return c.lam * c.F(s) - expr

    def edge(self, root, child):
        c = self.cert
        sp, sw = state(root), state(child)
        # All occupied states in the stated cutoffs are interior. Every parent
        # in an edge lies in the vertical-triple table's middle support.
        assert c.isinner(sw) and sp in c.hu and (not sw[1] or sw in c.hu)
        dp, bp = self.delta("G", (sp, sw), 0)
        dw, bw = self.delta("G", (sp, sw), 1)
        expr = (inv(child) - self.lookup("rch", sp)) * (c.F(up(sp)) - c.F(sp))
        expr += (inv(root) - self.lookup("rpar", sw)) * (c.F(up(sw)) - c.F(sw))
        expr += (sp[0] + inv(child) + self.lookup("xpar", sp)) * dp
        if sp in c.cu:
            expr += (sp[1] - 1) * self.lookup("xcc", sp) * dp
        elif sp[1] >= 2:
            # These terms occur when the apex's horizontal table is inactive.
            removed = False
            for sibling in root[1]:
                if not removed and sibling == child:
                    removed = True
                    continue
                expr += inv(sibling) * dp - c.A(c.kcc, sp, inv(sibling))
            expr += (sp[1] - 1) * c.A(c.kcc, sp, inv(child))
        expr += (sw[0] + inv(root) + sw[1] * self.lookup("xch", sw)) * dw
        if sp[0]:
            expr += sp[0] * inv(root) * (bp + c.HC(bud(sp), sw, NB))
        if sw[0]:
            expr += sw[0] * inv(child) * (bw + c.HV(sp, bud(sw), NB))
        return c.lam * c.G(sp, sw) - expr

    def triple(self, family, verts, ext_inputs):
        c = self.cert
        states = tuple(map(state, verts))
        expr = F()
        for i, x in enumerate(verts):
            gl, gb = self.delta(family, states, i)
            expr += (x[0] + ext_inputs[i] + sum(map(inv, x[1]), F())) * gl
            expr += x[0] * inv(x) * gb
        x, y, z = verts
        sx, sy, sz = states
        if family == "HV":
            expr += (inv(x) - self.lookup("xpar", sy)) * (c.G(up(sy), sz) - c.G(sy, sz))
            expr += (inv(z) - self.lookup("xch", sy)) * (c.G(sx, up(sy)) - c.G(sx, sy))
        else:
            offset = self.lookup("xcc", sx)
            expr += (inv(y) - offset) * (c.G(up(sx), sz) - c.G(sx, sz))
            expr += (inv(z) - offset) * (c.G(up(sx), sy) - c.G(sx, sy))
        return c.lam * getattr(c, family)(*states) - expr

    def score(self, tree, ext):
        c = self.cert
        s = state(tree)
        surplus = [self.vertex(s)]
        for y in tree[1]:
            surplus.append(self.edge(tree, y))
            if state(y) in c.hu:
                for z in y[1]:
                    surplus.append(self.triple("HV", (tree, y, z), (ext, inv(tree), inv(y))))
        if s in c.cu:
            for y, z in combinations(tree[1], 2):
                surplus.append(self.triple("HC", (tree, y, z), (ext, inv(tree), inv(tree))))
        sign = 1 if c.mode == 1 else -1
        assert min(sign * d for d in surplus) >= 0
        return sign * sum(surplus)

    def comparison(self):
        c = self.cert
        # One incoming edge and vertical path per core and finitely many
        # horizontal paths at each supported apex.
        raw = sum(max(0, max(t)) for t in (c.f, c.g, c.hV))
        raw += (c.KCU * (c.KCU - 1) // 2) * max(0, max(c.hC))
        return F(raw, c.S)


def solve(cert, cutoff, discount):
    levels = [{NEW}]
    for _ in range(2, cutoff):
        levels.append({y for x in levels[-1] for y in intrinsic(x)})
    apex = Apex(cert)
    next_value = Counter()
    for layer in reversed(levels):
        table = {}
        for x in layer:
            n = nverts(x)
            continuation = sum((r * next_value[y] for y, r in intrinsic(x).items()), F())
            trial = []
            for p in (F(), F(1, 2)):
                source = apex.score(x, p)
                rhs = UNIT * source + continuation + p * next_value[gain(x)]
                trial.append(rhs // (discount + n - inv(x) + p))
            table[x] = max(0, min(trial))
        next_value = Counter(table)
    return sum(map(len, levels)), F(next_value[NEW], UNIT), apex.comparison()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=9, choices=range(2, 17))
    parser.add_argument("--mode", choices=("both", "upper", "lower"), default="both")
    args = parser.parse_args()
    for mode, side in ((1, "upper"), (2, "lower")):
        if args.mode not in ("both", side):
            continue
        cert = load_reference(mode)
        assert cert.mode == mode
        n, reward, bound = solve(cert, args.size, U)
        adjustment = L * reward / bound
        endpoint = cert.lam + (-adjustment if mode == 1 else adjustment)
        unit = 10**10
        outward = -(-endpoint * unit // 1) if mode == 1 else endpoint * unit // 1
        print(f"{side} size <= {args.size}, states = {n}, lower reward = {reward}, "
              f"score/C <= {bound}, endpoint = {outward}/{unit}", flush=True)


if __name__ == "__main__":
    main()
