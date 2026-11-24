"""Check the gold2 drift after one actual discrete birth on all trees through size 12."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import file_digest, sha256
import json
from mmap import ACCESS_READ, mmap
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from size_epoch import core_state, enumerate_sizes


DIGEST = "5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284"


def check(reader, layers, events, record, direction):
    rate = Q(record["certified_exact"])
    assert reader.lam - direction*Q(record["err_per_core"],record["m_min"]) == rate
    values = []
    for layer in layers:
        current = {}
        for key, tree in layer.items():
            core, degree = core_state(tree)
            value = reader.eval_phi(core,degree)
            assert record["m_min"]*len(core) <= value <= record["m_max"]*len(core)
            current[key] = value
        values.append(current)

    seed, = layers[0]
    law = {seed: Q(1)}
    product = Q(1)
    faults = 0
    for step, transitions in enumerate(events):
        factor = 1+rate/(3+step)
        next_law = defaultdict(Q)
        for tree,(physical,uniform) in transitions.items():
            required = values[step][tree]*factor
            expected = sum(p*values[step+1][v] for v,p in physical.items())
            assert direction*(expected-required) >= 0
            incorrect = sum(p*values[step+1][v] for v,p in uniform.items())
            faults += direction*(incorrect-required) < 0
            for v,p in physical.items():
                next_law[v] += law[tree]*p
        law = next_law
        assert set(law)==set(layers[step+1]) and sum(law.values())==1
        product *= factor
        mean = sum(p*values[step+1][v] for v,p in law.items())
        assert direction*(mean-values[0][seed]*product)>=0
    return faults


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("directory",type=Path)
    work = ap.parse_args().directory
    source = work/"source/round2/neon2/quintet-cert"
    raw = (work/"gold.json").read_bytes()
    assert raw == (source/"certificate_gold2.json").read_bytes() and sha256(raw).hexdigest()==DIGEST
    gold = json.loads(raw)
    sys.path.insert(0,str(source.resolve()))
    from cert_check import PySnapshot

    layers,events = enumerate_sizes(12)
    failures = 0
    for side,direction in (("lower",1),("upper",-1)):
        with (work/f"quintet_{side}.bin").open("rb") as stream:
            assert file_digest(stream,"sha256").hexdigest()==gold[f"{side}_sha256"]
            with mmap(stream.fileno(),0,access=ACCESS_READ) as data:
                reader = PySnapshot(raw=data)
                faults = check(reader,layers,events,gold[side],direction)
                reader.buf.release()
                failures += faults
            print(side, "physical conditional rows", sum(map(len,events)),
                  "tree potential bounds", sum(map(len,layers)),
                  "product bounds", len(events), "wrong-law failures", faults,flush=True)
    assert failures > 0
    lo,hi = (Q(gold[side]["certified_exact"]) for side in ("lower","upper"))
    assert Q(142,251)-hi==Q(3697635659207,56786772132550000000)
    possible = {Q(num,den) for den in range(1,426)
                for num in range((lo*den).__ceil__(),(hi*den).__floor__()+1)}
    assert possible=={Q(99,175)}
    print("PASS gold2 size law through 12 and the sole rational through denominator 425")


if __name__=="__main__":
    main()
