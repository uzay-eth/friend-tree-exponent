#!/usr/bin/env python3
"""Cold cumulative physical replay, or a short check on existing snapshots."""

import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from arithmetic import (SIDES, TRIAL, S, CEILING, PRIOR, require, read_ledgers,
                        enclosure, ticks, crosscheck, plants)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = Path("round2/outputs")
MOM = Path("round2/moments")
FIVE = Path("round2/five-next")
MOTIF = Path("round2/motif-cert")
# The public Markdown is formatted for reading. Restore the authenticated
# source bytes only within the disposable tree before verifying the manifest.
ARCHIVAL_MARKDOWN = {Path("round2/occupation/physical-test.md"), FIVE/"split-audit.md"}
FLAGS = ("-std=c++20", "-O3", "-Wall", "-Wextra", "-fno-fast-math",
         "-fsanitize=undefined", "-fno-sanitize-recover=all")
PY = (sys.executable, "-B")
SHEET_HASH = {"lower": "eebfce2a3fa90826e20e36cb204e744c4fafb7eb955f9d2de505c676f858daca",
              "upper": "2de08960db1ed39b5eecdc65066af50be22cc3badd9c2241683c374a06cd8a96"}


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def command(root, args, log=None, live=False):
    if log is None:
        subprocess.run(args, cwd=root, check=True)
        return
    with log.open("w") as stream:
        if live:
            with subprocess.Popen(args, cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, bufsize=1) as proc:
                with proc.stdout as output:
                    for line in output:
                        stream.write(line)
                        print(line.replace(str(root.parent)+"/", "scratch/"), end="", flush=True)
                proc.wait()
        else:
            proc = subprocess.run(args, cwd=root, stdout=stream, stderr=subprocess.STDOUT)
    text = log.read_text()
    if proc.returncode and not live:
        print(text, file=sys.stderr, end="")
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, args)
    return text


def build_tree(tree):
    files = [Path(line) for line in (HERE/"inputs.txt").read_text().splitlines()
             if line and not line.startswith("#")]
    require(len(files) == len(set(files)), "duplicate replay input")
    for path in files:
        require(not path.is_absolute() and ".." not in path.parts and
                path.suffix not in (".bin", ".i64", ".xz"), "invalid static replay input")
        src = ROOT/path
        if path in ARCHIVAL_MARKDOWN:
            src = src.with_suffix(".source.txt")
        dest = tree/path
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
    return files


def sheet_info(tree, side, hashes=None):
    out = tree/OUT
    file = out/f"moments-v5-unsigned-fine-{side}-dual.json.gz"
    require(digest(file) == SHEET_HASH[side], f"{side}: frozen unsigned dual digest")
    with gzip.open(file, "rt") as stream:
        sheet = json.load(stream)
    require((sheet["side"], sheet["K"], sheet["M"]) == (side, 7, 144) and
            sheet["score_only"] is True and sheet["use_prior"] is False and
            sheet["score_overlap"] is True and sheet["edge_score"] is True and
            {"M-per-core", "M-per-birth"} == set(sheet["rows"]) and
            not any(sheet.get(key) for key in ("edge_donor", "lower_score", "tip_edges",
                    "star_children", "original_score", "singleton_gap", "tail_score")),
            f"{side}: expected the reviewed unsigned original-prior sheet")
    y = Q(sheet["duals"]["M-per-birth"]["bound"])
    require(y == CEILING[side] and sheet["duals"]["M-per-birth"]["direction"] == -1,
            f"{side}: per-newborn ceiling mismatch")
    f = f"moments-v5-{side}-k7.tsv.gz"
    with gzip.open(out/f, "rt") as stream:
        head = next(stream).split()
    require(int(head[1]) == (2 if side == "lower" else 1) and
            Q(int(head[2]), int(head[3])) == TRIAL[side] and int(head[4]) == S,
            f"{side}: profile normalization")
    if hashes is not None:
        pins = sheet["inputs"]
        require(set(pins["score_inputs"]) == {f"mc-{side}.bin", f"mc-mbar-{side}.i64",
                f"mc-mbar-{side}-tail.tsv", f}, f"{side}: incorrect unsigned input set")
        for name, expected in {**pins["score_inputs"],
                f"moments-v5-{side}-score-k7.tsv.gz": pins["score_sha256"],
                f"moments-v5-{side}-edge-score-k7.tsv.gz": pins["edge_score_sha256"]}.items():
            actual = hashes.get(OUT/name)
            require((actual if actual is not None else digest(out/name)) == expected,
                    f"{side}: input digest mismatch for {name}")
    return y, file


def frozen_hashes(tree, skip_raw=False):
    result = {}
    for line in (tree/FIVE/"split-inputs.sha256").read_text().splitlines():
        expected, name = line.split()
        rel = Path(name)
        if skip_raw and rel in {OUT/f"mc-{side}.bin" for side in SIDES}:
            continue
        require(not rel.is_absolute() and ".." not in rel.parts and rel not in result,
                "invalid upstream manifest path")
        require(digest(tree/rel) == expected, f"split manifest mismatch: {rel}")
        result[rel] = expected
    return result


def split_integers(tree, file, work):
    return command(tree, [*PY, FIVE/"split_integers.py", file], work/"split-integers.log")


def occupation_loader():
    # The occupation-tail script itself imports a different module named replay.
    return ("import runpy, sys; from pathlib import Path; from types import SimpleNamespace; "
            "driver=runpy.run_path('round2/occupation-tail/replay.py'); "
            "driver['run'](SimpleNamespace(snapshots=None),Path(sys.argv[1])) "
            "if len(sys.argv)>1 else None")


def imports_and_prior(tree):
    command(tree, [*PY, "-c", "import sys; from fractions import Fraction as F; "
        "sys.path.insert(0,'round2/moments'); import check_tables, v5; "
        "assert (v5.LOW,v5.HIGH)==(F(113,200),F(1133,2000))"])


def raw_import_probes(tree):
    probes = ((MOTIF, "low_oracle"), (MOTIF, "tree_check"),
              (Path("round2/motif-high"), "physical_check"), (MOTIF, "mutations"))
    for path, module in probes:
        command(tree, [*PY, "-c", "import importlib, sys; from pathlib import Path; "
            "sys.path.insert(0,sys.argv[1]); entry=importlib.import_module(sys.argv[2]); "
            "assert Path(entry.__file__).resolve() == Path(sys.argv[1],sys.argv[2]+'.py').resolve()",
            path, module])
    print("mandatory imports PASS: low_oracle, tree_check, physical_check, mutations", flush=True)


def run_rebuild(tree, work):
    # Calling the existing run() in its own process retains its raw work outputs
    # for downstream reuse; the original CLI would delete them at the end.
    base = work/"occupation"
    base.mkdir()
    print("cold: regenerate raw snapshots, comparison profiles and mask63 tables", flush=True)
    physical = command(tree, [*PY, "-c", occupation_loader(), base],
        work/"physical.log", live=True)
    for side in SIDES:
        for name in (f"mc-{side}.bin", f"mc-mbar-{side}.i64", f"mc-mbar-{side}-tail.tsv"):
            os.link(base/"base"/name, tree/OUT/name)

    print("cold: replay checked five-pass quartet replacement", flush=True)
    command(tree, ["bash", FIVE/"replay-split.sh", "cold"], work/"split-cold.log", live=True)
    split = split_integers(tree, tree/OUT/"five-next-split-cold.log", work)

    print("cold: regenerate K7/144 tables and replay exact unsigned price sheets", flush=True)
    for tag in ("v5_profile", "v5_score", "v5_edge_score"):
        command(tree, ["c++", *FLAGS, MOM/(tag+".cpp"), "-o", work/tag])
    check_profile_tables(tree, work)
    ceilings = {}
    for side in SIDES:
        y, sheet = sheet_info(tree, side)
        command(tree, [*PY, HERE/"check_tables_k7.py", tree, side])
        command(tree, [*PY, MOM/"check_unsigned.py", sheet])
        ceilings[side] = y
    return physical, split, ceilings


def check_profile_tables(tree, work, sides=SIDES, tags=("v5_profile", "v5_score", "v5_edge_score")):
    for side in sides:
        p, q = TRIAL[side].as_integer_ratio()
        # The profile header binds the unreduced upper trial's integer pair.
        if side == "upper":
            p, q = 56581165, 100000000
        snapshot = OUT/f"mc-{side}.bin"
        for tag, prefix, args in (
                ("v5_profile", "", [OUT/f"mc-mbar-{side}.i64", "7", "144", str(p), str(q)]),
                ("v5_score", "score-", ["7"]), ("v5_edge_score", "edge-score-", ["7"])):
            if tag not in tags:
                continue
            path = work/"profile.tsv"
            with path.open("wb") as stream:
                subprocess.run([work/tag, snapshot, *args], cwd=tree, stdout=stream, check=True)
            with gzip.open(tree/OUT/f"moments-v5-{side}-{prefix}k7.tsv.gz", "rb") as frozen:
                require(path.read_bytes() == frozen.read(), f"{side}: regenerated {tag} differs")
            print(f"{side}: regenerated {tag} matches frozen K7/144 contents", flush=True)


def run_quick(work):
    print("quick: snapshot and input digests; frozen transcripts; exact joint arithmetic", flush=True)
    tree = work/"tree"
    build_tree(tree)
    for side in SIDES:
        for name in (f"mc-{side}.bin", f"mc-mbar-{side}.i64"):
            (tree/OUT/name).symlink_to(ROOT/OUT/name)
    hashes = frozen_hashes(tree)
    ceilings = {side: sheet_info(tree, side, hashes)[0] for side in SIDES}
    split = split_integers(tree, tree/OUT/"five-next-split.log", work)
    transcript = tree/OUT/"occupation-affine-verified.txt"
    require(digest(transcript) == "2b5580fde088969ecc9d8da988ca48fb503ba78dddb494ffcc1cb2d0e4f00766",
            "digest mismatch for reviewed mask63 replay")
    physical = transcript.read_text()
    return physical, split, ceilings


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=("cold", "quick", "preflight"), nargs="?", default="cold")
    args = ap.parse_args()
    require(__debug__, "run with Python assertions enabled")
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    os.environ.update({"MOMENTS_CPU": str(cpu), "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                       "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
                       "PYTHONDONTWRITEBYTECODE": "1"})
    for key in ("PYTHONPATH", "PYTHONOPTIMIZE", "FT2_TAIL_DUMP", "FIVE_NEXT_SPLIT_PINNED"):
        os.environ.pop(key, None)
    with TemporaryDirectory(prefix="joint-") as scratch:
        work = Path(scratch)
        os.environ["TMPDIR"] = str(work)
        if args.mode == "quick":
            inputs = run_quick(work)
        else:
            tree = work/"tree"
            files = build_tree(tree)
            imports_and_prior(tree)
            raw_import_probes(tree)
            frozen_hashes(tree, skip_raw=True)
            if args.mode == "preflight":
                for side in SIDES:
                    sheet_info(tree, side)
                command(tree, [*PY, "-c", occupation_loader()])
                print(f"preflight: isolated {len(files)} static inputs; no raw, mbar or executable copied")
                return
            inputs = run_rebuild(tree, work)
        ledgers = read_ledgers(*inputs)
        ends = enclosure(ledgers)
        reported = ticks(ends)
        crosscheck(ledgers, ends, reported)
        plants(ledgers, ends, reported)
        for side, endpoint in zip(SIDES, ends, strict=True):
            row = ledgers[side]
            v = sum(value for _, value in row["sources"])
            print(f"{side}: v={v}, e={row['loss']}, y={row['ceiling']}, lambda={row['trial']}, "
                  f"prior={PRIOR[0]}..{PRIOR[1]}, endpoint={endpoint}")
        l, u = reported
        print(f"{args.mode} PASS: outward physical interval [0.{l:012d}, 0.{u:012d}]")


if __name__ == "__main__":
    main()
