#!/usr/bin/env python3
"""Rebuild the two witnesses and run the four independent physical-bracket audits."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROUND = HERE.parent
ROOT = ROUND.parent
CHECKERS = ("motif_base_hy4.py", "inside_hy4.py",
            "extra_sources_hy4.py", "occupation_tail_hy4.py")


def audit(base, sandbox):
    layout = sandbox / "round2"
    checks = layout / "checks"
    outputs = layout / "outputs"
    checks.mkdir(parents=True)
    outputs.mkdir()
    for folder in ("motif", "motif-cert", "motif-high", "occupation", "occupation-tail"):
        (layout / folder).symlink_to(ROUND / folder, target_is_directory=True)
    manifest = json.loads((ROUND / "motif-cert/certificate.json").read_text())
    for side in ("upper", "lower"):
        raw = base / f"mc-{side}.bin"
        with raw.open("rb") as source:
            digest = hashlib.file_digest(source, "sha256").hexdigest()
        if digest != manifest[side]["sha256"]:
            raise ValueError(f"incorrect {side} witness")
        (outputs / raw.name).symlink_to(raw)
        profile = base / f"mc-mbar-{side}.i64"
        (outputs / profile.name).symlink_to(profile)
    (outputs / "mc-replay-clean").symlink_to(base, target_is_directory=True)
    (outputs / "mh-check").symlink_to(base / "high")
    for filename in CHECKERS:
        target = checks / filename
        shutil.copyfile(HERE / filename, target)
        runs = [()] if filename == "inside_hy4.py" else [("upper",), ("lower",)]
        for arguments in runs:
            subprocess.run([sys.executable, "-B", target, *arguments], cwd=sandbox, check=True)
    print("PASS: all four independent physical-bracket audits", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path,
                    help="use an existing complete two-witness motif-cert/replay.py work directory")
    args = ap.parse_args()
    if not __debug__:
        raise RuntimeError("enable Python assertions")
    if hasattr(os, "sched_getaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    outputs = ROUND / "outputs"
    outputs.mkdir(exist_ok=True)
    with TemporaryDirectory(prefix="physical-reviews-", dir=outputs) as directory:
        work = Path(directory)
        base = args.base.resolve() if args.base else work / "base"
        if not args.base:
            subprocess.run([sys.executable, "-B", ROUND / "motif-cert/replay.py", "--work", base],
                           cwd=ROOT, check=True)
        audit(base, work / "audit")


if __name__ == "__main__":
    main()
