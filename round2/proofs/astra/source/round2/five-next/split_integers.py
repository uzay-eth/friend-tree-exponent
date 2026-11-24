"""Verify both signed five-pass transcripts, including each rank size."""
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
SCALE = 2**40
RANKS = (
    (2021486872, 56279691, 70491064, 32031570, 1061159),
    (7312312, 6449232, 2975049, 976140, 52278),
    *((1177302, 381038, 2975049, 976140, 52278),) * 3,
)


def check(log):
    lines = [s for s in log.splitlines() if s and not s.startswith("wall ")]
    it = iter(lines)
    reference = ["side\tarchived\tpass0\tpass1\tpass2\tpass3\tpass4\treplacement\tincrement"]
    for side, file in (("lower", "occupation-five-lower.txt"), ("upper", "occupation-five-6.txt")):
        old = (HERE.parent / "outputs" / file).read_text()
        archived = int(re.search(rf"newborn reserve (\d+)/{SCALE} at discount 1133/2000\b", old)[1])
        passes = []
        for n, counts in enumerate(RANKS):
            if next(it) != f"{side} split pass {n}":
                raise ValueError("incorrect signed pass")
            for rank, count in zip(range(5, 0, -1), counts, strict=True):
                if next(it) != f" rank {rank} rows {count}":
                    raise ValueError(f"incorrect mesh for {side} pass {n} rank {rank}")
            match = re.fullmatch(rf"{side} quartet pass {n} newborn (\d+)/{SCALE}", next(it))
            if not match:
                raise ValueError("incorrect newborn scale or signed pass")
            passes.append(int(match[1]))
        total = sum(passes)
        if next(it) != f"{side} quartet replacement {total}/{SCALE}":
            raise ValueError("incorrect five-pass total")
        reference.append("\t".join(map(str, (side, archived, *passes, total, total-archived))))
    if next(it, None) is not None:
        raise ValueError("extra transcript rows")
    result = "\n".join(reference) + "\n"
    if (HERE / "split-reference.tsv").read_text() != result:
        raise ValueError("transcript differs from frozen reference")
    return result


if __name__ == "__main__":
    read = Path(sys.argv[1]).read_text()
    result = check(read)
    for before, after in (("36290118/", "43692218/"), ("26942162/", "26942163/"),
                          ("2021486872", "2021486871")):
        try:
            check(read.replace(before, after, 1))
        except ValueError:
            continue
        raise ValueError("altered transcript passed")
    print(result, end="")
    print("all 50 rank counts agree; re-added baseline, perturbed pass and removed row rejected")
